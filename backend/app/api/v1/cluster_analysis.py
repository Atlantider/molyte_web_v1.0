"""
Cluster Analysis API - 统一的 Cluster 高级计算规划
支持 Binding、Desolvation、Redox、Reorganization 的统一规划和复用
"""
import logging
from typing import List, Optional, Dict, Any
from datetime import datetime
from collections import defaultdict

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from sqlalchemy import desc

from app.database import get_db
from app.models.user import User
from app.models.job import (
    MDJob, 
    AdvancedClusterJob as AdvancedClusterJobModel,
    AdvancedClusterJobStatus as AdvancedClusterJobStatusModel,
    ClusterCalcType as ClusterCalcTypeModel,
)
from app.models.result import SolvationStructure
from app.models.qc import QCJob, QCJobStatus, QCResult
from app.schemas.cluster_analysis import (
    ClusterCalcType,
    ClusterAnalysisPlanRequest,
    ClusterAnalysisPlanResponse,
    ClusterAnalysisSubmitRequest,
    AdvancedClusterJobResponse,
    PlannedQCTask,
    CalcTypeRequirements,
    AddCalcTypeRequest,
    AddCalcTypePlanResponse,
)
from app.dependencies import get_current_active_user

logger = logging.getLogger(__name__)
router = APIRouter()


def _get_md_job_or_404(db: Session, md_job_id: int, current_user: User) -> MDJob:
    """获取 MD 任务，检查权限"""
    md_job = db.query(MDJob).filter(MDJob.id == md_job_id).first()
    if not md_job:
        raise HTTPException(status_code=404, detail=f"MD 任务 {md_job_id} 不存在")
    if md_job.user_id != current_user.id and current_user.role.value != 'admin':
        raise HTTPException(status_code=403, detail="无权访问此 MD 任务")
    return md_job


def _get_selected_structures(
    db: Session, 
    md_job_id: int,
    solvation_structure_ids: Optional[List[int]] = None,
    composition_keys: Optional[List[str]] = None
) -> List[SolvationStructure]:
    """获取选中的溶剂化结构"""
    query = db.query(SolvationStructure).filter(
        SolvationStructure.md_job_id == md_job_id
    )
    
    if solvation_structure_ids:
        query = query.filter(SolvationStructure.id.in_(solvation_structure_ids))
    
    structures = query.all()
    
    # 如果指定了 composition_keys，进一步过滤
    if composition_keys:
        filtered = []
        for s in structures:
            comp = s.composition or {}
            key = "_".join(f"{k}_{v}" for k, v in sorted(comp.items()) if v > 0)
            if key in composition_keys or f"Li_{key}" in composition_keys:
                filtered.append(s)
        structures = filtered
    
    return structures


def _find_existing_qc_job(
    db: Session,
    smiles: str,
    charge: int,
    multiplicity: int,
    functional: str,
    basis_set: str,
    molecule_type: str = None
) -> Optional[QCJob]:
    """查找已有的、相同配置的 QC 任务"""
    query = db.query(QCJob).filter(
        QCJob.smiles == smiles,
        QCJob.charge == charge,
        QCJob.spin_multiplicity == multiplicity,
        QCJob.functional == functional,
        QCJob.basis_set == basis_set,
        QCJob.status == QCJobStatus.COMPLETED,
        QCJob.is_deleted == False
    )
    
    if molecule_type:
        query = query.filter(QCJob.molecule_type == molecule_type)
    
    # 获取最近完成的
    return query.order_by(desc(QCJob.finished_at)).first()


def _get_ligand_smiles_from_composition(composition: Dict[str, int]) -> List[Dict[str, Any]]:
    """从 composition 获取配体 SMILES 列表"""
    # 常见溶剂和阴离子的 SMILES 映射
    SMILES_MAP = {
        "EC": "C1COC(=O)O1",
        "DMC": "COC(=O)OC",
        "EMC": "CCOC(=O)OC",
        "DEC": "CCOC(=O)OCC",
        "PC": "CC1COC(=O)O1",
        "FEC": "FC1COC(=O)O1",
        "VC": "C=C1COC(=O)O1",
        "PF6": "F[P-](F)(F)(F)(F)F",
        "FSI": "FS(=O)(=O)[N-]S(=O)(=O)F",
        "TFSI": "FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F",
        "BF4": "F[B-](F)(F)F",
    }
    
    ligands = []
    for mol_name, count in composition.items():
        if count > 0 and mol_name != "Li":
            smiles = SMILES_MAP.get(mol_name, mol_name)
            ligands.append({
                "name": mol_name,
                "smiles": smiles,
                "count": count,
                "charge": -1 if mol_name in ["PF6", "FSI", "TFSI", "BF4"] else 0
            })
    return ligands


def _plan_qc_tasks_for_calc_type(
    db: Session,
    calc_type: ClusterCalcType,
    structures: List[SolvationStructure],
    qc_config: Dict[str, Any]
) -> CalcTypeRequirements:
    """为某个计算类型规划 QC 任务"""
    functional = qc_config.get("functional", "B3LYP")
    basis_set = qc_config.get("basis_set", "6-31G*")
    charge_ion = qc_config.get("charge_ion", 1)
    
    planned_tasks: List[PlannedQCTask] = []
    new_count = 0
    reused_count = 0
    
    # 收集所有需要的配体类型
    all_ligand_types = set()
    for s in structures:
        comp = s.composition or {}
        for mol_name, count in comp.items():
            if count > 0 and mol_name != "Li":
                all_ligand_types.add(mol_name)
    
    ligand_info = _get_ligand_smiles_from_composition(
        {name: 1 for name in all_ligand_types}
    )
    
    # 根据计算类型规划任务
    if calc_type in [ClusterCalcType.BINDING_TOTAL, ClusterCalcType.DESOLVATION_STEPWISE, 
                     ClusterCalcType.DESOLVATION_FULL]:
        # 这些计算需要 cluster 和 ligand 能量
        
        # 1. Cluster 任务 (每个结构一个)
        for s in structures:
            # Cluster 任务通常没有预先的 SMILES，需要从 xyz 生成
            task = PlannedQCTask(
                task_type="cluster",
                description=f"Cluster #{s.id} (CN={s.coordination_num})",
                structure_id=s.id,
                charge=charge_ion,  # Li+ 带 +1
                multiplicity=1,
                status="new"
            )
            planned_tasks.append(task)
            new_count += 1
        
        # 2. Ligand 任务 (每种配体一个，可复用)
        for lig in ligand_info:
            existing = _find_existing_qc_job(
                db, lig["smiles"], lig["charge"], 1, functional, basis_set, "ligand"
            )
            # task_type 使用 ligand_{name} 格式，以便结果计算时能识别配体
            ligand_task_type = f"ligand_{lig['name']}"
            if existing:
                task = PlannedQCTask(
                    task_type=ligand_task_type,
                    description=f"配体 {lig['name']}",
                    smiles=lig["smiles"],
                    charge=lig["charge"],
                    multiplicity=1,
                    status="reused",
                    existing_qc_job_id=existing.id,
                    existing_energy=existing.results[0].energy_au if existing.results else None
                )
                reused_count += 1
            else:
                task = PlannedQCTask(
                    task_type=ligand_task_type,
                    description=f"配体 {lig['name']}",
                    smiles=lig["smiles"],
                    charge=lig["charge"],
                    multiplicity=1,
                    status="new"
                )
                new_count += 1
            planned_tasks.append(task)
    
    if calc_type == ClusterCalcType.BINDING_TOTAL:
        # Binding 需要 ion 能量
        li_smiles = "[Li+]"
        existing = _find_existing_qc_job(db, li_smiles, 1, 1, functional, basis_set, "ion")
        if existing:
            task = PlannedQCTask(
                task_type="ion",
                description="Li+ 离子",
                smiles=li_smiles,
                charge=1,
                multiplicity=1,
                status="reused",
                existing_qc_job_id=existing.id,
                existing_energy=existing.results[0].energy_au if existing.results else None
            )
            reused_count += 1
        else:
            task = PlannedQCTask(
                task_type="ion",
                description="Li+ 离子",
                smiles=li_smiles,
                charge=1,
                multiplicity=1,
                status="new"
            )
            new_count += 1
        planned_tasks.append(task)

    if calc_type == ClusterCalcType.DESOLVATION_STEPWISE:
        # Stepwise 需要 cluster_minus 结构（每个结构每种配体一个）
        for s in structures:
            comp = s.composition or {}
            for mol_name, count in comp.items():
                if count > 0 and mol_name != "Li":
                    # task_type 使用 cluster_minus_{ligand_name} 格式
                    task = PlannedQCTask(
                        task_type=f"cluster_minus_{mol_name}",
                        description=f"Cluster #{s.id} 去除 {mol_name}",
                        structure_id=s.id,
                        charge=charge_ion,
                        multiplicity=1,
                        status="new"  # cluster_minus 很难复用，每个结构不同
                    )
                    planned_tasks.append(task)
                    new_count += 1

    if calc_type == ClusterCalcType.DESOLVATION_FULL:
        # Full desolvation 需要 ion 能量
        li_smiles = "[Li+]"
        existing = _find_existing_qc_job(db, li_smiles, 1, 1, functional, basis_set, "ion")
        if existing:
            task = PlannedQCTask(
                task_type="ion",
                description="Li+ 离子",
                smiles=li_smiles,
                charge=1,
                multiplicity=1,
                status="reused",
                existing_qc_job_id=existing.id,
                existing_energy=existing.results[0].energy_au if existing.results else None
            )
            reused_count += 1
        else:
            task = PlannedQCTask(
                task_type="ion",
                description="Li+ 离子",
                smiles=li_smiles,
                charge=1,
                multiplicity=1,
                status="new"
            )
            new_count += 1
        planned_tasks.append(task)

    if calc_type == ClusterCalcType.BINDING_PAIRWISE:
        # Pairwise binding: 需要 Li+、各配体、Li-配体二聚体
        li_smiles = "[Li+]"
        existing = _find_existing_qc_job(db, li_smiles, 1, 1, functional, basis_set, "ion")
        if existing:
            task = PlannedQCTask(
                task_type="ion",
                description="Li+ 离子",
                smiles=li_smiles,
                charge=1,
                multiplicity=1,
                status="reused",
                existing_qc_job_id=existing.id
            )
            reused_count += 1
        else:
            task = PlannedQCTask(
                task_type="ion", description="Li+ 离子", smiles=li_smiles,
                charge=1, multiplicity=1, status="new"
            )
            new_count += 1
        planned_tasks.append(task)

        # 各配体 (使用 ligand_{name} 格式)
        for lig in ligand_info:
            existing = _find_existing_qc_job(
                db, lig["smiles"], lig["charge"], 1, functional, basis_set, "ligand"
            )
            ligand_task_type = f"ligand_{lig['name']}"
            if existing:
                task = PlannedQCTask(
                    task_type=ligand_task_type, description=f"配体 {lig['name']}",
                    smiles=lig["smiles"], charge=lig["charge"], multiplicity=1,
                    status="reused", existing_qc_job_id=existing.id
                )
                reused_count += 1
            else:
                task = PlannedQCTask(
                    task_type=ligand_task_type, description=f"配体 {lig['name']}",
                    smiles=lig["smiles"], charge=lig["charge"], multiplicity=1, status="new"
                )
                new_count += 1
            planned_tasks.append(task)

        # Li-配体二聚体（需要新计算，使用 dimer_{name} 格式）
        for lig in ligand_info:
            task = PlannedQCTask(
                task_type=f"dimer_{lig['name']}",
                description=f"Li-{lig['name']} 二聚体",
                smiles=f"[Li+].{lig['smiles']}",
                charge=1 + lig["charge"],
                multiplicity=1,
                status="new"
            )
            planned_tasks.append(task)
            new_count += 1

    # 生成描述
    desc_map = {
        ClusterCalcType.BINDING_TOTAL: "总脱溶剂化能 (E_cluster - E_ion - ΣE_ligand)",
        ClusterCalcType.BINDING_PAIRWISE: "单分子-Li Binding (E_Li-X - E_Li - E_X)",
        ClusterCalcType.DESOLVATION_STEPWISE: "逐级去溶剂化能 (E_cluster - E_minus - E_ligand)",
        ClusterCalcType.DESOLVATION_FULL: "完全去溶剂化能 (E_cluster - E_ion - ΣE_ligand)",
        ClusterCalcType.REDOX: "氧化还原电位 (热力学循环)",
        ClusterCalcType.REORGANIZATION: "Marcus 重组能 (4点方案)",
    }

    return CalcTypeRequirements(
        calc_type=calc_type,
        description=desc_map.get(calc_type, str(calc_type)),
        required_qc_tasks=planned_tasks,
        new_tasks_count=new_count,
        reused_tasks_count=reused_count
    )


# ============================================================================
# API 端点
# ============================================================================

@router.post("/plan", response_model=ClusterAnalysisPlanResponse)
def plan_cluster_analysis(
    request: ClusterAnalysisPlanRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    规划 Cluster 高级计算

    根据选中的结构和计算类型，返回所需的 QC 任务列表和复用情况
    """
    # 检查 MD 任务
    md_job = _get_md_job_or_404(db, request.md_job_id, current_user)

    # 获取选中的结构
    structures = _get_selected_structures(
        db, request.md_job_id,
        request.solvation_structure_ids,
        request.composition_keys
    )

    if not structures:
        raise HTTPException(status_code=400, detail="未找到符合条件的溶剂化结构")

    # QC 配置
    qc_config = request.qc_config.model_dump() if request.qc_config else {}

    # 为每个计算类型规划 QC 任务
    calc_requirements = []
    total_new = 0
    total_reused = 0
    warnings = []

    for calc_type in request.calc_types:
        req = _plan_qc_tasks_for_calc_type(db, calc_type, structures, qc_config)
        calc_requirements.append(req)
        total_new += req.new_tasks_count
        total_reused += req.reused_tasks_count

        # 添加警告
        if calc_type == ClusterCalcType.REORGANIZATION:
            warnings.append("⚠️ 重组能计算量极大，每个物种需要多次优化，建议限制结构数量")
        if calc_type == ClusterCalcType.REDOX:
            warnings.append("⚠️ Redox 计算对方法/基组敏感，结果仅供参考")

    # 去重统计（同一个 QC 任务可能被多个计算类型复用）
    seen_tasks = set()
    unique_new = 0
    unique_reused = 0

    for req in calc_requirements:
        for task in req.required_qc_tasks:
            task_key = (task.task_type, task.smiles, task.structure_id, task.charge)
            if task_key not in seen_tasks:
                seen_tasks.add(task_key)
                if task.status == "reused":
                    unique_reused += 1
                else:
                    unique_new += 1

    # 预估计算时间（粗略估计）
    estimated_hours = unique_new * 0.5  # 每个新任务约 0.5 小时

    return ClusterAnalysisPlanResponse(
        md_job_id=request.md_job_id,
        selected_structures_count=len(structures),
        selected_structure_ids=[s.id for s in structures],
        calc_requirements=calc_requirements,
        total_new_qc_tasks=unique_new,
        total_reused_qc_tasks=unique_reused,
        estimated_compute_hours=estimated_hours,
        warnings=warnings
    )


@router.post("/submit", response_model=AdvancedClusterJobResponse)
def submit_cluster_analysis(
    request: ClusterAnalysisSubmitRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    提交 Cluster 高级计算任务

    创建 AdvancedClusterJob 并标记为 SUBMITTED，等待 Worker 处理
    """
    # 检查 MD 任务
    md_job = _get_md_job_or_404(db, request.md_job_id, current_user)

    # 获取选中的结构
    structures = _get_selected_structures(
        db, request.md_job_id,
        request.solvation_structure_ids,
        request.composition_keys
    )

    if not structures:
        raise HTTPException(status_code=400, detail="未找到符合条件的溶剂化结构")

    # QC 配置
    qc_config = request.qc_config.model_dump() if request.qc_config else {}

    # 规划 QC 任务
    all_planned_tasks = []
    reused_qc_jobs = []

    for calc_type in request.calc_types:
        req = _plan_qc_tasks_for_calc_type(db, calc_type, structures, qc_config)
        for task in req.required_qc_tasks:
            task_dict = task.model_dump()
            task_dict["calc_type"] = calc_type.value
            all_planned_tasks.append(task_dict)
            if task.existing_qc_job_id:
                reused_qc_jobs.append(task.existing_qc_job_id)

    # 创建任务
    job = AdvancedClusterJobModel(
        md_job_id=request.md_job_id,
        user_id=current_user.id,
        status=AdvancedClusterJobStatusModel.SUBMITTED,
        calc_types=[ct.value for ct in request.calc_types],
        selected_structures={
            "solvation_structure_ids": [s.id for s in structures],
            "count": len(structures)
        },
        qc_config=qc_config,
        qc_task_plan={
            "planned_qc_tasks": all_planned_tasks,
            "reused_qc_jobs": list(set(reused_qc_jobs)),
            "new_qc_jobs": [],
            "total_qc_tasks": len(all_planned_tasks),
            "completed_qc_tasks": len([t for t in all_planned_tasks if t["status"] == "reused"])
        },
        results={},
        progress=0.0
    )

    db.add(job)
    db.commit()
    db.refresh(job)

    logger.info(f"Created AdvancedClusterJob {job.id} for MD job {request.md_job_id}, "
                f"calc_types={request.calc_types}, structures={len(structures)}")

    return job


@router.get("/jobs", response_model=List[AdvancedClusterJobResponse])
def list_cluster_analysis_jobs(
    md_job_id: Optional[int] = Query(None, description="筛选指定 MD 任务"),
    skip: int = Query(0, ge=0),
    limit: int = Query(50, ge=1, le=100),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取 Cluster 高级计算任务列表"""
    query = db.query(AdvancedClusterJobModel)

    if current_user.role.value != 'admin':
        query = query.filter(AdvancedClusterJobModel.user_id == current_user.id)

    if md_job_id:
        query = query.filter(AdvancedClusterJobModel.md_job_id == md_job_id)

    jobs = query.order_by(desc(AdvancedClusterJobModel.created_at)).offset(skip).limit(limit).all()
    return jobs


@router.get("/jobs/{job_id}", response_model=AdvancedClusterJobResponse)
def get_cluster_analysis_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取单个 Cluster 高级计算任务详情"""
    job = db.query(AdvancedClusterJobModel).filter(AdvancedClusterJobModel.id == job_id).first()

    if not job:
        raise HTTPException(status_code=404, detail=f"任务 {job_id} 不存在")

    if job.user_id != current_user.id and current_user.role.value != 'admin':
        raise HTTPException(status_code=403, detail="无权访问此任务")

    return job


@router.post("/jobs/{job_id}/add-calc-types")
def add_calc_types_to_job(
    job_id: int,
    request: AddCalcTypeRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    追加计算类型到已有任务

    支持在已完成 Binding 后追加 Desolvation 等
    """
    job = db.query(AdvancedClusterJobModel).filter(AdvancedClusterJobModel.id == job_id).first()

    if not job:
        raise HTTPException(status_code=404, detail=f"任务 {job_id} 不存在")

    if job.user_id != current_user.id and current_user.role.value != 'admin':
        raise HTTPException(status_code=403, detail="无权访问此任务")

    # 检查状态
    if job.status not in [AdvancedClusterJobStatusModel.COMPLETED,
                          AdvancedClusterJobStatusModel.CREATED]:
        raise HTTPException(status_code=400, detail="只能在任务完成或创建状态下追加计算类型")

    # 获取已选结构
    structure_ids = job.selected_structures.get("solvation_structure_ids", [])
    structures = db.query(SolvationStructure).filter(
        SolvationStructure.id.in_(structure_ids)
    ).all()

    # 规划新的 QC 任务
    new_requirements = []
    for calc_type in request.additional_calc_types:
        if calc_type.value in job.calc_types:
            continue  # 跳过已有的计算类型
        req = _plan_qc_tasks_for_calc_type(db, calc_type, structures, job.qc_config)
        new_requirements.append(req)

    # 检查与已有 QC 任务的复用
    existing_plan = job.qc_task_plan or {}
    existing_qc_jobs = set(existing_plan.get("reused_qc_jobs", []) +
                          existing_plan.get("new_qc_jobs", []))

    new_tasks_count = 0
    reused_from_existing = 0

    for req in new_requirements:
        for task in req.required_qc_tasks:
            if task.existing_qc_job_id and task.existing_qc_job_id in existing_qc_jobs:
                reused_from_existing += 1
            elif task.status == "new":
                new_tasks_count += 1

    return AddCalcTypePlanResponse(
        job_id=job_id,
        existing_calc_types=job.calc_types,
        additional_calc_types=[ct.value for ct in request.additional_calc_types],
        new_qc_tasks_required=new_tasks_count,
        reused_from_existing=reused_from_existing,
        details=new_requirements
    )


# ============================================================================
# 结果计算和查询 API
# ============================================================================

@router.post("/jobs/{job_id}/calculate")
def calculate_cluster_analysis_results(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    计算 Cluster 高级计算结果

    当所有 QC 任务完成后，计算各种能量结果
    """
    job = db.query(AdvancedClusterJobModel).filter(AdvancedClusterJobModel.id == job_id).first()

    if not job:
        raise HTTPException(status_code=404, detail=f"任务 {job_id} 不存在")

    if job.user_id != current_user.id and current_user.role.value != 'admin':
        raise HTTPException(status_code=403, detail="无权访问此任务")

    # 获取所有相关的 QC 结果
    # 方法1: 通过 cluster_analysis_job_id 直接查询关联的 QC 任务
    from app.models.qc import QCJob as QCJobModel, QCResult

    qc_jobs = db.query(QCJobModel).filter(
        QCJobModel.cluster_analysis_job_id == job_id,
        QCJobModel.status == QCJobStatus.COMPLETED,
        QCJobModel.is_deleted == False
    ).all()

    # 构建 task_type -> energy 的映射
    qc_results_by_type = {}
    for qc_job in qc_jobs:
        if qc_job.results:
            result = qc_job.results[0] if qc_job.results else None
            if result and result.energy_au is not None:
                key = (qc_job.task_type, qc_job.solvation_structure_id)
                qc_results_by_type[key] = {
                    "qc_job_id": qc_job.id,
                    "task_type": qc_job.task_type,
                    "structure_id": qc_job.solvation_structure_id,
                    "energy_au": result.energy_au,
                    "energy_ev": result.energy_au * 27.2114,
                    "homo": result.homo,
                    "lumo": result.lumo,
                    "homo_lumo_gap": result.homo_lumo_gap,
                }

    # 方法2: 同时从 qc_task_plan 中获取复用的 QC 任务结果
    qc_task_plan = job.qc_task_plan or {}
    reused_qc_job_ids = qc_task_plan.get("reused_qc_jobs", [])
    if reused_qc_job_ids:
        reused_results = db.query(QCResult).filter(QCResult.qc_job_id.in_(reused_qc_job_ids)).all()
        for r in reused_results:
            # 从 planned_tasks 中找到对应的 task_type
            for task in qc_task_plan.get("planned_qc_tasks", []):
                if task.get("existing_qc_job_id") == r.qc_job_id:
                    key = (task.get("task_type"), task.get("structure_id"))
                    if key not in qc_results_by_type:
                        qc_results_by_type[key] = {
                            "qc_job_id": r.qc_job_id,
                            "task_type": task.get("task_type"),
                            "structure_id": task.get("structure_id"),
                            "energy_au": r.energy_au,
                            "energy_ev": r.energy_au * 27.2114 if r.energy_au else None,
                            "homo": r.homo,
                            "lumo": r.lumo,
                            "homo_lumo_gap": r.homo_lumo_gap,
                        }

    logger.info(f"Cluster Analysis {job_id}: 找到 {len(qc_results_by_type)} 个 QC 结果")

    # 计算各类型结果
    results = {}
    calc_types = job.calc_types or []

    for calc_type in calc_types:
        try:
            if calc_type == "BINDING_TOTAL":
                results["BINDING_TOTAL"] = _calculate_binding_total(job, qc_results_by_type, db)
            elif calc_type == "BINDING_PAIRWISE":
                results["BINDING_PAIRWISE"] = _calculate_binding_pairwise(job, qc_results_by_type, db)
            elif calc_type == "DESOLVATION_STEPWISE":
                results["DESOLVATION_STEPWISE"] = _calculate_desolvation_stepwise(job, qc_results_by_type, db)
            elif calc_type == "DESOLVATION_FULL":
                results["DESOLVATION_FULL"] = _calculate_desolvation_full(job, qc_results_by_type, db)
            elif calc_type == "REDOX":
                results["REDOX"] = _calculate_redox(job, qc_results_by_type, db)
            elif calc_type == "REORGANIZATION":
                results["REORGANIZATION"] = _calculate_reorganization(job, qc_results_by_type, db)
        except Exception as e:
            logger.error(f"计算 {calc_type} 结果失败: {e}")
            results[calc_type] = {"error": str(e)}

    # 更新任务结果
    job.results = results
    job.status = AdvancedClusterJobStatusModel.COMPLETED
    job.progress = 100.0
    job.finished_at = datetime.now()
    db.commit()

    return {"status": "ok", "results": results}


def _calculate_binding_total(job, qc_results_by_type: Dict, db: Session) -> Dict:
    """计算总 Binding Energy

    qc_results_by_type: Dict[(task_type, structure_id), result_dict]
    """
    # E_bind = E_cluster - (E_ion + Σ n_j × E_ligand_j)

    # 从 qc_results_by_type 中提取能量
    e_ion = None
    ligand_energies = {}
    cluster_results = []  # 可能有多个 cluster（多个结构）

    for (task_type, structure_id), result in qc_results_by_type.items():
        if task_type == "ion":
            e_ion = result["energy_au"]
        elif task_type and task_type.startswith("ligand_"):
            ligand_name = task_type.replace("ligand_", "")
            ligand_energies[ligand_name] = result["energy_au"]
        elif task_type == "cluster" and structure_id:
            cluster_results.append({
                "structure_id": structure_id,
                "energy_au": result["energy_au"]
            })

    if not cluster_results:
        return {"error": "缺少 cluster 能量"}
    if e_ion is None:
        return {"error": "缺少 ion 能量"}
    if not ligand_energies:
        return {"error": "缺少配体能量"}

    # 计算每个 cluster 的 binding energy
    binding_results = []
    total_ligand_energy = sum(ligand_energies.values())

    for cluster in cluster_results:
        e_cluster = cluster["energy_au"]
        e_bind_au = e_cluster - (e_ion + total_ligand_energy)
        binding_results.append({
            "structure_id": cluster["structure_id"],
            "e_cluster_au": e_cluster,
            "e_bind_au": e_bind_au,
            "e_bind_ev": e_bind_au * 27.2114,
            "e_bind_kcal_mol": e_bind_au * 627.509,
        })

    # 计算平均值
    avg_bind_au = sum(r["e_bind_au"] for r in binding_results) / len(binding_results)

    return {
        "e_ion_au": e_ion,
        "ligand_energies_au": ligand_energies,
        "cluster_binding_results": binding_results,
        "average_e_bind_au": avg_bind_au,
        "average_e_bind_ev": avg_bind_au * 27.2114,
        "average_e_bind_kcal_mol": avg_bind_au * 627.509,
    }


def _calculate_binding_pairwise(job, qc_results_by_type: Dict, db: Session) -> Dict:
    """计算分子-Li Pairwise Binding Energy

    qc_results_by_type: Dict[(task_type, structure_id), result_dict]
    """
    # E_bind = E(Li-X) - E(Li) - E(X)

    e_ion = None
    ligand_energies = {}
    dimer_energies = {}

    for (task_type, structure_id), result in qc_results_by_type.items():
        if task_type == "ion":
            e_ion = result["energy_au"]
        elif task_type and task_type.startswith("ligand_"):
            ligand_name = task_type.replace("ligand_", "")
            ligand_energies[ligand_name] = result["energy_au"]
        elif task_type and task_type.startswith("dimer_"):
            ligand_name = task_type.replace("dimer_", "")
            dimer_energies[ligand_name] = result["energy_au"]

    pairwise_results = []
    for ligand_name, e_dimer in dimer_energies.items():
        e_ligand = ligand_energies.get(ligand_name)
        if e_dimer is not None and e_ligand is not None and e_ion is not None:
            e_bind = e_dimer - e_ion - e_ligand
            pairwise_results.append({
                "ligand": ligand_name,
                "e_dimer_au": e_dimer,
                "e_ligand_au": e_ligand,
                "e_ion_au": e_ion,
                "e_bind_au": e_bind,
                "e_bind_ev": e_bind * 27.2114,
                "e_bind_kcal_mol": e_bind * 627.509,
            })

    return {"pairwise_bindings": pairwise_results}


def _calculate_desolvation_stepwise(job, qc_results_by_type: Dict, db: Session) -> Dict:
    """计算逐级去溶剂化能

    qc_results_by_type: Dict[(task_type, structure_id), result_dict]
    """
    # ΔE_i = E_cluster - (E_minus_i + E_ligand_i)

    # 按 structure_id 分组
    cluster_energies = {}  # structure_id -> energy
    ligand_energies = {}   # ligand_name -> energy
    minus_energies = {}    # (structure_id, ligand_name) -> energy

    for (task_type, structure_id), result in qc_results_by_type.items():
        if task_type == "cluster" and structure_id:
            cluster_energies[structure_id] = result["energy_au"]
        elif task_type and task_type.startswith("ligand_"):
            ligand_name = task_type.replace("ligand_", "")
            ligand_energies[ligand_name] = result["energy_au"]
        elif task_type and task_type.startswith("cluster_minus_") and structure_id:
            ligand_name = task_type.replace("cluster_minus_", "")
            minus_energies[(structure_id, ligand_name)] = result["energy_au"]

    stepwise_results = []
    for (structure_id, ligand_name), e_minus in minus_energies.items():
        e_cluster = cluster_energies.get(structure_id)
        e_ligand = ligand_energies.get(ligand_name)
        if e_cluster is not None and e_ligand is not None:
            delta_e = e_cluster - (e_minus + e_ligand)
            stepwise_results.append({
                "structure_id": structure_id,
                "ligand": ligand_name,
                "e_cluster_au": e_cluster,
                "e_minus_au": e_minus,
                "e_ligand_au": e_ligand,
                "delta_e_au": delta_e,
                "delta_e_ev": delta_e * 27.2114,
                "delta_e_kcal_mol": delta_e * 627.509,
            })

    return {"stepwise_desolvation": stepwise_results}


def _calculate_desolvation_full(job, qc_results_by_type: Dict, db: Session) -> Dict:
    """计算完全去溶剂化能"""
    # ΔE = E_cluster - (E_ion + Σ E_ligand_i)
    # 实际上与 BINDING_TOTAL 相同
    return _calculate_binding_total(job, qc_results_by_type, db)


def _calculate_redox(job, qc_results_by_type: Dict, db: Session) -> Dict:
    """计算氧化还原电位

    qc_results_by_type: Dict[(task_type, structure_id), result_dict]

    Redox 计算需要 4 个单点：neutral_gas, charged_gas, neutral_sol, charged_sol
    task_type 格式: redox_{smiles}_{state} 其中 state = neutral_gas/charged_gas/neutral_sol/charged_sol
    """
    # 按 SMILES 分组
    species_results = {}

    for (task_type, structure_id), result in qc_results_by_type.items():
        if task_type and task_type.startswith("redox_"):
            # 解析 task_type: redox_{smiles}_{state}
            parts = task_type.split("_", 2)  # ["redox", smiles, state]
            if len(parts) >= 3:
                smiles = parts[1]
                state = parts[2]  # neutral_gas, charged_gas, etc.
                if smiles not in species_results:
                    species_results[smiles] = {}
                species_results[smiles][state] = result["energy_au"]

    # 计算电位
    redox_results = []
    for smiles, energies in species_results.items():
        e_neutral_gas = energies.get("neutral_gas")
        e_charged_gas = energies.get("charged_gas")
        e_neutral_sol = energies.get("neutral_sol")
        e_charged_sol = energies.get("charged_sol")

        if all([e_neutral_gas, e_charged_gas, e_neutral_sol, e_charged_sol]):
            # 简化计算（实际需要更复杂的热力学校正）
            delta_g_sol = (e_charged_sol - e_neutral_sol) * 27.2114  # eV
            # 氧化电位 (vs SHE, 假设 SHE = 4.44 V)
            ox_potential = delta_g_sol - 4.44

            redox_results.append({
                "smiles": smiles,
                "e_neutral_gas_au": e_neutral_gas,
                "e_charged_gas_au": e_charged_gas,
                "e_neutral_sol_au": e_neutral_sol,
                "e_charged_sol_au": e_charged_sol,
                "delta_g_sol_ev": delta_g_sol,
                "oxidation_potential_v": ox_potential,
            })

    return {"redox_potentials": redox_results}


def _calculate_reorganization(job, qc_results_by_type: Dict, db: Session) -> Dict:
    """计算 Marcus 重组能

    qc_results_by_type: Dict[(task_type, structure_id), result_dict]
    """
    # λ = (λ_ox + λ_red) / 2
    # λ_ox = E(N|N+) - E(N+|N+)
    # λ_red = E(N+|N) - E(N|N)

    # 这里简化处理，实际需要更复杂的几何优化和单点计算
    return {
        "message": "重组能计算需要多步优化，请使用专门的重组能计算功能",
        "status": "not_implemented"
    }


@router.get("/jobs/{job_id}/results")
def get_cluster_analysis_results(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取 Cluster 高级计算结果"""
    job = db.query(AdvancedClusterJobModel).filter(AdvancedClusterJobModel.id == job_id).first()

    if not job:
        raise HTTPException(status_code=404, detail=f"任务 {job_id} 不存在")

    if job.user_id != current_user.id and current_user.role.value != 'admin':
        raise HTTPException(status_code=403, detail="无权访问此任务")

    return {
        "job_id": job_id,
        "status": job.status.value,
        "progress": job.progress,
        "calc_types": job.calc_types,
        "results": job.results,
        "qc_task_plan": job.qc_task_plan,
    }


@router.get("/jobs/{job_id}/qc-status")
def get_cluster_analysis_qc_status(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取 Cluster 高级计算的 QC 任务状态"""
    job = db.query(AdvancedClusterJobModel).filter(AdvancedClusterJobModel.id == job_id).first()

    if not job:
        raise HTTPException(status_code=404, detail=f"任务 {job_id} 不存在")

    if job.user_id != current_user.id and current_user.role.value != 'admin':
        raise HTTPException(status_code=403, detail="无权访问此任务")

    # 方法 1: 从 qc_task_plan 获取 QC 任务 ID（旧数据兼容）
    qc_task_plan = job.qc_task_plan or {}
    plan_qc_job_ids = set(
        qc_task_plan.get("reused_qc_jobs", []) +
        qc_task_plan.get("new_qc_jobs", [])
    )

    # 方法 2: 直接通过 cluster_analysis_job_id 查询关联的 QC 任务（新数据）
    linked_qc_jobs = db.query(QCJob).filter(
        QCJob.cluster_analysis_job_id == job_id
    ).all()
    linked_qc_job_ids = set(qc.id for qc in linked_qc_jobs)

    # 合并两种方式的 ID
    all_qc_job_ids = list(plan_qc_job_ids | linked_qc_job_ids)

    if not all_qc_job_ids:
        return {
            "job_id": job_id,
            "total_qc_jobs": 0,
            "completed": 0,
            "running": 0,
            "pending": 0,
            "failed": 0,
            "all_completed": True,
        }

    # 查询所有 QC 任务状态
    qc_jobs = db.query(QCJob).filter(QCJob.id.in_(all_qc_job_ids)).all()

    status_counts = {"COMPLETED": 0, "RUNNING": 0, "SUBMITTED": 0, "CREATED": 0, "FAILED": 0}
    for qc_job in qc_jobs:
        status = qc_job.status.value if hasattr(qc_job.status, 'value') else str(qc_job.status)
        if status in status_counts:
            status_counts[status] += 1

    all_completed = status_counts["COMPLETED"] == len(qc_jobs) and status_counts["FAILED"] == 0

    return {
        "job_id": job_id,
        "total_qc_jobs": len(qc_jobs),
        "completed": status_counts["COMPLETED"],
        "running": status_counts["RUNNING"],
        "pending": status_counts["SUBMITTED"] + status_counts["CREATED"],
        "failed": status_counts["FAILED"],
        "all_completed": all_completed,
        "qc_jobs": [
            {
                "id": qc.id,
                "status": qc.status.value if hasattr(qc.status, 'value') else str(qc.status),
                "molecule_name": qc.molecule_name,
                "task_type": qc.task_type,
            }
            for qc in qc_jobs
        ]
    }
