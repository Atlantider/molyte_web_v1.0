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
    # 权限检查：任务所有者或管理员可以访问
    is_owner = md_job.user_id == current_user.id
    is_admin = current_user.role.value == 'ADMIN'  # UserRole.ADMIN 的值是 'ADMIN'
    if not (is_owner or is_admin):
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
    solvent_model: str = "gas",
    solvent_name: str = None,
    molecule_type: str = None,
    molecule_name: str = None
) -> Optional[QCJob]:
    """
    查找已有的、相同配置的 QC 任务（全平台复用）

    复用条件：
    - 相同 SMILES（分子结构）或相同基础 molecule_name（忽略 #数字 后缀）
    - 相同电荷和自旋多重度
    - 相同计算参数（泛函、基组）
    - 相同溶剂模型配置
    - 任务已完成且未删除

    注意：不限制用户，实现全平台复用
    """
    import re
    from sqlalchemy.orm import joinedload

    # 策略 1：先尝试 SMILES 精确匹配
    query = db.query(QCJob).options(
        joinedload(QCJob.results)
    ).filter(
        QCJob.smiles == smiles,
        QCJob.charge == charge,
        QCJob.spin_multiplicity == multiplicity,
        QCJob.functional == functional,
        QCJob.basis_set == basis_set,
        QCJob.solvent_model == solvent_model,
        QCJob.status == QCJobStatus.COMPLETED,
        QCJob.is_deleted == False
    )

    if solvent_model != "gas" and solvent_name:
        query = query.filter(QCJob.solvent_name == solvent_name)

    if molecule_type:
        query = query.filter(QCJob.molecule_type == molecule_type)

    job = query.order_by(desc(QCJob.finished_at)).first()

    # 策略 2：如果 SMILES 匹配不到，尝试用基础 molecule_name 匹配
    # 这对于 EC#1, EC#2 等同类型分子很有用
    if not job and molecule_name:
        base_name = re.sub(r'#\d+$', '', molecule_name)
        name_pattern = f"^{re.escape(base_name)}(#[0-9]+)?$"

        name_query = db.query(QCJob).options(
            joinedload(QCJob.results)
        ).filter(
            QCJob.molecule_name.op('~')(name_pattern),
            QCJob.charge == charge,
            QCJob.spin_multiplicity == multiplicity,
            QCJob.functional == functional,
            QCJob.basis_set == basis_set,
            QCJob.solvent_model == solvent_model,
            QCJob.status == QCJobStatus.COMPLETED,
            QCJob.is_deleted == False
        )

        if solvent_model != "gas" and solvent_name:
            name_query = name_query.filter(QCJob.solvent_name == solvent_name)

        job = name_query.order_by(desc(QCJob.finished_at)).first()

        if job:
            logger.debug(f"[全局复用-名称匹配] 找到已有 QC 任务 {job.id}: "
                        f"'{molecule_name}' matches '{job.molecule_name}'")

    if job:
        logger.debug(f"[全局复用] 找到已有 QC 任务 {job.id}: {smiles}, charge={charge}, "
                     f"solvent={solvent_model}/{solvent_name}, energy={job.results[0].energy_au if job.results else 'N/A'}")

    return job


# 常见分子的 SMILES 和电荷信息（用作后备）
MOLECULE_INFO_MAP = {
    # 阳离子
    "Li": {"smiles": "[Li+]", "charge": 1},
    "Na": {"smiles": "[Na+]", "charge": 1},
    "K": {"smiles": "[K+]", "charge": 1},
    "Mg": {"smiles": "[Mg+2]", "charge": 2},
    "Ca": {"smiles": "[Ca+2]", "charge": 2},
    "Zn": {"smiles": "[Zn+2]", "charge": 2},
    # 阴离子
    "PF6": {"smiles": "F[P-](F)(F)(F)(F)F", "charge": -1},
    "BF4": {"smiles": "F[B-](F)(F)F", "charge": -1},
    "TFSI": {"smiles": "FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F", "charge": -1},
    "FSI": {"smiles": "FS(=O)(=O)[N-]S(=O)(=O)F", "charge": -1},
    "DFOB": {"smiles": "FB1OC(=O)C(=O)O[B-]1F", "charge": -1},
    "ClO4": {"smiles": "[O-]Cl(=O)(=O)=O", "charge": -1},
    "NO3": {"smiles": "[O-][N+](=O)[O-]", "charge": -1},
    "OTf": {"smiles": "[O-]S(=O)(=O)C(F)(F)F", "charge": -1},
    "F": {"smiles": "[F-]", "charge": -1},
    "Cl": {"smiles": "[Cl-]", "charge": -1},
    # 碳酸酯溶剂
    "EC": {"smiles": "C1COC(=O)O1", "charge": 0},
    "PC": {"smiles": "CC1COC(=O)O1", "charge": 0},
    "DMC": {"smiles": "COC(=O)OC", "charge": 0},
    "DEC": {"smiles": "CCOC(=O)OCC", "charge": 0},
    "EMC": {"smiles": "CCOC(=O)OC", "charge": 0},
    "FEC": {"smiles": "FC1COC(=O)O1", "charge": 0},
    "VC": {"smiles": "C=C1COC(=O)O1", "charge": 0},
    # 醚类溶剂
    "DME": {"smiles": "COCCOC", "charge": 0},
    "DOL": {"smiles": "C1COCO1", "charge": 0},
    "TEGDME": {"smiles": "COCCOCCOCCOCCOC", "charge": 0},
    "DEGDME": {"smiles": "COCCOCCOC", "charge": 0},
    # 其他溶剂
    "MPN": {"smiles": "CCCCC#N", "charge": 0},  # 3-Methoxypropionitrile
    "AN": {"smiles": "CC#N", "charge": 0},      # Acetonitrile
    "THF": {"smiles": "C1CCOC1", "charge": 0},
    "DMSO": {"smiles": "CS(=O)C", "charge": 0},
    "TTE": {"smiles": "FC(F)(F)C(F)(F)OCC(F)(F)F", "charge": 0},  # 1,1,2,2-tetrafluoroethyl-2,2,3,3-tetrafluoropropyl ether
}


def _get_molecule_info_from_md_job(db: Session, md_job: MDJob) -> Dict[str, Dict[str, Any]]:
    """
    从 MD 任务获取分子信息（SMILES 和电荷）

    优先级：
    1. ResultSummary.molecule_structures（来自 MD 模拟结果）
    2. ElectrolyteSystem 的阴离子和溶剂配置
    3. 后备映射表
    """
    from app.models.result import ResultSummary
    from app.models.electrolyte import ElectrolyteSystem

    mol_info: Dict[str, Dict[str, Any]] = {}

    # 1. 尝试从 ResultSummary 获取
    result_summary = db.query(ResultSummary).filter(
        ResultSummary.md_job_id == md_job.id
    ).first()

    if result_summary and result_summary.molecule_structures:
        for mol in result_summary.molecule_structures:
            name = mol.get('name')
            if name:
                smiles = mol.get('smiles')
                # total_charge 是浮点数，需要四舍五入
                charge = round(mol.get('total_charge', 0)) if mol.get('total_charge') is not None else None
                if smiles:
                    mol_info[name] = {"smiles": smiles, "charge": charge}
                    logger.debug(f"Got molecule info from ResultSummary: {name} -> smiles={smiles}, charge={charge}")

    # 2. 尝试从 ElectrolyteSystem 获取（补充 ResultSummary 中缺失的）
    if md_job.system_id:
        electrolyte = db.query(ElectrolyteSystem).filter(
            ElectrolyteSystem.id == md_job.system_id
        ).first()

        if electrolyte:
            # 从阴离子配置获取
            for anion in (electrolyte.anions or []):
                name = anion.get('name')
                if name and name not in mol_info:
                    smiles = anion.get('smiles')
                    charge = anion.get('charge', -1)
                    if smiles:
                        mol_info[name] = {"smiles": smiles, "charge": charge}
                        logger.debug(f"Got anion info from ElectrolyteSystem: {name} -> smiles={smiles}, charge={charge}")

            # 从溶剂配置获取
            for solvent in (electrolyte.solvents or []):
                name = solvent.get('name')
                if name and name not in mol_info:
                    smiles = solvent.get('smiles')
                    charge = solvent.get('charge', 0)
                    if smiles:
                        mol_info[name] = {"smiles": smiles, "charge": charge}
                        logger.debug(f"Got solvent info from ElectrolyteSystem: {name} -> smiles={smiles}, charge={charge}")

    return mol_info


def _get_ligand_smiles_from_composition(
    composition: Dict[str, int],
    mol_info_from_db: Optional[Dict[str, Dict[str, Any]]] = None
) -> List[Dict[str, Any]]:
    """
    从 composition 获取配体 SMILES 列表

    Args:
        composition: 配位组成，如 {"EC": 2, "PF6": 1}
        mol_info_from_db: 从数据库获取的分子信息（优先使用）
    """
    mol_info_from_db = mol_info_from_db or {}

    ligands = []
    for mol_name, count in composition.items():
        if count > 0 and mol_name not in ["Li", "Na", "K", "Mg", "Ca", "Zn"]:  # 排除阳离子
            # 优先使用数据库中的信息
            if mol_name in mol_info_from_db:
                info = mol_info_from_db[mol_name]
                smiles = info.get("smiles", mol_name)
                charge = info.get("charge")
                # 如果 charge 是 None，尝试从后备映射获取
                if charge is None and mol_name in MOLECULE_INFO_MAP:
                    charge = MOLECULE_INFO_MAP[mol_name]["charge"]
                elif charge is None:
                    charge = 0  # 默认中性
            # 后备：使用映射表
            elif mol_name in MOLECULE_INFO_MAP:
                info = MOLECULE_INFO_MAP[mol_name]
                smiles = info["smiles"]
                charge = info["charge"]
            else:
                # 未知分子，使用名称作为占位符
                smiles = mol_name
                charge = 0
                logger.warning(f"Unknown molecule {mol_name}, using name as SMILES placeholder")

            ligands.append({
                "name": mol_name,
                "smiles": smiles,
                "count": count,
                "charge": charge
            })
    return ligands


def _plan_qc_tasks_for_calc_type(
    db: Session,
    calc_type: ClusterCalcType,
    structures: List[SolvationStructure],
    qc_config: Dict[str, Any],
    mol_info_from_db: Optional[Dict[str, Dict[str, Any]]] = None,
    redox_options: Optional[Dict[str, bool]] = None,
    reorganization_options: Optional[Dict[str, bool]] = None
) -> CalcTypeRequirements:
    """为某个计算类型规划 QC 任务（支持全局复用）

    Args:
        redox_options: {"include_molecule": bool, "include_dimer": bool}
        reorganization_options: {"include_molecule": bool, "include_cluster": bool}
    """
    functional = qc_config.get("functional", "B3LYP")
    basis_set = qc_config.get("basis_set", "6-31G*")
    charge_ion = qc_config.get("charge_ion", 1)
    # 溶剂配置：支持 solvent_model/solvent 或 solvent_model/solvent_name
    solvent_model = qc_config.get("solvent_model") or "gas"
    solvent_name = qc_config.get("solvent") or qc_config.get("solvent_name")

    logger.info(f"[规划] calc_type={calc_type}, functional={functional}, basis_set={basis_set}, "
                f"solvent_model={solvent_model}, solvent_name={solvent_name}")

    planned_tasks: List[PlannedQCTask] = []
    new_count = 0
    reused_count = 0

    # 收集所有需要的配体类型
    all_ligand_types = set()
    for s in structures:
        comp = s.composition or {}
        for mol_name, count in comp.items():
            if count > 0 and mol_name not in ["Li", "Na", "K", "Mg", "Ca", "Zn"]:
                all_ligand_types.add(mol_name)

    # 使用从数据库获取的分子信息
    ligand_info = _get_ligand_smiles_from_composition(
        {name: 1 for name in all_ligand_types},
        mol_info_from_db
    )
    
    # 根据计算类型规划任务
    if calc_type in [ClusterCalcType.BINDING_TOTAL, ClusterCalcType.DESOLVATION_STEPWISE,
                     ClusterCalcType.DESOLVATION_FULL]:
        # 这些计算需要 cluster 和 ligand 能量

        # 1. Cluster 任务 (每个结构一个)
        first_structure_id = None  # 记录第一个结构 ID，用于配体提取
        for s in structures:
            if first_structure_id is None:
                first_structure_id = s.id
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

        # 2. Ligand 任务 (每种配体一个，从 cluster 中提取几何结构)
        # 注意：配体几何结构从 cluster 中提取，不使用 RDKit 生成
        # 这样可以和 BINDING_PAIRWISE 等共享相同的配体任务
        for lig in ligand_info:
            # task_type 使用 ligand_{name} 格式，以便结果计算时能识别配体
            ligand_task_type = f"ligand_{lig['name']}"
            # 配体从 cluster 提取，需要新计算（几何结构依赖源 cluster）
            task = PlannedQCTask(
                task_type=ligand_task_type,
                description=f"配体 {lig['name']}（从 cluster 提取）",
                smiles=lig["smiles"],
                charge=lig["charge"],
                multiplicity=1,
                status="new",
                structure_id=first_structure_id  # 关联到源 cluster
            )
            planned_tasks.append(task)
            new_count += 1
    
    if calc_type == ClusterCalcType.BINDING_TOTAL:
        # Binding 需要 ion 能量（全局可复用）
        li_smiles = "[Li+]"
        existing = _find_existing_qc_job(
            db, li_smiles, 1, 1, functional, basis_set,
            solvent_model, solvent_name, "ion"
        )
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
        # Stepwise 需要所有中间态的 cluster_minus 结构
        # 例如 Li·EC₂·DMC 需要计算：
        #   - Li·EC₂·DMC (完整 cluster，已在上面添加)
        #   - Li·EC·DMC (移除1个EC)
        #   - Li·DMC (移除2个EC)
        #   - Li·EC₂ (移除1个DMC)
        #   - Li·EC (移除1个EC和1个DMC)
        #   - Li (裸离子，需要 ion 任务)

        from itertools import product

        for s in structures:
            comp = s.composition or {}
            # 提取配体及其数量（排除中心离子）
            ligand_counts = {k: v for k, v in comp.items()
                           if v > 0 and k not in ["Li", "Na", "K", "Mg", "Ca", "Zn"]}

            if not ligand_counts:
                continue

            # 生成所有可能的中间态组合
            # 对于每个配体，可以保留 0 到 count 个
            ligand_names = list(ligand_counts.keys())
            ranges = [range(ligand_counts[name] + 1) for name in ligand_names]

            for combo in product(*ranges):
                # combo 是每个配体保留的数量
                remaining = dict(zip(ligand_names, combo))

                # 跳过完整 cluster（已在上面添加）
                if remaining == ligand_counts:
                    continue

                # 跳过裸离子（需要单独的 ion 任务）
                if all(c == 0 for c in combo):
                    continue

                # 生成描述：Li·EC₁·DMC₂ 格式
                parts = [f"{name}_{remaining[name]}" for name in ligand_names if remaining[name] > 0]
                desc_parts = [f"{name}×{remaining[name]}" for name in ligand_names if remaining[name] > 0]

                task_type = f"cluster_minus_{'_'.join(parts)}"
                description = f"中间态 Li·{'·'.join(desc_parts)}"

                task = PlannedQCTask(
                    task_type=task_type,
                    description=description,
                    structure_id=s.id,
                    charge=charge_ion,
                    multiplicity=1,
                    status="new"
                )
                planned_tasks.append(task)
                new_count += 1

        # Stepwise 还需要 ion 能量（最终态）
        li_smiles = "[Li+]"
        existing = _find_existing_qc_job(
            db, li_smiles, 1, 1, functional, basis_set,
            solvent_model, solvent_name, "ion"
        )
        if existing:
            task = PlannedQCTask(
                task_type="ion",
                description="Li+ 离子（最终态）",
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
                description="Li+ 离子（最终态）",
                smiles=li_smiles,
                charge=1,
                multiplicity=1,
                status="new"
            )
            new_count += 1
        planned_tasks.append(task)

    if calc_type == ClusterCalcType.DESOLVATION_FULL:
        # Full desolvation 需要 ion 能量（全局可复用）
        li_smiles = "[Li+]"
        existing = _find_existing_qc_job(
            db, li_smiles, 1, 1, functional, basis_set,
            solvent_model, solvent_name, "ion"
        )
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
        # Pairwise binding:
        # 1. Li+ 离子能量
        # 2. 各配体能量（从 cluster 中提取的实际几何结构）
        # 3. Li-配体二聚体（从 cluster 中提取 Li+ 和某个配体的组合）

        # 1. Li+ 离子（全局可复用）
        li_smiles = "[Li+]"
        existing = _find_existing_qc_job(
            db, li_smiles, 1, 1, functional, basis_set,
            solvent_model, solvent_name, "ion"
        )
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
                task_type="ion", description="Li+ 离子", smiles=li_smiles,
                charge=1, multiplicity=1, status="new"
            )
            new_count += 1
        planned_tasks.append(task)

        # 2 & 3. 从每个 cluster 结构中提取 dimer（Li+ 加上某种配体）
        # 每种配体只需要一个代表性的 dimer
        # 同时也记录配体任务（使用 cluster 中提取的实际几何结构）
        processed_ligand_types = set()  # 记录已处理的配体类型，避免重复

        for s in structures:
            comp = s.composition or {}
            for mol_name, count in comp.items():
                if count > 0 and mol_name not in ["Li", "Na", "K", "Mg", "Ca", "Zn"]:
                    if mol_name not in processed_ligand_types:
                        processed_ligand_types.add(mol_name)

                        # 配体任务（从 cluster 中提取的几何结构，需要新计算）
                        lig_info = mol_info_from_db.get(mol_name, {}) if mol_info_from_db else {}
                        lig_smiles = lig_info.get("smiles") or MOLECULE_INFO_MAP.get(mol_name, {}).get("smiles", mol_name)
                        lig_charge = lig_info.get("charge")
                        if lig_charge is None:
                            lig_charge = MOLECULE_INFO_MAP.get(mol_name, {}).get("charge", 0)

                        # 配体任务 - 需要从 cluster 中提取真实几何结构
                        task = PlannedQCTask(
                            task_type=f"ligand_{mol_name}",
                            description=f"配体 {mol_name}（从 cluster 提取）",
                            smiles=lig_smiles,
                            charge=lig_charge,
                            multiplicity=1,
                            status="new",
                            structure_id=s.id  # 关联到源 cluster，用于提取几何结构
                        )
                        planned_tasks.append(task)
                        new_count += 1

                        # Li-配体 dimer 任务 - 从 cluster 中提取 Li+ 和该配体
                        dimer_charge = charge_ion + lig_charge
                        task = PlannedQCTask(
                            task_type=f"dimer_{mol_name}",
                            description=f"Li-{mol_name} 二聚体（从 cluster 提取）",
                            smiles=f"[Li+].{lig_smiles}",  # 仅用于标识，实际几何从 cluster 提取
                            charge=dimer_charge,
                            multiplicity=1,
                            status="new",
                            structure_id=s.id  # 关联到源 cluster
                        )
                        planned_tasks.append(task)
                        new_count += 1

    if calc_type == ClusterCalcType.REDOX:
        # Redox 计算包含两部分（可通过 redox_options 控制）：
        # A. 单独配体分子的 Redox（4 个单点）
        # B. Li-配体 Dimer 的 Redox（4 个单点，从 cluster 提取几何）
        #
        # 每个物种需要 4 个计算（热力学循环）：
        # 1. 中性态-气相 (neutral_gas)
        # 2. 氧化态-气相 (charged_gas)
        # 3. 中性态-溶液相 (neutral_sol)
        # 4. 氧化态-溶液相 (charged_sol)

        include_molecule = redox_options.get("include_molecule", True) if redox_options else True
        include_dimer = redox_options.get("include_dimer", True) if redox_options else True

        processed_species = set()
        first_structure_id = structures[0].id if structures else None

        for s in structures:
            comp = s.composition or {}
            for mol_name, count in comp.items():
                if count > 0 and mol_name not in ["Li", "Na", "K", "Mg", "Ca", "Zn"]:
                    if mol_name in processed_species:
                        continue
                    processed_species.add(mol_name)

                    lig_info = mol_info_from_db.get(mol_name, {}) if mol_info_from_db else {}
                    lig_smiles = lig_info.get("smiles") or MOLECULE_INFO_MAP.get(mol_name, {}).get("smiles", mol_name)
                    lig_charge = lig_info.get("charge")
                    if lig_charge is None:
                        lig_charge = MOLECULE_INFO_MAP.get(mol_name, {}).get("charge", 0)

                    # A. 单独配体分子的 Redox（4 个状态）
                    if include_molecule:
                        mol_states = [
                            ("neutral_gas", f"[分子] {mol_name} 中性-气相", lig_charge, "gas", None),
                            ("charged_gas", f"[分子] {mol_name} 氧化态-气相", lig_charge + 1, "gas", None),
                            ("neutral_sol", f"[分子] {mol_name} 中性-溶液相", lig_charge, solvent_model, solvent_name),
                            ("charged_sol", f"[分子] {mol_name} 氧化态-溶液相", lig_charge + 1, solvent_model, solvent_name),
                        ]

                        for state_suffix, desc, charge, sol_model, sol_name in mol_states:
                            task_type = f"redox_mol_{mol_name}_{state_suffix}"
                            mult = 1 if (charge - lig_charge) % 2 == 0 else 2

                            existing = _find_existing_qc_job(
                                db, lig_smiles, charge, mult, functional, basis_set, sol_model, sol_name, task_type
                            )

                            if existing:
                                task = PlannedQCTask(
                                    task_type=task_type, description=desc, smiles=lig_smiles,
                                    charge=charge, multiplicity=mult, status="reused",
                                    existing_qc_job_id=existing.id,
                                    existing_energy=existing.results[0].energy_au if existing.results else None
                                )
                                reused_count += 1
                            else:
                                task = PlannedQCTask(
                                    task_type=task_type, description=desc, smiles=lig_smiles,
                                    charge=charge, multiplicity=mult, status="new"
                                )
                                new_count += 1
                            planned_tasks.append(task)

                    # B. Li-配体 Dimer 的 Redox（4 个状态，从 cluster 提取几何）
                    if include_dimer:
                        dimer_base_charge = charge_ion + lig_charge  # Li+ + 配体电荷
                        dimer_smiles = f"[Li+].{lig_smiles}"

                        dimer_states = [
                            ("neutral_gas", f"[Dimer] Li-{mol_name} 中性-气相", dimer_base_charge, "gas", None),
                            ("charged_gas", f"[Dimer] Li-{mol_name} 氧化态-气相", dimer_base_charge + 1, "gas", None),
                            ("neutral_sol", f"[Dimer] Li-{mol_name} 中性-溶液相", dimer_base_charge, solvent_model, solvent_name),
                            ("charged_sol", f"[Dimer] Li-{mol_name} 氧化态-溶液相", dimer_base_charge + 1, solvent_model, solvent_name),
                        ]

                        for state_suffix, desc, charge, sol_model, sol_name in dimer_states:
                            task_type = f"redox_dimer_{mol_name}_{state_suffix}"
                            mult = 1 if charge == dimer_base_charge else 2

                            # Dimer 也尝试全局复用（相同 SMILES + 电荷 + 多重度 + 溶剂配置）
                            existing = _find_existing_qc_job(
                                db, dimer_smiles, charge, mult, functional, basis_set, sol_model, sol_name, task_type
                            )

                            if existing:
                                task = PlannedQCTask(
                                    task_type=task_type, description=desc, smiles=dimer_smiles,
                                    charge=charge, multiplicity=mult, status="reused",
                                    existing_qc_job_id=existing.id,
                                    existing_energy=existing.results[0].energy_au if existing.results else None,
                                    structure_id=first_structure_id
                                )
                                reused_count += 1
                            else:
                                task = PlannedQCTask(
                                    task_type=task_type, description=desc, smiles=dimer_smiles,
                                    charge=charge, multiplicity=mult, status="new",
                                    structure_id=first_structure_id
                                )
                                new_count += 1
                            planned_tasks.append(task)

    if calc_type == ClusterCalcType.REORGANIZATION:
        # Marcus 重组能包含两部分（可通过 reorganization_options 控制）：
        # A. 单独配体分子的重组能（4 点计算）
        # B. 整个 Cluster 的重组能（4 点计算，从 MD 提取几何）
        #
        # λ = (λ_ox + λ_red) / 2
        # λ_ox = E(N@N+) - E(N+@N+)  # 中性几何下氧化态能量 - 氧化几何下氧化态能量
        # λ_red = E(N+@N) - E(N@N)   # 氧化几何下中性态能量 - 中性几何下中性态能量
        #
        # 需要 4 个计算：
        # 1. 中性态优化几何 (opt_neutral)
        # 2. 氧化态优化几何 (opt_charged)
        # 3. 中性几何+氧化态单点 (sp_charged_at_neutral)
        # 4. 氧化几何+中性态单点 (sp_neutral_at_charged)

        include_molecule = reorganization_options.get("include_molecule", True) if reorganization_options else True
        include_cluster = reorganization_options.get("include_cluster", True) if reorganization_options else True

        processed_species = set()

        # A. 单独配体分子的重组能
        if include_molecule:
            for s in structures:
                comp = s.composition or {}
                for mol_name, count in comp.items():
                    if count > 0 and mol_name not in ["Li", "Na", "K", "Mg", "Ca", "Zn"]:
                        if mol_name in processed_species:
                            continue
                        processed_species.add(mol_name)

                        lig_info = mol_info_from_db.get(mol_name, {}) if mol_info_from_db else {}
                        lig_smiles = lig_info.get("smiles") or MOLECULE_INFO_MAP.get(mol_name, {}).get("smiles", mol_name)
                        lig_charge = lig_info.get("charge")
                        if lig_charge is None:
                            lig_charge = MOLECULE_INFO_MAP.get(mol_name, {}).get("charge", 0)

                        # 分子的 4 个计算任务
                        mol_tasks = [
                            (f"reorg_mol_{mol_name}_opt_neutral", f"[分子] {mol_name} 中性态优化", lig_charge, 1, "opt"),
                            (f"reorg_mol_{mol_name}_opt_charged", f"[分子] {mol_name} 氧化态优化", lig_charge + 1, 2, "opt"),
                            (f"reorg_mol_{mol_name}_sp_charged_at_neutral", f"[分子] {mol_name} 氧化态@中性几何", lig_charge + 1, 2, "sp"),
                            (f"reorg_mol_{mol_name}_sp_neutral_at_charged", f"[分子] {mol_name} 中性态@氧化几何", lig_charge, 1, "sp"),
                        ]

                        for task_type, desc, charge, mult, calc_mode in mol_tasks:
                            # 尝试全局复用（优化任务也可以复用，如果配置相同）
                            existing = _find_existing_qc_job(
                                db, lig_smiles, charge, mult, functional, basis_set,
                                solvent_model, solvent_name, task_type
                            )

                            if existing:
                                task = PlannedQCTask(
                                    task_type=task_type, description=desc, smiles=lig_smiles,
                                    charge=charge, multiplicity=mult, status="reused",
                                    existing_qc_job_id=existing.id,
                                    existing_energy=existing.results[0].energy_au if existing.results else None
                                )
                                reused_count += 1
                            else:
                                task = PlannedQCTask(
                                    task_type=task_type, description=desc, smiles=lig_smiles,
                                    charge=charge, multiplicity=mult, status="new"
                                )
                                new_count += 1
                            planned_tasks.append(task)

        # B. 整个 Cluster 的重组能（每个 cluster 结构）
        if include_cluster:
            for s in structures:
                cluster_charge = charge_ion  # cluster 电荷（通常是 Li+ 的 +1）
                cluster_desc = f"Cluster #{s.id} (CN={s.coordination_num})"

                cluster_tasks = [
                    (f"reorg_cluster_{s.id}_opt_neutral", f"[Cluster] {cluster_desc} 中性态优化", cluster_charge, 1, "opt"),
                    (f"reorg_cluster_{s.id}_opt_charged", f"[Cluster] {cluster_desc} 氧化态优化", cluster_charge + 1, 2, "opt"),
                    (f"reorg_cluster_{s.id}_sp_charged_at_neutral", f"[Cluster] {cluster_desc} 氧化态@中性几何", cluster_charge + 1, 2, "sp"),
                    (f"reorg_cluster_{s.id}_sp_neutral_at_charged", f"[Cluster] {cluster_desc} 中性态@氧化几何", cluster_charge, 1, "sp"),
                ]

                for task_type, desc, charge, mult, calc_mode in cluster_tasks:
                    # Cluster 任务从 MD 提取几何，不尝试全局复用（每个 cluster 几何不同）
                    task = PlannedQCTask(
                        task_type=task_type, description=desc, smiles=None,
                        structure_id=s.id, charge=charge, multiplicity=mult, status="new"
                    )
                    new_count += 1
                    planned_tasks.append(task)

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

    # 从 MD 任务获取分子信息
    mol_info_from_db = _get_molecule_info_from_md_job(db, md_job)
    logger.info(f"Got molecule info from MD job {md_job.id}: {list(mol_info_from_db.keys())}")

    # QC 配置
    qc_config = request.qc_config.model_dump() if request.qc_config else {}

    # 子选项
    redox_opts = request.redox_options.model_dump() if request.redox_options else None
    reorg_opts = request.reorganization_options.model_dump() if request.reorganization_options else None

    # 为每个计算类型规划 QC 任务
    calc_requirements = []
    warnings = []

    # 跨计算类型复用追踪：记录本次请求中已规划的任务
    # key = (task_type_base, smiles, charge) -> 第一次出现的任务信息
    planned_in_request: Dict[tuple, PlannedQCTask] = {}

    for calc_type in request.calc_types:
        req = _plan_qc_tasks_for_calc_type(
            db, calc_type, structures, qc_config, mol_info_from_db,
            redox_options=redox_opts,
            reorganization_options=reorg_opts
        )

        # 更新任务状态：如果任务已经在之前的计算类型中规划过，标记为"跨类型复用"
        updated_tasks = []
        new_count = 0
        reused_count = 0

        for task in req.required_qc_tasks:
            # 生成唯一键用于跨计算类型复用检测
            # 策略：
            # - cluster/intermediate 类型：用 structure_id 区分（不同 cluster 不复用）
            # - ligand/dimer/ion 等类型：用 task_type + smiles + charge 匹配（可跨类型复用）
            if task.task_type in ["cluster", "intermediate"]:
                # 这些任务依赖具体的 cluster 结构，用 structure_id 区分
                task_key = (task.task_type, None, task.structure_id, task.charge)
            else:
                # ligand_XX, dimer_XX, ion 等可以跨类型复用
                task_key = (task.task_type, task.smiles, None, task.charge)

            if task.status == "reused":
                # 已经在数据库中找到复用（全局复用）
                # 更新描述添加标识
                task.description = f"{task.description}【全局复用】"
                updated_tasks.append(task)
                reused_count += 1
            elif task_key in planned_in_request:
                # 本次请求中已经规划过（局部复用 - 跨计算类型）
                first_task = planned_in_request[task_key]
                updated_task = PlannedQCTask(
                    task_type=task.task_type,
                    description=f"{task.description}【局部复用】",
                    smiles=task.smiles,
                    structure_id=task.structure_id,
                    charge=task.charge,
                    multiplicity=task.multiplicity,
                    status="local_reused",  # 局部复用状态
                    existing_qc_job_id=None,  # 还没有 job
                    existing_energy=None
                )
                updated_tasks.append(updated_task)
                reused_count += 1
            else:
                # 新任务，记录到追踪表
                planned_in_request[task_key] = task
                updated_tasks.append(task)
                new_count += 1

        # 更新 requirements
        req.required_qc_tasks = updated_tasks
        req.new_tasks_count = new_count
        req.reused_tasks_count = reused_count

        calc_requirements.append(req)

        # 添加警告
        if calc_type == ClusterCalcType.REORGANIZATION:
            warnings.append("⚠️ 重组能计算量极大，每个物种需要多次优化，建议限制结构数量")
        if calc_type == ClusterCalcType.REDOX:
            warnings.append("⚠️ Redox 计算对方法/基组敏感，结果仅供参考")

    # 统计唯一任务数
    unique_new = len([t for t in planned_in_request.values() if t.status == "new"])
    unique_reused_from_db = sum(1 for req in calc_requirements
                                 for task in req.required_qc_tasks
                                 if task.status == "reused" and task.existing_qc_job_id)
    # 去重计算跨类型复用
    seen_reused = set()
    for req in calc_requirements:
        for task in req.required_qc_tasks:
            if task.status == "reused" and task.existing_qc_job_id:
                seen_reused.add(task.existing_qc_job_id)
    unique_reused = len(seen_reused)

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

    # 从 MD 任务获取分子信息
    mol_info_from_db = _get_molecule_info_from_md_job(db, md_job)

    # QC 配置
    qc_config = request.qc_config.model_dump() if request.qc_config else {}

    # 规划 QC 任务
    all_planned_tasks = []
    reused_qc_jobs = []

    for calc_type in request.calc_types:
        req = _plan_qc_tasks_for_calc_type(db, calc_type, structures, qc_config, mol_info_from_db)
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

    # 获取 MD 任务和分子信息
    md_job = db.query(MDJob).filter(MDJob.id == job.md_job_id).first()
    mol_info_from_db = _get_molecule_info_from_md_job(db, md_job) if md_job else None

    # 规划新的 QC 任务
    new_requirements = []
    for calc_type in request.additional_calc_types:
        if calc_type.value in job.calc_types:
            continue  # 跳过已有的计算类型
        req = _plan_qc_tasks_for_calc_type(db, calc_type, structures, job.qc_config, mol_info_from_db)
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

    E_bind = E_cluster - (E_ion + Σ n_j × E_ligand_j)
    其中 n_j 是每种配体的数量，来自 SolvationStructure.composition
    """
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

    # 获取每个 cluster 结构的 composition 信息
    structure_ids = [c["structure_id"] for c in cluster_results]
    structures = db.query(SolvationStructure).filter(
        SolvationStructure.id.in_(structure_ids)
    ).all()
    structure_composition_map = {s.id: s.composition or {} for s in structures}

    # 计算每个 cluster 的 binding energy
    binding_results = []

    for cluster in cluster_results:
        structure_id = cluster["structure_id"]
        e_cluster = cluster["energy_au"]

        # 根据 composition 计算配体总能量
        composition = structure_composition_map.get(structure_id, {})
        total_ligand_energy = 0.0
        ligand_details = {}

        for ligand_name, e_ligand in ligand_energies.items():
            count = composition.get(ligand_name, 0)
            if count > 0:
                total_ligand_energy += count * e_ligand
                ligand_details[ligand_name] = {"count": count, "energy_au": e_ligand}

        e_bind_au = e_cluster - (e_ion + total_ligand_energy)
        binding_results.append({
            "structure_id": structure_id,
            "composition": composition,
            "e_cluster_au": e_cluster,
            "e_ion_au": e_ion,
            "total_ligand_energy_au": total_ligand_energy,
            "ligand_details": ligand_details,
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
