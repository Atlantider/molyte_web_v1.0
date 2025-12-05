"""
Desolvation energy calculation API endpoints
"""
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
from sqlalchemy import func
from typing import List, Dict, Any
from datetime import datetime
from collections import defaultdict

from app.database import get_db
from app.models.user import User
from app.models.job import PostprocessJob, PostprocessType, JobStatus, MDJob
from app.models.electrolyte import ElectrolyteSystem
from app.models.result import SolvationStructure, DesolvationEnergyResult
from app.models.qc import QCJob
from app.schemas.desolvation import (
    DesolvationJobCreate,
    DesolvationJobResponse,
    DesolvationEnergyResultSchema,
    LigandDesolvationResult,
    TypeSummary,
    BatchDesolvationJobCreate,
    BatchDesolvationJobResponse,
    DesolvationOverviewResponse
)
from app.dependencies import get_current_active_user

router = APIRouter()


@router.post("/jobs", response_model=DesolvationJobResponse)
async def create_desolvation_job(
    job_data: DesolvationJobCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    创建去溶剂化能计算任务
    
    如果该 cluster + method_level 已有 COMPLETED 的结果，直接返回
    """
    # 1. 检查溶剂化结构是否存在
    solvation_structure = db.query(SolvationStructure).filter(
        SolvationStructure.id == job_data.solvation_structure_id
    ).first()
    
    if not solvation_structure:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Solvation structure {job_data.solvation_structure_id} not found"
        )
    
    # 2. 检查是否已有完成的任务
    existing_job = db.query(PostprocessJob).join(
        DesolvationEnergyResult
    ).filter(
        PostprocessJob.md_job_id == job_data.md_job_id,
        PostprocessJob.job_type == PostprocessType.DESOLVATION_ENERGY,
        PostprocessJob.status == JobStatus.COMPLETED,
        DesolvationEnergyResult.solvation_structure_id == job_data.solvation_structure_id,
        DesolvationEnergyResult.method_level == job_data.method_level
    ).first()
    
    if existing_job:
        # 返回已有任务
        return _build_job_response(existing_job, db)
    
    # 3. 创建新任务
    new_job = PostprocessJob(
        md_job_id=job_data.md_job_id,
        job_type=PostprocessType.DESOLVATION_ENERGY,
        status=JobStatus.SUBMITTED,  # 等待 Worker 拉取
        config={
            "solvation_structure_id": job_data.solvation_structure_id,
            "method_level": job_data.method_level,
            "desolvation_mode": job_data.desolvation_mode  # stepwise or full
        }
    )

    db.add(new_job)
    db.commit()
    db.refresh(new_job)

    return _build_job_response(new_job, db)


@router.get("/jobs/{job_id}", response_model=DesolvationJobResponse)
async def get_desolvation_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取去溶剂化能任务详情"""
    job = db.query(PostprocessJob).filter(
        PostprocessJob.id == job_id,
        PostprocessJob.job_type == PostprocessType.DESOLVATION_ENERGY
    ).first()
    
    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Desolvation job {job_id} not found"
        )
    
    return _build_job_response(job, db)


@router.get("/cluster/{cluster_id}/jobs", response_model=List[DesolvationJobResponse])
async def list_cluster_desolvation_jobs(
    cluster_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取某个 cluster 的所有去溶剂化能任务"""
    jobs = db.query(PostprocessJob).join(
        DesolvationEnergyResult,
        PostprocessJob.id == DesolvationEnergyResult.postprocess_job_id,
        isouter=True
    ).filter(
        PostprocessJob.job_type == PostprocessType.DESOLVATION_ENERGY
    ).filter(
        (DesolvationEnergyResult.solvation_structure_id == cluster_id) |
        (PostprocessJob.config['solvation_structure_id'].astext == str(cluster_id))
    ).order_by(PostprocessJob.created_at.desc()).all()

    return [_build_job_response(job, db) for job in jobs]


@router.post("/batch", response_model=BatchDesolvationJobResponse)
async def batch_create_desolvation_jobs(
    batch_data: BatchDesolvationJobCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    批量创建去溶剂化能计算任务
    自动跳过已完成的相同任务
    """
    created_jobs = []
    skipped_count = 0

    for structure_id in batch_data.structure_ids:
        # 检查溶剂化结构是否存在
        solvation_structure = db.query(SolvationStructure).filter(
            SolvationStructure.id == structure_id
        ).first()

        if not solvation_structure:
            continue

        # 检查是否已有完成或正在进行的任务
        existing_job = db.query(PostprocessJob).filter(
            PostprocessJob.job_type == PostprocessType.DESOLVATION_ENERGY,
            PostprocessJob.config['solvation_structure_id'].astext == str(structure_id),
            PostprocessJob.config['method_level'].astext == batch_data.method_level,
            PostprocessJob.status.in_([JobStatus.COMPLETED, JobStatus.RUNNING, JobStatus.QUEUED, JobStatus.SUBMITTED])
        ).first()

        if existing_job:
            skipped_count += 1
            created_jobs.append(_build_job_response(existing_job, db))
            continue

        # 创建新任务
        config = {
            "solvation_structure_id": structure_id,
            "method_level": batch_data.method_level,
            "desolvation_mode": batch_data.desolvation_mode,
            # Slurm 资源配置
            "slurm_partition": batch_data.slurm_partition or "cpu",
            "slurm_cpus": batch_data.slurm_cpus or 16,
            "slurm_time": batch_data.slurm_time or 7200,
        }
        if batch_data.solvent_config:
            config["solvent_config"] = batch_data.solvent_config.model_dump()

        new_job = PostprocessJob(
            md_job_id=batch_data.md_job_id,
            job_type=PostprocessType.DESOLVATION_ENERGY,
            status=JobStatus.SUBMITTED,
            config=config
        )

        db.add(new_job)
        db.commit()
        db.refresh(new_job)
        created_jobs.append(_build_job_response(new_job, db))

    return BatchDesolvationJobResponse(
        created_count=len(batch_data.structure_ids) - skipped_count,
        skipped_count=skipped_count,
        jobs=created_jobs
    )


@router.get("/md/{md_job_id}/overview", response_model=DesolvationOverviewResponse)
async def get_md_desolvation_overview(
    md_job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取某个 MD 任务下所有去溶剂化能计算的总览
    用于监控面板显示
    """
    # 获取 MD 任务和电解液信息
    md_job = db.query(MDJob).filter(MDJob.id == md_job_id).first()
    electrolyte_name = None
    if md_job and md_job.system_id:
        electrolyte = db.query(ElectrolyteSystem).filter(
            ElectrolyteSystem.id == md_job.system_id
        ).first()
        if electrolyte:
            electrolyte_name = electrolyte.name

    # 获取所有去溶剂化任务
    jobs = db.query(PostprocessJob).filter(
        PostprocessJob.md_job_id == md_job_id,
        PostprocessJob.job_type == PostprocessType.DESOLVATION_ENERGY
    ).order_by(PostprocessJob.created_at.desc()).all()

    # 统计状态
    status_summary = defaultdict(int)
    job_responses = []
    for job in jobs:
        status_summary[job.status.value] += 1
        job_responses.append(_build_job_response(job, db))

    return DesolvationOverviewResponse(
        md_job_id=md_job_id,
        electrolyte_name=electrolyte_name,
        total_jobs=len(jobs),
        status_summary=dict(status_summary),
        jobs=job_responses
    )


@router.get("/jobs/{job_id}/qc-tasks")
async def get_desolvation_qc_tasks(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取某个去溶剂化任务的 QC 子任务列表
    用于展示多级计算结构
    """
    # 检查任务是否存在
    job = db.query(PostprocessJob).filter(
        PostprocessJob.id == job_id,
        PostprocessJob.job_type == PostprocessType.DESOLVATION_ENERGY
    ).first()

    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="去溶剂化任务不存在"
        )

    # 从 config 中获取 qc_job_ids
    qc_job_ids = job.config.get("qc_job_ids", [])

    # 如果没有 qc_job_ids，尝试从其他字段获取
    if not qc_job_ids:
        # 尝试获取 cluster_qc_job_id
        cluster_qc_job_id = job.config.get("cluster_qc_job_id")
        if cluster_qc_job_id:
            qc_job_ids.append(cluster_qc_job_id)

        # 获取 cluster_minus_job_ids
        cluster_minus_ids = job.config.get("cluster_minus_job_ids", [])
        qc_job_ids.extend(cluster_minus_ids)

        # 获取 ligand_qc_jobs
        ligand_qc_jobs = job.config.get("ligand_qc_jobs", {})
        qc_job_ids.extend(ligand_qc_jobs.values())

        # 获取 center_ion_job_id
        center_ion_job_id = job.config.get("center_ion_job_id")
        if center_ion_job_id:
            qc_job_ids.append(center_ion_job_id)

    # 查询关联的 QC 任务
    qc_jobs = []
    if qc_job_ids:
        qc_jobs = db.query(QCJob).filter(
            QCJob.id.in_(qc_job_ids)
        ).order_by(QCJob.created_at.asc()).all()

    # 如果还是没有，尝试用 desolvation_postprocess_job_id 查询
    if not qc_jobs:
        qc_jobs = db.query(QCJob).filter(
            QCJob.desolvation_postprocess_job_id == job_id
        ).order_by(QCJob.created_at.asc()).all()

    # 构建响应
    qc_tasks = []
    for qc in qc_jobs:
        # 从 molecule_name 提取类型信息
        # 格式如: Cluster_995_minus_PF6_2 或 EC 或 Li
        mol_name = qc.molecule_name
        task_type = "cluster"  # 默认
        if mol_name.startswith("Cluster_") and "_minus_" in mol_name:
            task_type = "cluster_minus"
        elif not mol_name.startswith("Cluster_"):
            task_type = "ligand"

        qc_tasks.append({
            "id": qc.id,
            "molecule_name": mol_name,
            "task_type": task_type,  # cluster, cluster_minus, ligand
            "status": qc.status.value,
            "progress": qc.progress,
            "charge": qc.charge,
            "spin_multiplicity": qc.spin_multiplicity,
            "basis_set": qc.basis_set,
            "functional": qc.functional,
            "is_reused": qc.is_reused,
            "reused_from_job_id": qc.reused_from_job_id,
            "slurm_job_id": qc.slurm_job_id,
            "error_message": qc.error_message,
            "created_at": qc.created_at.isoformat() if qc.created_at else None,
            "started_at": qc.started_at.isoformat() if qc.started_at else None,
            "finished_at": qc.finished_at.isoformat() if qc.finished_at else None,
        })

    # 统计
    total = len(qc_tasks)
    completed = sum(1 for t in qc_tasks if t["status"] == "COMPLETED")
    running = sum(1 for t in qc_tasks if t["status"] == "RUNNING")
    failed = sum(1 for t in qc_tasks if t["status"] == "FAILED")
    queued = sum(1 for t in qc_tasks if t["status"] in ["QUEUED", "SUBMITTED", "CREATED"])
    reused = sum(1 for t in qc_tasks if t["is_reused"])

    return {
        "job_id": job_id,
        "composition_key": job.config.get("solvation_structure_id"),
        "total": total,
        "completed": completed,
        "running": running,
        "failed": failed,
        "queued": queued,
        "reused": reused,
        "qc_tasks": qc_tasks
    }


def _build_job_response(job: PostprocessJob, db: Session) -> DesolvationJobResponse:
    """构建任务响应，包含溯源信息和 QC 进度"""
    elapsed_seconds = None
    if job.started_at and job.finished_at:
        elapsed_seconds = (job.finished_at - job.started_at).total_seconds()

    # 获取溯源信息
    solvation_structure_id = job.config.get("solvation_structure_id")
    composition_key = None
    electrolyte_name = None

    if solvation_structure_id:
        solvation = db.query(SolvationStructure).filter(
            SolvationStructure.id == solvation_structure_id
        ).first()
        if solvation:
            # 构建 composition_key，如 "Li-EC2-DMC1-PF6_1"
            composition = solvation.composition or {}
            parts = [solvation.center_ion] if solvation.center_ion else []
            for mol, count in sorted(composition.items()):
                if count > 0:
                    parts.append(f"{mol}{count}")
            composition_key = "-".join(parts)

    # 获取电解液名称
    if job.md_job_id:
        md_job = db.query(MDJob).filter(MDJob.id == job.md_job_id).first()
        if md_job and md_job.system_id:
            electrolyte = db.query(ElectrolyteSystem).filter(
                ElectrolyteSystem.id == md_job.system_id
            ).first()
            if electrolyte:
                electrolyte_name = electrolyte.name

    # 获取 QC 任务进度（从 config 中读取 qc_job_ids）
    qc_progress = None
    if job.status in [JobStatus.RUNNING, JobStatus.QUEUED, JobStatus.SUBMITTED, JobStatus.POSTPROCESSING]:
        # 从 config 获取 qc_job_ids
        qc_job_ids = job.config.get("qc_job_ids", [])
        if not qc_job_ids:
            # 尝试从其他字段获取
            cluster_qc_job_id = job.config.get("cluster_qc_job_id")
            if cluster_qc_job_id:
                qc_job_ids.append(cluster_qc_job_id)
            qc_job_ids.extend(job.config.get("cluster_minus_job_ids", []))
            ligand_qc_jobs = job.config.get("ligand_qc_jobs", {})
            qc_job_ids.extend(ligand_qc_jobs.values())
            center_ion_job_id = job.config.get("center_ion_job_id")
            if center_ion_job_id:
                qc_job_ids.append(center_ion_job_id)

        # 查询关联的 QC 任务
        qc_jobs = []
        if qc_job_ids:
            qc_jobs = db.query(QCJob).filter(QCJob.id.in_(qc_job_ids)).all()

        if qc_jobs:
            total = len(qc_jobs)
            completed = sum(1 for qc in qc_jobs if qc.status.value == 'COMPLETED')
            running = sum(1 for qc in qc_jobs if qc.status.value == 'RUNNING')
            failed = sum(1 for qc in qc_jobs if qc.status.value == 'FAILED')
            qc_progress = {
                "total": total,
                "completed": completed,
                "running": running,
                "failed": failed,
                "progress_percent": round(completed / total * 100) if total > 0 else 0
            }

    result = None
    if job.status == JobStatus.COMPLETED:
        result_obj = db.query(DesolvationEnergyResult).filter(
            DesolvationEnergyResult.postprocess_job_id == job.id
        ).first()

        if result_obj:
            # 转换 per_ligand_results
            per_ligand_results = []
            if result_obj.per_ligand_results:
                for item in result_obj.per_ligand_results:
                    per_ligand_results.append(LigandDesolvationResult(**item))

            # 转换 per_type_summary
            per_type_summary = []
            if result_obj.per_type_summary:
                for item in result_obj.per_type_summary:
                    per_type_summary.append(TypeSummary(**item))

            result = DesolvationEnergyResultSchema(
                id=result_obj.id,
                postprocess_job_id=result_obj.postprocess_job_id,
                solvation_structure_id=result_obj.solvation_structure_id,
                method_level=result_obj.method_level,
                basis_set=result_obj.basis_set,
                functional=result_obj.functional,
                e_cluster=result_obj.e_cluster,
                per_ligand_results=per_ligand_results,
                per_type_summary=per_type_summary,
                created_at=result_obj.created_at
            )

    return DesolvationJobResponse(
        job_id=job.id,
        status=job.status.value,
        method_level=job.config.get("method_level", "fast_xtb"),
        desolvation_mode=job.config.get("desolvation_mode", "stepwise"),
        solvent_config=job.config.get("solvent_config"),
        created_at=job.created_at,
        started_at=job.started_at,
        finished_at=job.finished_at,
        elapsed_seconds=elapsed_seconds,
        error_message=job.error_message,
        result=result,
        # 溯源信息
        solvation_structure_id=solvation_structure_id,
        composition_key=composition_key,
        md_job_id=job.md_job_id,
        electrolyte_name=electrolyte_name,
        qc_progress=qc_progress
    )


@router.get("/preview/{structure_id}")
async def preview_desolvation_structures(
    structure_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
) -> Dict[str, Any]:
    """
    预览去溶剂化结构（不创建任务）

    返回完整 cluster 和每个 cluster_minus 的 XYZ 结构，
    用于用户在提交前检查结构是否正确。
    """
    from app.tasks.desolvation import (
        parse_solvation_cluster,
        generate_cluster_minus_xyz,
        generate_center_ion_xyz,
        generate_dimer_xyz
    )
    from app.models.result import ResultSummary

    # 获取溶剂化结构
    solvation_structure = db.query(SolvationStructure).filter(
        SolvationStructure.id == structure_id
    ).first()

    if not solvation_structure:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Solvation structure {structure_id} not found"
        )

    # 获取 molecule_structures
    result_summary = db.query(ResultSummary).filter(
        ResultSummary.md_job_id == solvation_structure.md_job_id
    ).first()

    molecule_structures = []
    if result_summary and result_summary.molecule_structures:
        molecule_structures = result_summary.molecule_structures

    # 解析 cluster
    try:
        cluster_data = parse_solvation_cluster(solvation_structure, molecule_structures)
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to parse cluster structure: {str(e)}"
        )

    # 生成配位组成名称
    composition = solvation_structure.composition or {}
    composition_str = "".join([
        f"{name}{count}" for name, count in sorted(composition.items()) if count > 0
    ])
    cluster_name = f"{cluster_data['center_ion']}_{composition_str}_{structure_id}"

    # 构建响应
    response = {
        "structure_id": structure_id,
        "cluster_name": cluster_name,
        "center_ion": cluster_data['center_ion'],
        "total_charge": cluster_data['total_charge'],
        "composition": composition,
        # 完整 cluster
        "cluster": {
            "name": cluster_name,
            "xyz_content": cluster_data['xyz_content'],
            "atom_count": len(cluster_data['all_atoms']),
            "charge": cluster_data['total_charge']
        },
        # 每个配体
        "ligands": [],
        # 每个 cluster_minus
        "cluster_minus_structures": [],
        # 每种配体类型的 dimer（Li + 配体）结构（用于 pairwise binding）
        "dimer_structures": [],
        # 中心离子（用于 full 模式）
        "center_ion_structure": {
            "name": cluster_data['center_ion'],
            "xyz_content": generate_center_ion_xyz(cluster_data),
            "atom_count": 1,
            "charge": cluster_data['center_ion_charge']
        }
    }

    # 添加配体信息
    for ligand in cluster_data['ligands']:
        response["ligands"].append({
            "ligand_id": ligand['ligand_id'],
            "ligand_type": ligand['ligand_type'],
            "ligand_label": ligand['ligand_label'],
            "xyz_content": ligand['xyz_content'],
            "atom_count": len(ligand['atom_indices']),
            "charge": ligand['charge']
        })

    # 统计每种配体类型的数量，用于标记等价配体
    ligand_type_counts = {}
    for ligand in cluster_data['ligands']:
        ltype = ligand['ligand_type']
        ligand_type_counts[ltype] = ligand_type_counts.get(ltype, 0) + 1

    # 记录已处理的配体类型（用于标记等价配体）
    processed_types = set()

    # 生成每个 cluster_minus 结构
    for ligand in cluster_data['ligands']:
        cluster_minus_xyz = generate_cluster_minus_xyz(cluster_data, ligand)
        cluster_minus_name = f"{cluster_name}_minus_{ligand['ligand_label']}"
        cluster_minus_charge = cluster_data['total_charge'] - ligand['charge']

        # 计算原子数
        xyz_lines = cluster_minus_xyz.strip().split('\n')
        atom_count = int(xyz_lines[0]) if xyz_lines else 0

        # 判断是否是等价配体（同类型配体超过1个）
        ltype = ligand['ligand_type']
        is_equivalent = ligand_type_counts.get(ltype, 0) > 1
        is_representative = ltype not in processed_types  # 第一个作为代表
        processed_types.add(ltype)

        response["cluster_minus_structures"].append({
            "name": cluster_minus_name,
            "removed_ligand": ligand['ligand_label'],
            "removed_ligand_type": ligand['ligand_type'],
            "xyz_content": cluster_minus_xyz,
            "atom_count": atom_count,
            "charge": cluster_minus_charge,
            "is_equivalent": is_equivalent,  # 是否有等价配体
            "is_representative": is_representative,  # 是否是该类型的代表（第一个）
            "equivalent_count": ligand_type_counts.get(ltype, 1)  # 该类型的配体数量
        })

    # 添加配体类型统计到响应
    response["ligand_type_summary"] = ligand_type_counts

    # 生成每种配体类型的 dimer 结构（Li + 配体）
    # 每种类型只生成一个代表性的 dimer
    dimer_processed_types = set()
    for ligand in cluster_data['ligands']:
        ltype = ligand['ligand_type']
        if ltype not in dimer_processed_types:
            dimer_processed_types.add(ltype)
            try:
                dimer_xyz = generate_dimer_xyz(cluster_data, ligand)
                dimer_name = f"{cluster_data['center_ion']}-{ltype}"
                dimer_charge = cluster_data['center_ion_charge'] + ligand['charge']

                xyz_lines = dimer_xyz.strip().split('\n')
                atom_count = int(xyz_lines[0]) if xyz_lines else 0

                response["dimer_structures"].append({
                    "name": dimer_name,
                    "ligand_type": ltype,
                    "xyz_content": dimer_xyz,
                    "atom_count": atom_count,
                    "charge": dimer_charge,
                    "source_ligand_id": ligand['ligand_id']
                })
            except Exception as e:
                logger.warning(f"Failed to generate dimer for {ltype}: {e}")

    return response


@router.get("/jobs/{job_id}/binding-summary")
async def get_binding_energy_summary(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取去溶剂化任务的 Binding Energy 汇总
    从已有的 per_ligand_results 派生，不需要额外计算

    返回:
    - per_ligand_binding: 每个配体的 binding energy
    - per_type_stats: 按类型汇总的统计 (均值、标准差、最小、最大)
    - last_layer_binding: 最后一级去溶剂化能（最内层 binding）
    """
    import numpy as np

    # 获取任务
    job = db.query(PostprocessJob).filter(
        PostprocessJob.id == job_id,
        PostprocessJob.job_type == PostprocessType.DESOLVATION_ENERGY
    ).first()

    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="去溶剂化任务不存在"
        )

    # 获取结果
    result = db.query(DesolvationEnergyResult).filter(
        DesolvationEnergyResult.postprocess_job_id == job_id
    ).first()

    if not result:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="去溶剂化结果不存在，任务可能未完成"
        )

    per_ligand_results = result.per_ligand_results or []

    # 1. 提取每个配体的 binding energy
    per_ligand_binding = []
    for ligand in per_ligand_results:
        per_ligand_binding.append({
            "ligand_id": ligand.get("ligand_id"),
            "ligand_type": ligand.get("ligand_type"),
            "ligand_label": ligand.get("ligand_label"),
            "binding_energy_kcal": ligand.get("delta_e"),  # kcal/mol
            "e_ligand_au": ligand.get("e_ligand"),
            "e_cluster_minus_au": ligand.get("e_cluster_minus")
        })

    # 2. 按类型统计
    type_energies = defaultdict(list)
    for ligand in per_ligand_binding:
        if ligand["binding_energy_kcal"] is not None:
            type_energies[ligand["ligand_type"]].append(ligand["binding_energy_kcal"])

    per_type_stats = {}
    for ltype, energies in type_energies.items():
        arr = np.array(energies)
        per_type_stats[ltype] = {
            "count": len(energies),
            "mean": float(np.mean(arr)),
            "std": float(np.std(arr)) if len(arr) > 1 else 0.0,
            "min": float(np.min(arr)),
            "max": float(np.max(arr)),
            "values": energies
        }

    # 3. 找出"最后一级去溶剂化"（最内层 binding）
    # 定义：最大的 binding energy（最难脱除的配体）
    last_layer_binding = None
    if per_ligand_binding:
        valid_bindings = [l for l in per_ligand_binding if l["binding_energy_kcal"] is not None]
        if valid_bindings:
            # 按 binding energy 绝对值排序，取最大的
            last_layer = max(valid_bindings, key=lambda x: abs(x["binding_energy_kcal"]))
            last_layer_binding = {
                "ligand_type": last_layer["ligand_type"],
                "ligand_label": last_layer["ligand_label"],
                "binding_energy_kcal": last_layer["binding_energy_kcal"]
            }

    # 4. 获取溶剂化结构信息
    solvation_info = None
    if result.solvation_structure_id:
        solvation = db.query(SolvationStructure).filter(
            SolvationStructure.id == result.solvation_structure_id
        ).first()
        if solvation:
            solvation_info = {
                "center_ion": solvation.center_ion,
                "coordination_num": solvation.coordination_num,
                "composition_key": solvation.composition_key
            }

    return {
        "job_id": job_id,
        "method_level": result.method_level,
        "e_cluster_au": result.e_cluster,
        "solvation_info": solvation_info,
        "per_ligand_binding": per_ligand_binding,
        "per_type_stats": per_type_stats,
        "last_layer_binding": last_layer_binding,
        "total_ligands": len(per_ligand_binding)
    }


@router.get("/overview/{md_job_id}/binding-statistics")
async def get_binding_statistics_overview(
    md_job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取某个 MD 任务下所有已完成去溶剂化任务的 Binding Energy 统计汇总
    用于总览分析面板
    """
    import numpy as np

    # 获取所有已完成的去溶剂化任务
    completed_jobs = db.query(PostprocessJob).filter(
        PostprocessJob.md_job_id == md_job_id,
        PostprocessJob.job_type == PostprocessType.DESOLVATION_ENERGY,
        PostprocessJob.status == JobStatus.COMPLETED
    ).all()

    if not completed_jobs:
        return {
            "md_job_id": md_job_id,
            "total_completed": 0,
            "message": "没有已完成的去溶剂化任务"
        }

    # 收集所有配体的 binding energy
    all_bindings = []  # 所有配体
    type_bindings = defaultdict(list)  # 按类型分组
    last_layer_bindings = []  # 最后一级

    for job in completed_jobs:
        result = db.query(DesolvationEnergyResult).filter(
            DesolvationEnergyResult.postprocess_job_id == job.id
        ).first()

        if not result or not result.per_ligand_results:
            continue

        # 获取溶剂化结构信息
        composition_key = None
        if result.solvation_structure_id:
            solvation = db.query(SolvationStructure).filter(
                SolvationStructure.id == result.solvation_structure_id
            ).first()
            if solvation:
                composition_key = solvation.composition_key

        # 提取每个配体的 binding
        valid_ligands = []
        for ligand in result.per_ligand_results:
            delta_e = ligand.get("delta_e")
            if delta_e is not None:
                binding_data = {
                    "job_id": job.id,
                    "composition_key": composition_key,
                    "ligand_type": ligand.get("ligand_type"),
                    "binding_energy_kcal": delta_e
                }
                all_bindings.append(binding_data)
                type_bindings[ligand.get("ligand_type")].append(delta_e)
                valid_ligands.append(binding_data)

        # 找出该 cluster 的最后一级
        if valid_ligands:
            last_layer = max(valid_ligands, key=lambda x: abs(x["binding_energy_kcal"]))
            last_layer_bindings.append(last_layer)

    # 统计
    def calc_stats(values):
        if not values:
            return None
        arr = np.array(values)
        return {
            "count": len(values),
            "mean": float(np.mean(arr)),
            "std": float(np.std(arr)) if len(arr) > 1 else 0.0,
            "min": float(np.min(arr)),
            "max": float(np.max(arr)),
            "percentile_25": float(np.percentile(arr, 25)),
            "percentile_75": float(np.percentile(arr, 75))
        }

    # 按类型统计
    per_type_overall = {}
    for ltype, values in type_bindings.items():
        per_type_overall[ltype] = calc_stats(values)
        per_type_overall[ltype]["values"] = values  # 用于前端画分布图

    # 最后一级统计
    last_layer_stats = None
    if last_layer_bindings:
        last_values = [l["binding_energy_kcal"] for l in last_layer_bindings]
        last_layer_stats = calc_stats(last_values)
        last_layer_stats["details"] = last_layer_bindings

    return {
        "md_job_id": md_job_id,
        "total_completed": len(completed_jobs),
        "total_ligand_bindings": len(all_bindings),
        "per_type_statistics": per_type_overall,
        "last_layer_statistics": last_layer_stats
    }

