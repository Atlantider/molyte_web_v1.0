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
            "desolvation_mode": batch_data.desolvation_mode
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

    # 查询关联的 QC 任务
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

    # 获取 QC 任务进度（通过 desolvation_postprocess_job_id 查询）
    qc_progress = None
    if job.status in [JobStatus.RUNNING, JobStatus.QUEUED, JobStatus.SUBMITTED]:
        # 查询关联的 QC 任务
        qc_jobs = db.query(QCJob).filter(
            QCJob.desolvation_postprocess_job_id == job.id
        ).all()

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

