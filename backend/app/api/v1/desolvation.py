"""
Desolvation energy calculation API endpoints
"""
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
from typing import List
from datetime import datetime

from app.database import get_db
from app.models.user import User
from app.models.job import PostprocessJob, PostprocessType, JobStatus
from app.models.result import SolvationStructure, DesolvationEnergyResult
from app.schemas.desolvation import (
    DesolvationJobCreate,
    DesolvationJobResponse,
    DesolvationEnergyResultSchema,
    LigandDesolvationResult,
    TypeSummary
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
            "method_level": job_data.method_level
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


def _build_job_response(job: PostprocessJob, db: Session) -> DesolvationJobResponse:
    """构建任务响应"""
    elapsed_seconds = None
    if job.started_at and job.finished_at:
        elapsed_seconds = (job.finished_at - job.started_at).total_seconds()
    
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
        created_at=job.created_at,
        started_at=job.started_at,
        finished_at=job.finished_at,
        elapsed_seconds=elapsed_seconds,
        error_message=job.error_message,
        result=result
    )

