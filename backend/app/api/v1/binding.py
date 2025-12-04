"""
Binding Analysis API routes
Li-配体 Binding Energy 分析相关的 API 接口
"""
import logging
from typing import List, Optional
from datetime import datetime

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from sqlalchemy import desc

from app.database import get_db
from app.models.user import User
from app.models.job import MDJob, BindingAnalysisJob as BindingAnalysisJobModel, BindingAnalysisStatus as BindingStatusModel
from app.schemas.binding import (
    BindingAnalysisCreate,
    BindingAnalysisJob,
    BindingAnalysisJobList,
    BindingAnalysisStatus,
)
from app.dependencies import get_current_active_user

logger = logging.getLogger(__name__)
router = APIRouter()


@router.post("/jobs", response_model=BindingAnalysisJob)
def create_binding_analysis_job(
    request: BindingAnalysisCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    创建 Binding 分析任务
    
    简化版 Li-配体 binding energy 计算：
    E_bind_shell = E_cluster - (E_center_ion + Σ n_j × E_ligand_j)
    """
    # 检查 MD 任务是否存在
    md_job = db.query(MDJob).filter(MDJob.id == request.md_job_id).first()
    if not md_job:
        raise HTTPException(status_code=404, detail=f"MD 任务 {request.md_job_id} 不存在")
    
    # 检查权限
    if md_job.user_id != current_user.id and current_user.role.value != 'admin':
        raise HTTPException(status_code=403, detail="无权访问此 MD 任务")
    
    # 创建任务
    job = BindingAnalysisJobModel(
        md_job_id=request.md_job_id,
        user_id=current_user.id,
        status=BindingStatusModel.CREATED,
        config=request.config.model_dump(),
        progress=0.0,
    )
    
    db.add(job)
    db.commit()
    db.refresh(job)
    
    logger.info(f"Created binding analysis job {job.id} for MD job {request.md_job_id}")
    
    return job


@router.get("/jobs", response_model=BindingAnalysisJobList)
def list_binding_analysis_jobs(
    md_job_id: Optional[int] = Query(None, description="筛选指定 MD 任务下的 binding 分析"),
    status: Optional[BindingAnalysisStatus] = Query(None, description="筛选状态"),
    skip: int = Query(0, ge=0),
    limit: int = Query(50, ge=1, le=100),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取 Binding 分析任务列表"""
    query = db.query(BindingAnalysisJobModel)
    
    # 权限过滤
    if current_user.role.value != 'admin':
        query = query.filter(BindingAnalysisJobModel.user_id == current_user.id)
    
    if md_job_id:
        query = query.filter(BindingAnalysisJobModel.md_job_id == md_job_id)
    
    if status:
        query = query.filter(BindingAnalysisJobModel.status == status.value)
    
    total = query.count()
    items = query.order_by(desc(BindingAnalysisJobModel.created_at)).offset(skip).limit(limit).all()
    
    return BindingAnalysisJobList(items=items, total=total)


@router.get("/jobs/{job_id}", response_model=BindingAnalysisJob)
def get_binding_analysis_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取 Binding 分析任务详情"""
    job = db.query(BindingAnalysisJobModel).filter(BindingAnalysisJobModel.id == job_id).first()
    
    if not job:
        raise HTTPException(status_code=404, detail=f"Binding 分析任务 {job_id} 不存在")
    
    # 权限检查
    if job.user_id != current_user.id and current_user.role.value != 'admin':
        raise HTTPException(status_code=403, detail="无权访问此任务")
    
    return job


@router.post("/jobs/{job_id}/submit")
def submit_binding_analysis_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """提交 Binding 分析任务"""
    job = db.query(BindingAnalysisJobModel).filter(BindingAnalysisJobModel.id == job_id).first()

    if not job:
        raise HTTPException(status_code=404, detail=f"Binding 分析任务 {job_id} 不存在")

    if job.user_id != current_user.id and current_user.role.value != 'admin':
        raise HTTPException(status_code=403, detail="无权操作此任务")

    if job.status != BindingStatusModel.CREATED:
        raise HTTPException(status_code=400, detail=f"任务状态为 {job.status}，无法提交")

    job.status = BindingStatusModel.SUBMITTED
    job.updated_at = datetime.utcnow()
    db.commit()

    return {"message": "任务已提交", "job_id": job_id, "status": "SUBMITTED"}


@router.delete("/jobs/{job_id}")
def delete_binding_analysis_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """删除 Binding 分析任务"""
    job = db.query(BindingAnalysisJobModel).filter(BindingAnalysisJobModel.id == job_id).first()

    if not job:
        raise HTTPException(status_code=404, detail=f"Binding 分析任务 {job_id} 不存在")

    if job.user_id != current_user.id and current_user.role.value != 'admin':
        raise HTTPException(status_code=403, detail="无权操作此任务")

    # 只允许删除 CREATED 或 FAILED 状态的任务
    if job.status not in [BindingStatusModel.CREATED, BindingStatusModel.FAILED, BindingStatusModel.COMPLETED]:
        raise HTTPException(status_code=400, detail=f"任务状态为 {job.status}，无法删除")

    db.delete(job)
    db.commit()

    return {"message": "任务已删除", "job_id": job_id}


@router.get("/available-clusters/{md_job_id}")
def get_available_clusters_for_binding(
    md_job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取可用于 Binding 分析的 Cluster 列表

    返回该 MD 任务下已有溶剂化结构（SolvationStructure）和已完成的 QC 计算的汇总
    """
    from app.models.result import SolvationStructure
    from app.models.qc import QCJob, QCJobStatus as QCJobStatusModel
    from collections import defaultdict

    # 检查 MD 任务
    md_job = db.query(MDJob).filter(MDJob.id == md_job_id).first()
    if not md_job:
        raise HTTPException(status_code=404, detail=f"MD 任务 {md_job_id} 不存在")

    if md_job.user_id != current_user.id and current_user.role.value != 'admin':
        raise HTTPException(status_code=403, detail="无权访问此 MD 任务")

    # 获取溶剂化结构
    solvation_structures = db.query(SolvationStructure).filter(
        SolvationStructure.md_job_id == md_job_id
    ).all()

    # 按 composition_key 分组
    clusters_by_composition = defaultdict(list)
    for struct in solvation_structures:
        key = struct.composition_key or f"cluster_{struct.id}"
        clusters_by_composition[key].append({
            "id": struct.id,
            "frame": struct.frame_number,
            "center_ion": struct.center_ion,
            "ligand_types": struct.ligand_types,
            "has_xyz": bool(struct.xyz_content),
        })

    # 获取已有的 QC 计算结果（按分子类型统计）
    from app.models.job import PostprocessJob, PostprocessType

    postprocess_jobs = db.query(PostprocessJob).filter(
        PostprocessJob.md_job_id == md_job_id,
        PostprocessJob.job_type == PostprocessType.DESOLVATION_ENERGY
    ).all()

    postprocess_ids = [pj.id for pj in postprocess_jobs]

    existing_qc = {}
    if postprocess_ids:
        qc_jobs = db.query(QCJob).filter(
            QCJob.desolvation_postprocess_job_id.in_(postprocess_ids),
            QCJob.status == QCJobStatusModel.COMPLETED
        ).all()

        for qc in qc_jobs:
            mol_type = qc.molecule_type or 'unknown'
            if mol_type not in existing_qc:
                existing_qc[mol_type] = {"count": 0, "examples": []}
            existing_qc[mol_type]["count"] += 1
            if len(existing_qc[mol_type]["examples"]) < 3:
                existing_qc[mol_type]["examples"].append(qc.molecule_name)

    return {
        "md_job_id": md_job_id,
        "total_solvation_structures": len(solvation_structures),
        "composition_keys": list(clusters_by_composition.keys()),
        "clusters_by_composition": dict(clusters_by_composition),
        "existing_qc_by_type": existing_qc,
        "note": "可选择 composition_keys 中的项目进行 Binding 分析"
    }

