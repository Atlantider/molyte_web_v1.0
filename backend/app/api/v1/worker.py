"""
Worker API 端点

用于轮询 Worker 与阿里云后端通信
"""

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
from typing import List, Optional
from datetime import datetime
import logging

from app.database import get_db
from app.models.job import MDJob, JobStatus
from app.models.qc import QCJob, QCJobStatus
from app.models.user import User
from app.dependencies import get_current_user
from pydantic import BaseModel

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/workers", tags=["workers"])


# ==================== Schemas ====================

class WorkerHeartbeat(BaseModel):
    """Worker 心跳"""
    worker_name: str
    status: str
    running_jobs: int
    timestamp: str


class JobStatusUpdate(BaseModel):
    """任务状态更新"""
    status: str
    job_type: str  # MD or QC
    worker_name: str
    slurm_job_id: Optional[str] = None
    work_dir: Optional[str] = None
    error_message: Optional[str] = None
    result_files: Optional[List[str]] = None


class PendingJobResponse(BaseModel):
    """待处理任务响应"""
    id: int
    type: str  # MD or QC
    config: dict
    created_at: datetime


# ==================== API Endpoints ====================

@router.post("/heartbeat")
async def worker_heartbeat(
    heartbeat: WorkerHeartbeat,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Worker 心跳接口
    
    Worker 定期发送心跳，表明自己在线
    """
    # 验证是否是 Worker 用户（可以添加特殊的 Worker 角色）
    if current_user.user_type != "admin":
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only worker users can send heartbeat"
        )
    
    logger.info(
        f"Received heartbeat from {heartbeat.worker_name}: "
        f"{heartbeat.running_jobs} jobs running"
    )
    
    # TODO: 可以将心跳信息存储到数据库或 Redis
    # 用于监控 Worker 状态
    
    return {"status": "ok", "timestamp": datetime.now().isoformat()}


@router.get("/jobs/pending", response_model=List[PendingJobResponse])
async def get_pending_jobs(
    job_type: str = "MD",  # MD or QC
    limit: int = 10,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    获取待处理的任务

    Worker 轮询此接口获取新任务
    """
    # 验证是否是管理员用户（Worker 用户应该是 ADMIN 角色）
    from app.models.user import UserRole
    if current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only admin/worker users can fetch pending jobs"
        )
    
    job_type = job_type.upper()
    
    if job_type == "MD":
        # 获取 PENDING 状态的 MD 任务
        jobs = db.query(MDJob).filter(
            MDJob.status == JobStatus.PENDING
        ).order_by(MDJob.created_at).limit(limit).all()
        
        return [
            PendingJobResponse(
                id=job.id,
                type="MD",
                config=job.config,
                created_at=job.created_at
            )
            for job in jobs
        ]
    
    elif job_type == "QC":
        # 获取 PENDING 状态的 QC 任务
        jobs = db.query(QCJob).filter(
            QCJob.status == QCJobStatus.PENDING
        ).order_by(QCJob.created_at).limit(limit).all()
        
        return [
            PendingJobResponse(
                id=job.id,
                type="QC",
                config={
                    "molecule_name": job.molecule_name,
                    "smiles": job.smiles,
                    "basis_set": job.basis_set,
                    "functional": job.functional,
                    "charge": job.charge,
                    "spin_multiplicity": job.spin_multiplicity,
                    "solvent_model": job.solvent_model,
                    "solvent_name": job.solvent_name,
                },
                created_at=job.created_at
            )
            for job in jobs
        ]
    
    else:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid job type: {job_type}"
        )


@router.put("/jobs/{job_id}/status")
async def update_job_status(
    job_id: int,
    status_update: JobStatusUpdate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    更新任务状态

    Worker 在任务状态变化时调用此接口
    """
    # 验证是否是管理员用户（Worker 用户应该是 ADMIN 角色）
    from app.models.user import UserRole
    if current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only admin/worker users can update job status"
        )
    
    job_type = status_update.job_type.upper()
    
    if job_type == "MD":
        job = db.query(MDJob).filter(MDJob.id == job_id).first()
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"MD Job {job_id} not found"
            )
        
        # 更新状态
        job.status = JobStatus[status_update.status]
        
        if status_update.slurm_job_id:
            job.slurm_job_id = status_update.slurm_job_id
        
        if status_update.work_dir:
            job.work_dir = status_update.work_dir
        
        if status_update.error_message:
            job.error_message = status_update.error_message
        
        if status_update.status == "RUNNING" and not job.started_at:
            job.started_at = datetime.now()
        
        if status_update.status in ["COMPLETED", "FAILED"]:
            job.finished_at = datetime.now()
        
        # 如果有结果文件，可以存储到 config 中
        if status_update.result_files:
            if not job.config:
                job.config = {}
            job.config['result_files'] = status_update.result_files
        
        db.commit()
        
        logger.info(
            f"MD Job {job_id} status updated to {status_update.status} "
            f"by worker {status_update.worker_name}"
        )
        
        return {"status": "ok", "job_id": job_id, "new_status": status_update.status}
    
    elif job_type == "QC":
        job = db.query(QCJob).filter(QCJob.id == job_id).first()
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"QC Job {job_id} not found"
            )
        
        # 更新状态
        job.status = QCJobStatus[status_update.status]
        
        if status_update.slurm_job_id:
            job.slurm_job_id = status_update.slurm_job_id
        
        if status_update.work_dir:
            job.work_dir = status_update.work_dir
        
        if status_update.error_message:
            job.error_message = status_update.error_message
        
        if status_update.status == "RUNNING" and not job.started_at:
            job.started_at = datetime.now()
        
        if status_update.status in ["COMPLETED", "FAILED"]:
            job.finished_at = datetime.now()
        
        # 如果有结果文件，可以存储到 config 中
        if status_update.result_files:
            if not job.config:
                job.config = {}
            job.config['result_files'] = status_update.result_files
        
        db.commit()
        
        logger.info(
            f"QC Job {job_id} status updated to {status_update.status} "
            f"by worker {status_update.worker_name}"
        )
        
        return {"status": "ok", "job_id": job_id, "new_status": status_update.status}
    
    else:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid job type: {job_type}"
        )


@router.get("/jobs/{job_id}/input")
async def get_job_input_data(
    job_id: int,
    job_type: str = "MD",
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    获取任务输入数据
    
    Worker 可以通过此接口下载任务所需的输入数据
    """
    # 验证是否是 Worker 用户
    if current_user.user_type != "admin":
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only worker users can fetch job input data"
        )
    
    job_type = job_type.upper()
    
    if job_type == "MD":
        job = db.query(MDJob).filter(MDJob.id == job_id).first()
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"MD Job {job_id} not found"
            )
        
        return {
            "job_id": job_id,
            "type": "MD",
            "config": job.config,
            "system_id": job.system_id
        }
    
    elif job_type == "QC":
        job = db.query(QCJob).filter(QCJob.id == job_id).first()
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"QC Job {job_id} not found"
            )
        
        return {
            "job_id": job_id,
            "type": "QC",
            "molecule_name": job.molecule_name,
            "smiles": job.smiles,
            "basis_set": job.basis_set,
            "functional": job.functional,
            "charge": job.charge,
            "spin_multiplicity": job.spin_multiplicity,
            "solvent_model": job.solvent_model,
            "solvent_name": job.solvent_name,
        }
    
    else:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid job type: {job_type}"
        )

