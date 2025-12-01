"""
Worker API 端点

用于轮询 Worker 与云端后端通信
"""

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
from typing import List, Optional, Dict, Any
from datetime import datetime
import logging

from app.database import get_db
from app.models.job import MDJob, JobStatus
from app.models.qc import QCJob, QCJobStatus
from app.models.user import User, UserRole
from app.dependencies import get_current_user
from pydantic import BaseModel


def is_worker_user(user: User) -> bool:
    """检查用户是否是 Worker 用户（ADMIN 角色）"""
    return user.role == UserRole.ADMIN

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/workers", tags=["workers"])


# ==================== 全局缓存 ====================
# 存储 Worker 上报的分区信息
_partition_cache: Dict[str, Any] = {
    "partitions": [],
    "last_updated": None,
    "worker_name": None,
}


def get_cached_partitions() -> List[dict]:
    """获取缓存的分区信息"""
    return _partition_cache["partitions"]


def get_partition_cache_info() -> Dict[str, Any]:
    """获取分区缓存的完整信息"""
    return _partition_cache.copy()


# ==================== Schemas ====================

class PartitionReport(BaseModel):
    """分区信息上报"""
    name: str
    state: str
    total_nodes: int
    available_nodes: int
    total_cpus: int
    available_cpus: int
    max_time: Optional[str] = None


class WorkerHeartbeat(BaseModel):
    """Worker 心跳"""
    worker_name: str
    status: str
    running_jobs: int
    timestamp: str
    partitions: Optional[List[PartitionReport]] = None  # 可选：同时上报分区信息


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

    Worker 定期发送心跳，表明自己在线。
    可以同时上报分区信息。
    """
    # 验证是否是 Worker 用户（ADMIN 角色）
    if not is_worker_user(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only worker users can send heartbeat"
        )

    logger.info(
        f"Received heartbeat from {heartbeat.worker_name}: "
        f"{heartbeat.running_jobs} jobs running"
    )

    # 如果心跳中包含分区信息，更新缓存
    if heartbeat.partitions:
        _partition_cache["partitions"] = [p.dict() for p in heartbeat.partitions]
        _partition_cache["last_updated"] = datetime.now().isoformat()
        _partition_cache["worker_name"] = heartbeat.worker_name
        logger.info(f"Updated partition cache with {len(heartbeat.partitions)} partitions from {heartbeat.worker_name}")

    return {"status": "ok", "timestamp": datetime.now().isoformat()}


@router.post("/partitions")
async def report_partitions(
    partitions: List[PartitionReport],
    current_user: User = Depends(get_current_user),
):
    """
    Worker 上报分区信息

    Worker 定期调用此接口上报校园网集群的分区状态。
    云端会缓存这些信息，供前端获取。
    """
    # 验证是否是 Worker 用户（ADMIN 角色）
    if not is_worker_user(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only worker users can report partitions"
        )

    # 更新分区缓存
    _partition_cache["partitions"] = [p.dict() for p in partitions]
    _partition_cache["last_updated"] = datetime.now().isoformat()
    _partition_cache["worker_name"] = current_user.username

    logger.info(f"Received {len(partitions)} partitions from worker {current_user.username}")

    return {
        "status": "ok",
        "partitions_count": len(partitions),
        "timestamp": datetime.now().isoformat()
    }


@router.get("/partitions")
async def get_worker_partitions(
    current_user: User = Depends(get_current_user),
):
    """
    获取 Worker 上报的分区信息

    返回校园网 Worker 最近上报的分区状态。
    """
    cache_info = get_partition_cache_info()

    return {
        "partitions": cache_info["partitions"],
        "last_updated": cache_info["last_updated"],
        "worker_name": cache_info["worker_name"],
    }


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
    # 验证是否是 Worker 用户（ADMIN 角色）
    if not is_worker_user(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only worker users can fetch pending jobs"
        )
    
    job_type = job_type.upper()
    
    if job_type == "MD":
        # 获取 CREATED 状态的 MD 任务（待处理）
        from app.models.electrolyte import ElectrolyteSystem

        jobs = db.query(MDJob).filter(
            MDJob.status == JobStatus.CREATED
        ).order_by(MDJob.created_at).limit(limit).all()

        result = []
        for job in jobs:
            # 获取电解质系统数据
            electrolyte = db.query(ElectrolyteSystem).filter(
                ElectrolyteSystem.id == job.system_id
            ).first()

            if electrolyte:
                # 合并 config 和电解质数据
                job_config = {
                    **(job.config or {}),
                    "name": job.config.get("job_name") if job.config else f"MD-{job.id}",
                    "cations": electrolyte.cations,
                    "anions": electrolyte.anions,
                    "solvents": electrolyte.solvents,
                    "additives": getattr(electrolyte, 'additives', None),
                    "box_size": electrolyte.box_size,
                    "temperature": electrolyte.temperature,
                    "pressure": electrolyte.pressure,
                }
            else:
                job_config = job.config

            result.append(PendingJobResponse(
                id=job.id,
                type="MD",
                config=job_config,
                created_at=job.created_at
            ))

        return result
    
    elif job_type == "QC":
        # 获取 CREATED 状态的 QC 任务（待处理）
        jobs = db.query(QCJob).filter(
            QCJob.status == QCJobStatus.CREATED
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
                    "slurm_partition": job.slurm_partition or "cpu",
                    "slurm_cpus": job.slurm_cpus or 16,
                    "slurm_time": job.slurm_time or 7200,
                    **(job.config or {}),  # 包含额外的配置信息
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
    # 验证是否是 Worker 用户（ADMIN 角色）
    if not is_worker_user(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only worker users can update job status"
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


@router.get("/jobs/{job_id}/check_cancelled")
async def check_job_cancelled(
    job_id: int,
    job_type: str = "MD",
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    检查任务是否已被用户取消

    Worker 定期调用此接口检查任务是否需要取消
    """
    # 验证是否是 Worker 用户（ADMIN 角色）
    if not is_worker_user(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only worker users can check job cancellation"
        )

    job_type = job_type.upper()

    if job_type == "MD":
        job = db.query(MDJob).filter(MDJob.id == job_id).first()
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"MD Job {job_id} not found"
            )

        # 检查是否已被取消
        cancelled = job.status == JobStatus.CANCELLED

        return {
            "job_id": job_id,
            "type": "MD",
            "cancelled": cancelled,
            "status": job.status.value
        }

    elif job_type == "QC":
        job = db.query(QCJob).filter(QCJob.id == job_id).first()
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"QC Job {job_id} not found"
            )

        # 检查是否已被取消
        cancelled = job.status == QCJobStatus.CANCELLED

        return {
            "job_id": job_id,
            "type": "QC",
            "cancelled": cancelled,
            "status": job.status.value
        }

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
    # 验证是否是 Worker 用户（ADMIN 角色）
    if not is_worker_user(current_user):
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

        # 获取电解质系统数据
        from app.models.electrolyte import ElectrolyteSystem
        electrolyte = db.query(ElectrolyteSystem).filter(
            ElectrolyteSystem.id == job.system_id
        ).first()

        if not electrolyte:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Electrolyte system {job.system_id} not found"
            )

        # 合并 config 和电解质数据
        job_data = {
            **(job.config or {}),
            "name": job.config.get("job_name") if job.config else f"MD-{job_id}",
            "cations": electrolyte.cations,
            "anions": electrolyte.anions,
            "solvents": electrolyte.solvents,
            "additives": getattr(electrolyte, 'additives', None),
            "box_size": electrolyte.box_size,
            "temperature": electrolyte.temperature,
            "pressure": electrolyte.pressure,
        }

        return {
            "job_id": job_id,
            "type": "MD",
            "config": job_data,
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

