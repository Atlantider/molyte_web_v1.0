"""
Slurm 任务状态监控 Celery Worker

定时检查所有 QUEUED/RUNNING 任务的 Slurm 状态并更新数据库
"""

import logging
import subprocess
from datetime import datetime
from typing import Dict, List, Optional

from sqlalchemy.orm import Session

from app.celery_app import celery_app
from app.database import SessionLocal
from app.models.job import MDJob, JobStatus

logger = logging.getLogger(__name__)


def get_slurm_job_status(slurm_job_id: str) -> Optional[Dict[str, str]]:
    """
    查询 Slurm 任务状态
    
    Args:
        slurm_job_id: Slurm 任务 ID
        
    Returns:
        Dict 包含任务状态信息，如果查询失败返回 None
    """
    try:
        # 使用 sacct 查询任务状态
        cmd = [
            "sacct",
            "-j", slurm_job_id,
            "--format=JobID,State,ExitCode,Start,End",
            "--parsable2",
            "--noheader",
        ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=10,
        )
        
        if result.returncode != 0:
            logger.warning(f"sacct failed for job {slurm_job_id}: {result.stderr}")
            return None
        
        # 解析输出（取第一行，即主任务）
        lines = result.stdout.strip().split("\n")
        if not lines or not lines[0]:
            logger.warning(f"No sacct output for job {slurm_job_id}")
            return None
        
        fields = lines[0].split("|")
        if len(fields) < 5:
            logger.warning(f"Invalid sacct output for job {slurm_job_id}: {lines[0]}")
            return None
        
        return {
            "job_id": fields[0],
            "state": fields[1],
            "exit_code": fields[2],
            "start_time": fields[3],
            "end_time": fields[4],
        }
        
    except subprocess.TimeoutExpired:
        logger.error(f"sacct timeout for job {slurm_job_id}")
        return None
    except Exception as e:
        logger.exception(f"Error querying Slurm status for job {slurm_job_id}: {e}")
        return None


def map_slurm_status_to_job_status(slurm_state: str) -> JobStatus:
    """
    将 Slurm 状态映射到 JobStatus

    Slurm 状态参考：
    - PENDING: 等待资源
    - RUNNING: 正在运行
    - COMPLETED: 成功完成
    - FAILED: 失败
    - CANCELLED: 被取消
    - TIMEOUT: 超时
    - OUT_OF_MEMORY: 内存不足
    - NODE_FAIL: 节点故障
    """
    slurm_state = slurm_state.upper()

    # 处理 "CANCELLED BY xxx" 格式
    if slurm_state.startswith("CANCELLED"):
        return JobStatus.CANCELLED

    if slurm_state in ["PENDING", "CONFIGURING", "RESIZING"]:
        return JobStatus.QUEUED
    elif slurm_state in ["RUNNING", "COMPLETING"]:
        return JobStatus.RUNNING
    elif slurm_state == "COMPLETED":
        return JobStatus.COMPLETED
    elif slurm_state in ["FAILED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED"]:
        return JobStatus.FAILED
    elif slurm_state in ["BOOT_FAIL", "DEADLINE", "REVOKED"]:
        return JobStatus.CANCELLED
    else:
        logger.warning(f"Unknown Slurm state: {slurm_state}, treating as FAILED")
        return JobStatus.FAILED


@celery_app.task(
    name="app.tasks.job_monitor.monitor_slurm_jobs_task",
    bind=True,
)
def monitor_slurm_jobs_task(self) -> Dict[str, int]:
    """
    定时监控所有 QUEUED/RUNNING 任务的 Slurm 状态
    
    Returns:
        Dict 包含监控统计信息
    """
    db = SessionLocal()
    
    try:
        logger.info(f"[Task {self.request.id}] Starting Slurm job monitoring...")
        
        # 查询所有 QUEUED 或 RUNNING 状态的任务
        jobs = db.query(MDJob).filter(
            MDJob.status.in_([JobStatus.QUEUED, JobStatus.RUNNING])
        ).all()
        
        if not jobs:
            logger.info(f"[Task {self.request.id}] No jobs to monitor")
            return {"total": 0, "updated": 0, "failed": 0}
        
        logger.info(f"[Task {self.request.id}] Monitoring {len(jobs)} jobs")
        
        updated_count = 0
        failed_count = 0
        
        for job in jobs:
            if not job.slurm_job_id:
                logger.warning(f"[Task {self.request.id}] Job {job.id} has no slurm_job_id, skipping")
                continue
            
            # 查询 Slurm 状态
            slurm_status = get_slurm_job_status(job.slurm_job_id)
            
            if not slurm_status:
                logger.warning(f"[Task {self.request.id}] Failed to get status for job {job.id} (slurm_id={job.slurm_job_id})")
                failed_count += 1
                continue
            
            # 映射状态
            new_status = map_slurm_status_to_job_status(slurm_status["state"])
            
            # 如果状态发生变化，更新数据库
            if new_status != job.status:
                old_status = job.status
                job.status = new_status
                
                # 记录状态变化
                if "status_history" not in job.config:
                    job.config["status_history"] = []
                job.config["status_history"].append({
                    "from": old_status.value,
                    "to": new_status.value,
                    "timestamp": datetime.now().isoformat(),
                    "slurm_state": slurm_status["state"],
                })
                
                # 如果任务失败，记录错误信息
                if new_status == JobStatus.FAILED:
                    job.error_message = f"Slurm job failed with state: {slurm_status['state']}, exit_code: {slurm_status['exit_code']}"

                # 如果任务完成，触发后处理
                if new_status == JobStatus.COMPLETED:
                    from app.tasks.postprocess import postprocess_md_job_task
                    logger.info(f"[Task {self.request.id}] Triggering postprocessing for job {job.id}")
                    postprocess_md_job_task.delay(job.id)

                db.commit()
                updated_count += 1
                
                logger.info(
                    f"[Task {self.request.id}] Job {job.id} status changed: "
                    f"{old_status.value} → {new_status.value} (Slurm: {slurm_status['state']})"
                )
        
        logger.info(
            f"[Task {self.request.id}] Monitoring complete: "
            f"total={len(jobs)}, updated={updated_count}, failed={failed_count}"
        )
        
        return {
            "total": len(jobs),
            "updated": updated_count,
            "failed": failed_count,
        }
        
    except Exception as exc:
        logger.exception(f"[Task {self.request.id}] Exception during job monitoring: {exc}")
        raise
    finally:
        db.close()

