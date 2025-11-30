"""
QC 任务监控 Celery Worker

负责监控 QC 任务在 Slurm 集群上的运行状态
"""

import logging
import subprocess
from datetime import datetime
from typing import Dict, Any, List

from celery import Task
from sqlalchemy.orm import Session
from sqlalchemy import or_

from app.celery_app import celery_app
from app.database import SessionLocal
from app.models.qc import QCJob, QCJobStatus

logger = logging.getLogger(__name__)


class DatabaseTask(Task):
    """自动管理数据库会话的任务基类"""
    _db: Session = None

    def after_return(self, *args, **kwargs):
        if self._db is not None:
            self._db.close()
            self._db = None

    @property
    def db(self) -> Session:
        if self._db is None:
            self._db = SessionLocal()
        return self._db


def get_slurm_job_status(slurm_job_id: str) -> Dict[str, Any]:
    """
    查询Slurm任务状态
    
    Returns:
        Dict with keys: state, exit_code, reason
    """
    try:
        result = subprocess.run(
            ["squeue", "-j", slurm_job_id, "-h", "-o", "%T|%r"],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode == 0 and result.stdout.strip():
            parts = result.stdout.strip().split("|")
            state = parts[0] if parts else "UNKNOWN"
            reason = parts[1] if len(parts) > 1 else ""
            return {"state": state, "reason": reason, "running": True}
        
        # 任务不在队列中，检查sacct
        result = subprocess.run(
            ["sacct", "-j", slurm_job_id, "-n", "-o", "State,ExitCode", "-P"],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split("\n")
            for line in lines:
                if line.strip():
                    parts = line.split("|")
                    state = parts[0] if parts else "UNKNOWN"
                    exit_code = parts[1] if len(parts) > 1 else "0:0"
                    return {"state": state, "exit_code": exit_code, "running": False}
        
        return {"state": "UNKNOWN", "running": False}
        
    except subprocess.TimeoutExpired:
        logger.warning(f"Timeout querying Slurm job {slurm_job_id}")
        return {"state": "TIMEOUT", "running": True}
    except Exception as e:
        logger.error(f"Error querying Slurm job {slurm_job_id}: {e}")
        return {"state": "ERROR", "error": str(e), "running": True}


@celery_app.task(
    bind=True,
    base=DatabaseTask,
    name="app.tasks.qc_monitor.monitor_qc_jobs",
)
def monitor_qc_jobs(self) -> Dict[str, Any]:
    """
    监控所有运行中的QC任务
    
    这个任务应该被定期调度（如每分钟一次）
    """
    db = self.db
    
    try:
        # 查询所有需要监控的任务
        jobs = db.query(QCJob).filter(
            or_(
                QCJob.status == QCJobStatus.QUEUED,
                QCJob.status == QCJobStatus.RUNNING
            )
        ).all()
        
        if not jobs:
            logger.debug("No QC jobs to monitor")
            return {"monitored": 0, "updated": 0}
        
        logger.info(f"Monitoring {len(jobs)} QC jobs")
        
        updated_count = 0
        completed_jobs = []
        
        for job in jobs:
            if not job.slurm_job_id:
                continue
            
            status = get_slurm_job_status(job.slurm_job_id)
            slurm_state = status.get("state", "UNKNOWN")
            
            logger.debug(f"QC Job {job.id} (Slurm {job.slurm_job_id}): {slurm_state}")
            
            # 更新任务状态
            if slurm_state in ["RUNNING"]:
                if job.status != QCJobStatus.RUNNING:
                    job.status = QCJobStatus.RUNNING
                    job.progress = 50.0  # 运行中设为50%
                    updated_count += 1
                    
            elif slurm_state in ["PENDING", "CONFIGURING"]:
                if job.status != QCJobStatus.QUEUED:
                    job.status = QCJobStatus.QUEUED
                    updated_count += 1
                    
            elif slurm_state in ["COMPLETED"]:
                job.status = QCJobStatus.POSTPROCESSING
                job.progress = 80.0
                job.finished_at = datetime.now()
                updated_count += 1
                completed_jobs.append(job.id)
                
            elif slurm_state in ["FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL"]:
                job.status = QCJobStatus.FAILED
                job.error_message = f"Slurm job {slurm_state}: {status.get('reason', '')}"
                job.finished_at = datetime.now()
                updated_count += 1
        
        db.commit()
        
        # 触发后处理任务
        for job_id in completed_jobs:
            from app.tasks.qc_postprocess import postprocess_qc_job
            postprocess_qc_job.delay(job_id)
            logger.info(f"Triggered postprocessing for QC job {job_id}")
        
        return {
            "monitored": len(jobs),
            "updated": updated_count,
            "completed": len(completed_jobs)
        }
        
    except Exception as exc:
        logger.exception(f"Error monitoring QC jobs: {exc}")
        return {"error": str(exc)}

