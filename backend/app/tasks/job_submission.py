"""
MD 任务提交 Celery Worker

负责异步提交 MD 任务到 Slurm 集群
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Any

from celery import Task
from sqlalchemy.orm import Session

from app.celery_app import celery_app
from app.database import SessionLocal
from app.models.job import MDJob, JobStatus
from app.models.electrolyte import ElectrolyteSystem
from app.workers.molyte_wrapper import MolyteWrapper
from app.core.config import settings

logger = logging.getLogger(__name__)


class DatabaseTask(Task):
    """自动管理数据库会话的任务基类"""
    _db: Session = None

    def after_return(self, *args, **kwargs):
        if self._db is not None:
            self._db.close()

    @property
    def db(self) -> Session:
        if self._db is None:
            self._db = SessionLocal()
        return self._db


@celery_app.task(
    bind=True,
    base=DatabaseTask,
    name="app.tasks.job_submission.submit_md_job_task",
    max_retries=3,
    default_retry_delay=60,  # 60秒后重试
)
def submit_md_job_task(self, job_id: int) -> Dict[str, Any]:
    """
    异步提交 MD 任务到 Slurm 集群
    
    Args:
        job_id: MD 任务 ID
        
    Returns:
        Dict 包含任务提交结果
    """
    db = self.db
    
    try:
        logger.info(f"[Task {self.request.id}] Starting MD job submission for job_id={job_id}")
        
        # 1. 从数据库读取任务信息
        job = db.query(MDJob).filter(MDJob.id == job_id).first()
        if not job:
            error_msg = f"Job {job_id} not found in database"
            logger.error(f"[Task {self.request.id}] {error_msg}")
            return {"success": False, "error": error_msg}
        
        # 2. 读取电解液系统信息
        electrolyte = db.query(ElectrolyteSystem).filter(
            ElectrolyteSystem.id == job.system_id
        ).first()
        if not electrolyte:
            error_msg = f"Electrolyte system {job.system_id} not found"
            logger.error(f"[Task {self.request.id}] {error_msg}")
            job.status = JobStatus.FAILED
            job.error_message = error_msg
            db.commit()
            return {"success": False, "error": error_msg}
        
        logger.info(f"[Task {self.request.id}] Job: {job.config.get('job_name')}, Electrolyte: {electrolyte.name}")
        
        # 3. 准备任务数据
        job_name = job.config.get("job_name")
        job_data = {
            "name": job_name,  # molyte_wrapper 期望的键名是 "name"
            "job_name": job_name,  # 保留兼容性
            "cations": electrolyte.cations,
            "anions": electrolyte.anions,
            "solvents": electrolyte.solvents,
            "temperature": job.config.get("temperature", 298.15),
            "pressure": job.config.get("pressure", 1.0),
            "box_size": electrolyte.box_size if electrolyte.box_size else 40.0,  # 从电解质系统获取盒子大小，默认40
            "nsteps_npt": job.config.get("nsteps_npt", 5000000),
            "nsteps_nvt": job.config.get("nsteps_nvt", 10000000),
            "dump_freq": job.config.get("dump_freq", 1000),
            "thermo_freq": job.config.get("thermo_freq", 100),
            "restart_freq": job.config.get("restart_freq", 10000),
            "charge_method": job.config.get("charge_method", "cm1a"),
            "timestep": job.config.get("timestep", 2.0),
            # Slurm 资源配置
            "slurm_partition": job.config.get("slurm_partition", "cpu"),
            "slurm_nodes": job.config.get("slurm_nodes", 1),
            "slurm_ntasks": job.config.get("slurm_ntasks", 8),
            "slurm_cpus_per_task": job.config.get("slurm_cpus_per_task", 8),
            "slurm_time": job.config.get("slurm_time", 7200),
        }
        
        # 4. 初始化 MolyteWrapper
        wrapper = MolyteWrapper(
            work_base_path=settings.MOLYTE_WORK_BASE_PATH,
            initial_salts_path=settings.MOLYTE_INITIAL_SALTS_PATH,
            ligpargen_path=settings.MOLYTE_LIGPARGEN_PATH,
            packmol_path=settings.MOLYTE_PACKMOL_PATH,
            ltemplify_path=settings.MOLYTE_LTEMPLIFY_PATH,
            moltemplate_path=settings.MOLYTE_MOLTEMPLATE_PATH,
            charge_save_path=settings.MOLYTE_CHARGE_SAVE_PATH,
        )
        
        # 5. 生成 LAMMPS 输入文件
        logger.info(f"[Task {self.request.id}] Generating LAMMPS input files...")
        result = wrapper.generate_lammps_input(job_data)
        
        if not result.get("success"):
            error_msg = result.get("error", "Unknown error during LAMMPS input generation")
            logger.error(f"[Task {self.request.id}] {error_msg}")
            job.status = JobStatus.FAILED
            job.error_message = error_msg
            db.commit()
            return {"success": False, "error": error_msg}
        
        work_dir = result.get("work_dir")
        logger.info(f"[Task {self.request.id}] LAMMPS input files generated in {work_dir}")
        
        # 6. 提交到 Slurm
        logger.info(f"[Task {self.request.id}] Submitting to Slurm...")
        slurm_result = wrapper.submit_to_slurm(work_dir)
        
        if not slurm_result.get("success"):
            error_msg = slurm_result.get("error", "Unknown error during Slurm submission")
            logger.error(f"[Task {self.request.id}] {error_msg}")
            job.status = JobStatus.FAILED
            job.error_message = error_msg
            db.commit()
            return {"success": False, "error": error_msg}
        
        slurm_job_id = slurm_result.get("slurm_job_id")
        logger.info(f"[Task {self.request.id}] Submitted to Slurm with job_id={slurm_job_id}")
        
        # 7. 更新数据库
        job.status = JobStatus.QUEUED
        job.slurm_job_id = str(slurm_job_id)
        job.work_dir = str(work_dir)
        job.error_message = None
        job.config["work_dir"] = str(work_dir)
        job.config["slurm_job_id"] = str(slurm_job_id)
        job.config["submitted_at"] = datetime.now().isoformat()
        db.commit()
        
        logger.info(f"[Task {self.request.id}] Job {job_id} submitted successfully")
        
        return {
            "success": True,
            "job_id": job_id,
            "slurm_job_id": slurm_job_id,
            "work_dir": str(work_dir),
        }
        
    except Exception as exc:
        logger.exception(f"[Task {self.request.id}] Exception during job submission: {exc}")
        
        # 更新任务状态为失败
        try:
            job = db.query(MDJob).filter(MDJob.id == job_id).first()
            if job:
                job.status = JobStatus.FAILED
                job.error_message = f"Task exception: {str(exc)}"
                db.commit()
        except Exception as db_exc:
            logger.error(f"[Task {self.request.id}] Failed to update job status: {db_exc}")
        
        # 重试任务
        raise self.retry(exc=exc, countdown=60)

