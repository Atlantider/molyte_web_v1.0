"""
Celery 应用配置

用于异步任务处理：
- MD 任务提交到 Slurm
- Slurm 状态监控
- 后处理任务
"""

from celery import Celery
from celery.schedules import crontab
import os
from app.core.config import settings

# 创建 Celery 应用
celery_app = Celery(
    "molyte_tasks",
    broker=os.getenv("CELERY_BROKER_URL", "redis://localhost:6379/0"),
    backend=os.getenv("CELERY_RESULT_BACKEND", "redis://localhost:6379/1"),
)

# Celery 配置
celery_app.conf.update(
    # 序列化
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    
    # 时区
    timezone="Asia/Shanghai",
    enable_utc=True,
    
    # 任务结果过期时间（24小时）
    result_expires=86400,
    
    # 任务超时（30分钟）
    task_time_limit=1800,
    task_soft_time_limit=1700,
    
    # 任务重试
    task_acks_late=True,
    task_reject_on_worker_lost=True,
    
    # Worker 配置
    worker_prefetch_multiplier=1,
    worker_max_tasks_per_child=1000,
    
    # 队列配置 - 暂时都使用默认队列
    # task_routes={
    #     "app.tasks.job_submission.submit_md_job_task": {"queue": "submit"},
    #     "app.tasks.job_monitor.monitor_slurm_jobs_task": {"queue": "monitor"},
    #     "app.tasks.postprocess.postprocess_md_job_task": {"queue": "postprocess"},
    # },
    
    # 定时任务配置（Celery Beat）
    beat_schedule={
        # 每 30 秒监控一次 MD 任务状态
        "monitor-slurm-jobs": {
            "task": "app.tasks.job_monitor.monitor_slurm_jobs_task",
            "schedule": 30.0,  # 30 秒
        },
        # 每 30 秒监控一次 QC 任务状态
        "monitor-qc-jobs": {
            "task": "app.tasks.qc_monitor.monitor_qc_jobs",
            "schedule": 30.0,  # 30 秒
        },
    },
)

# 自动发现任务模块
celery_app.autodiscover_tasks(["app.tasks"])

# 日志配置
celery_app.conf.update(
    worker_log_format="[%(asctime)s: %(levelname)s/%(processName)s] %(message)s",
    worker_task_log_format="[%(asctime)s: %(levelname)s/%(processName)s] [%(task_name)s(%(task_id)s)] %(message)s",
)


@celery_app.task(bind=True)
def debug_task(self):
    """调试任务，用于测试 Celery 是否正常工作"""
    print(f"Request: {self.request!r}")
    return {"status": "ok", "message": "Celery is working!"}

