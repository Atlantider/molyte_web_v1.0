"""
User quota checking and management utilities
"""
from datetime import date, datetime
from sqlalchemy.orm import Session
from sqlalchemy import func
from app.models.user import User, UserRole
from app.models.job import MDJob, JobStatus
from app.models.user_stats import UserUsageStats
from app.utils.timezone import ensure_timezone_aware


def calculate_job_cpu_hours(job: MDJob) -> float:
    """
    Calculate CPU hours used by a job

    Args:
        job: MD job instance

    Returns:
        float: CPU hours used
    """
    if not job.started_at:
        return 0.0

    # Calculate elapsed time - 确保时区一致
    end_time = ensure_timezone_aware(job.finished_at) if job.finished_at else ensure_timezone_aware(datetime.now())
    start_time = ensure_timezone_aware(job.started_at)
    elapsed_hours = (end_time - start_time).total_seconds() / 3600.0

    # Get number of cores
    cores = 1
    if job.config:
        ntasks = job.config.get('slurm_ntasks', 1)
        cpus_per_task = job.config.get('slurm_cpus_per_task', 1)
        cores = ntasks * cpus_per_task

    return elapsed_hours * cores


def calculate_user_cpu_hours(user_id: int, db: Session) -> float:
    """
    Calculate total CPU hours used by a user
    
    Args:
        user_id: User ID
        db: Database session
        
    Returns:
        float: Total CPU hours used
    """
    jobs = db.query(MDJob).filter(
        MDJob.user_id == user_id,
        MDJob.status.in_([JobStatus.COMPLETED, JobStatus.RUNNING, JobStatus.FAILED])
    ).all()
    
    total_hours = sum(calculate_job_cpu_hours(job) for job in jobs)
    return total_hours


def get_today_job_count(user_id: int, db: Session) -> int:
    """
    Get number of jobs submitted today by a user
    
    Args:
        user_id: User ID
        db: Database session
        
    Returns:
        int: Number of jobs submitted today
    """
    today = date.today()
    count = db.query(MDJob).filter(
        MDJob.user_id == user_id,
        func.date(MDJob.created_at) == today
    ).count()
    
    return count


def get_running_job_count(user_id: int, db: Session) -> int:
    """
    Get number of currently running/queued jobs by a user
    
    Args:
        user_id: User ID
        db: Database session
        
    Returns:
        int: Number of running/queued jobs
    """
    count = db.query(MDJob).filter(
        MDJob.user_id == user_id,
        MDJob.status.in_([JobStatus.RUNNING, JobStatus.QUEUED, JobStatus.POSTPROCESSING])
    ).count()
    
    return count


def update_user_usage_stats(user_id: int, job: MDJob, db: Session):
    """
    Update user usage statistics when a job changes status

    Args:
        user_id: User ID
        job: MD job instance
        db: Database session
    """
    today = date.today()

    # Get or create today's stats
    stats = db.query(UserUsageStats).filter(
        UserUsageStats.user_id == user_id,
        UserUsageStats.date == today
    ).first()

    if not stats:
        stats = UserUsageStats(
            user_id=user_id,
            date=today
        )
        db.add(stats)

    # Update based on job status
    if job.status == JobStatus.COMPLETED:
        stats.jobs_completed += 1
        stats.cpu_hours_used += calculate_job_cpu_hours(job)
    elif job.status == JobStatus.FAILED:
        stats.jobs_failed += 1
        stats.cpu_hours_used += calculate_job_cpu_hours(job)
    elif job.status == JobStatus.CANCELLED:
        stats.jobs_cancelled += 1
        stats.cpu_hours_used += calculate_job_cpu_hours(job)

    db.commit()


def check_user_quota(user: User, db: Session) -> dict:
    """
    Check if user has quota to submit a new job
    
    Args:
        user: User instance
        db: Database session
        
    Returns:
        dict: {
            "allowed": bool,
            "reason": str (if not allowed),
            "details": dict (quota details)
        }
    """
    # Admin has unlimited quota
    if user.role == UserRole.ADMIN:
        return {
            "allowed": True,
            "details": {
                "role": "admin",
                "unlimited": True
            }
        }
    
    # Check if user is active
    if not user.is_active:
        return {
            "allowed": False,
            "reason": "User account is disabled",
            "details": {}
        }
    
    # 1. Check CPU hours quota
    used_hours = calculate_user_cpu_hours(user.id, db)
    if used_hours >= user.total_cpu_hours:
        return {
            "allowed": False,
            "reason": f"CPU 核时超出配额，请充值 ({used_hours:.2f}/{user.total_cpu_hours:.2f} hours)",
            "details": {
                "used_hours": used_hours,
                "total_hours": user.total_cpu_hours,
                "remaining_hours": 0
            }
        }
    
    # 2. Check daily job limit
    today_jobs = get_today_job_count(user.id, db)
    if today_jobs >= user.daily_job_limit:
        return {
            "allowed": False,
            "reason": f"Daily job limit exceeded ({today_jobs}/{user.daily_job_limit} jobs today)",
            "details": {
                "today_jobs": today_jobs,
                "daily_limit": user.daily_job_limit
            }
        }
    
    # 3. Check concurrent job limit
    running_jobs = get_running_job_count(user.id, db)
    if running_jobs >= user.concurrent_job_limit:
        return {
            "allowed": False,
            "reason": f"任务书超上限，请联系管理员升级 ({running_jobs}/{user.concurrent_job_limit} running jobs)",
            "details": {
                "running_jobs": running_jobs,
                "concurrent_limit": user.concurrent_job_limit
            }
        }
    
    # All checks passed
    return {
        "allowed": True,
        "details": {
            "used_hours": used_hours,
            "total_hours": user.total_cpu_hours,
            "remaining_hours": user.total_cpu_hours - used_hours,
            "today_jobs": today_jobs,
            "daily_limit": user.daily_job_limit,
            "running_jobs": running_jobs,
            "concurrent_limit": user.concurrent_job_limit
        }
    }

