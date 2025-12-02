"""
Admin API endpoints
"""
from typing import List, Optional
from datetime import datetime, date, timedelta
from fastapi import APIRouter, Depends, HTTPException, status, Request
from sqlalchemy.orm import Session
from sqlalchemy import func, desc
import logging

from app.database import get_db
from app.dependencies import get_current_admin_user
from app.models.user import User, UserRole
from app.models.job import MDJob, JobStatus
from app.models.user_stats import UserUsageStats, AuditLog
from app.schemas.admin import (
    UserListItem, UserDetail, UserUpdate, UserCreate,
    GlobalStats, UserUsageStatsItem, UserRanking,
    TrendDataPoint, StatisticsSummary, AuditLogItem,
    QuotaCheckResponse
)
from app.utils.quota import (
    calculate_user_cpu_hours, get_today_job_count,
    get_running_job_count, check_user_quota
)
from app.utils.audit import (
    log_user_update, log_quota_update, log_user_status_change,
    log_job_cancel, create_audit_log
)
from app.core.security import get_password_hash
from app.services.slurm import list_partitions

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/admin", tags=["admin"])


# ============ User Management Endpoints ============

@router.get("/users", response_model=List[UserListItem])
async def get_all_users(
    skip: int = 0,
    limit: int = 100,
    role: Optional[UserRole] = None,
    user_type: Optional[str] = None,  # UserType枚举值
    is_active: Optional[bool] = None,
    search: Optional[str] = None,  # 搜索用户名、邮箱、组织
    organization: Optional[str] = None,  # 按组织筛选
    sort_by: str = "created_at",  # created_at, username, email, last_login_at
    sort_order: str = "desc",  # asc, desc
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get all users with enhanced filtering (admin only)

    Args:
        skip: Number of records to skip
        limit: Maximum number of records to return
        role: Filter by role (ADMIN, PREMIUM, USER, GUEST)
        user_type: Filter by user type (STUDENT, RESEARCHER, COMPANY)
        is_active: Filter by active status
        search: Search in username, email, organization
        organization: Filter by organization name
        sort_by: Sort field (created_at, username, email, last_login_at)
        sort_order: Sort order (asc, desc)
        db: Database session
        admin: Current admin user

    Returns:
        List of users
    """
    from app.models.user import UserType

    query = db.query(User)

    # 角色筛选
    if role:
        query = query.filter(User.role == role)

    # 用户类型筛选
    if user_type:
        try:
            user_type_enum = UserType(user_type)
            query = query.filter(User.user_type == user_type_enum)
        except ValueError:
            pass  # 忽略无效的user_type值

    # 激活状态筛选
    if is_active is not None:
        query = query.filter(User.is_active == is_active)

    # 组织筛选
    if organization:
        query = query.filter(User.organization.ilike(f"%{organization}%"))

    # 搜索功能（用户名、邮箱、组织）
    if search:
        search_pattern = f"%{search}%"
        query = query.filter(
            (User.username.ilike(search_pattern)) |
            (User.email.ilike(search_pattern)) |
            (User.organization.ilike(search_pattern))
        )

    # 排序
    sort_column = {
        "created_at": User.created_at,
        "username": User.username,
        "email": User.email,
        "last_login_at": User.last_login_at,
    }.get(sort_by, User.created_at)

    if sort_order == "asc":
        query = query.order_by(sort_column.asc())
    else:
        query = query.order_by(sort_column.desc())

    users = query.offset(skip).limit(limit).all()
    return users


@router.get("/users/count/total")
async def get_users_count(
    role: Optional[UserRole] = None,
    user_type: Optional[str] = None,
    is_active: Optional[bool] = None,
    search: Optional[str] = None,
    organization: Optional[str] = None,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get total count of users with filters (admin only)

    Args:
        role: Filter by role
        user_type: Filter by user type
        is_active: Filter by active status
        search: Search in username, email, organization
        organization: Filter by organization name
        db: Database session
        admin: Current admin user

    Returns:
        Total count of users matching filters
    """
    from app.models.user import UserType

    query = db.query(User)

    if role:
        query = query.filter(User.role == role)

    if user_type:
        try:
            user_type_enum = UserType(user_type)
            query = query.filter(User.user_type == user_type_enum)
        except ValueError:
            pass

    if is_active is not None:
        query = query.filter(User.is_active == is_active)

    if organization:
        query = query.filter(User.organization.ilike(f"%{organization}%"))

    if search:
        search_pattern = f"%{search}%"
        query = query.filter(
            (User.username.ilike(search_pattern)) |
            (User.email.ilike(search_pattern)) |
            (User.organization.ilike(search_pattern))
        )

    total = query.count()
    return {"total": total}


@router.get("/users/{user_id}", response_model=UserDetail)
async def get_user_detail(
    user_id: int,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get detailed user information (admin only)
    
    Args:
        user_id: User ID
        db: Database session
        admin: Current admin user
        
    Returns:
        Detailed user information
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")
    
    # Calculate usage statistics
    used_cpu_hours = calculate_user_cpu_hours(user_id, db)
    today_jobs = get_today_job_count(user_id, db)
    running_jobs = get_running_job_count(user_id, db)
    
    total_jobs = db.query(MDJob).filter(MDJob.user_id == user_id).count()
    completed_jobs = db.query(MDJob).filter(
        MDJob.user_id == user_id,
        MDJob.status == JobStatus.COMPLETED
    ).count()
    failed_jobs = db.query(MDJob).filter(
        MDJob.user_id == user_id,
        MDJob.status == JobStatus.FAILED
    ).count()
    
    # Convert to dict and add statistics
    user_dict = {
        "id": user.id,
        "username": user.username,
        "email": user.email,
        "role": user.role,
        "is_active": user.is_active,
        "total_cpu_hours": user.total_cpu_hours,
        "daily_job_limit": user.daily_job_limit,
        "concurrent_job_limit": user.concurrent_job_limit,
        "storage_quota_gb": user.storage_quota_gb,
        "last_login_at": user.last_login_at,
        "created_at": user.created_at,
        "updated_at": user.updated_at,
        "used_cpu_hours": used_cpu_hours,
        "today_jobs": today_jobs,
        "running_jobs": running_jobs,
        "total_jobs": total_jobs,
        "completed_jobs": completed_jobs,
        "failed_jobs": failed_jobs
    }
    
    return UserDetail(**user_dict)


@router.put("/users/{user_id}", response_model=UserDetail)
async def update_user(
    user_id: int,
    user_update: UserUpdate,
    request: Request,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Update user information and quotas (admin only)
    
    Args:
        user_id: User ID
        user_update: User update data
        request: FastAPI request
        db: Database session
        admin: Current admin user
        
    Returns:
        Updated user information
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")
    
    # Track changes for audit log
    changes = {}
    old_quotas = {}
    new_quotas = {}
    
    # Update fields
    update_data = user_update.model_dump(exclude_unset=True)
    for field, value in update_data.items():
        if hasattr(user, field):
            old_value = getattr(user, field)
            if old_value != value:
                changes[field] = {"old": str(old_value), "new": str(value)}
                setattr(user, field, value)
                
                # Track quota changes
                if field in ["total_cpu_hours", "daily_job_limit", "concurrent_job_limit", "storage_quota_gb"]:
                    old_quotas[field] = old_value
                    new_quotas[field] = value
    
    db.commit()
    db.refresh(user)
    
    # Log the update
    if changes:
        log_user_update(db, admin, user, changes, request)
        
        if old_quotas:
            log_quota_update(db, admin, user, old_quotas, new_quotas, request)
    
    # Return detailed user info
    return await get_user_detail(user_id, db, admin)


@router.post("/users", response_model=UserDetail, status_code=status.HTTP_201_CREATED)
async def create_user(
    user_create: UserCreate,
    request: Request,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Create a new user (admin only)

    Args:
        user_create: User creation data
        request: FastAPI request
        db: Database session
        admin: Current admin user

    Returns:
        Created user information
    """
    # Check if username already exists
    existing_user = db.query(User).filter(User.username == user_create.username).first()
    if existing_user:
        raise HTTPException(status_code=400, detail="Username already exists")

    # Check if email already exists
    existing_email = db.query(User).filter(User.email == user_create.email).first()
    if existing_email:
        raise HTTPException(status_code=400, detail="Email already exists")

    # Create new user
    new_user = User(
        username=user_create.username,
        email=user_create.email,
        password_hash=get_password_hash(user_create.password),
        role=user_create.role,
        is_active=True,
        total_cpu_hours=user_create.total_cpu_hours,
        daily_job_limit=user_create.daily_job_limit,
        concurrent_job_limit=user_create.concurrent_job_limit,
        storage_quota_gb=user_create.storage_quota_gb,
        allowed_partitions=user_create.allowed_partitions
    )

    db.add(new_user)
    db.commit()
    db.refresh(new_user)

    # Log the creation
    create_audit_log(
        db=db,
        user=admin,
        action="create_user",
        resource_type="user",
        resource_id=new_user.id,
        details={
            "username": new_user.username,
            "email": new_user.email,
            "role": new_user.role.value
        },
        request=request
    )

    return await get_user_detail(new_user.id, db, admin)


@router.delete("/users/{user_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_user(
    user_id: int,
    request: Request,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Delete a user (admin only)

    Args:
        user_id: User ID
        request: FastAPI request
        db: Database session
        admin: Current admin user
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")

    # Prevent deleting yourself
    if user.id == admin.id:
        raise HTTPException(status_code=400, detail="Cannot delete yourself")

    # Log the deletion
    create_audit_log(
        db=db,
        user=admin,
        action="delete_user",
        resource_type="user",
        resource_id=user.id,
        details={
            "username": user.username,
            "email": user.email
        },
        request=request
    )

    db.delete(user)
    db.commit()


@router.put("/users/{user_id}/status")
async def update_user_status(
    user_id: int,
    is_active: bool,
    request: Request,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Enable or disable a user (admin only)

    Args:
        user_id: User ID
        is_active: New active status
        request: FastAPI request
        db: Database session
        admin: Current admin user

    Returns:
        Success message
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")

    # Prevent disabling yourself
    if user.id == admin.id and not is_active:
        raise HTTPException(status_code=400, detail="Cannot disable yourself")

    old_status = user.is_active
    user.is_active = is_active
    db.commit()

    # Log the status change
    log_user_status_change(db, admin, user, old_status, is_active, request)

    return {
        "message": f"User {'enabled' if is_active else 'disabled'} successfully",
        "user_id": user_id,
        "is_active": is_active
    }


@router.get("/users/{user_id}/quota", response_model=QuotaCheckResponse)
async def check_user_quota_endpoint(
    user_id: int,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Check user's quota status (admin only)

    Args:
        user_id: User ID
        db: Database session
        admin: Current admin user

    Returns:
        Quota check result
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")

    result = check_user_quota(user, db)
    return QuotaCheckResponse(**result)


# ============ Statistics and Monitoring Endpoints ============

@router.get("/stats/overview", response_model=GlobalStats)
async def get_global_stats(
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get global system statistics (admin only)

    Args:
        db: Database session
        admin: Current admin user

    Returns:
        Global statistics
    """
    # User statistics
    total_users = db.query(User).count()
    active_users = db.query(User).filter(User.is_active == True).count()

    # Job statistics
    total_jobs = db.query(MDJob).count()
    running_jobs = db.query(MDJob).filter(MDJob.status == JobStatus.RUNNING).count()
    queued_jobs = db.query(MDJob).filter(MDJob.status == JobStatus.QUEUED).count()
    completed_jobs = db.query(MDJob).filter(MDJob.status == JobStatus.COMPLETED).count()
    failed_jobs = db.query(MDJob).filter(MDJob.status == JobStatus.FAILED).count()

    # Resource statistics
    total_cpu_hours_allocated = db.query(func.sum(User.total_cpu_hours)).scalar() or 0.0
    total_storage_allocated_gb = db.query(func.sum(User.storage_quota_gb)).scalar() or 0.0

    # Calculate total CPU hours used - 优化版本：使用交易记录统计
    from app.models.billing import QuotaTransaction
    total_consumed_from_transactions = db.query(func.sum(QuotaTransaction.amount)).filter(
        QuotaTransaction.type == 'consume'
    ).scalar() or 0.0
    total_cpu_hours_used = abs(total_consumed_from_transactions)  # consume是负数

    # 如果没有交易记录（旧系统），回退到慢速计算（但只计算前100个用户）
    if total_cpu_hours_used == 0:
        logger.warning("No billing transactions found, using slow calculation method (limited to 100 users)")
        sample_users = db.query(User).limit(100).all()
        total_cpu_hours_used = sum(calculate_user_cpu_hours(u.id, db) for u in sample_users)

    # Storage used (placeholder - would need actual implementation)
    total_storage_used_gb = 0.0

    return GlobalStats(
        total_users=total_users,
        active_users=active_users,
        total_jobs=total_jobs,
        running_jobs=running_jobs,
        queued_jobs=queued_jobs,
        completed_jobs=completed_jobs,
        failed_jobs=failed_jobs,
        total_cpu_hours_used=total_cpu_hours_used,
        total_cpu_hours_allocated=total_cpu_hours_allocated,
        total_storage_used_gb=total_storage_used_gb,
        total_storage_allocated_gb=total_storage_allocated_gb
    )


@router.get("/stats/users", response_model=List[UserUsageStatsItem])
async def get_user_usage_stats(
    skip: int = 0,
    limit: int = 100,
    sort_by: str = "used_cpu_hours",  # used_cpu_hours, total_jobs, usage_percentage
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get usage statistics for all users (admin only)

    Args:
        skip: Number of records to skip
        limit: Maximum number of records to return
        sort_by: Sort field
        db: Database session
        admin: Current admin user

    Returns:
        List of user usage statistics
    """
    users = db.query(User).offset(skip).limit(limit).all()

    stats_list = []
    for user in users:
        used_cpu_hours = calculate_user_cpu_hours(user.id, db)
        total_jobs = db.query(MDJob).filter(MDJob.user_id == user.id).count()
        running_jobs = get_running_job_count(user.id, db)
        completed_jobs = db.query(MDJob).filter(
            MDJob.user_id == user.id,
            MDJob.status == JobStatus.COMPLETED
        ).count()
        failed_jobs = db.query(MDJob).filter(
            MDJob.user_id == user.id,
            MDJob.status == JobStatus.FAILED
        ).count()

        # Get last job time
        last_job = db.query(MDJob).filter(
            MDJob.user_id == user.id
        ).order_by(desc(MDJob.created_at)).first()
        last_job_at = last_job.created_at if last_job else None

        usage_percentage = (used_cpu_hours / user.total_cpu_hours * 100) if user.total_cpu_hours > 0 else 0

        stats_list.append(UserUsageStatsItem(
            user_id=user.id,
            username=user.username,
            email=user.email,
            role=user.role,
            used_cpu_hours=used_cpu_hours,
            total_cpu_hours=user.total_cpu_hours,
            usage_percentage=usage_percentage,
            total_jobs=total_jobs,
            running_jobs=running_jobs,
            completed_jobs=completed_jobs,
            failed_jobs=failed_jobs,
            last_job_at=last_job_at
        ))

    # Sort the results
    if sort_by == "used_cpu_hours":
        stats_list.sort(key=lambda x: x.used_cpu_hours, reverse=True)
    elif sort_by == "total_jobs":
        stats_list.sort(key=lambda x: x.total_jobs, reverse=True)
    elif sort_by == "usage_percentage":
        stats_list.sort(key=lambda x: x.usage_percentage, reverse=True)

    return stats_list


@router.get("/stats/ranking/cpu", response_model=List[UserRanking])
async def get_cpu_usage_ranking(
    limit: int = 10,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get top users by CPU hours usage (admin only)

    Args:
        limit: Number of top users to return
        db: Database session
        admin: Current admin user

    Returns:
        List of top users by CPU usage
    """
    from app.models.billing import QuotaTransaction

    users = db.query(User).all()

    user_cpu_hours = []
    for user in users:
        # 优先从交易记录计算总消耗
        total_consumed = db.query(func.sum(QuotaTransaction.amount)).filter(
            QuotaTransaction.user_id == user.id,
            QuotaTransaction.type == 'consume'
        ).scalar() or 0.0
        total_consumed = abs(total_consumed)  # consume是负数

        # 如果没有交易记录，使用旧的计算方法
        if total_consumed == 0:
            total_consumed = calculate_user_cpu_hours(user.id, db)

        user_cpu_hours.append({
            "user_id": user.id,
            "username": user.username,
            "email": user.email,
            "cpu_hours": total_consumed
        })

    # Sort by CPU hours
    user_cpu_hours.sort(key=lambda x: x["cpu_hours"], reverse=True)

    # Create ranking
    ranking = []
    for rank, item in enumerate(user_cpu_hours[:limit], start=1):
        ranking.append(UserRanking(
            user_id=item["user_id"],
            username=item["username"],
            email=item["email"],
            metric_value=item["cpu_hours"],
            rank=rank
        ))

    return ranking


@router.get("/stats/ranking/jobs", response_model=List[UserRanking])
async def get_job_count_ranking(
    limit: int = 10,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get top users by job count (admin only)

    Args:
        limit: Number of top users to return
        db: Database session
        admin: Current admin user

    Returns:
        List of top users by job count
    """
    # Query job counts grouped by user
    job_counts = db.query(
        MDJob.user_id,
        func.count(MDJob.id).label("job_count")
    ).group_by(MDJob.user_id).order_by(desc("job_count")).limit(limit).all()

    ranking = []
    for rank, (user_id, job_count) in enumerate(job_counts, start=1):
        user = db.query(User).filter(User.id == user_id).first()
        if user:
            ranking.append(UserRanking(
                user_id=user.id,
                username=user.username,
                email=user.email,
                metric_value=float(job_count),
                rank=rank
            ))

    return ranking


# ============ Job Management Endpoints ============

@router.get("/jobs/all")
async def get_all_jobs(
    skip: int = 0,
    limit: int = 100,
    status: Optional[JobStatus] = None,
    user_id: Optional[int] = None,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get all jobs from all users (admin only)

    Args:
        skip: Number of records to skip
        limit: Maximum number of records to return
        status: Filter by job status
        user_id: Filter by user ID
        db: Database session
        admin: Current admin user

    Returns:
        List of jobs
    """
    query = db.query(MDJob)

    if status:
        query = query.filter(MDJob.status == status)
    if user_id:
        query = query.filter(MDJob.user_id == user_id)

    jobs = query.order_by(desc(MDJob.created_at)).offset(skip).limit(limit).all()

    # Convert to dict and add user info
    result = []
    for job in jobs:
        job_dict = {
            "id": job.id,
            "user_id": job.user_id,
            "system_id": job.system_id,
            "status": job.status,
            "slurm_job_id": job.slurm_job_id,
            "created_at": job.created_at,
            "started_at": job.started_at,
            "finished_at": job.finished_at,
            "config": job.config
        }

        # Add user info
        user = db.query(User).filter(User.id == job.user_id).first()
        if user:
            job_dict["username"] = user.username
            job_dict["user_email"] = user.email

        result.append(job_dict)

    return result


@router.post("/jobs/{job_id}/cancel")
async def admin_cancel_job(
    job_id: int,
    reason: Optional[str] = None,
    request: Request = None,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Cancel a job (admin only)

    Args:
        job_id: Job ID
        reason: Reason for cancellation
        request: FastAPI request
        db: Database session
        admin: Current admin user

    Returns:
        Success message
    """
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # Check if job can be cancelled
    if job.status not in [JobStatus.QUEUED, JobStatus.RUNNING, JobStatus.POSTPROCESSING]:
        raise HTTPException(
            status_code=400,
            detail=f"Cannot cancel job with status {job.status}"
        )

    # Update job status
    job.status = JobStatus.CANCELLED
    job.finished_at = datetime.now()
    db.commit()

    # Log the cancellation
    log_job_cancel(db, admin, job_id, reason, request)

    return {
        "message": "Job cancelled successfully",
        "job_id": job_id,
        "reason": reason
    }


# ============ Audit Log Endpoints ============

@router.get("/logs", response_model=List[AuditLogItem])
async def get_audit_logs(
    skip: int = 0,
    limit: int = 100,
    action: Optional[str] = None,
    user_id: Optional[int] = None,
    resource_type: Optional[str] = None,
    start_date: Optional[datetime] = None,
    end_date: Optional[datetime] = None,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get audit logs (admin only)

    Args:
        skip: Number of records to skip
        limit: Maximum number of records to return
        action: Filter by action
        user_id: Filter by user ID
        resource_type: Filter by resource type
        start_date: Filter by start date
        end_date: Filter by end date
        db: Database session
        admin: Current admin user

    Returns:
        List of audit logs
    """
    query = db.query(AuditLog)

    if action:
        query = query.filter(AuditLog.action == action)
    if user_id:
        query = query.filter(AuditLog.user_id == user_id)
    if resource_type:
        query = query.filter(AuditLog.resource_type == resource_type)
    if start_date:
        query = query.filter(AuditLog.created_at >= start_date)
    if end_date:
        query = query.filter(AuditLog.created_at <= end_date)

    logs = query.order_by(desc(AuditLog.created_at)).offset(skip).limit(limit).all()

    # Add username to each log
    result = []
    for log in logs:
        log_dict = {
            "id": log.id,
            "user_id": log.user_id,
            "action": log.action,
            "resource_type": log.resource_type,
            "resource_id": log.resource_id,
            "details": log.details,
            "ip_address": log.ip_address,
            "created_at": log.created_at,
            "username": None
        }

        if log.user_id:
            user = db.query(User).filter(User.id == log.user_id).first()
            if user:
                log_dict["username"] = user.username

        result.append(AuditLogItem(**log_dict))

    return result


# ============ Statistics Report Endpoints ============

@router.get("/statistics/summary", response_model=StatisticsSummary)
async def get_statistics_summary(
    period: str = "7days",  # today, 7days, 30days
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """
    Get statistics summary for a period (admin only)

    Args:
        period: Time period (today, 7days, 30days)
        db: Database session
        admin: Current admin user

    Returns:
        Statistics summary
    """
    # Calculate date range
    end_date = datetime.now()
    if period == "today":
        start_date = datetime.now().replace(hour=0, minute=0, second=0, microsecond=0)
    elif period == "7days":
        start_date = end_date - timedelta(days=7)
    elif period == "30days":
        start_date = end_date - timedelta(days=30)
    else:
        raise HTTPException(status_code=400, detail="Invalid period")

    # Query jobs in the period
    jobs = db.query(MDJob).filter(
        MDJob.created_at >= start_date,
        MDJob.created_at <= end_date
    ).all()

    jobs_submitted = len(jobs)
    jobs_completed = sum(1 for j in jobs if j.status == JobStatus.COMPLETED)
    jobs_failed = sum(1 for j in jobs if j.status == JobStatus.FAILED)

    # Calculate CPU hours and average duration
    total_cpu_hours = sum(calculate_user_cpu_hours(j.user_id, db) for j in jobs if j.status == JobStatus.COMPLETED)

    completed_durations = []
    for job in jobs:
        if job.status == JobStatus.COMPLETED and job.started_at and job.finished_at:
            duration = (job.finished_at - job.started_at).total_seconds() / 3600.0
            completed_durations.append(duration)

    avg_duration = sum(completed_durations) / len(completed_durations) if completed_durations else 0.0

    # Peak concurrent jobs (simplified - would need time-series data for accurate calculation)
    peak_concurrent = db.query(MDJob).filter(
        MDJob.status.in_([JobStatus.RUNNING, JobStatus.QUEUED])
    ).count()

    return StatisticsSummary(
        period=period,
        jobs_submitted=jobs_submitted,
        jobs_completed=jobs_completed,
        jobs_failed=jobs_failed,
        cpu_hours_used=total_cpu_hours,
        avg_job_duration_hours=avg_duration,
        peak_concurrent_jobs=peak_concurrent
    )

# ============ Slurm Partition Management ============

@router.get("/partitions")
async def get_all_partitions(
    admin: User = Depends(get_current_admin_user)
):
    """
    Get all available Slurm partitions (admin only)

    This endpoint returns ALL partitions without filtering,
    used for admin to assign partition permissions to users.
    """
    try:
        partitions = list_partitions()

        return [
            {
                "name": p.name,
                "state": p.state,
                "total_nodes": p.total_nodes,
                "available_nodes": p.available_nodes,
                "total_cpus": p.total_cpus,
                "available_cpus": p.available_cpus,
                "max_time": p.max_time,
            }
            for p in partitions
        ]
    except Exception as e:
        # 如果 Slurm 命令失败，返回默认分区列表
        return [
            {
                "name": "cpu",
                "state": "up",
                "total_nodes": 10,
                "available_nodes": 8,
                "total_cpus": 320,
                "available_cpus": 256,
                "max_time": "7-00:00:00",
            },
            {
                "name": "gpu",
                "state": "up",
                "total_nodes": 4,
                "available_nodes": 3,
                "total_cpus": 128,
                "available_cpus": 96,
                "max_time": "3-00:00:00",
            },
            {
                "name": "debug",
                "state": "up",
                "total_nodes": 2,
                "available_nodes": 2,
                "total_cpus": 64,
                "available_cpus": 64,
                "max_time": "01:00:00",
            },
        ]


