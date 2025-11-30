"""
数据可见性控制 API
"""
from fastapi import APIRouter, Depends, HTTPException, status, Query
from sqlalchemy.orm import Session
from sqlalchemy import func, or_, and_
from datetime import datetime, timedelta
from typing import Optional, List
from pydantic import BaseModel, Field

from app.database import get_db
from app.dependencies import get_current_user, get_current_admin
from app.models.user import User, UserRole, UserType, USER_TYPE_QUOTAS
from app.models.job import MDJob, JobStatus, DataVisibility
from app.core.logger import logger

router = APIRouter()


# ============ Schemas ============

class VisibilityUpdate(BaseModel):
    """更新可见性请求"""
    visibility: DataVisibility
    delay_days: Optional[int] = Field(None, ge=0, le=1095, description="延期天数（最多3年）")
    anonymous_public: Optional[bool] = Field(False, description="匿名公开")
    allow_download: Optional[bool] = Field(True, description="允许下载")
    reason: Optional[str] = Field(None, max_length=500, description="修改原因")


class VisibilityBatchUpdate(BaseModel):
    """批量更新可见性请求"""
    job_ids: List[int]
    visibility: DataVisibility
    delay_days: Optional[int] = None
    reason: Optional[str] = None


class VisibilityStats(BaseModel):
    """可见性统计"""
    total: int
    public: int
    private: int
    delayed: int
    admin_only: int


# ============ 用户端 API ============

@router.get("/my-jobs")
def get_my_jobs_visibility(
    visibility: Optional[DataVisibility] = None,
    page: int = Query(1, ge=1),
    page_size: int = Query(20, ge=1, le=100),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """获取当前用户的任务可见性列表"""
    query = db.query(MDJob).filter(
        MDJob.user_id == current_user.id,
        MDJob.status == JobStatus.COMPLETED
    )
    
    if visibility:
        query = query.filter(MDJob.visibility == visibility)
    
    total = query.count()
    jobs = query.order_by(MDJob.created_at.desc()).offset((page - 1) * page_size).limit(page_size).all()
    
    return {
        "total": total,
        "page": page,
        "page_size": page_size,
        "items": [
            {
                "id": job.id,
                "name": job.config.get("job_name", f"任务 #{job.id}") if job.config else f"任务 #{job.id}",
                "visibility": job.visibility.value if job.visibility else "DELAYED",
                "visibility_delay_until": job.visibility_delay_until.isoformat() if job.visibility_delay_until else None,
                "anonymous_public": job.anonymous_public,
                "allow_download": job.allow_download,
                "view_count": job.view_count or 0,
                "download_count": job.download_count or 0,
                "reward_claimed": job.reward_claimed,
                "is_free_quota": job.is_free_quota,
                "created_at": job.created_at.isoformat() if job.created_at else None,
            }
            for job in jobs
        ]
    }


@router.put("/job/{job_id}")
def update_job_visibility(
    job_id: int,
    update: VisibilityUpdate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """更新任务可见性"""
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="任务不存在")
    
    # 检查权限
    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="无权修改此任务")
    
    # 检查是否使用免费核时（免费核时不能设为永久私有）
    if job.is_free_quota and update.visibility == DataVisibility.PRIVATE:
        raise HTTPException(
            status_code=400, 
            detail="使用免费核时的任务不能设为永久私有，请选择延期公开或立即公开"
        )
    
    # 检查私有配额
    user_quotas = USER_TYPE_QUOTAS.get(current_user.user_type, USER_TYPE_QUOTAS[UserType.STUDENT])
    if update.visibility == DataVisibility.PRIVATE:
        private_count = db.query(func.count(MDJob.id)).filter(
            MDJob.user_id == current_user.id,
            MDJob.visibility == DataVisibility.PRIVATE
        ).scalar() or 0
        if private_count >= user_quotas["private_quota"]:
            raise HTTPException(
                status_code=400,
                detail=f"私有配额已用完（{private_count}/{user_quotas['private_quota']}），请升级账户或选择其他可见性"
            )
    
    # 检查是否已领取奖励且试图改为非公开状态
    old_visibility = job.visibility
    if job.reward_claimed and update.visibility != DataVisibility.PUBLIC:
        # 已领取奖励，改为非公开需要扣除奖励
        reward_hours = 10.0
        if current_user.balance_cpu_hours >= reward_hours:
            current_user.balance_cpu_hours = (current_user.balance_cpu_hours or 0) - reward_hours
            current_user.contribution_points = max(0, (current_user.contribution_points or 0) - 10)
            current_user.public_data_count = max(0, (current_user.public_data_count or 0) - 1)
            job.reward_claimed = False
            logger.info(f"User {current_user.id} revoked public status for job {job_id}, deducted {reward_hours} CPU hours")
        else:
            raise HTTPException(
                status_code=400,
                detail=f"余额不足！取消公开需要扣除 {reward_hours} 核时奖励，当前余额：{current_user.balance_cpu_hours:.2f} 核时"
            )

    # 更新可见性
    job.visibility = update.visibility
    job.anonymous_public = update.anonymous_public
    job.allow_download = update.allow_download
    job.visibility_changed_by = current_user.id
    job.visibility_changed_at = datetime.utcnow()
    job.visibility_reason = update.reason

    # 设置延期时间
    if update.visibility == DataVisibility.DELAYED and update.delay_days:
        max_delay_years = user_quotas["max_delay_years"]
        max_delay_days = max_delay_years * 365
        actual_delay = min(update.delay_days, max_delay_days)
        job.visibility_delay_until = datetime.utcnow() + timedelta(days=actual_delay)
    elif update.visibility == DataVisibility.PUBLIC:
        job.visibility_delay_until = None

    db.commit()

    logger.info(f"Job {job_id} visibility changed: {old_visibility} -> {update.visibility} by user {current_user.id}")

    return {"message": "可见性已更新", "visibility": job.visibility.value}


@router.post("/job/{job_id}/claim-reward")
def claim_public_reward(
    job_id: int,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """领取公开数据奖励"""
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="任务不存在")

    if job.user_id != current_user.id:
        raise HTTPException(status_code=403, detail="无权操作此任务")

    if job.visibility != DataVisibility.PUBLIC:
        raise HTTPException(status_code=400, detail="只有公开数据才能领取奖励")

    if job.reward_claimed:
        raise HTTPException(status_code=400, detail="奖励已领取")

    # 发放奖励：+10 核时
    reward_hours = 10.0
    job.reward_claimed = True
    current_user.balance_cpu_hours = (current_user.balance_cpu_hours or 0) + reward_hours
    current_user.contribution_points = (current_user.contribution_points or 0) + 10
    current_user.public_data_count = (current_user.public_data_count or 0) + 1

    db.commit()

    logger.info(f"User {current_user.id} claimed reward for job {job_id}: +{reward_hours} CPU hours")

    return {
        "message": f"奖励已领取：+{reward_hours} 核时",
        "reward_hours": reward_hours,
        "new_balance": current_user.balance_cpu_hours
    }


@router.get("/my-stats")
def get_my_visibility_stats(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """获取当前用户的可见性统计"""
    base_query = db.query(MDJob).filter(
        MDJob.user_id == current_user.id,
        MDJob.status == JobStatus.COMPLETED
    )

    total = base_query.count()
    public = base_query.filter(MDJob.visibility == DataVisibility.PUBLIC).count()
    private = base_query.filter(MDJob.visibility == DataVisibility.PRIVATE).count()
    delayed = base_query.filter(MDJob.visibility == DataVisibility.DELAYED).count()

    # 获取用户配额信息
    user_quotas = USER_TYPE_QUOTAS.get(current_user.user_type, USER_TYPE_QUOTAS[UserType.STUDENT])

    return {
        "total": total,
        "public": public,
        "private": private,
        "delayed": delayed,
        "private_quota_used": private,
        "private_quota_limit": user_quotas["private_quota"],
        "max_delay_years": user_quotas["max_delay_years"],
        "contribution_points": current_user.contribution_points or 0,
        "public_data_count": current_user.public_data_count or 0,
        "balance_cpu_hours": current_user.balance_cpu_hours or 0,
    }


# 积分兑换比例: 10积分 = 1核时
POINTS_TO_CPU_HOURS_RATIO = 10.0


class PointsExchangeRequest(BaseModel):
    """积分兑换请求"""
    points: float = Field(..., gt=0, description="要兑换的积分数量")


@router.post("/exchange-points")
def exchange_points_for_cpu_hours(
    request: PointsExchangeRequest,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    积分兑换核时
    兑换比例: 10积分 = 1核时
    """
    current_points = current_user.contribution_points or 0

    if request.points > current_points:
        raise HTTPException(
            status_code=400,
            detail=f"积分不足！当前积分：{current_points:.1f}，需要：{request.points:.1f}"
        )

    # 计算可兑换的核时
    cpu_hours = request.points / POINTS_TO_CPU_HOURS_RATIO

    # 扣除积分，增加核时
    current_user.contribution_points = current_points - request.points
    current_user.balance_cpu_hours = (current_user.balance_cpu_hours or 0) + cpu_hours

    db.commit()

    logger.info(f"User {current_user.id} exchanged {request.points} points for {cpu_hours} CPU hours")

    return {
        "message": f"兑换成功！{request.points:.1f} 积分 → {cpu_hours:.2f} 核时",
        "points_used": request.points,
        "cpu_hours_gained": cpu_hours,
        "remaining_points": current_user.contribution_points,
        "new_balance": current_user.balance_cpu_hours,
    }


@router.get("/exchange-rate")
def get_exchange_rate(current_user: User = Depends(get_current_user)):
    """获取积分兑换比例信息"""
    current_points = current_user.contribution_points or 0
    max_cpu_hours = current_points / POINTS_TO_CPU_HOURS_RATIO

    return {
        "ratio": POINTS_TO_CPU_HOURS_RATIO,
        "description": f"{int(POINTS_TO_CPU_HOURS_RATIO)} 积分 = 1 核时",
        "current_points": current_points,
        "max_exchangeable_cpu_hours": max_cpu_hours,
        "current_balance": current_user.balance_cpu_hours or 0,
    }


# ============ 管理员 API ============

@router.get("/admin/all-jobs")
def admin_get_all_jobs_visibility(
    visibility: Optional[DataVisibility] = None,
    user_id: Optional[int] = None,
    search: Optional[str] = None,
    page: int = Query(1, ge=1),
    page_size: int = Query(20, ge=1, le=100),
    current_user: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """管理员获取所有任务可见性列表"""
    query = db.query(MDJob).filter(MDJob.status == JobStatus.COMPLETED)

    if visibility:
        query = query.filter(MDJob.visibility == visibility)
    if user_id:
        query = query.filter(MDJob.user_id == user_id)

    total = query.count()
    jobs = query.order_by(MDJob.created_at.desc()).offset((page - 1) * page_size).limit(page_size).all()

    # 获取用户信息
    user_ids = list(set(job.user_id for job in jobs))
    users = {u.id: u for u in db.query(User).filter(User.id.in_(user_ids)).all()}

    return {
        "total": total,
        "page": page,
        "page_size": page_size,
        "items": [
            {
                "id": job.id,
                "name": job.config.get("job_name", f"任务 #{job.id}") if job.config else f"任务 #{job.id}",
                "user_id": job.user_id,
                "username": users.get(job.user_id, {}).username if users.get(job.user_id) else "Unknown",
                "visibility": job.visibility.value if job.visibility else "DELAYED",
                "visibility_delay_until": job.visibility_delay_until.isoformat() if job.visibility_delay_until else None,
                "anonymous_public": job.anonymous_public,
                "allow_download": job.allow_download,
                "view_count": job.view_count or 0,
                "download_count": job.download_count or 0,
                "is_free_quota": job.is_free_quota,
                "created_at": job.created_at.isoformat() if job.created_at else None,
            }
            for job in jobs
        ]
    }


@router.put("/admin/job/{job_id}")
def admin_update_job_visibility(
    job_id: int,
    update: VisibilityUpdate,
    current_user: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """管理员更新任务可见性"""
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="任务不存在")

    old_visibility = job.visibility
    job.visibility = update.visibility
    job.anonymous_public = update.anonymous_public
    job.allow_download = update.allow_download
    job.visibility_changed_by = current_user.id
    job.visibility_changed_at = datetime.utcnow()
    job.visibility_reason = update.reason

    if update.visibility == DataVisibility.DELAYED and update.delay_days:
        job.visibility_delay_until = datetime.utcnow() + timedelta(days=update.delay_days)
    elif update.visibility == DataVisibility.PUBLIC:
        job.visibility_delay_until = None

    db.commit()

    logger.info(f"Admin {current_user.id} changed job {job_id} visibility: {old_visibility} -> {update.visibility}")

    return {"message": "可见性已更新", "visibility": job.visibility.value}


@router.post("/admin/batch-update")
def admin_batch_update_visibility(
    update: VisibilityBatchUpdate,
    current_user: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """管理员批量更新可见性"""
    jobs = db.query(MDJob).filter(MDJob.id.in_(update.job_ids)).all()

    updated_count = 0
    for job in jobs:
        job.visibility = update.visibility
        job.visibility_changed_by = current_user.id
        job.visibility_changed_at = datetime.utcnow()
        job.visibility_reason = update.reason

        if update.visibility == DataVisibility.DELAYED and update.delay_days:
            job.visibility_delay_until = datetime.utcnow() + timedelta(days=update.delay_days)
        elif update.visibility == DataVisibility.PUBLIC:
            job.visibility_delay_until = None

        updated_count += 1

    db.commit()

    logger.info(f"Admin {current_user.id} batch updated {updated_count} jobs visibility to {update.visibility}")

    return {"message": f"已更新 {updated_count} 个任务", "updated_count": updated_count}


@router.get("/admin/stats")
def admin_get_visibility_stats(
    current_user: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """管理员获取全局可见性统计"""
    base_query = db.query(MDJob).filter(MDJob.status == JobStatus.COMPLETED)

    total = base_query.count()
    public = base_query.filter(MDJob.visibility == DataVisibility.PUBLIC).count()
    private = base_query.filter(MDJob.visibility == DataVisibility.PRIVATE).count()
    delayed = base_query.filter(MDJob.visibility == DataVisibility.DELAYED).count()
    admin_only = base_query.filter(MDJob.visibility == DataVisibility.ADMIN_ONLY).count()

    # 即将到期的延期数据（30天内）
    soon_public = base_query.filter(
        MDJob.visibility == DataVisibility.DELAYED,
        MDJob.visibility_delay_until <= datetime.utcnow() + timedelta(days=30)
    ).count()

    return {
        "total": total,
        "public": public,
        "private": private,
        "delayed": delayed,
        "admin_only": admin_only,
        "soon_public": soon_public,
    }

