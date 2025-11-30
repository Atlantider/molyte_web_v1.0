"""
Authentication API endpoints
"""
from fastapi import APIRouter, Depends, HTTPException, status
from fastapi.security import OAuth2PasswordRequestForm
from sqlalchemy.orm import Session
from sqlalchemy import func
from datetime import timedelta
from app.database import get_db
from app.dependencies import get_current_user
from app.models.user import User, UserRole, UserType, USER_TYPE_QUOTAS
from app.models.job import MDJob, JobStatus
from app.schemas.user import UserCreate, User as UserSchema, ChangePassword
from app.schemas.token import Token
from app.core.security import verify_password, get_password_hash, create_access_token
from app.config import settings
from app.core.logger import logger
from app.utils.quota import calculate_user_cpu_hours, get_today_job_count, get_running_job_count

router = APIRouter()


@router.post("/register", response_model=UserSchema, status_code=status.HTTP_201_CREATED)
def register(user_data: UserCreate, db: Session = Depends(get_db)):
    """
    Register a new user

    Args:
        user_data: User registration data including user_type and organization
        db: Database session

    Returns:
        User: Created user with 100 free CPU hours

    Raises:
        HTTPException: If email or username already exists
    """
    # Check if email exists
    if db.query(User).filter(User.email == user_data.email).first():
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="该邮箱已被注册"
        )

    # Check if username exists
    if db.query(User).filter(User.username == user_data.username).first():
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="该用户名已被使用"
        )

    # 验证邮箱类型与用户类型是否匹配
    is_free_email = UserCreate.is_free_email(user_data.email)
    is_academic_email = UserCreate.is_academic_email(user_data.email)

    # 对于学生和研究者，建议使用学术邮箱（但不强制）
    email_warning = None
    if user_data.user_type in [UserType.STUDENT, UserType.RESEARCHER]:
        if is_free_email:
            email_warning = "建议使用单位邮箱注册，可获得更好的服务体验"

    # 获取用户类型对应的默认配额
    quotas = USER_TYPE_QUOTAS.get(user_data.user_type, USER_TYPE_QUOTAS[UserType.STUDENT])

    # Create new user with default quotas based on user_type
    hashed_password = get_password_hash(user_data.password)
    db_user = User(
        email=user_data.email,
        username=user_data.username,
        password_hash=hashed_password,
        role=user_data.role,
        user_type=user_data.user_type,
        organization=user_data.organization,
        department=user_data.department,
        # 根据用户类型设置配额
        daily_job_limit=quotas["daily_job_limit"],
        concurrent_job_limit=quotas["concurrent_job_limit"],
        # 免费赠送100核时
        balance_cpu_hours=quotas["free_cpu_hours"],
        free_cpu_hours_granted=quotas["free_cpu_hours"],
        total_cpu_hours=quotas["free_cpu_hours"],
    )

    db.add(db_user)
    db.commit()
    db.refresh(db_user)

    logger.info(f"New user registered: {db_user.username}, type={db_user.user_type}, org={db_user.organization}")
    return db_user


@router.post("/login", response_model=Token)
def login(
    form_data: OAuth2PasswordRequestForm = Depends(),
    db: Session = Depends(get_db)
):
    """
    Login and get access token

    Args:
        form_data: Login form data (username/email and password)
        db: Database session

    Returns:
        Token: Access token

    Raises:
        HTTPException: If credentials are invalid
    """
    # Try to get user by username first, then by email
    user = db.query(User).filter(User.username == form_data.username).first()

    # If not found by username, try email
    if not user:
        user = db.query(User).filter(User.email == form_data.username).first()

    # Verify user and password
    if not user or not verify_password(form_data.password, user.password_hash):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect username/email or password",
            headers={"WWW-Authenticate": "Bearer"},
        )

    # Create access token
    access_token_expires = timedelta(minutes=settings.ACCESS_TOKEN_EXPIRE_MINUTES)
    access_token = create_access_token(
        data={"sub": user.username},
        expires_delta=access_token_expires
    )

    logger.info(f"User logged in: {user.username}")
    return {"access_token": access_token, "token_type": "bearer"}


@router.get("/me", response_model=UserSchema)
def get_current_user_info(current_user: User = Depends(get_current_user)):
    """
    Get current user information

    Args:
        current_user: Current authenticated user

    Returns:
        User: Current user data
    """
    return current_user


@router.post("/change-password", status_code=status.HTTP_200_OK)
def change_password(
    password_data: ChangePassword,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Change user password

    Args:
        password_data: Old and new password
        current_user: Current authenticated user
        db: Database session

    Returns:
        dict: Success message

    Raises:
        HTTPException: If old password is incorrect
    """
    # Verify old password
    if not verify_password(password_data.old_password, current_user.password_hash):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Incorrect old password"
        )

    # Update password
    current_user.password_hash = get_password_hash(password_data.new_password)
    db.commit()

    logger.info(f"Password changed for user: {current_user.username}")
    return {"message": "Password changed successfully"}


@router.get("/me/profile")
def get_user_profile(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    获取当前用户完整资料（包含配额使用情况和贡献统计）

    Returns:
        dict: 用户完整资料
    """
    # 计算任务统计
    total_jobs = db.query(func.count(MDJob.id)).filter(MDJob.user_id == current_user.id).scalar() or 0
    completed_jobs = db.query(func.count(MDJob.id)).filter(
        MDJob.user_id == current_user.id,
        MDJob.status == JobStatus.COMPLETED
    ).scalar() or 0
    running_jobs = get_running_job_count(current_user.id, db)
    failed_jobs = db.query(func.count(MDJob.id)).filter(
        MDJob.user_id == current_user.id,
        MDJob.status == JobStatus.FAILED
    ).scalar() or 0

    # 计算配额使用情况
    used_cpu_hours = calculate_user_cpu_hours(current_user.id, db)
    today_jobs = get_today_job_count(current_user.id, db)

    # 获取余额系统数据
    balance = current_user.balance_cpu_hours or 0.0
    frozen = current_user.frozen_cpu_hours or 0.0
    debt = current_user.debt_cpu_hours or 0.0
    free_granted = current_user.free_cpu_hours_granted or 100.0

    # 计算总充值核时（从交易记录统计）
    from app.models.billing import QuotaTransaction
    total_recharged = db.query(func.sum(QuotaTransaction.amount)).filter(
        QuotaTransaction.user_id == current_user.id,
        QuotaTransaction.type.in_(['recharge', 'points_exchange', 'admin_adjust', 'bonus'])
    ).scalar() or 0.0

    # 计算总消耗核时
    total_consumed = db.query(func.sum(QuotaTransaction.amount)).filter(
        QuotaTransaction.user_id == current_user.id,
        QuotaTransaction.type == 'consume'
    ).scalar() or 0.0
    total_consumed = abs(total_consumed)  # consume是负数

    # 可用余额 = 账户余额 - 冻结
    available_balance = max(0, balance - frozen)

    # 总核时 = 初始赠送 + 总充值
    total_cpu_hours = free_granted + total_recharged

    return {
        # 基本信息
        "id": current_user.id,
        "username": current_user.username,
        "email": current_user.email,
        "role": current_user.role.value if current_user.role else "USER",
        "user_type": current_user.user_type.value if current_user.user_type else "STUDENT",
        "organization": current_user.organization,
        "department": current_user.department,
        "email_verified": getattr(current_user, 'email_verified', False),
        "is_active": current_user.is_active,
        "created_at": current_user.created_at.isoformat() if current_user.created_at else None,
        "last_login_at": current_user.last_login_at.isoformat() if current_user.last_login_at else None,

        # 配额配置
        "daily_job_limit": current_user.daily_job_limit,
        "concurrent_job_limit": current_user.concurrent_job_limit,
        "storage_quota_gb": current_user.storage_quota_gb,
        "allowed_partitions": current_user.allowed_partitions,

        # 核时余额系统
        "balance_cpu_hours": round(balance, 2),           # 账户余额
        "frozen_cpu_hours": round(frozen, 2),             # 冻结核时（运行中的任务）
        "available_cpu_hours": round(available_balance, 2), # 可用余额
        "debt_cpu_hours": round(debt, 2),                 # 欠费核时

        # 核时统计
        "free_cpu_hours_granted": round(free_granted, 2),   # 初始赠送
        "total_recharged": round(total_recharged, 2),       # 总充值
        "total_consumed": round(total_consumed, 2),         # 总消耗
        "total_cpu_hours": round(total_cpu_hours, 2),       # 总核时额度

        # 兼容旧字段（给前端使用）
        "used_cpu_hours": round(total_consumed, 2),
        "remaining_cpu_hours": round(available_balance, 2),

        # 贡献统计
        "public_data_count": current_user.public_data_count or 0,
        "contribution_points": round(current_user.contribution_points or 0.0, 1),
        "private_quota_used": current_user.private_quota_used or 0,
        "private_quota_limit": getattr(current_user, 'private_quota_limit', 0),

        # 使用情况
        "today_jobs": today_jobs,
        "running_jobs": running_jobs,

        # 任务统计
        "total_jobs": total_jobs,
        "completed_jobs": completed_jobs,
        "failed_jobs": failed_jobs,

        # 计算使用率（基于余额系统）
        "cpu_hours_usage_percent": round(min(100, total_consumed / total_cpu_hours * 100), 1) if total_cpu_hours > 0 else 0,
        "daily_jobs_usage_percent": round(min(100, today_jobs / current_user.daily_job_limit * 100), 1) if current_user.daily_job_limit > 0 else 0,
    }

