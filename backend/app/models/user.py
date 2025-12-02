"""
User model
"""
from sqlalchemy import Column, Integer, String, DateTime, Enum, Boolean, Float, JSON, Text
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
from app.database import Base
import enum


class UserRole(str, enum.Enum):
    """User role enumeration"""
    ADMIN = "ADMIN"
    PREMIUM = "PREMIUM"
    USER = "USER"
    GUEST = "GUEST"


class UserType(str, enum.Enum):
    """User type enumeration - 用户类型"""
    STUDENT = "STUDENT"        # 学生
    RESEARCHER = "RESEARCHER"  # 研究者
    COMPANY = "COMPANY"        # 企业用户


# 不同用户类型的默认配额配置
USER_TYPE_QUOTAS = {
    UserType.STUDENT: {
        "free_cpu_hours": 100.0,      # 免费核时
        "price_per_hour": 0.05,       # 单价
        "daily_job_limit": 5,         # 日任务限制
        "concurrent_job_limit": 2,    # 并发限制
        "private_quota": 0,           # 私有配额（0表示不能永久私有）
        "max_delay_years": 1,         # 最长延期年数
    },
    UserType.RESEARCHER: {
        "free_cpu_hours": 100.0,
        "price_per_hour": 0.08,
        "daily_job_limit": 10,
        "concurrent_job_limit": 3,
        "private_quota": 50,
        "max_delay_years": 2,
    },
    UserType.COMPANY: {
        "free_cpu_hours": 100.0,
        "price_per_hour": 0.15,
        "daily_job_limit": 20,
        "concurrent_job_limit": 5,
        "private_quota": 200,
        "max_delay_years": 3,
    },
}


class User(Base):
    """User model"""
    __tablename__ = "users"

    id = Column(Integer, primary_key=True, index=True)
    email = Column(String, unique=True, index=True, nullable=False)
    username = Column(String, unique=True, index=True, nullable=False)
    password_hash = Column(String, nullable=False)
    role = Column(Enum(UserRole), default=UserRole.USER, nullable=False)

    # 用户类型和单位信息
    user_type = Column(Enum(UserType), default=UserType.STUDENT, nullable=False)  # 用户类型
    organization = Column(String(200), nullable=True)  # 单位名称
    department = Column(String(100), nullable=True)    # 部门（可选）

    # 邮箱验证
    email_verified = Column(Boolean, default=False, nullable=False)  # 邮箱是否验证
    verification_token = Column(String(100), nullable=True)          # 验证令牌
    verification_expires = Column(DateTime(timezone=True), nullable=True)  # 令牌过期时间

    # Status
    is_active = Column(Boolean, default=True, nullable=False)
    last_login_at = Column(DateTime(timezone=True), nullable=True)

    # Resource Quotas (旧版，保留兼容)
    total_cpu_hours = Column(Float, default=100.0, nullable=False)  # 总机时配额
    daily_job_limit = Column(Integer, default=10, nullable=False)   # 每日任务限制
    concurrent_job_limit = Column(Integer, default=3, nullable=False)  # 并发任务限制
    storage_quota_gb = Column(Float, default=10.0, nullable=False)  # 存储配额（GB）

    # 余额系统（新版计费）
    balance_cpu_hours = Column(Float, default=100.0, nullable=False)  # 可用余额（机时）- 默认100核时
    frozen_cpu_hours = Column(Float, default=0.0, nullable=False)    # 冻结机时（运行中任务）
    debt_cpu_hours = Column(Float, default=0.0, nullable=False)      # 欠费机时
    free_cpu_hours_granted = Column(Float, default=100.0, nullable=False)  # 已赠送的免费核时

    # 公开贡献统计
    public_data_count = Column(Integer, default=0, nullable=False)    # 公开数据数量
    contribution_points = Column(Float, default=0.0, nullable=False)  # 贡献积分
    private_quota_used = Column(Integer, default=0, nullable=False)   # 已用私有配额

    # Queue/Partition permissions (JSON array of allowed partition names)
    # Example: ["cpu", "gpu"] or null for admin (all partitions)
    allowed_partitions = Column(JSON, nullable=True)

    # Timestamps
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)

    # Relationships
    projects = relationship("Project", back_populates="owner", cascade="all, delete-orphan")
    md_jobs = relationship("MDJob", back_populates="user", cascade="all, delete-orphan",
                          foreign_keys="MDJob.user_id")
    qc_jobs = relationship("QCJob", back_populates="user", cascade="all, delete-orphan",
                          foreign_keys="QCJob.user_id")
    resp_jobs = relationship("RESPJob", back_populates="user", cascade="all, delete-orphan",
                            foreign_keys="RESPJob.user_id")
    usage_stats = relationship("UserUsageStats", back_populates="user", cascade="all, delete-orphan")
    recharge_orders = relationship("RechargeOrder", back_populates="user", cascade="all, delete-orphan")
    quota_transactions = relationship("QuotaTransaction", back_populates="user", cascade="all, delete-orphan")

    # User preferences
    solvent_combinations = relationship("UserSolventCombination", back_populates="user", cascade="all, delete-orphan")
    ion_combinations = relationship("UserIonCombination", back_populates="user", cascade="all, delete-orphan")

    @property
    def private_quota_limit(self) -> int:
        """获取用户的私有配额限制"""
        return USER_TYPE_QUOTAS.get(self.user_type, {}).get("private_quota", 0)

    @property
    def max_delay_years(self) -> int:
        """获取用户的最长延期年数"""
        return USER_TYPE_QUOTAS.get(self.user_type, {}).get("max_delay_years", 1)

    @property
    def can_set_private(self) -> bool:
        """检查用户是否可以设置永久私有"""
        if self.role == UserRole.ADMIN:
            return True
        return self.private_quota_used < self.private_quota_limit

    def __repr__(self):
        return f"<User(id={self.id}, username={self.username}, email={self.email}, role={self.role}, type={self.user_type})>"

