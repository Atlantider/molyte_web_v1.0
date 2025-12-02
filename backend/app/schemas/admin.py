"""
Admin API schemas
"""
from typing import Optional, List, Dict, Any
from datetime import datetime, date
from pydantic import BaseModel, EmailStr, Field
from app.models.user import UserRole, UserType


# ============ User Management Schemas ============

class UserListItem(BaseModel):
    """User list item for admin"""
    id: int
    username: str
    email: EmailStr
    role: UserRole
    user_type: UserType
    organization: Optional[str] = None
    department: Optional[str] = None
    is_active: bool
    total_cpu_hours: float
    balance_cpu_hours: float = 0.0
    daily_job_limit: int
    concurrent_job_limit: int
    storage_quota_gb: float
    allowed_partitions: Optional[List[str]] = None
    last_login_at: Optional[datetime] = None
    created_at: datetime

    class Config:
        from_attributes = True


class UserDetail(BaseModel):
    """Detailed user information for admin"""
    id: int
    username: str
    email: EmailStr
    role: UserRole
    user_type: UserType
    organization: Optional[str] = None
    department: Optional[str] = None
    is_active: bool
    total_cpu_hours: float
    balance_cpu_hours: float = 0.0
    daily_job_limit: int
    concurrent_job_limit: int
    storage_quota_gb: float
    allowed_partitions: Optional[List[str]] = None
    last_login_at: Optional[datetime] = None
    created_at: datetime
    updated_at: datetime

    # Usage statistics
    used_cpu_hours: float = 0.0
    today_jobs: int = 0
    running_jobs: int = 0
    total_jobs: int = 0
    completed_jobs: int = 0
    failed_jobs: int = 0

    class Config:
        from_attributes = True


class UserUpdate(BaseModel):
    """Update user information"""
    email: Optional[EmailStr] = None
    role: Optional[UserRole] = None
    is_active: Optional[bool] = None
    total_cpu_hours: Optional[float] = Field(None, ge=0)
    daily_job_limit: Optional[int] = Field(None, ge=0)
    concurrent_job_limit: Optional[int] = Field(None, ge=0)
    storage_quota_gb: Optional[float] = Field(None, ge=0)
    allowed_partitions: Optional[List[str]] = None  # None means all partitions (admin only)


class UserCreate(BaseModel):
    """Create new user"""
    username: str = Field(..., min_length=3, max_length=50)
    email: EmailStr
    password: str = Field(..., min_length=6)
    role: UserRole = Field(default=UserRole.USER, description="用户角色，默认为普通用户")
    total_cpu_hours: float = Field(default=100.0, ge=0, description="总CPU核时配额")
    daily_job_limit: int = Field(default=10, ge=0, description="每日任务限制")
    concurrent_job_limit: int = Field(default=3, ge=0, description="并发任务限制")
    storage_quota_gb: float = Field(default=10.0, ge=0, description="存储配额(GB)")
    allowed_partitions: Optional[List[str]] = Field(default_factory=lambda: ["cpu"], description="允许使用的队列，默认仅cpu队列")


# ============ Statistics Schemas ============

class GlobalStats(BaseModel):
    """Global system statistics"""
    total_users: int
    active_users: int
    total_jobs: int
    running_jobs: int
    queued_jobs: int
    completed_jobs: int
    failed_jobs: int
    total_cpu_hours_used: float
    total_cpu_hours_allocated: float
    total_storage_used_gb: float
    total_storage_allocated_gb: float


class UserUsageStatsItem(BaseModel):
    """User usage statistics item"""
    user_id: int
    username: str
    email: str
    role: UserRole
    used_cpu_hours: float
    total_cpu_hours: float
    usage_percentage: float
    total_jobs: int
    running_jobs: int
    completed_jobs: int
    failed_jobs: int
    last_job_at: Optional[datetime] = None


class UserRanking(BaseModel):
    """User ranking by resource usage"""
    user_id: int
    username: str
    email: str
    metric_value: float
    rank: int


class TrendDataPoint(BaseModel):
    """Trend data point"""
    date: date
    value: float
    label: Optional[str] = None


class StatisticsSummary(BaseModel):
    """Statistics summary"""
    period: str  # "today", "7days", "30days"
    jobs_submitted: int
    jobs_completed: int
    jobs_failed: int
    cpu_hours_used: float
    avg_job_duration_hours: float
    peak_concurrent_jobs: int


# ============ Audit Log Schemas ============

class AuditLogItem(BaseModel):
    """Audit log item"""
    id: int
    user_id: Optional[int] = None
    username: Optional[str] = None
    action: str
    resource_type: Optional[str] = None
    resource_id: Optional[int] = None
    details: Optional[Dict[str, Any]] = None
    ip_address: Optional[str] = None
    created_at: datetime
    
    class Config:
        from_attributes = True


# ============ Response Schemas ============

class QuotaCheckResponse(BaseModel):
    """Quota check response"""
    allowed: bool
    reason: Optional[str] = None
    details: Dict[str, Any]

