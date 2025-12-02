"""
User schemas
"""
from pydantic import BaseModel, Field, field_validator
from typing import Optional, List
from datetime import datetime
from app.models.user import UserRole, UserType
import re


# 免费邮箱域名列表（需要人工审核）
FREE_EMAIL_DOMAINS = [
    'gmail.com', 'qq.com', '163.com', '126.com', 'hotmail.com',
    'outlook.com', 'yahoo.com', 'sina.com', 'sohu.com', 'foxmail.com',
    'icloud.com', 'mail.com', 'protonmail.com', 'yandex.com'
]

# 学术邮箱域名后缀
ACADEMIC_EMAIL_SUFFIXES = [
    '.edu', '.edu.cn', '.ac.cn', '.cas.cn', '.org.cn'
]


class UserBase(BaseModel):
    """Base user schema"""
    email: str = Field(..., pattern=r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$')
    username: str = Field(..., min_length=3, max_length=50)

    @field_validator('email')
    @classmethod
    def validate_email(cls, v: str) -> str:
        """Validate email format"""
        if not re.match(r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$', v):
            raise ValueError('Invalid email format')
        return v.lower()


class UserCreate(UserBase):
    """Schema for creating a user"""
    password: str = Field(..., min_length=6, max_length=100)
    role: UserRole = Field(default=UserRole.USER, description="用户角色，默认为普通用户")
    # 新增字段
    user_type: UserType = Field(default=UserType.STUDENT, description="用户类型：学生/研究者/企业")
    organization: str = Field(..., min_length=2, max_length=200, description="单位名称")
    department: Optional[str] = Field(None, max_length=100, description="部门（可选）")

    @field_validator('organization')
    @classmethod
    def validate_organization(cls, v: str) -> str:
        """Validate organization name"""
        if not v or len(v.strip()) < 2:
            raise ValueError('单位名称至少2个字符')
        return v.strip()

    @classmethod
    def is_free_email(cls, email: str) -> bool:
        """检查是否为免费邮箱"""
        domain = email.split('@')[-1].lower()
        return domain in FREE_EMAIL_DOMAINS

    @classmethod
    def is_academic_email(cls, email: str) -> bool:
        """检查是否为学术邮箱"""
        domain = email.split('@')[-1].lower()
        return any(domain.endswith(suffix) for suffix in ACADEMIC_EMAIL_SUFFIXES)


class UserUpdate(BaseModel):
    """Schema for updating a user"""
    email: Optional[str] = Field(None, pattern=r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$')
    username: Optional[str] = Field(None, min_length=3, max_length=50)
    password: Optional[str] = Field(None, min_length=6, max_length=100)
    role: Optional[UserRole] = None
    user_type: Optional[UserType] = None
    organization: Optional[str] = Field(None, max_length=200)
    department: Optional[str] = Field(None, max_length=100)

    @field_validator('email')
    @classmethod
    def validate_email(cls, v: Optional[str]) -> Optional[str]:
        """Validate email format"""
        if v is not None:
            if not re.match(r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$', v):
                raise ValueError('Invalid email format')
            return v.lower()
        return v


class UserInDB(UserBase):
    """User schema with database fields"""
    id: int
    role: UserRole
    user_type: UserType
    organization: Optional[str] = None
    department: Optional[str] = None
    email_verified: bool = False
    is_active: bool
    last_login_at: Optional[datetime] = None
    total_cpu_hours: float
    daily_job_limit: int
    concurrent_job_limit: int
    storage_quota_gb: float
    balance_cpu_hours: float = 100.0
    frozen_cpu_hours: float = 0.0
    public_data_count: int = 0
    contribution_points: float = 0.0
    allowed_partitions: Optional[List[str]] = None
    created_at: datetime
    updated_at: datetime

    class Config:
        from_attributes = True


class User(UserInDB):
    """User response schema"""
    pass


class UserLogin(BaseModel):
    """Schema for user login"""
    username: str
    password: str


class ChangePassword(BaseModel):
    """Schema for changing password"""
    old_password: str = Field(..., min_length=6, max_length=100)
    new_password: str = Field(..., min_length=6, max_length=100)


class UserProfile(BaseModel):
    """用户个人资料响应"""
    id: int
    username: str
    email: str
    role: UserRole
    user_type: UserType
    organization: Optional[str] = None
    department: Optional[str] = None
    email_verified: bool
    balance_cpu_hours: float
    frozen_cpu_hours: float
    public_data_count: int
    contribution_points: float
    created_at: datetime

    class Config:
        from_attributes = True

