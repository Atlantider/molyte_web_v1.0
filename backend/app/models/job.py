"""
Job models (MD and Postprocess)
"""
from sqlalchemy import Column, Integer, String, Text, Float, DateTime, ForeignKey, Enum, Index, Boolean
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
from sqlalchemy.dialects.postgresql import JSONB
from app.database import Base
import enum
from datetime import datetime, timedelta


class JobStatus(str, enum.Enum):
    """Job status enumeration

    状态流程：
    CREATED -> SUBMITTED -> QUEUED -> RUNNING -> COMPLETED/FAILED

    - CREATED: 任务创建，可修改配置参数
    - SUBMITTED: 用户提交，等待 Worker 拉取
    - QUEUED: Worker 已拉取，Slurm 排队等资源中（对应 Slurm PENDING）
    - RUNNING: Slurm 正在执行（对应 Slurm RUNNING）
    - POSTPROCESSING: 后处理中
    - COMPLETED: 完成
    - FAILED: 失败
    - CANCELLED: 取消
    """
    CREATED = "CREATED"
    SUBMITTED = "SUBMITTED"  # 新增：用户已提交，等待 Worker 拉取
    QUEUED = "QUEUED"
    RUNNING = "RUNNING"
    POSTPROCESSING = "POSTPROCESSING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"


class PostprocessType(str, enum.Enum):
    """Postprocess job type enumeration"""
    RDF = "RDF"
    MSD = "MSD"
    SOLVATION = "SOLVATION"


class RESPJobStatus(str, enum.Enum):
    """RESP charge calculation job status"""
    CREATED = "CREATED"
    QUEUED = "QUEUED"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"


class DataVisibility(str, enum.Enum):
    """Data visibility enumeration - 数据展示状态"""
    PRIVATE = "PRIVATE"          # 私有 - 仅自己和管理员可见
    DELAYED = "DELAYED"          # 延期公开 - 在指定日期后自动公开
    PUBLIC = "PUBLIC"            # 公开 - 所有人可见
    ADMIN_ONLY = "ADMIN_ONLY"    # 仅管理员 - 管理员强制设置（低质量数据等）


class MDJob(Base):
    """Molecular dynamics job model"""
    __tablename__ = "md_jobs"

    id = Column(Integer, primary_key=True, index=True)
    system_id = Column(Integer, ForeignKey("electrolyte_systems.id", ondelete="CASCADE"), nullable=False, index=True)
    user_id = Column(Integer, ForeignKey("users.id", ondelete="CASCADE"), nullable=False, index=True)
    status = Column(Enum(JobStatus), default=JobStatus.CREATED, nullable=False, index=True)
    slurm_job_id = Column(String, index=True)
    progress = Column(Float, default=0.0)
    work_dir = Column(String)
    log_file = Column(String)
    error_message = Column(Text)
    config = Column(JSONB)  # 存储计算参数配置

    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False, index=True)
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)
    started_at = Column(DateTime(timezone=True))
    finished_at = Column(DateTime(timezone=True))

    # 计费相关字段
    cpu_cores = Column(Integer, default=1)  # 使用的CPU核数
    estimated_cpu_hours = Column(Float, default=0.0)  # 预估机时
    actual_cpu_hours = Column(Float, default=0.0)     # 实际消耗机时（MD 计算）
    resp_cpu_hours = Column(Float, default=0.0)       # RESP 电荷计算消耗机时
    result_locked = Column(Boolean, default=False)    # 结果是否锁定（因欠费）
    locked_reason = Column(String(200))               # 锁定原因
    billed = Column(Boolean, default=False)           # 是否已结算
    is_free_quota = Column(Boolean, default=True)     # 是否使用免费核时计算

    # 数据展示控制字段
    visibility = Column(Enum(DataVisibility), default=DataVisibility.DELAYED, nullable=False, index=True)
    visibility_delay_until = Column(DateTime(timezone=True), nullable=True)  # 延期公开日期
    anonymous_public = Column(Boolean, default=False)    # 匿名公开（隐藏用户名和单位）
    allow_download = Column(Boolean, default=True)       # 允许他人下载数据
    visibility_changed_by = Column(Integer, ForeignKey("users.id", ondelete="SET NULL"), nullable=True)  # 修改人
    visibility_changed_at = Column(DateTime(timezone=True), nullable=True)  # 修改时间
    visibility_reason = Column(String(500), nullable=True)  # 修改原因（管理员填写）

    # 贡献奖励追踪
    view_count = Column(Integer, default=0, nullable=False)       # 被查看次数
    download_count = Column(Integer, default=0, nullable=False)   # 被下载次数
    reward_claimed = Column(Boolean, default=False)               # 公开奖励是否已领取

    # 软删除字段
    is_deleted = Column(Boolean, default=False, index=True)  # 是否已删除
    deleted_at = Column(DateTime(timezone=True))  # 删除时间
    deleted_by = Column(Integer, ForeignKey("users.id", ondelete="SET NULL"))  # 删除操作者
    delete_reason = Column(String(500))  # 删除原因

    # Relationships
    system = relationship("ElectrolyteSystem", back_populates="md_jobs")
    user = relationship("User", back_populates="md_jobs", foreign_keys=[user_id])
    visibility_changed_by_user = relationship("User", foreign_keys=[visibility_changed_by])
    deleted_by_user = relationship("User", foreign_keys=[deleted_by])
    postprocess_jobs = relationship("PostprocessJob", back_populates="md_job", cascade="all, delete-orphan")
    resp_jobs = relationship("RESPJob", back_populates="md_job", cascade="all, delete-orphan")
    result_summary = relationship("ResultSummary", back_populates="md_job", uselist=False, cascade="all, delete-orphan")
    rdf_results = relationship("RDFResult", back_populates="md_job", cascade="all, delete-orphan")
    msd_results = relationship("MSDResult", back_populates="md_job", cascade="all, delete-orphan")
    solvation_structures = relationship("SolvationStructure", back_populates="md_job", cascade="all, delete-orphan")
    qc_jobs = relationship("QCJob", back_populates="md_job", cascade="all, delete-orphan")

    # Indexes
    __table_args__ = (
        Index('idx_jobs_user_id', 'user_id'),
        Index('idx_jobs_system_id', 'system_id'),
        Index('idx_jobs_status', 'status'),
        Index('idx_jobs_slurm_job_id', 'slurm_job_id'),
        Index('idx_jobs_created_at', 'created_at'),
        Index('idx_jobs_visibility', 'visibility'),
        Index('idx_jobs_visibility_delay', 'visibility_delay_until'),
    )

    @property
    def is_publicly_visible(self) -> bool:
        """检查数据当前是否公开可见"""
        if self.visibility == DataVisibility.PUBLIC:
            return True
        if self.visibility == DataVisibility.DELAYED:
            if self.visibility_delay_until and datetime.now(self.visibility_delay_until.tzinfo) >= self.visibility_delay_until:
                return True
        return False

    def set_default_delay(self, years: int = 1):
        """设置默认延期公开时间"""
        self.visibility = DataVisibility.DELAYED
        self.visibility_delay_until = datetime.now() + timedelta(days=365 * years)

    def __repr__(self):
        return f"<MDJob(id={self.id}, status={self.status}, visibility={self.visibility})>"


class PostprocessJob(Base):
    """Postprocess job model"""
    __tablename__ = "postprocess_jobs"
    
    id = Column(Integer, primary_key=True, index=True)
    md_job_id = Column(Integer, ForeignKey("md_jobs.id", ondelete="CASCADE"), nullable=False, index=True)
    job_type = Column(Enum(PostprocessType), nullable=False)
    status = Column(Enum(JobStatus), default=JobStatus.CREATED, nullable=False, index=True)
    config = Column(JSONB)
    output_file = Column(String)
    error_message = Column(Text)
    
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)
    started_at = Column(DateTime(timezone=True))
    finished_at = Column(DateTime(timezone=True))
    
    # Relationships
    md_job = relationship("MDJob", back_populates="postprocess_jobs")
    
    # Indexes
    __table_args__ = (
        Index('idx_postprocess_md_job_id', 'md_job_id'),
        Index('idx_postprocess_status', 'status'),
    )

    def __repr__(self):
        return f"<PostprocessJob(id={self.id}, type={self.job_type}, status={self.status})>"


class RESPJob(Base):
    """RESP charge calculation job model"""
    __tablename__ = "resp_jobs"

    id = Column(Integer, primary_key=True, index=True)
    md_job_id = Column(Integer, ForeignKey("md_jobs.id", ondelete="CASCADE"), nullable=False, index=True)
    user_id = Column(Integer, ForeignKey("users.id", ondelete="CASCADE"), nullable=False, index=True)

    # 分子信息
    molecule_name = Column(String(255), nullable=False)
    smiles = Column(Text)

    # 任务状态
    status = Column(Enum(RESPJobStatus), default=RESPJobStatus.CREATED, nullable=False, index=True)
    slurm_job_id = Column(String(50), index=True)

    # 工作目录和文件
    work_dir = Column(Text)
    charge_file = Column(Text)  # 生成的电荷文件路径
    log_file = Column(Text)
    error_message = Column(Text)

    # 核时数统计
    cpu_hours = Column(Float, default=0.0)  # 实际消耗的核时数
    estimated_cpu_hours = Column(Float)  # 预估核时数

    # 配置
    config = Column(JSONB)  # 存储计算参数

    # 时间戳
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)
    started_at = Column(DateTime(timezone=True))
    finished_at = Column(DateTime(timezone=True))

    # Relationships
    md_job = relationship("MDJob", back_populates="resp_jobs", foreign_keys=[md_job_id])
    user = relationship("User", back_populates="resp_jobs", foreign_keys=[user_id])

    # Indexes
    __table_args__ = (
        Index('idx_resp_md_job_id', 'md_job_id'),
        Index('idx_resp_user_id', 'user_id'),
        Index('idx_resp_status', 'status'),
        Index('idx_resp_created_at', 'created_at'),
    )

    def __repr__(self):
        return f"<RESPJob(id={self.id}, molecule={self.molecule_name}, status={self.status})>"

