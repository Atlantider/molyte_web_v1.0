"""
Quantum Chemistry (QC) models
量子化学计算相关的数据库模型
"""
from enum import Enum as PyEnum
from sqlalchemy import Column, Integer, String, Float, DateTime, ForeignKey, Text, Index, Enum, Boolean
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
from sqlalchemy.dialects.postgresql import JSONB

from app.database import Base


class QCJobStatus(str, PyEnum):
    """QC任务状态枚举

    状态流程：
    CREATED -> SUBMITTED -> QUEUED -> RUNNING -> COMPLETED/FAILED

    - CREATED: 任务创建，可修改配置参数
    - SUBMITTED: 用户提交，等待 Worker 拉取
    - QUEUED: Worker 已拉取，Slurm 排队等资源中
    - RUNNING: Slurm 正在执行
    - POSTPROCESSING: 后处理中
    - COMPLETED: 完成
    - FAILED: 失败
    - CANCELLED: 取消
    """
    CREATED = "CREATED"
    SUBMITTED = "SUBMITTED"
    QUEUED = "QUEUED"
    RUNNING = "RUNNING"
    POSTPROCESSING = "POSTPROCESSING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"


class MoleculeType(str, PyEnum):
    """分子类型枚举"""
    SOLVENT = "solvent"
    CATION = "cation"
    ANION = "anion"
    CUSTOM = "custom"


class QCJob(Base):
    """量子化学任务模型"""
    __tablename__ = "qc_jobs"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id", ondelete="CASCADE"), nullable=False, index=True)
    
    # 关联MD任务（可选）
    md_job_id = Column(Integer, ForeignKey("md_jobs.id", ondelete="SET NULL"), nullable=True, index=True)

    # 关联去溶剂化后处理任务（可选）
    desolvation_postprocess_job_id = Column(Integer, ForeignKey("postprocess_jobs.id", ondelete="SET NULL"), nullable=True, index=True)
    
    # 分子信息
    molecule_name = Column(String(255), nullable=False)
    smiles = Column(Text, nullable=False, index=True)
    molecule_type = Column(String(20), default="custom")
    
    # 计算参数
    basis_set = Column(String(50), nullable=False, default="6-31++g(d,p)")
    functional = Column(String(50), default="B3LYP")
    charge = Column(Integer, default=0)
    spin_multiplicity = Column(Integer, default=1)
    solvent_model = Column(String(20), default="gas")  # gas, pcm, smd
    solvent_name = Column(String(50))  # 隐式溶剂名称
    accuracy_level = Column(String(20), default="standard")  # fast, standard, accurate, custom

    # 复用已有计算结果
    is_reused = Column(Boolean, default=False)  # 是否复用已有计算结果
    reused_from_job_id = Column(Integer, ForeignKey("qc_jobs.id", ondelete="SET NULL"), nullable=True)  # 原始任务ID

    # 额外配置
    config = Column(JSONB, default={})
    
    # 任务状态
    status = Column(Enum(QCJobStatus), default=QCJobStatus.CREATED, nullable=False, index=True)
    slurm_job_id = Column(String(50), index=True)
    progress = Column(Float, default=0.0)
    work_dir = Column(Text)
    log_file = Column(Text)
    error_message = Column(Text)
    
    # 时间戳
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())
    started_at = Column(DateTime(timezone=True))
    finished_at = Column(DateTime(timezone=True))

    # 软删除字段
    is_deleted = Column(Boolean, default=False, index=True)  # 是否已删除
    deleted_at = Column(DateTime(timezone=True))  # 删除时间
    deleted_by = Column(Integer, ForeignKey("users.id", ondelete="SET NULL"))  # 删除操作者
    delete_reason = Column(String(500))  # 删除原因

    # 可见性管理字段（QC数据默认公开）
    visibility = Column(String(20), default="PUBLIC", index=True)  # PUBLIC, DELAYED, PRIVATE
    visibility_delay_until = Column(DateTime(timezone=True))  # 延期公开日期

    # 统计信息
    view_count = Column(Integer, default=0)  # 查看次数
    download_count = Column(Integer, default=0)  # 下载次数

    # 关系
    user = relationship("User", back_populates="qc_jobs", foreign_keys=[user_id])
    md_job = relationship("MDJob", back_populates="qc_jobs")
    results = relationship("QCResult", back_populates="qc_job", cascade="all, delete-orphan")
    deleted_by_user = relationship("User", foreign_keys=[deleted_by])

    # Slurm资源配置
    slurm_partition = Column(String(50), default="cpu")
    slurm_cpus = Column(Integer, default=16)
    slurm_time = Column(Integer, default=7200)  # 分钟

    __table_args__ = (
        Index('idx_qc_jobs_user_status', 'user_id', 'status'),
        Index('idx_qc_jobs_created_at', 'created_at'),
    )

    @property
    def solvent_config(self):
        """返回溶剂配置对象"""
        # 优先使用数据库字段
        if self.solvent_model and self.solvent_model != 'gas':
            return {
                'model': self.solvent_model,
                'solvent_name': self.solvent_name
            }
        # 兼容旧数据：从config中读取
        if self.config and 'solvent_config' in self.config:
            sc = self.config['solvent_config']
            if sc and sc.get('model') and sc.get('model') != 'gas':
                return sc
        return None

    def __repr__(self):
        return f"<QCJob(id={self.id}, molecule={self.molecule_name}, status={self.status})>"


class QCResult(Base):
    """量子化学计算结果模型"""
    __tablename__ = "qc_results"

    id = Column(Integer, primary_key=True, index=True)
    qc_job_id = Column(Integer, ForeignKey("qc_jobs.id", ondelete="CASCADE"), nullable=False, index=True)
    smiles = Column(Text, nullable=False, index=True)
    
    # 核心结果
    energy_au = Column(Float)  # 能量 (A.U.)
    homo = Column(Float)  # HOMO 能量 (Hartree)
    lumo = Column(Float)  # LUMO 能量 (Hartree)
    homo_lumo_gap = Column(Float)  # HOMO-LUMO 能隙 (eV)
    
    # ESP 相关
    esp_min_kcal = Column(Float)  # ESP最小值 (kcal/mol)
    esp_max_kcal = Column(Float)  # ESP最大值 (kcal/mol)
    
    # 文件路径
    esp_image_path = Column(Text)
    homo_image_path = Column(Text)  # HOMO轨道图片
    lumo_image_path = Column(Text)  # LUMO轨道图片
    fchk_file_path = Column(Text)
    log_file_path = Column(Text)
    cube_density_path = Column(Text)
    cube_esp_path = Column(Text)
    cube_homo_path = Column(Text)  # HOMO轨道cube文件
    cube_lumo_path = Column(Text)  # LUMO轨道cube文件

    # 图片内容（base64编码，用于混合云架构）
    homo_image_content = Column(Text)  # HOMO轨道图片base64
    lumo_image_content = Column(Text)  # LUMO轨道图片base64
    esp_image_content = Column(Text)  # ESP图片base64
    
    # 其他属性
    dipole_moment = Column(Float)
    polarizability = Column(Float)
    additional_properties = Column(JSONB, default={})

    # VIP/VEA 相关字段（用于电化学窗口估计）
    vip_ev = Column(Float)  # 垂直电离势 VIP (eV)
    vea_ev = Column(Float)  # 垂直电子亲和能 VEA (eV)
    oxidation_potential_v = Column(Float)  # 氧化电位 vs Li/Li+ (V)
    reduction_potential_v = Column(Float)  # 还原电位 vs Li/Li+ (V)

    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)

    # 关系
    qc_job = relationship("QCJob", back_populates="results")

    def __repr__(self):
        return f"<QCResult(id={self.id}, smiles={self.smiles[:20]}...)>"

    @property
    def homo_ev(self) -> float:
        """HOMO能量转换为eV"""
        if self.homo is None:
            return None
        return self.homo * 27.2114  # Hartree to eV

    @property
    def lumo_ev(self) -> float:
        """LUMO能量转换为eV"""
        if self.lumo is None:
            return None
        return self.lumo * 27.2114  # Hartree to eV


class MoleculeQCCache(Base):
    """分子QC结果缓存表"""
    __tablename__ = "molecule_qc_cache"

    id = Column(Integer, primary_key=True, index=True)
    smiles = Column(Text, unique=True, nullable=False, index=True)
    molecule_name = Column(String(255))
    
    # 关联的QC结果
    preferred_qc_result_id = Column(Integer, ForeignKey("qc_results.id", ondelete="SET NULL"))
    basis_set = Column(String(50))
    functional = Column(String(50))
    
    # 缓存的核心数据 (eV单位)
    energy_au = Column(Float)
    homo_ev = Column(Float)
    lumo_ev = Column(Float)
    homo_lumo_gap_ev = Column(Float)
    esp_min_kcal = Column(Float)
    esp_max_kcal = Column(Float)
    
    # 图像路径
    esp_image_path = Column(Text)
    homo_image_path = Column(Text)  # HOMO轨道图片
    lumo_image_path = Column(Text)  # LUMO轨道图片

    # 统计
    calculation_count = Column(Integer, default=1)
    
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())
    
    def __repr__(self):
        return f"<MoleculeQCCache(smiles={self.smiles[:20]}...)>"

