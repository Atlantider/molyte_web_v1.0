"""
Job schemas
"""
from pydantic import BaseModel, Field
from typing import Optional, Dict, Any, List
from datetime import datetime
from app.models.job import JobStatus, PostprocessType
from app.schemas.accuracy_level import AccuracyLevel


class MDJobBase(BaseModel):
    """Base MD job schema"""
    system_id: int


class QCResultSummary(BaseModel):
    """QC结果简要信息（用于MD任务关联展示）"""
    energy_au: Optional[float] = None  # 能量 (A.U.)
    homo_ev: Optional[float] = None  # HOMO能量 (eV)
    lumo_ev: Optional[float] = None  # LUMO能量 (eV)
    homo_lumo_gap: Optional[float] = None  # HOMO-LUMO能隙 (eV)
    esp_min_kcal: Optional[float] = None  # ESP最小值 (kcal/mol)
    esp_max_kcal: Optional[float] = None  # ESP最大值 (kcal/mol)
    dipole_moment: Optional[float] = None  # 偶极矩
    has_esp_image: bool = False  # 是否有ESP图片
    has_homo_image: bool = False  # 是否有HOMO图片
    has_lumo_image: bool = False  # 是否有LUMO图片

    class Config:
        from_attributes = True


class QCJobSummary(BaseModel):
    """QC任务简要信息（用于MD任务关联展示）"""
    id: int
    molecule_name: str
    smiles: str
    molecule_type: str
    status: str
    progress: float
    basis_set: str
    functional: str
    charge: int = 0
    spin_multiplicity: int = 1
    solvent_model: Optional[str] = None  # gas, pcm, smd
    solvent_name: Optional[str] = None  # 隐式溶剂名称
    accuracy_level: Optional[str] = None  # fast, standard, accurate, custom
    is_reused: bool = False  # 是否复用已有计算结果
    reused_from_job_id: Optional[int] = None  # 如果是复用，原始任务ID
    slurm_job_id: Optional[str] = None
    work_dir: Optional[str] = None
    created_at: datetime
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None
    error_message: Optional[str] = None
    # 计算结果（仅已完成任务有）
    result: Optional[QCResultSummary] = None

    class Config:
        from_attributes = True


class QCJobsStatusSummary(BaseModel):
    """QC任务状态汇总"""
    total: int = 0
    created: int = 0
    queued: int = 0
    running: int = 0
    postprocessing: int = 0
    completed: int = 0
    failed: int = 0
    cancelled: int = 0

    @property
    def all_completed(self) -> bool:
        """所有任务是否都已完成或失败"""
        return self.total > 0 and (self.completed + self.failed + self.cancelled) == self.total

    @property
    def has_running(self) -> bool:
        """是否有正在运行的任务"""
        return (self.queued + self.running + self.postprocessing) > 0


class MDJobQCOptions(BaseModel):
    """QC计算选项（用于MD任务附带QC计算）"""
    enabled: bool = False
    molecules: Optional[List[str]] = None  # SMILES列表，为空则自动从电解质配方提取
    accuracy_level: Optional[str] = "standard"  # fast, standard, accurate, custom

    # 计算参数（支持多选）
    basis_sets: Optional[List[str]] = ["6-31++g(d,p)"]  # 基组列表
    functionals: Optional[List[str]] = ["B3LYP"]  # 泛函列表
    solvent_models: Optional[List[str]] = ["pcm"]  # 溶剂模型列表: gas, pcm, smd
    solvents: Optional[List[str]] = ["Water"]  # 溶剂列表

    # 兼容旧版单选字段
    basis_set: Optional[str] = None  # 单个基组（兼容旧版）
    functional: Optional[str] = None  # 单个泛函（兼容旧版）
    solvent_model: Optional[str] = None  # 单个溶剂模型（兼容旧版）
    solvent_name: Optional[str] = None  # 单个溶剂名称（兼容旧版）

    # 高级选项
    use_recommended_params: bool = True  # 是否对不同分子类型使用推荐参数


class MDJobCreate(MDJobBase):
    """Schema for creating an MD job"""
    job_name: Optional[str] = None
    accuracy_level: Optional[AccuracyLevel] = AccuracyLevel.STANDARD  # 默认标准模式
    nsteps_npt: Optional[int] = None  # 如果为 None，使用 accuracy_level 的默认值
    nsteps_nvt: Optional[int] = None
    timestep: Optional[float] = 1.0
    temperature: Optional[float] = None
    pressure: Optional[float] = None
    freq_trj_npt: Optional[int] = None
    freq_trj_nvt: Optional[int] = None
    thermo_freq: Optional[int] = None
    submit_to_cluster: Optional[bool] = False

    # Slurm 资源配置
    slurm_partition: Optional[str] = "cpu"  # 队列/分区
    slurm_nodes: Optional[int] = 1  # 节点数
    slurm_ntasks: Optional[int] = 8  # 任务数
    slurm_cpus_per_task: Optional[int] = 8  # 每个任务的 CPU 核心数
    slurm_time: Optional[int] = 7200  # 最大运行时间（分钟）

    # QC计算选项
    qc_options: Optional[MDJobQCOptions] = None


class BatchMDJobCreate(BaseModel):
    """Schema for batch creating MD jobs"""
    system_ids: List[int] = Field(..., description="配方系统ID列表")
    job_name: Optional[str] = None
    accuracy_level: Optional[AccuracyLevel] = AccuracyLevel.STANDARD
    nsteps_npt: Optional[int] = None
    nsteps_nvt: Optional[int] = None
    timestep: Optional[float] = 1.0
    temperature: Optional[float] = None
    pressure: Optional[float] = None
    freq_trj_npt: Optional[int] = None
    freq_trj_nvt: Optional[int] = None
    thermo_freq: Optional[int] = None
    submit_to_cluster: Optional[bool] = False
    slurm_partition: Optional[str] = "cpu"
    slurm_nodes: Optional[int] = 1
    slurm_ntasks: Optional[int] = 8
    slurm_cpus_per_task: Optional[int] = 8
    slurm_time: Optional[int] = 7200
    qc_options: Optional[MDJobQCOptions] = None


class MDJobUpdate(BaseModel):
    """Schema for updating an MD job"""
    status: Optional[JobStatus] = None
    progress: Optional[float] = Field(None, ge=0, le=100)
    slurm_job_id: Optional[str] = None
    work_dir: Optional[str] = None
    log_file: Optional[str] = None
    error_message: Optional[str] = None


class MDJobInDB(MDJobBase):
    """MD job schema with database fields"""
    id: int
    user_id: int
    status: JobStatus
    slurm_job_id: Optional[str] = None
    progress: float
    work_dir: Optional[str] = None
    log_file: Optional[str] = None
    error_message: Optional[str] = None
    config: Optional[Dict[str, Any]] = None
    created_at: datetime
    updated_at: datetime
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None

    class Config:
        from_attributes = True


class MDJob(MDJobInDB):
    """MD job response schema"""
    username: Optional[str] = None  # 用户名，用于管理端显示
    user_email: Optional[str] = None  # 用户邮箱，用于管理端显示


class MDJobWithQC(MDJob):
    """MD任务响应（包含关联的QC任务信息）"""
    qc_jobs: List[QCJobSummary] = []
    qc_status_summary: Optional[QCJobsStatusSummary] = None

    class Config:
        from_attributes = True


class PostprocessJobBase(BaseModel):
    """Base postprocess job schema"""
    md_job_id: int
    job_type: PostprocessType
    config: Optional[Dict[str, Any]] = None


class PostprocessJobCreate(PostprocessJobBase):
    """Schema for creating a postprocess job"""
    pass


class PostprocessJobUpdate(BaseModel):
    """Schema for updating a postprocess job"""
    status: Optional[JobStatus] = None
    output_file: Optional[str] = None
    error_message: Optional[str] = None


class PostprocessJobInDB(PostprocessJobBase):
    """Postprocess job schema with database fields"""
    id: int
    status: JobStatus
    output_file: Optional[str] = None
    error_message: Optional[str] = None
    created_at: datetime
    updated_at: datetime
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None
    
    class Config:
        from_attributes = True


class PostprocessJob(PostprocessJobInDB):
    """Postprocess job response schema"""
    pass

