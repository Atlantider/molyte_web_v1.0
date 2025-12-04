"""
Binding Analysis schemas
Li-配体 Binding Energy 分析相关的 Pydantic 模式
"""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from datetime import datetime
from enum import Enum


class BindingAnalysisStatus(str, Enum):
    """Binding 分析任务状态"""
    CREATED = "CREATED"
    SUBMITTED = "SUBMITTED"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"


class BindingAnalysisConfig(BaseModel):
    """Binding 分析任务配置"""
    composition_keys: List[str] = Field(default=[], description="要分析的 cluster 组成列表")
    functional: str = Field(default="B3LYP", description="DFT 泛函")
    basis_set: str = Field(default="6-31G(d)", description="基组")
    solvent_model: str = Field(default="gas", description="溶剂模型: gas, pcm, smd")
    solvent_name: Optional[str] = Field(default=None, description="溶剂名称（如使用 PCM/SMD）")
    reuse_existing_qc: bool = Field(default=True, description="是否复用已有 QC 结果")


class BindingAnalysisCreate(BaseModel):
    """创建 Binding 分析任务"""
    md_job_id: int = Field(..., description="关联的 MD 任务 ID")
    config: BindingAnalysisConfig = Field(default_factory=BindingAnalysisConfig)


class ClusterBindingResult(BaseModel):
    """单个 Cluster 的 Binding 结果"""
    composition_key: str = Field(..., description="配位组成 key，如 'Li_EC3_PF6'")
    cluster_energy_au: Optional[float] = Field(None, description="Cluster 能量 (Hartree)")
    center_ion_energy_au: Optional[float] = Field(None, description="中心离子能量 (Hartree)")
    ligand_energies_au: Dict[str, float] = Field(default_factory=dict, description="各配体类型能量 {type: energy}")
    ligand_counts: Dict[str, int] = Field(default_factory=dict, description="各配体数量 {type: count}")
    
    # 计算得到的 Binding
    binding_energy_au: Optional[float] = Field(None, description="总 binding energy (Hartree)")
    binding_energy_kcal: Optional[float] = Field(None, description="总 binding energy (kcal/mol)")
    per_ligand_binding_kcal: Dict[str, float] = Field(default_factory=dict, description="平均每个配体的 binding")
    
    # QC 任务关联
    cluster_qc_job_id: Optional[int] = None
    center_ion_qc_job_id: Optional[int] = None
    ligand_qc_job_ids: Dict[str, int] = Field(default_factory=dict)
    
    # 状态
    converged: bool = Field(default=False)
    warnings: List[str] = Field(default_factory=list)


class TypeBindingStatistics(BaseModel):
    """按配体类型的 Binding 统计"""
    ligand_type: str
    count: int
    mean_binding_kcal: float
    std_binding_kcal: float
    min_binding_kcal: float
    max_binding_kcal: float


class BindingAnalysisSummary(BaseModel):
    """Binding 分析汇总"""
    total_clusters: int = 0
    completed_clusters: int = 0
    failed_clusters: int = 0
    
    # 整体统计
    mean_total_binding_kcal: Optional[float] = None
    std_total_binding_kcal: Optional[float] = None
    
    # 按类型统计
    per_type_statistics: List[TypeBindingStatistics] = Field(default_factory=list)
    
    # 警告
    warnings: List[str] = Field(default_factory=list)


class BindingAnalysisResult(BaseModel):
    """Binding 分析结果"""
    per_cluster_results: List[ClusterBindingResult] = Field(default_factory=list)
    summary: BindingAnalysisSummary = Field(default_factory=BindingAnalysisSummary)


class BindingAnalysisJob(BaseModel):
    """Binding 分析任务响应"""
    id: int
    md_job_id: int
    user_id: int
    status: BindingAnalysisStatus
    progress: float = 0.0
    error_message: Optional[str] = None
    config: Optional[Dict[str, Any]] = None
    result: Optional[Dict[str, Any]] = None
    qc_job_ids: Optional[List[int]] = None
    created_at: datetime
    updated_at: datetime
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None

    class Config:
        from_attributes = True


class BindingAnalysisJobList(BaseModel):
    """Binding 分析任务列表响应"""
    items: List[BindingAnalysisJob]
    total: int

