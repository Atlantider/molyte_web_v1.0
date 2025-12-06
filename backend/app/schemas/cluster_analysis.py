"""
Cluster Analysis Schemas - 统一的 Cluster 高级计算规划
"""
from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
from enum import Enum
from datetime import datetime


class ClusterCalcType(str, Enum):
    """Cluster 计算类型"""
    BINDING_TOTAL = "BINDING_TOTAL"           # Binding A: 总脱溶剂化能
    BINDING_PAIRWISE = "BINDING_PAIRWISE"     # Binding B: 单分子-Li Binding
    DESOLVATION_STEPWISE = "DESOLVATION_STEPWISE"  # 逐级去溶剂化能
    DESOLVATION_FULL = "DESOLVATION_FULL"     # 完全去溶剂化能
    REDOX = "REDOX"                           # 氧化还原电位
    REORGANIZATION = "REORGANIZATION"         # Marcus 重组能


class AdvancedClusterJobStatus(str, Enum):
    """高级 Cluster 计算任务状态"""
    CREATED = "CREATED"
    SUBMITTED = "SUBMITTED"
    RUNNING = "RUNNING"
    WAITING_QC = "WAITING_QC"
    CALCULATING = "CALCULATING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"


# ============================================================================
# 规划请求/响应
# ============================================================================

class QCConfig(BaseModel):
    """QC 计算配置"""
    functional: str = Field("B3LYP", description="DFT 泛函")
    basis_set: str = Field("6-31G*", description="基组")
    solvent_model: Optional[str] = Field(None, description="溶剂模型 (SMD/PCM/None)")
    solvent: Optional[str] = Field(None, description="溶剂名称")
    use_dispersion: bool = Field(True, description="是否使用色散校正 (D3BJ)")
    charge_cluster: int = Field(1, description="Cluster 电荷")
    charge_ion: int = Field(1, description="离子电荷 (Li+ = 1)")
    # Slurm 资源配置
    slurm_partition: Optional[str] = Field("cpu", description="Slurm 队列/分区")
    slurm_cpus: Optional[int] = Field(16, ge=1, le=128, description="CPU 核心数")
    slurm_time: Optional[int] = Field(7200, ge=10, description="最大运行时间（分钟）")


class RedoxOptions(BaseModel):
    """REDOX 子选项"""
    include_molecule: bool = Field(True, description="计算单独配体分子的氧化还原电位")
    include_dimer: bool = Field(True, description="计算 Li-配体 Dimer 的氧化还原电位")
    include_cluster: bool = Field(False, description="计算整个 Cluster 的氧化还原电位（计算量大）")


class ReorganizationOptions(BaseModel):
    """REORGANIZATION 子选项"""
    include_molecule: bool = Field(True, description="计算单独配体分子的重组能")
    include_cluster: bool = Field(True, description="计算整个 Cluster 的重组能")


class ClusterAnalysisPlanRequest(BaseModel):
    """规划请求 - 用户选择结构和计算类型"""
    md_job_id: int = Field(..., description="MD 任务 ID")

    # 选中的结构（二选一）
    solvation_structure_ids: Optional[List[int]] = Field(None, description="选中的溶剂化结构 ID 列表")
    composition_keys: Optional[List[str]] = Field(None, description="选中的 composition_key 列表")

    # 选中的计算类型
    calc_types: List[ClusterCalcType] = Field(..., description="选中的计算类型")

    # 子选项
    redox_options: Optional[RedoxOptions] = Field(default_factory=RedoxOptions, description="REDOX 子选项")
    reorganization_options: Optional[ReorganizationOptions] = Field(default_factory=ReorganizationOptions, description="REORGANIZATION 子选项")

    # QC 配置
    qc_config: Optional[QCConfig] = Field(default_factory=QCConfig, description="QC 计算配置")


class PlannedQCTask(BaseModel):
    """规划的 QC 任务"""
    task_type: str = Field(..., description="任务类型: cluster/ligand/ion/cluster_minus/dimer/charged")
    description: str = Field(..., description="任务描述")
    smiles: Optional[str] = Field(None, description="分子 SMILES")
    structure_id: Optional[int] = Field(None, description="关联的结构 ID")
    charge: int = Field(0, description="电荷")
    multiplicity: int = Field(1, description="多重度")
    calc_mode: str = Field("opt", description="计算模式: opt (几何优化) / sp (单点能量)")

    # 复用状态
    status: str = Field("new", description="状态: new/reused/pending")
    existing_qc_job_id: Optional[int] = Field(None, description="已有的 QC 任务 ID（复用时）")
    existing_energy: Optional[float] = Field(None, description="已有的能量值（Hartree）")


class CalcTypeRequirements(BaseModel):
    """某个计算类型的 QC 需求"""
    calc_type: ClusterCalcType
    description: str
    required_qc_tasks: List[PlannedQCTask]
    new_tasks_count: int = Field(0, description="需要新建的 QC 任务数")
    reused_tasks_count: int = Field(0, description="可复用的 QC 任务数")


class ClusterAnalysisPlanResponse(BaseModel):
    """规划响应 - 返回 QC 任务规划详情"""
    md_job_id: int
    
    # 选中的结构信息
    selected_structures_count: int
    selected_structure_ids: List[int]
    
    # 各计算类型的 QC 需求
    calc_requirements: List[CalcTypeRequirements]
    
    # 汇总
    total_new_qc_tasks: int = Field(..., description="总共需要新建的 QC 任务数")
    total_reused_qc_tasks: int = Field(..., description="总共可复用的 QC 任务数")
    estimated_compute_hours: float = Field(0.0, description="预估计算时间（小时）")
    
    # 警告
    warnings: List[str] = Field(default_factory=list, description="警告信息")


# ============================================================================
# 提交/创建任务
# ============================================================================

class ClusterAnalysisSubmitRequest(BaseModel):
    """提交计算任务请求"""
    md_job_id: int = Field(..., description="MD 任务 ID")
    solvation_structure_ids: Optional[List[int]] = Field(None, description="选中的溶剂化结构 ID 列表")
    composition_keys: Optional[List[str]] = Field(None, description="选中的 composition_key 列表")
    calc_types: List[ClusterCalcType] = Field(..., description="选中的计算类型")
    qc_config: Optional[QCConfig] = Field(default_factory=QCConfig, description="QC 计算配置")
    # 子选项（与 PlanRequest 保持一致）
    redox_options: Optional[RedoxOptions] = Field(default_factory=RedoxOptions, description="REDOX 子选项")
    reorganization_options: Optional[ReorganizationOptions] = Field(default_factory=ReorganizationOptions, description="REORGANIZATION 子选项")


class AdvancedClusterJobResponse(BaseModel):
    """任务响应"""
    id: int
    md_job_id: int
    user_id: int
    status: AdvancedClusterJobStatus
    progress: float
    calc_types: List[str]
    selected_structures: Dict[str, Any]
    qc_config: Dict[str, Any]
    qc_task_plan: Dict[str, Any]
    results: Dict[str, Any]
    error_message: Optional[str]
    created_at: datetime
    updated_at: datetime
    started_at: Optional[datetime]
    finished_at: Optional[datetime]
    
    class Config:
        from_attributes = True


# ============================================================================
# 追加计算
# ============================================================================

class AddCalcTypeRequest(BaseModel):
    """追加计算类型请求"""
    job_id: int = Field(..., description="已有任务 ID")
    additional_calc_types: List[ClusterCalcType] = Field(..., description="追加的计算类型")


class AddCalcTypePlanResponse(BaseModel):
    """追加计算规划响应"""
    job_id: int
    existing_calc_types: List[str]
    additional_calc_types: List[str]
    new_qc_tasks_required: int
    reused_from_existing: int
    details: List[CalcTypeRequirements]

