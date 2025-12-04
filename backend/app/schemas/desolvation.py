"""
Desolvation energy calculation schemas
"""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any, Literal
from datetime import datetime


class SolventConfigSchema(BaseModel):
    """溶剂配置（用于去溶剂化能计算的隐式溶剂模型）"""
    model: Literal["gas", "pcm", "smd", "custom"] = Field(default="gas", description="溶剂模型")
    solvent_name: Optional[str] = Field(default=None, description="预定义溶剂名称")
    # 自定义溶剂参数（SMD模型需要）
    eps: Optional[float] = Field(default=None, description="介电常数 ε")
    eps_inf: Optional[float] = Field(default=None, description="光学介电常数 n²")
    hbond_acidity: Optional[float] = Field(default=None, description="Abraham氢键酸度 α")
    hbond_basicity: Optional[float] = Field(default=None, description="Abraham氢键碱度 β")
    surface_tension: Optional[float] = Field(default=None, description="表面张力 γ (cal/mol·Å²)")
    carbon_aromaticity: Optional[float] = Field(default=None, description="芳香碳原子比例 φ")
    halogenicity: Optional[float] = Field(default=None, description="卤素原子比例 ψ")


class DesolvationJobCreate(BaseModel):
    """创建去溶剂化能任务"""
    md_job_id: int = Field(..., description="MD job ID")
    solvation_structure_id: int = Field(..., description="Solvation structure ID")
    method_level: str = Field(default="standard", description="Calculation method level: fast (6-31G(d)/B3LYP), standard (6-31++G(d,p)/B3LYP), accurate (6-311++G(2d,2p)/wB97XD)")
    desolvation_mode: Literal["stepwise", "full"] = Field(
        default="stepwise",
        description="Desolvation mode: stepwise (remove one ligand at a time) or full (remove all ligands at once)"
    )
    # 溶剂配置（可选，默认气相计算）
    solvent_config: Optional[SolventConfigSchema] = Field(
        default=None,
        description="Solvent configuration for implicit solvent model. If not provided, gas phase calculation is used."
    )


class LigandDesolvationResult(BaseModel):
    """单个配体的去溶剂化能结果"""
    ligand_id: int = Field(..., description="Ligand ID in the cluster")
    ligand_type: str = Field(..., description="Ligand type (e.g., FSI, EC, DMC)")
    ligand_label: str = Field(..., description="Ligand label (e.g., FSI_1, EC_2)")
    e_ligand: float = Field(..., description="Isolated ligand energy in A.U.")
    e_cluster_minus: float = Field(..., description="Cluster energy without this ligand in A.U.")
    delta_e: float = Field(..., description="Desolvation energy in kcal/mol")


class TypeSummary(BaseModel):
    """按类型汇总的统计"""
    ligand_type: str = Field(..., description="Ligand type")
    avg_delta_e: float = Field(..., description="Average desolvation energy in kcal/mol")
    std_delta_e: float = Field(..., description="Standard deviation in kcal/mol")
    count: int = Field(..., description="Number of ligands of this type")
    min_delta_e: float = Field(..., description="Minimum desolvation energy in kcal/mol")
    max_delta_e: float = Field(..., description="Maximum desolvation energy in kcal/mol")


class DesolvationEnergyResultSchema(BaseModel):
    """去溶剂化能结果"""
    id: int
    postprocess_job_id: int
    solvation_structure_id: int
    method_level: str
    basis_set: Optional[str] = None
    functional: Optional[str] = None
    e_cluster: float = Field(..., description="Complete cluster energy in A.U.")
    per_ligand_results: List[LigandDesolvationResult] = Field(default_factory=list)
    per_type_summary: List[TypeSummary] = Field(default_factory=list)
    created_at: datetime

    class Config:
        from_attributes = True


class DesolvationJobResponse(BaseModel):
    """去溶剂化能任务响应"""
    job_id: int
    status: str
    method_level: str
    desolvation_mode: str = "stepwise"  # stepwise or full
    solvent_config: Optional[SolventConfigSchema] = None  # 溶剂模型配置
    created_at: datetime
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None
    elapsed_seconds: Optional[float] = None
    error_message: Optional[str] = None
    result: Optional[DesolvationEnergyResultSchema] = None
    # 溯源信息
    solvation_structure_id: Optional[int] = None
    composition_key: Optional[str] = None  # 如 "Li-EC2-DMC1-PF6_1"
    md_job_id: Optional[int] = None
    electrolyte_name: Optional[str] = None  # 电解液配方名称
    qc_progress: Optional[Dict[str, Any]] = None  # QC任务进度 {"total": 10, "completed": 5, "running": 2}

    class Config:
        from_attributes = True


class BatchDesolvationJobCreate(BaseModel):
    """批量创建去溶剂化能任务"""
    md_job_id: int = Field(..., description="MD job ID")
    structure_ids: List[int] = Field(..., description="要计算的溶剂化结构ID列表")
    method_level: str = Field(default="standard", description="计算方法级别")
    desolvation_mode: Literal["stepwise", "full"] = Field(default="stepwise")
    solvent_config: Optional[SolventConfigSchema] = None
    # Slurm 资源配置
    slurm_partition: Optional[str] = Field(default="cpu", description="Slurm队列/分区")
    slurm_cpus: Optional[int] = Field(default=16, ge=1, le=64, description="CPU核心数")
    slurm_time: Optional[int] = Field(default=7200, ge=10, description="最大运行时间（分钟）")


class BatchDesolvationJobResponse(BaseModel):
    """批量创建响应"""
    created_count: int
    skipped_count: int  # 已存在的任务数
    jobs: List[DesolvationJobResponse]


class DesolvationOverviewResponse(BaseModel):
    """去溶剂化任务总览（用于监控面板）"""
    md_job_id: int
    electrolyte_name: Optional[str] = None
    total_jobs: int
    status_summary: Dict[str, int]  # {"COMPLETED": 3, "RUNNING": 2, "QUEUED": 1}
    jobs: List[DesolvationJobResponse]

