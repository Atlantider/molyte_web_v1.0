"""
热力学循环计算氧化还原电位 - Pydantic Schemas

⚠️ 高风险警告：
- 结果对方法/基组/溶剂模型/构型高度敏感
- 计算量大，经常不收敛
- 数值可能存在数百 mV 的系统性偏差
- 仅供研究参考，不应作为定量预测
"""
from datetime import datetime
from enum import Enum
from typing import Optional, List, Dict, Any
from pydantic import BaseModel, Field


# ============================================================================
# 物理常数
# ============================================================================

class PhysicalConstants:
    """物理常数（用于电位计算）"""
    FARADAY_C_MOL = 96485.33  # C/mol
    FARADAY_KCAL_MOL_V = 23.061  # kcal/(mol·V)
    HARTREE_TO_KCAL = 627.509  # kcal/mol
    HARTREE_TO_EV = 27.2114  # eV

    # Li+/Li 绝对电位参考值（vs SHE），文献有差异
    # 常用值: -3.04 V (Trasatti), -3.04 V (IUPAC推荐)
    LI_ABSOLUTE_POTENTIAL_VS_SHE = -3.04  # V

    # SHE 绝对电位（vs 真空），文献有差异
    # 常用值: 4.44 V (Trasatti), 4.28 V (Kelly 2006)
    SHE_ABSOLUTE_POTENTIAL = 4.44  # V


# ============================================================================
# 枚举类型
# ============================================================================

class RedoxCalculationMode(str, Enum):
    """计算模式"""
    CHEAP = "cheap"        # 垂直近似 + 单点溶剂化（低机时，仅供趋势参考）
    STANDARD = "standard"  # gas优化+freq，solvent单点（中等机时，半定量）
    HEAVY = "heavy"        # gas+solvent都优化+freq（极高机时，较准确）


class RedoxJobStatus(str, Enum):
    """任务状态"""
    CREATED = "CREATED"
    SUBMITTED = "SUBMITTED"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"


class RedoxType(str, Enum):
    """氧化还原类型"""
    OXIDATION = "oxidation"    # 氧化（失电子）
    REDUCTION = "reduction"    # 还原（得电子）


# ============================================================================
# 物种配置
# ============================================================================

class SpeciesConfig(BaseModel):
    """单个物种的配置（支持两种方式：基于已有 QC 任务或手动输入）"""
    name: str = Field(..., description="物种名称/cluster 类型名")

    # 方式1：基于已有 QC 任务（推荐）
    qc_job_id: Optional[int] = Field(None, description="已完成的 QC 任务 ID（优先使用此方式）")

    # 方式2：手动输入（不推荐，仅用于单分子测试）
    smiles: Optional[str] = Field(None, description="SMILES 字符串")
    xyz_content: Optional[str] = Field(None, description="XYZ 坐标（如果不从 SMILES 生成）")

    charge: int = Field(0, description="电荷")
    multiplicity: int = Field(1, description="自旋多重度")
    redox_type: RedoxType = Field(RedoxType.OXIDATION, description="氧化还原类型")


# ============================================================================
# 任务配置
# ============================================================================

class RedoxJobConfig(BaseModel):
    """热力学循环任务配置"""
    species_list: List[SpeciesConfig] = Field(
        ...,
        description="物种列表（建议不超过10个）",
        max_length=20
    )
    mode: RedoxCalculationMode = Field(
        RedoxCalculationMode.CHEAP,
        description="计算模式"
    )
    functional: str = Field("B3LYP", description="DFT 泛函")
    basis_set: str = Field("6-31G*", description="基组")
    solvent_model: str = Field("SMD", description="溶剂模型 (SMD/PCM/CPCM)")
    solvent: str = Field("water", description="溶剂名称")

    # 高级选项
    use_dispersion: bool = Field(True, description="是否使用色散校正 (D3BJ)")
    scf_max_cycles: int = Field(200, description="SCF 最大迭代次数")
    opt_max_cycles: int = Field(100, description="优化最大迭代次数")

    # QC 复用选项
    reuse_existing_qc: bool = Field(
        True,
        description="是否复用已有 QC 结果（全局复用，节省计算资源）"
    )

    # Li+/Li 参考电位（可自定义）
    li_reference_potential: float = Field(
        PhysicalConstants.LI_ABSOLUTE_POTENTIAL_VS_SHE,
        description="Li+/Li 绝对电位 vs SHE (V)，文献有差异，默认 -3.04 V"
    )


class RedoxJobCreate(BaseModel):
    """创建热力学循环任务"""
    md_job_id: Optional[int] = Field(None, description="关联的 MD 任务 ID（可选）")
    config: RedoxJobConfig


# ============================================================================
# 单个物种的计算结果
# ============================================================================

class SpeciesRedoxResult(BaseModel):
    """单个物种的氧化还原结果"""
    name: str
    redox_type: RedoxType

    # 能量数据 (Hartree)
    e_neutral_gas: Optional[float] = Field(None, description="中性态气相能量")
    e_charged_gas: Optional[float] = Field(None, description="带电态气相能量")
    e_neutral_sol: Optional[float] = Field(None, description="中性态溶液能量")
    e_charged_sol: Optional[float] = Field(None, description="带电态溶液能量")

    # 热力学量 (kcal/mol)
    dg_gas_kcal: Optional[float] = Field(None, description="气相反应自由能")
    dg_solv_neutral_kcal: Optional[float] = Field(None, description="中性态溶剂化能")
    dg_solv_charged_kcal: Optional[float] = Field(None, description="带电态溶剂化能")
    dg_sol_kcal: Optional[float] = Field(None, description="溶液中反应自由能")

    # 电位 (V)
    e_abs_v: Optional[float] = Field(None, description="绝对电位 vs SHE")
    e_vs_li_v: Optional[float] = Field(None, description="电位 vs Li+/Li")

    # 状态
    converged: bool = Field(True, description="是否收敛")
    warnings: List[str] = Field(default_factory=list, description="警告信息")
    qc_job_ids: List[int] = Field(default_factory=list, description="关联的 QC 任务 ID")


# ============================================================================
# 任务结果
# ============================================================================

class RedoxJobResult(BaseModel):
    """热力学循环任务结果"""
    species_results: List[SpeciesRedoxResult] = Field(
        default_factory=list,
        description="各物种的计算结果"
    )

    # 汇总统计
    oxidation_potentials_v: List[float] = Field(
        default_factory=list,
        description="所有氧化电位 vs Li+/Li"
    )
    reduction_potentials_v: List[float] = Field(
        default_factory=list,
        description="所有还原电位 vs Li+/Li"
    )

    # 电化学窗口估计
    oxidation_limit_v: Optional[float] = Field(
        None, description="氧化极限 (5% 分位数) vs Li+/Li"
    )
    reduction_limit_v: Optional[float] = Field(
        None, description="还原极限 (95% 分位数) vs Li+/Li"
    )
    electrochemical_window_v: Optional[float] = Field(
        None, description="电化学窗口宽度"
    )

    # 警告和元数据
    global_warnings: List[str] = Field(
        default_factory=list,
        description="全局警告信息"
    )
    calculation_mode: str = Field("cheap", description="使用的计算模式")
    reference_note: str = Field(
        "电位参考: Li+/Li, 计算方法存在系统性偏差，仅供研究参考",
        description="参考说明"
    )


# ============================================================================
# API 响应
# ============================================================================

class RedoxJobResponse(BaseModel):
    """热力学循环任务响应"""
    id: int
    md_job_id: Optional[int]
    status: RedoxJobStatus
    progress: float = Field(0.0, ge=0.0, le=100.0)
    error_message: Optional[str] = None

    config: Optional[Dict[str, Any]] = None
    result: Optional[RedoxJobResult] = None

    created_at: datetime
    updated_at: Optional[datetime] = None
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None

    class Config:
        from_attributes = True


class RedoxJobListResponse(BaseModel):
    """热力学循环任务列表响应"""
    total: int
    jobs: List[RedoxJobResponse]


# ============================================================================
# 重组能计算 (Marcus 理论)
# ============================================================================

class ReorgEnergyJobStatus(str, Enum):
    """重组能任务状态"""
    CREATED = "CREATED"
    SUBMITTED = "SUBMITTED"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"


class ReorgSpeciesConfig(BaseModel):
    """重组能物种配置（支持两种方式：基于已有 QC 任务或手动输入）"""
    name: str = Field(..., description="物种名称/cluster 类型名")

    # 方式1：基于已有 QC 任务（推荐）
    qc_job_id: Optional[int] = Field(None, description="已完成的 QC 任务 ID（优先使用此方式）")

    # 方式2：手动输入（不推荐）
    smiles: Optional[str] = Field(None, description="SMILES 字符串")
    xyz_content: Optional[str] = Field(None, description="XYZ 坐标")

    charge_neutral: int = Field(0, description="中性态电荷")
    charge_oxidized: int = Field(1, description="氧化态电荷")
    multiplicity_neutral: int = Field(1, description="中性态自旋多重度")
    multiplicity_oxidized: int = Field(2, description="氧化态自旋多重度")


class ReorgEnergyJobConfig(BaseModel):
    """重组能任务配置

    ⚠️ 极高风险警告：
    - 每个物种至少 2 次优化 + 4 次单点
    - Cluster 体系极易不收敛
    - 构型依赖极强
    - 默认限制：最多 5 个物种
    """
    species_list: List[ReorgSpeciesConfig] = Field(
        ...,
        description="物种列表（强烈建议不超过5个）",
        max_length=10
    )
    functional: str = Field("B3LYP", description="DFT 泛函")
    basis_set: str = Field("6-31G*", description="基组")
    use_dispersion: bool = Field(True, description="是否使用色散校正")
    scf_max_cycles: int = Field(200, description="SCF 最大迭代次数")
    opt_max_cycles: int = Field(150, description="优化最大迭代次数")

    # QC 复用选项
    reuse_existing_qc: bool = Field(
        True,
        description="是否复用已有 QC 结果（全局复用，节省计算资源）"
    )


class ReorgEnergyJobCreate(BaseModel):
    """创建重组能任务"""
    md_job_id: Optional[int] = Field(None, description="关联的 MD 任务 ID")
    config: ReorgEnergyJobConfig


class SpeciesReorgResult(BaseModel):
    """单个物种的重组能结果"""
    name: str

    # 4 点能量 (Hartree)
    e_neutral_at_neutral_geom: Optional[float] = Field(
        None, description="E(q, R_q): 中性态在中性几何"
    )
    e_oxidized_at_neutral_geom: Optional[float] = Field(
        None, description="E(q+1, R_q): 氧化态在中性几何"
    )
    e_neutral_at_oxidized_geom: Optional[float] = Field(
        None, description="E(q, R_{q+1}): 中性态在氧化几何"
    )
    e_oxidized_at_oxidized_geom: Optional[float] = Field(
        None, description="E(q+1, R_{q+1}): 氧化态在氧化几何"
    )

    # 重组能 (eV)
    lambda_ox_ev: Optional[float] = Field(
        None, description="氧化重组能 λ_ox (eV)"
    )
    lambda_red_ev: Optional[float] = Field(
        None, description="还原重组能 λ_red (eV)"
    )
    lambda_total_ev: Optional[float] = Field(
        None, description="总重组能 λ = (λ_ox + λ_red) / 2 (eV)"
    )

    converged: bool = Field(True, description="是否收敛")
    warnings: List[str] = Field(default_factory=list)
    qc_job_ids: List[int] = Field(default_factory=list)


class ReorgEnergyJobResult(BaseModel):
    """重组能任务结果"""
    species_results: List[SpeciesReorgResult] = Field(default_factory=list)

    # 汇总
    lambda_ox_mean_ev: Optional[float] = None
    lambda_red_mean_ev: Optional[float] = None
    lambda_total_mean_ev: Optional[float] = None

    global_warnings: List[str] = Field(default_factory=list)
    reference_note: str = Field(
        "重组能计算基于 Marcus 理论 4 点方案，结果对构型极其敏感，仅供研究参考",
        description="参考说明"
    )


class ReorgEnergyJobResponse(BaseModel):
    """重组能任务响应"""
    id: int
    md_job_id: Optional[int]
    status: ReorgEnergyJobStatus
    progress: float = Field(0.0, ge=0.0, le=100.0)
    error_message: Optional[str] = None

    config: Optional[Dict[str, Any]] = None
    result: Optional[ReorgEnergyJobResult] = None

    created_at: datetime
    updated_at: Optional[datetime] = None
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None

    class Config:
        from_attributes = True
