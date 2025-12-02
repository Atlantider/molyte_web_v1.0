"""
Quantum Chemistry (QC) schemas
量子化学计算相关的Pydantic模式
"""
from pydantic import BaseModel, Field, model_validator
from typing import Optional, List, Dict, Any
from datetime import datetime
from enum import Enum


class QCJobStatus(str, Enum):
    """QC任务状态

    状态流程：CREATED -> SUBMITTED -> QUEUED -> RUNNING -> COMPLETED/FAILED
    """
    CREATED = "CREATED"
    SUBMITTED = "SUBMITTED"
    QUEUED = "QUEUED"
    RUNNING = "RUNNING"
    POSTPROCESSING = "POSTPROCESSING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"


class MoleculeType(str, Enum):
    """分子类型"""
    SOLVENT = "solvent"
    CATION = "cation"
    ANION = "anion"
    CUSTOM = "custom"


class QCAccuracyLevel(str, Enum):
    """QC计算精度等级"""
    FAST = "fast"           # 快速：HF/STO-3G
    STANDARD = "standard"   # 标准：B3LYP/6-31G(d)
    ACCURATE = "accurate"   # 精确：B3LYP/6-311++G(d,p)
    CUSTOM = "custom"       # 自定义


class SolventModel(str, Enum):
    """溶剂模型"""
    GAS = "gas"             # 气相（无溶剂）
    PCM = "pcm"             # PCM隐式溶剂
    SMD = "smd"             # SMD隐式溶剂
    CUSTOM = "custom"       # 自定义溶剂参数


class BasisSet(str, Enum):
    """常用基组"""
    STO3G = "STO-3G"
    B321G = "3-21G"
    B631G = "6-31G"
    B631GD = "6-31G(d)"
    B631GDP = "6-31G(d,p)"
    B631_PLUSPLUS_GDP = "6-31++G(d,p)"
    B6311_GDP = "6-311G(d,p)"
    B6311_PLUSPLUS_GDP = "6-311++G(d,p)"
    DEF2SVP = "Def2SVP"
    DEF2TZVP = "Def2TZVP"
    DEF2QZVP = "Def2QZVP"
    CCPVDZ = "cc-pVDZ"
    CCPVTZ = "cc-pVTZ"
    AUGCCPVDZ = "aug-cc-pVDZ"


class Functional(str, Enum):
    """常用泛函"""
    HF = "HF"
    B3LYP = "B3LYP"
    M062X = "M062X"
    WB97XD = "wB97XD"
    PBE0 = "PBE0"
    CAM_B3LYP = "CAM-B3LYP"
    B3PW91 = "B3PW91"
    BLYP = "BLYP"
    PBE = "PBE"


# 预定义的Gaussian溶剂列表
GAUSSIAN_SOLVENTS = {
    "Water": {"eps": 78.3553, "description": "水"},
    "Acetonitrile": {"eps": 35.688, "description": "乙腈"},
    "Methanol": {"eps": 32.613, "description": "甲醇"},
    "Ethanol": {"eps": 24.852, "description": "乙醇"},
    "Acetone": {"eps": 20.493, "description": "丙酮"},
    "DiMethylSulfoxide": {"eps": 46.826, "description": "二甲亚砜(DMSO)"},
    "Dichloromethane": {"eps": 8.93, "description": "二氯甲烷"},
    "Chloroform": {"eps": 4.7113, "description": "氯仿"},
    "CarbonTetraChloride": {"eps": 2.228, "description": "四氯化碳"},
    "Benzene": {"eps": 2.2706, "description": "苯"},
    "Toluene": {"eps": 2.3741, "description": "甲苯"},
    "DiethylEther": {"eps": 4.24, "description": "乙醚"},
    "TetraHydroFuran": {"eps": 7.4257, "description": "四氢呋喃(THF)"},
    "n-Hexane": {"eps": 1.8819, "description": "正己烷"},
    "CycloHexane": {"eps": 2.0165, "description": "环己烷"},
    "Heptane": {"eps": 1.9113, "description": "正庚烷"},
    "n-Octanol": {"eps": 9.8629, "description": "正辛醇"},
    "1,2-EthaneDiol": {"eps": 40.245, "description": "乙二醇"},
    "1,4-Dioxane": {"eps": 2.2099, "description": "1,4-二噁烷"},
    "NitroMethane": {"eps": 36.562, "description": "硝基甲烷"},
    "Pyridine": {"eps": 12.978, "description": "吡啶"},
    "Aniline": {"eps": 6.8882, "description": "苯胺"},
    "FormicAcid": {"eps": 51.1, "description": "甲酸"},
    "AceticAcid": {"eps": 6.2528, "description": "乙酸"},
    "n,n-DiMethylFormamide": {"eps": 37.219, "description": "N,N-二甲基甲酰胺(DMF)"},
}

# 精度等级对应的参数
QC_ACCURACY_PRESETS = {
    QCAccuracyLevel.FAST: {
        "functional": "HF",
        "basis_set": "STO-3G",
        "description": "快速计算，适合大分子初步结构优化",
        "estimated_time": "1-5分钟"
    },
    QCAccuracyLevel.STANDARD: {
        "functional": "B3LYP",
        "basis_set": "6-31G(d)",
        "description": "标准精度，适合大多数计算",
        "estimated_time": "5-30分钟"
    },
    QCAccuracyLevel.ACCURATE: {
        "functional": "B3LYP",
        "basis_set": "6-311++G(d,p)",
        "description": "高精度计算，适合需要精确结果的场景",
        "estimated_time": "30分钟-数小时"
    }
}


# ============================================================================
# 溶剂环境配置
# ============================================================================

class SolventConfig(BaseModel):
    """溶剂环境配置"""
    model: SolventModel = Field(default=SolventModel.GAS, description="溶剂模型")
    solvent_name: Optional[str] = Field(default=None, description="预定义溶剂名称")
    # 自定义溶剂参数（SMD模型需要）
    eps: Optional[float] = Field(default=None, description="介电常数 ε")
    eps_inf: Optional[float] = Field(default=None, description="光学介电常数 n²")
    hbond_acidity: Optional[float] = Field(default=None, description="Abraham氢键酸度 α")
    hbond_basicity: Optional[float] = Field(default=None, description="Abraham氢键碱度 β")
    surface_tension: Optional[float] = Field(default=None, description="表面张力 γ (cal/mol·Å²)")
    carbon_aromaticity: Optional[float] = Field(default=None, description="芳香碳原子比例 φ")
    halogenicity: Optional[float] = Field(default=None, description="卤素原子比例 ψ")


# ============================================================================
# QC Job Schemas
# ============================================================================

class QCJobBase(BaseModel):
    """QC任务基础Schema"""
    molecule_name: str = Field(
        ...,
        min_length=1,
        max_length=255,
        description="分子名称（仅支持英文、数字和符号）",
        pattern=r'^[A-Za-z0-9+\-_\s,()]+$'
    )
    smiles: str = Field(..., min_length=1, description="SMILES表达式")
    molecule_type: MoleculeType = Field(default=MoleculeType.CUSTOM, description="分子类型")
    basis_set: str = Field(default="6-31G(d)", description="基组")
    functional: str = Field(default="B3LYP", description="泛函")
    charge: int = Field(default=0, description="分子电荷")
    spin_multiplicity: int = Field(default=1, ge=1, description="自旋多重度")


class QCJobCreate(QCJobBase):
    """创建QC任务的Schema"""
    md_job_id: Optional[int] = Field(default=None, description="关联的MD任务ID（可选）")
    accuracy_level: QCAccuracyLevel = Field(default=QCAccuracyLevel.STANDARD, description="计算精度等级")
    solvent_config: Optional[SolventConfig] = Field(default=None, description="溶剂环境配置")
    auto_spin: bool = Field(default=True, description="自动计算自旋多重度")
    config: Optional[Dict[str, Any]] = Field(default={}, description="额外配置")

    # Slurm 资源配置
    slurm_partition: Optional[str] = Field(default="cpu", description="Slurm队列/分区")
    slurm_cpus: Optional[int] = Field(default=16, ge=1, le=64, description="CPU核心数")
    slurm_time: Optional[int] = Field(default=7200, ge=10, description="最大运行时间（分钟）")


class QCJobBatchCreate(BaseModel):
    """批量创建QC任务"""
    molecules: List[QCJobCreate] = Field(..., min_length=1, description="分子列表")
    md_job_id: Optional[int] = Field(default=None, description="关联的MD任务ID")
    basis_set: str = Field(default="6-31++g(d,p)", description="统一基组")
    functional: str = Field(default="B3LYP", description="统一泛函")


class QCJobUpdate(BaseModel):
    """更新QC任务的Schema"""
    status: Optional[QCJobStatus] = None
    progress: Optional[float] = Field(default=None, ge=0, le=100)
    error_message: Optional[str] = None


class QCJobEdit(BaseModel):
    """编辑QC任务参数的Schema（仅CREATED状态可编辑）"""
    molecule_name: Optional[str] = Field(None, min_length=1, max_length=255, description="分子名称")
    smiles: Optional[str] = Field(None, min_length=1, description="SMILES表达式")
    molecule_type: Optional[MoleculeType] = Field(None, description="分子类型")
    basis_set: Optional[str] = Field(None, description="基组")
    functional: Optional[str] = Field(None, description="泛函")
    charge: Optional[int] = Field(None, description="分子电荷")
    spin_multiplicity: Optional[int] = Field(None, ge=1, description="自旋多重度")
    accuracy_level: Optional[QCAccuracyLevel] = Field(None, description="计算精度等级")


class QCJobRecalculate(BaseModel):
    """重新计算QC任务的Schema"""
    functional: Optional[str] = Field(None, description="新的泛函（不指定则使用原值）")
    basis_set: Optional[str] = Field(None, description="新的基组（不指定则使用原值）")
    solvent_config: Optional[SolventConfig] = Field(None, description="新的溶剂配置（不指定则使用原值）")
    accuracy_level: Optional[QCAccuracyLevel] = Field(None, description="计算精度等级")
    slurm_partition: Optional[str] = Field(None, description="Slurm队列/分区")
    slurm_cpus: Optional[int] = Field(None, ge=1, le=64, description="CPU核心数")
    slurm_time: Optional[int] = Field(None, ge=10, description="最大运行时间（分钟）")
    solvent_config: Optional[SolventConfig] = Field(None, description="溶剂环境配置")
    slurm_partition: Optional[str] = Field(None, description="Slurm队列/分区")
    slurm_cpus: Optional[int] = Field(None, ge=1, le=64, description="CPU核心数")
    slurm_time: Optional[int] = Field(None, ge=10, description="最大运行时间（分钟）")


class QCJobInDB(QCJobBase):
    """数据库中的QC任务"""
    id: int
    user_id: int
    md_job_id: Optional[int] = None
    status: QCJobStatus
    slurm_job_id: Optional[str] = None
    progress: float = 0.0
    work_dir: Optional[str] = None
    log_file: Optional[str] = None
    error_message: Optional[str] = None
    config: Dict[str, Any] = {}
    solvent_config: Optional[SolventConfig] = None  # 溶剂配置
    accuracy_level: Optional[QCAccuracyLevel] = None  # 精度等级
    slurm_partition: Optional[str] = None  # Slurm队列
    slurm_cpus: Optional[int] = None  # CPU核心数
    slurm_time: Optional[int] = None  # 最大运行时间

    # 可见性管理字段
    visibility: Optional[str] = "PRIVATE"  # PUBLIC, DELAYED, PRIVATE
    visibility_delay_until: Optional[datetime] = None  # 延期公开日期

    # 统计信息
    view_count: Optional[int] = 0  # 查看次数
    download_count: Optional[int] = 0  # 下载次数

    created_at: datetime
    updated_at: Optional[datetime] = None
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None

    @model_validator(mode='before')
    @classmethod
    def build_solvent_config(cls, data: Any) -> Any:
        """从数据库的solvent_model和solvent_name字段构建solvent_config对象"""
        if isinstance(data, dict):
            # 如果已经有solvent_config，直接返回
            if data.get('solvent_config'):
                return data

            # 从solvent_model和solvent_name构建solvent_config
            solvent_model = data.get('solvent_model')
            solvent_name = data.get('solvent_name')

            if solvent_model:
                data['solvent_config'] = {
                    'model': solvent_model,
                    'solvent_name': solvent_name
                }
        else:
            # 处理ORM对象
            if hasattr(data, 'solvent_model') and data.solvent_model:
                if not hasattr(data, 'solvent_config') or data.solvent_config is None:
                    # 创建一个临时字典来存储solvent_config
                    # 注意：这里我们需要将ORM对象转换为字典
                    from sqlalchemy.inspection import inspect
                    mapper = inspect(data.__class__)
                    result = {}
                    for column in mapper.columns:
                        result[column.name] = getattr(data, column.name)

                    result['solvent_config'] = {
                        'model': data.solvent_model,
                        'solvent_name': data.solvent_name if hasattr(data, 'solvent_name') else None
                    }
                    return result

        return data

    class Config:
        from_attributes = True


class QCJob(QCJobInDB):
    """QC任务响应Schema"""
    pass


# ============================================================================
# QC Result Schemas
# ============================================================================

class QCResultBase(BaseModel):
    """QC结果基础Schema"""
    smiles: str
    energy_au: Optional[float] = Field(default=None, description="能量 (A.U.)")
    homo: Optional[float] = Field(default=None, description="HOMO能量 (Hartree)")
    lumo: Optional[float] = Field(default=None, description="LUMO能量 (Hartree)")
    homo_lumo_gap: Optional[float] = Field(default=None, description="HOMO-LUMO能隙 (eV)")
    esp_min_kcal: Optional[float] = Field(default=None, description="ESP最小值 (kcal/mol)")
    esp_max_kcal: Optional[float] = Field(default=None, description="ESP最大值 (kcal/mol)")
    dipole_moment: Optional[float] = Field(default=None, description="偶极矩")
    polarizability: Optional[float] = Field(default=None, description="极化率")


class QCResultCreate(QCResultBase):
    """创建QC结果的Schema"""
    qc_job_id: int
    esp_image_path: Optional[str] = None
    homo_image_path: Optional[str] = None
    lumo_image_path: Optional[str] = None
    fchk_file_path: Optional[str] = None
    log_file_path: Optional[str] = None
    cube_homo_path: Optional[str] = None
    cube_lumo_path: Optional[str] = None
    additional_properties: Optional[Dict[str, Any]] = {}


class QCResultInDB(QCResultBase):
    """数据库中的QC结果"""
    id: int
    qc_job_id: int
    esp_image_path: Optional[str] = None
    homo_image_path: Optional[str] = None
    lumo_image_path: Optional[str] = None
    fchk_file_path: Optional[str] = None
    log_file_path: Optional[str] = None
    cube_density_path: Optional[str] = None
    cube_esp_path: Optional[str] = None
    cube_homo_path: Optional[str] = None
    cube_lumo_path: Optional[str] = None
    additional_properties: Dict[str, Any] = {}
    created_at: datetime

    # 计算属性
    homo_ev: Optional[float] = None
    lumo_ev: Optional[float] = None

    class Config:
        from_attributes = True


class QCResult(QCResultInDB):
    """QC结果响应Schema"""

    def __init__(self, **data):
        super().__init__(**data)
        # 计算 eV 单位的 HOMO 和 LUMO 能量
        if self.homo is not None and self.homo_ev is None:
            self.homo_ev = self.homo * 27.2114  # Hartree to eV
        if self.lumo is not None and self.lumo_ev is None:
            self.lumo_ev = self.lumo * 27.2114  # Hartree to eV


# ============================================================================
# Molecule QC Cache Schemas
# ============================================================================

class MoleculeQCCacheBase(BaseModel):
    """分子QC缓存基础Schema"""
    smiles: str
    molecule_name: Optional[str] = None
    basis_set: Optional[str] = None
    functional: Optional[str] = None
    energy_au: Optional[float] = None
    homo_ev: Optional[float] = None
    lumo_ev: Optional[float] = None
    homo_lumo_gap_ev: Optional[float] = None
    esp_min_kcal: Optional[float] = None
    esp_max_kcal: Optional[float] = None
    esp_image_path: Optional[str] = None
    homo_image_path: Optional[str] = None
    lumo_image_path: Optional[str] = None


class MoleculeQCCacheInDB(MoleculeQCCacheBase):
    """数据库中的分子QC缓存"""
    id: int
    preferred_qc_result_id: Optional[int] = None
    calculation_count: int = 1
    created_at: datetime
    updated_at: Optional[datetime] = None

    class Config:
        from_attributes = True


class MoleculeQCCache(MoleculeQCCacheInDB):
    """分子QC缓存响应Schema"""
    pass


# ============================================================================
# Combined Response Schemas
# ============================================================================

class QCJobWithResults(QCJob):
    """带结果的QC任务"""
    results: List[QCResult] = []


class QCJobListResponse(BaseModel):
    """QC任务列表响应"""
    total: int
    jobs: List[QCJobWithResults]


class QCSearchParams(BaseModel):
    """QC搜索参数"""
    smiles: Optional[str] = None
    molecule_name: Optional[str] = None
    basis_set: Optional[str] = None
    lumo_min: Optional[float] = None
    lumo_max: Optional[float] = None
    homo_min: Optional[float] = None
    homo_max: Optional[float] = None
    has_esp: Optional[bool] = None
    skip: int = Field(default=0, ge=0)
    limit: int = Field(default=20, ge=1, le=100)


# ============================================================================
# 重复计算检查
# ============================================================================

class MoleculeCheckRequest(BaseModel):
    """单个分子检查请求"""
    smiles: str
    molecule_name: Optional[str] = None
    functional: str = "B3LYP"
    basis_set: str = "6-31G(d)"
    solvent_model: str = "gas"
    solvent_name: Optional[str] = None
    charge: int = 0
    spin_multiplicity: int = 1


class MoleculeCheckResult(BaseModel):
    """单个分子检查结果"""
    smiles: str
    molecule_name: Optional[str] = None
    has_existing_result: bool = False
    existing_qc_job_id: Optional[int] = None
    existing_result_id: Optional[int] = None
    functional: Optional[str] = None
    basis_set: Optional[str] = None
    solvent_model: Optional[str] = None
    solvent_name: Optional[str] = None
    # 如果有结果，返回基本信息
    energy_au: Optional[float] = None
    homo_ev: Optional[float] = None
    lumo_ev: Optional[float] = None
    homo_lumo_gap_ev: Optional[float] = None
    completed_at: Optional[datetime] = None


class DuplicateCheckRequest(BaseModel):
    """重复计算检查请求"""
    molecules: List[MoleculeCheckRequest]


class DuplicateCheckResponse(BaseModel):
    """重复计算检查响应"""
    total_molecules: int
    existing_count: int
    new_count: int
    results: List[MoleculeCheckResult]


# ============================================================================
# MD Job with QC Options
# ============================================================================

class MDJobQCOptions(BaseModel):
    """MD任务的QC计算选项"""
    enabled: bool = False
    molecules: List[str] = Field(default=[], description="需要计算QC的分子SMILES列表")
    basis_set: str = Field(default="6-31++g(d,p)", description="基组")
    functional: str = Field(default="B3LYP", description="泛函")

