"""
批量导入相关的Schema定义
"""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from datetime import datetime


class BatchElectrolyteRow(BaseModel):
    """批量导入配方的单行数据"""
    # 配方基本信息
    formulation_name: str = Field(..., description="配方名称")
    
    # 阳离子信息（支持多个，用分号分隔）
    cation_names: str = Field(..., description="阳离子名称，多个用分号分隔，如: Li+;Na+")
    cation_smiles: str = Field(..., description="阳离子SMILES，多个用分号分隔")
    cation_numbers: str = Field(..., description="阳离子数量，多个用分号分隔，如: 50;50")
    
    # 阴离子信息（支持多个，用分号分隔）
    anion_names: str = Field(..., description="阴离子名称，多个用分号分隔，如: PF6-;TFSI-")
    anion_smiles: str = Field(..., description="阴离子SMILES，多个用分号分隔")
    anion_numbers: str = Field(..., description="阴离子数量，多个用分号分隔，如: 50;50")
    
    # 溶剂信息（支持多个，用分号分隔）
    solvent_names: Optional[str] = Field(None, description="溶剂名称，多个用分号分隔，如: EC;DMC")
    solvent_smiles: Optional[str] = Field(None, description="溶剂SMILES，多个用分号分隔")
    solvent_numbers: Optional[str] = Field(None, description="溶剂数量，多个用分号分隔，如: 100;100")
    
    # 模拟参数
    temperature: Optional[float] = Field(298.15, description="温度(K)")
    pressure: Optional[float] = Field(1.0, description="压强(atm)")
    box_size: Optional[float] = Field(None, description="盒子尺寸(Å)")
    
    # MD计算参数
    nsteps_npt: Optional[int] = Field(5000000, description="NPT步数")
    nsteps_nvt: Optional[int] = Field(10000000, description="NVT步数")
    timestep: Optional[float] = Field(1.0, description="时间步长(fs)")
    
    # 是否提交MD计算
    submit_md: Optional[bool] = Field(True, description="是否提交MD计算")
    
    # Slurm资源配置
    slurm_partition: Optional[str] = Field("cpu", description="队列/分区")
    slurm_nodes: Optional[int] = Field(1, description="节点数")
    slurm_ntasks: Optional[int] = Field(8, description="任务数")
    slurm_cpus_per_task: Optional[int] = Field(8, description="每任务CPU数")
    slurm_time: Optional[int] = Field(7200, description="最大运行时间(分钟)")


class BatchQCRow(BaseModel):
    """批量导入QC计算的单行数据"""
    # 分子信息
    molecule_name: str = Field(..., description="分子名称")
    smiles: str = Field(..., description="分子SMILES")
    molecule_type: str = Field(..., description="分子类型: solvent/cation/anion/custom")
    
    # QC计算参数
    functional: str = Field("B3LYP", description="泛函")
    basis_set: str = Field("6-31++G(d,p)", description="基组")
    charge: int = Field(0, description="电荷")
    spin_multiplicity: int = Field(1, description="自旋多重度")
    
    # 溶剂模型
    solvent_model: Optional[str] = Field("pcm", description="溶剂模型: gas/pcm/smd")
    solvent_name: Optional[str] = Field("Water", description="溶剂名称")
    
    # 是否提交计算
    submit_qc: Optional[bool] = Field(True, description="是否提交QC计算")
    
    # Slurm资源配置
    slurm_partition: Optional[str] = Field("cpu", description="队列/分区")
    slurm_nodes: Optional[int] = Field(1, description="节点数")
    slurm_ntasks: Optional[int] = Field(1, description="任务数")
    slurm_cpus_per_task: Optional[int] = Field(8, description="每任务CPU数")
    slurm_time: Optional[int] = Field(1440, description="最大运行时间(分钟)")


class BatchImportRequest(BaseModel):
    """批量导入请求"""
    project_id: Optional[int] = Field(None, description="项目ID，如果为空则创建新项目")
    project_name: Optional[str] = Field(None, description="新项目名称（当project_id为空时必填）")
    project_description: Optional[str] = Field(None, description="新项目描述")
    
    # 配方数据
    electrolytes: Optional[List[BatchElectrolyteRow]] = Field(None, description="配方列表")
    
    # QC计算数据
    qc_jobs: Optional[List[BatchQCRow]] = Field(None, description="QC任务列表")


class BatchImportResult(BaseModel):
    """批量导入结果"""
    success: bool
    message: str

    # 项目信息
    project_id: Optional[int] = None
    project_name: Optional[str] = None

    # 导入统计
    total_electrolytes: int = 0
    success_electrolytes: int = 0
    failed_electrolytes: int = 0

    total_md_jobs: int = 0
    success_md_jobs: int = 0
    failed_md_jobs: int = 0

    total_qc_jobs: int = 0
    success_qc_jobs: int = 0
    failed_qc_jobs: int = 0

    # 详细结果
    electrolyte_results: List[Dict[str, Any]] = []  # 成功导入的配方详情
    md_job_results: List[Dict[str, Any]] = []  # 成功创建的MD任务详情
    qc_job_results: List[Dict[str, Any]] = []  # 成功创建的QC任务详情

    # 错误信息
    errors: List[Dict[str, Any]] = []

    # 成功导入的ID列表（用于批量操作）
    success_electrolyte_ids: List[int] = []
    success_md_job_ids: List[int] = []
    success_qc_job_ids: List[int] = []


class TemplateDownloadRequest(BaseModel):
    """模板下载请求"""
    template_type: str = Field(..., description="模板类型: electrolyte/qc/combined")
    include_examples: bool = Field(True, description="是否包含示例数据")

