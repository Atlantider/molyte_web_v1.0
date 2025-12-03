"""
Desolvation energy calculation schemas
"""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from datetime import datetime


class DesolvationJobCreate(BaseModel):
    """创建去溶剂化能任务"""
    md_job_id: int = Field(..., description="MD job ID")
    solvation_structure_id: int = Field(..., description="Solvation structure ID")
    method_level: str = Field(default="standard", description="Calculation method level: fast (6-31G(d)/B3LYP), standard (6-31++G(d,p)/B3LYP), accurate (6-311++G(2d,2p)/wB97XD)")


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
    created_at: datetime
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None
    elapsed_seconds: Optional[float] = None
    error_message: Optional[str] = None
    result: Optional[DesolvationEnergyResultSchema] = None

    class Config:
        from_attributes = True

