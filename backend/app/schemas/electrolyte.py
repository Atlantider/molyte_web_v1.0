"""
Electrolyte system schemas
"""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from datetime import datetime


class MoleculeSpec(BaseModel):
    """Molecule specification (old format - using number)"""
    name: str
    smiles: str
    number: int = Field(..., gt=0)
    concentration: Optional[float] = None  # Store original concentration for editing
    charge: Optional[int] = None  # Store charge for editing

    class Config:
        extra = "allow"  # Allow extra fields to be stored


class IonSpec(BaseModel):
    """Ion specification (new format - using concentration)"""
    name: str
    charge: int
    concentration: float = Field(..., gt=0)  # mol/L
    count: Optional[int] = None  # Calculated from concentration and volume


class SolventSpec(BaseModel):
    """Solvent specification (new format - using molar ratio)"""
    name: str
    smiles: str
    molar_ratio: float = Field(..., gt=0)  # Relative to first cation
    count: Optional[int] = None  # Calculated from molar ratio


class BoxConfig(BaseModel):
    """Simulation box configuration"""
    type: str = Field(..., pattern="^(cubic|rectangular)$")  # cubic or rectangular
    dimensions: List[float] = Field(..., min_items=1, max_items=3)  # [size] or [x, y, z] in Angstroms


class ElectrolyteBase(BaseModel):
    """Base electrolyte system schema"""
    name: str = Field(..., min_length=1, max_length=200)
    cations: List[MoleculeSpec]
    anions: List[MoleculeSpec]
    solvents: Optional[List[MoleculeSpec]] = None
    temperature: Optional[float] = Field(default=298.15, gt=0)
    pressure: Optional[float] = Field(default=1.0, gt=0)
    density: Optional[float] = Field(default=None, gt=0)
    concentration: Optional[float] = Field(default=None, gt=0)
    box_size: Optional[float] = Field(default=None, gt=0)
    nsteps_npt: Optional[int] = Field(default=5000000, gt=0)
    nsteps_nvt: Optional[int] = Field(default=10000000, gt=0)
    timestep: Optional[float] = Field(default=1.0, gt=0)
    force_field: Optional[str] = Field(default="OPLS")


class ElectrolyteCreate(ElectrolyteBase):
    """Schema for creating an electrolyte system (old format)"""
    project_id: int


class ElectrolyteCreateNew(BaseModel):
    """Schema for creating an electrolyte system (new format with concentration)"""
    project_id: int
    name: Optional[str] = Field(default=None, max_length=200)
    description: Optional[str] = None

    # Ions with concentration
    cations: List[IonSpec]
    anions: List[IonSpec]

    # Solvents with molar ratio (optional)
    solvents: Optional[List[SolventSpec]] = Field(default=None)

    # Box configuration
    box: BoxConfig

    # Simulation parameters
    temperature: Optional[float] = Field(default=298.15, gt=0)
    pressure: Optional[float] = Field(default=1.0, gt=0)
    nsteps_npt: Optional[int] = Field(default=5000000, gt=0)
    nsteps_nvt: Optional[int] = Field(default=10000000, gt=0)
    timestep: Optional[float] = Field(default=1.0, gt=0)
    force_field: Optional[str] = Field(default="OPLS")


class ElectrolyteUpdate(BaseModel):
    """Schema for updating an electrolyte system"""
    name: Optional[str] = Field(None, min_length=1, max_length=200)
    temperature: Optional[float] = Field(None, gt=0)
    pressure: Optional[float] = Field(None, gt=0)
    density: Optional[float] = Field(None, gt=0)
    concentration: Optional[float] = Field(None, gt=0)
    box_size: Optional[float] = Field(None, gt=0)
    nsteps_npt: Optional[int] = Field(None, gt=0)
    nsteps_nvt: Optional[int] = Field(None, gt=0)
    timestep: Optional[float] = Field(None, gt=0)
    force_field: Optional[str] = None


class ElectrolyteInDB(ElectrolyteBase):
    """Electrolyte system schema with database fields"""
    id: int
    project_id: int
    hash_key: str
    created_at: datetime
    
    class Config:
        from_attributes = True


class Electrolyte(ElectrolyteInDB):
    """Electrolyte system response schema"""
    username: Optional[str] = None  # 用户名，用于管理端显示
    user_email: Optional[str] = None  # 用户邮箱，用于管理端显示

