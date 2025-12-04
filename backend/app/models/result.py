"""
Result models (Summary, RDF, MSD, Solvation)
"""
from sqlalchemy import Column, Integer, String, Float, DateTime, ForeignKey, Index
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
from sqlalchemy.dialects.postgresql import JSONB
from app.database import Base


class ResultSummary(Base):
    """Result summary model"""
    __tablename__ = "result_summary"

    id = Column(Integer, primary_key=True, index=True)
    md_job_id = Column(Integer, ForeignKey("md_jobs.id", ondelete="CASCADE"), unique=True, nullable=False, index=True)

    # System information
    total_atoms = Column(Integer)
    total_molecules = Column(Integer)

    # Box dimensions (Å)
    box_x = Column(Float)  # 盒子 X 尺寸
    box_y = Column(Float)  # 盒子 Y 尺寸
    box_z = Column(Float)  # 盒子 Z 尺寸
    initial_box_x = Column(Float)  # 初始盒子 X
    initial_box_y = Column(Float)  # 初始盒子 Y
    initial_box_z = Column(Float)  # 初始盒子 Z

    # Concentration (mol/L)
    concentration = Column(Float)  # 计算浓度
    initial_concentration = Column(Float)  # 初始浓度

    # Final state
    final_density = Column(Float)
    initial_density = Column(Float)  # 初始密度
    final_temperature = Column(Float)
    final_pressure = Column(Float)

    # Energy
    total_energy = Column(Float)
    potential_energy = Column(Float)
    kinetic_energy = Column(Float)

    # System structure (最后一帧 XYZ 内容)
    system_xyz_content = Column(String)  # XYZ 格式的系统结构

    # Molecule structures (分子结构信息 JSON)
    molecule_structures = Column(JSONB)  # [{name, type, pdb_content, smiles, total_charge, atoms}]

    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)

    # Relationships
    md_job = relationship("MDJob", back_populates="result_summary")

    def __repr__(self):
        return f"<ResultSummary(id={self.id}, md_job_id={self.md_job_id})>"


class RDFResult(Base):
    """Radial distribution function result model"""
    __tablename__ = "rdf_results"

    id = Column(Integer, primary_key=True, index=True)
    md_job_id = Column(Integer, ForeignKey("md_jobs.id", ondelete="CASCADE"), nullable=False, index=True)
    center_species = Column(String, nullable=False)
    shell_species = Column(String, nullable=False)

    # RDF data
    r_values = Column(JSONB)
    g_r_values = Column(JSONB)
    coordination_number_values = Column(JSONB)  # 配位数数组

    # Analysis results
    first_peak_position = Column(Float)
    first_peak_height = Column(Float)
    coordination_number = Column(Float)  # 最终配位数（向后兼容）

    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)

    # Relationships
    md_job = relationship("MDJob", back_populates="rdf_results")

    # Indexes
    __table_args__ = (
        Index('idx_rdf_md_job_id', 'md_job_id'),
        Index('idx_rdf_species', 'center_species', 'shell_species'),
    )

    def __repr__(self):
        return f"<RDFResult(id={self.id}, {self.center_species}-{self.shell_species})>"


class MSDResult(Base):
    """Mean square displacement result model"""
    __tablename__ = "msd_results"

    id = Column(Integer, primary_key=True, index=True)
    md_job_id = Column(Integer, ForeignKey("md_jobs.id", ondelete="CASCADE"), nullable=False, index=True)
    species = Column(String, nullable=False, index=True)

    # MSD data (arrays stored as JSON)
    t_values = Column(JSONB)  # Time array (fs) - 使用数据库中的列名
    msd_x_values = Column(JSONB)  # MSD in x direction (Å²)
    msd_y_values = Column(JSONB)  # MSD in y direction (Å²)
    msd_z_values = Column(JSONB)  # MSD in z direction (Å²)
    msd_total_values = Column(JSONB)  # Total MSD (Å²)

    # Labels for plotting
    labels = Column(JSONB)  # {'time': 'fs', 'x': 'Li_x', 'y': 'Li_y', 'z': 'Li_z', 'total': 'Li_total'}

    # Analysis results
    diffusion_coefficient = Column('diffusion_coeff', Float)  # Diffusion coefficient (cm²/s) - 使用数据库中的列名

    # Transport properties (新增)
    ionic_conductivity = Column(Float)  # 离子电导率 (S/cm)
    mobility = Column(Float)  # 离子迁移率 (cm²/(V·s))
    charge = Column(Integer)  # 离子电荷数

    # Legacy columns (keep for compatibility)
    fit_range = Column(JSONB)
    file_path = Column(String)

    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)

    # Relationships
    md_job = relationship("MDJob", back_populates="msd_results")

    # Indexes
    __table_args__ = (
        Index('idx_msd_md_job_id', 'md_job_id'),
        Index('idx_msd_species', 'species'),
    )

    def __repr__(self):
        return f"<MSDResult(id={self.id}, species={self.species})>"


class SolvationStructure(Base):
    """Solvation structure model - 溶剂化结构模型"""
    __tablename__ = "solvation_structures"

    id = Column(Integer, primary_key=True, index=True)
    md_job_id = Column(Integer, ForeignKey("md_jobs.id", ondelete="CASCADE"), nullable=False, index=True)
    center_ion = Column(String, nullable=False)  # 中心离子，如 'Li'
    structure_type = Column(String)  # 结构类型，如 'first_shell', 'typical_cluster'

    # 配位信息
    coordination_num = Column(Integer)  # 配位数
    composition = Column(JSONB)  # 溶剂壳组成，如 {"EC": 3, "DMC": 1, "FSI": 0}
    mol_order = Column(JSONB)  # XYZ中分子的实际顺序，如 [{"mol_name": "EC", "atom_count": 10}, ...]

    # 文件信息
    file_path = Column(String)  # xyz/pdb 文件路径（本地路径，可选）
    xyz_content = Column(String)  # XYZ 文件内容（直接存储，用于 3D 显示）
    snapshot_frame = Column(Integer)  # 快照帧号
    description = Column(String)  # 描述信息

    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)

    # Relationships
    md_job = relationship("MDJob", back_populates="solvation_structures")

    # Indexes
    __table_args__ = (
        Index('idx_solvation_md_job_id', 'md_job_id'),
    )

    def __repr__(self):
        return f"<SolvationStructure(id={self.id}, center_ion={self.center_ion}, cn={self.coordination_num})>"


class DesolvationEnergyResult(Base):
    """Desolvation energy calculation result model - 去溶剂化能计算结果模型"""
    __tablename__ = "desolvation_energy_results"

    id = Column(Integer, primary_key=True, index=True)
    postprocess_job_id = Column(Integer, ForeignKey("postprocess_jobs.id", ondelete="CASCADE"), nullable=False, index=True)
    solvation_structure_id = Column(Integer, ForeignKey("solvation_structures.id", ondelete="CASCADE"), nullable=False, index=True)

    # 计算参数
    method_level = Column(String(50), default="standard")  # fast, standard, accurate
    basis_set = Column(String(50))
    functional = Column(String(50))

    # 完整簇能量
    e_cluster = Column(Float)  # A.U.

    # 每个配体的结果 (JSON)
    per_ligand_results = Column(JSONB)  # [{ligand_id, ligand_type, ligand_label, e_ligand, e_cluster_minus, delta_e}]

    # 按类型汇总 (JSON)
    per_type_summary = Column(JSONB)  # {ligand_type: {avg_delta_e, std_delta_e, count}}

    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)

    # Relationships
    postprocess_job = relationship("PostprocessJob", back_populates="desolvation_energy_result")
    solvation_structure = relationship("SolvationStructure")

    # Indexes
    __table_args__ = (
        Index('idx_desolvation_postprocess_job_id', 'postprocess_job_id'),
        Index('idx_desolvation_solvation_structure_id', 'solvation_structure_id'),
    )

    def __repr__(self):
        return f"<DesolvationEnergyResult(id={self.id}, method={self.method_level}, e_cluster={self.e_cluster})>"


class SystemStructure(Base):
    """System structure model - 系统结构模型（用于3D可视化）"""
    __tablename__ = "system_structures"

    id = Column(Integer, primary_key=True, index=True)
    md_job_id = Column(Integer, ForeignKey("md_jobs.id", ondelete="CASCADE"), unique=True, nullable=False, index=True)

    # 帧信息
    frame_index = Column(Integer)  # 帧索引
    total_frames = Column(Integer)  # 总帧数
    atom_count = Column(Integer)  # 原子数

    # 盒子信息
    box = Column(JSONB)  # [lx, ly, lz]

    # XYZ 内容
    xyz_content = Column(String)  # XYZ 格式的系统结构（用于3D显示）

    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)

    # Relationships
    md_job = relationship("MDJob", back_populates="system_structure", uselist=False)

    # Indexes
    __table_args__ = (
        Index('idx_system_structure_md_job_id', 'md_job_id'),
    )

    def __repr__(self):
        return f"<SystemStructure(id={self.id}, md_job_id={self.md_job_id}, atoms={self.atom_count})>"

