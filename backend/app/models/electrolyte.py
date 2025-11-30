"""
Electrolyte system model
"""
from sqlalchemy import Column, Integer, String, Float, DateTime, ForeignKey, Index, Boolean
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
from sqlalchemy.dialects.postgresql import JSONB
from app.database import Base


class ElectrolyteSystem(Base):
    """Electrolyte system model"""
    __tablename__ = "electrolyte_systems"

    id = Column(Integer, primary_key=True, index=True)
    project_id = Column(Integer, ForeignKey("projects.id", ondelete="CASCADE"), nullable=False, index=True)
    hash_key = Column(String(64), unique=True, nullable=False, index=True)
    name = Column(String, nullable=False)

    # JSONB fields for flexible composition
    cations = Column(JSONB, nullable=False)
    anions = Column(JSONB, nullable=False)
    solvents = Column(JSONB)

    # Simulation parameters
    temperature = Column(Float, default=298.15)
    pressure = Column(Float, default=1.0)
    density = Column(Float)
    concentration = Column(Float)
    box_size = Column(Float)

    # MD parameters
    nsteps_npt = Column(Integer, default=5000000)
    nsteps_nvt = Column(Integer, default=10000000)
    timestep = Column(Float, default=1.0)
    force_field = Column(String, default="OPLS")

    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)

    # 软删除字段
    is_deleted = Column(Boolean, default=False, index=True)  # 是否已删除
    deleted_at = Column(DateTime(timezone=True))  # 删除时间
    deleted_by = Column(Integer, ForeignKey("users.id", ondelete="SET NULL"))  # 删除操作者
    delete_reason = Column(String(500))  # 删除原因

    # Relationships
    project = relationship("Project", back_populates="electrolyte_systems")
    md_jobs = relationship("MDJob", back_populates="system", cascade="all, delete-orphan")
    deleted_by_user = relationship("User", foreign_keys=[deleted_by])
    
    # Indexes
    __table_args__ = (
        Index('idx_systems_project_id', 'project_id'),
        Index('idx_systems_hash_key', 'hash_key'),
    )
    
    def __repr__(self):
        return f"<ElectrolyteSystem(id={self.id}, name={self.name})>"

