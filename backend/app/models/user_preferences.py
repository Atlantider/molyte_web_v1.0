"""
用户偏好设置模型
存储用户自定义的常用溶剂组合、离子组合等
"""
from sqlalchemy import Column, Integer, String, Text, ForeignKey, TIMESTAMP
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship
from datetime import datetime

from app.database import Base


class UserSolventCombination(Base):
    """用户自定义溶剂组合"""
    __tablename__ = "user_solvent_combinations"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id", ondelete="CASCADE"), nullable=False)
    name = Column(String(100), nullable=False)  # 组合名称
    description = Column(Text)  # 描述
    solvents = Column(JSONB, nullable=False)  # 溶剂列表 [{"name": "EC", "smiles": "...", "molar_ratio": 1.0}]
    created_at = Column(TIMESTAMP(timezone=True), default=datetime.utcnow)
    updated_at = Column(TIMESTAMP(timezone=True), default=datetime.utcnow, onupdate=datetime.utcnow)

    # 关系
    user = relationship("User", back_populates="solvent_combinations")

    def __repr__(self):
        return f"<UserSolventCombination(id={self.id}, user_id={self.user_id}, name='{self.name}')>"


class UserIonCombination(Base):
    """用户自定义离子组合"""
    __tablename__ = "user_ion_combinations"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id", ondelete="CASCADE"), nullable=False)
    name = Column(String(100), nullable=False)  # 组合名称
    description = Column(Text)  # 描述
    cations = Column(JSONB, nullable=False)  # 阳离子列表 [{"name": "Li", "charge": 1, "concentration": 1.0}]
    anions = Column(JSONB, nullable=False)  # 阴离子列表
    created_at = Column(TIMESTAMP(timezone=True), default=datetime.utcnow)
    updated_at = Column(TIMESTAMP(timezone=True), default=datetime.utcnow, onupdate=datetime.utcnow)

    # 关系
    user = relationship("User", back_populates="ion_combinations")

    def __repr__(self):
        return f"<UserIonCombination(id={self.id}, user_id={self.user_id}, name='{self.name}')>"

