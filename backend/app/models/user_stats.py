"""
User usage statistics and audit log models
"""
from sqlalchemy import Column, Integer, String, Float, DateTime, Date, ForeignKey, JSON, UniqueConstraint
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
from app.database import Base


class UserUsageStats(Base):
    """User daily usage statistics"""
    __tablename__ = "user_usage_stats"
    
    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id", ondelete="CASCADE"), nullable=False)
    date = Column(Date, nullable=False, index=True)
    
    # Job statistics
    jobs_submitted = Column(Integer, default=0, nullable=False)
    jobs_completed = Column(Integer, default=0, nullable=False)
    jobs_failed = Column(Integer, default=0, nullable=False)
    jobs_cancelled = Column(Integer, default=0, nullable=False)
    
    # Resource usage
    cpu_hours_used = Column(Float, default=0.0, nullable=False)
    storage_used_gb = Column(Float, default=0.0, nullable=False)
    
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)
    
    # Relationships
    user = relationship("User", back_populates="usage_stats")
    
    # Unique constraint: one record per user per day
    __table_args__ = (
        UniqueConstraint('user_id', 'date', name='uq_user_date'),
    )
    
    def __repr__(self):
        return f"<UserUsageStats(user_id={self.user_id}, date={self.date}, jobs={self.jobs_submitted})>"


class AuditLog(Base):
    """Audit log for sensitive operations"""
    __tablename__ = "audit_logs"
    
    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id", ondelete="SET NULL"), nullable=True)
    
    # Action details
    action = Column(String(100), nullable=False, index=True)  # e.g., "update_user_quota", "cancel_job"
    resource_type = Column(String(50), nullable=True)  # e.g., "user", "job", "project"
    resource_id = Column(Integer, nullable=True)
    
    # Additional details
    details = Column(JSON, nullable=True)  # Store additional context as JSON
    ip_address = Column(String(50), nullable=True)
    user_agent = Column(String(500), nullable=True)
    
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False, index=True)
    
    def __repr__(self):
        return f"<AuditLog(id={self.id}, action={self.action}, user_id={self.user_id})>"

