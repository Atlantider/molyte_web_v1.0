"""
Audit logging utilities
"""
from typing import Optional, Dict, Any
from datetime import datetime
from sqlalchemy.orm import Session
from fastapi import Request
from app.models.user_stats import AuditLog
from app.models.user import User


def create_audit_log(
    db: Session,
    user: Optional[User],
    action: str,
    resource_type: Optional[str] = None,
    resource_id: Optional[int] = None,
    details: Optional[Dict[str, Any]] = None,
    request: Optional[Request] = None
) -> AuditLog:
    """
    Create an audit log entry
    
    Args:
        db: Database session
        user: User performing the action (None for system actions)
        action: Action name (e.g., "update_user_quota", "cancel_job")
        resource_type: Type of resource (e.g., "user", "job", "project")
        resource_id: ID of the resource
        details: Additional details as JSON
        request: FastAPI request object (to extract IP and user agent)
        
    Returns:
        AuditLog: Created audit log entry
    """
    ip_address = None
    user_agent = None
    
    if request:
        # Extract IP address
        ip_address = request.client.host if request.client else None
        
        # Extract user agent
        user_agent = request.headers.get("user-agent")
    
    log = AuditLog(
        user_id=user.id if user else None,
        action=action,
        resource_type=resource_type,
        resource_id=resource_id,
        details=details,
        ip_address=ip_address,
        user_agent=user_agent
    )
    
    db.add(log)
    db.commit()
    db.refresh(log)
    
    return log


def log_user_update(
    db: Session,
    admin: User,
    target_user: User,
    changes: Dict[str, Any],
    request: Optional[Request] = None
):
    """
    Log user update action
    
    Args:
        db: Database session
        admin: Admin user performing the update
        target_user: User being updated
        changes: Dictionary of changed fields
        request: FastAPI request object
    """
    return create_audit_log(
        db=db,
        user=admin,
        action="update_user",
        resource_type="user",
        resource_id=target_user.id,
        details={
            "target_username": target_user.username,
            "changes": changes
        },
        request=request
    )


def log_quota_update(
    db: Session,
    admin: User,
    target_user: User,
    old_quotas: Dict[str, Any],
    new_quotas: Dict[str, Any],
    request: Optional[Request] = None
):
    """
    Log quota update action
    
    Args:
        db: Database session
        admin: Admin user performing the update
        target_user: User whose quota is being updated
        old_quotas: Old quota values
        new_quotas: New quota values
        request: FastAPI request object
    """
    return create_audit_log(
        db=db,
        user=admin,
        action="update_user_quota",
        resource_type="user",
        resource_id=target_user.id,
        details={
            "target_username": target_user.username,
            "old_quotas": old_quotas,
            "new_quotas": new_quotas
        },
        request=request
    )


def log_job_cancel(
    db: Session,
    user: User,
    job_id: int,
    reason: Optional[str] = None,
    request: Optional[Request] = None
):
    """
    Log job cancellation action
    
    Args:
        db: Database session
        user: User cancelling the job
        job_id: ID of the job being cancelled
        reason: Reason for cancellation
        request: FastAPI request object
    """
    return create_audit_log(
        db=db,
        user=user,
        action="cancel_job",
        resource_type="job",
        resource_id=job_id,
        details={
            "reason": reason
        },
        request=request
    )


def log_user_status_change(
    db: Session,
    admin: User,
    target_user: User,
    old_status: bool,
    new_status: bool,
    request: Optional[Request] = None
):
    """
    Log user status change (enable/disable)
    
    Args:
        db: Database session
        admin: Admin user performing the change
        target_user: User whose status is being changed
        old_status: Old is_active status
        new_status: New is_active status
        request: FastAPI request object
    """
    action = "enable_user" if new_status else "disable_user"
    
    return create_audit_log(
        db=db,
        user=admin,
        action=action,
        resource_type="user",
        resource_id=target_user.id,
        details={
            "target_username": target_user.username,
            "old_status": old_status,
            "new_status": new_status
        },
        request=request
    )

