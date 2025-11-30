"""
Electrolyte system API endpoints
"""
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
from typing import List, Dict
import os
from pathlib import Path
from app.database import get_db
from app.models.user import User, UserRole
from app.models.project import Project
from app.models.electrolyte import ElectrolyteSystem
from app.schemas.electrolyte import (
    Electrolyte as ElectrolyteSchema,
    ElectrolyteCreate,
    ElectrolyteCreateNew,
    ElectrolyteUpdate
)
from app.dependencies import get_current_active_user
from app.utils.hash import calculate_system_hash
from app.utils.ion_parser import scan_available_ions, get_cations_and_anions
from app.utils.smiles_validator import validate_smiles
from app.utils.electrolyte_converter import convert_new_to_old_format, convert_old_to_new_format
from app.core.logger import logger
from pydantic import BaseModel

router = APIRouter()

# Path to initial salts directory
SALTS_DIR = Path("/public/home/xiaoji/molyte_web/data/initial_salts")

# Cache for ion information (loaded on startup)
_ions_cache = None


def get_ions_info():
    """Get cached ion information, or scan if not cached"""
    global _ions_cache
    if _ions_cache is None:
        _ions_cache = scan_available_ions(SALTS_DIR)
    return _ions_cache


@router.get("/available-ions")
def get_available_ions_endpoint(current_user: User = Depends(get_current_active_user)):
    """
    Get list of available cations and anions with their charges
    Automatically parsed from .lt files in the inital_salts directory

    Args:
        current_user: Current authenticated user

    Returns:
        Dict with 'cations' and 'anions' lists, each containing {name, charge}
    """
    ions_info = get_ions_info()
    cations, anions = get_cations_and_anions(ions_info)

    logger.info(f"Available ions requested by {current_user.username}: {len(cations)} cations, {len(anions)} anions")

    return {
        "cations": cations,
        "anions": anions
    }


@router.post("/refresh-ions")
def refresh_ions_cache(current_user: User = Depends(get_current_active_user)):
    """
    Refresh the ion cache by re-scanning the salts directory
    Useful after adding new .lt files

    Args:
        current_user: Current authenticated user

    Returns:
        Updated ion lists
    """
    global _ions_cache
    _ions_cache = None  # Clear cache

    ions_info = get_ions_info()  # Re-scan
    cations, anions = get_cations_and_anions(ions_info)

    logger.info(f"Ion cache refreshed by {current_user.username}: {len(cations)} cations, {len(anions)} anions")

    return {
        "message": "Ion cache refreshed successfully",
        "cations": cations,
        "anions": anions
    }


class ValidateSmilesRequest(BaseModel):
    smiles: str


@router.post("/validate-smiles")
def validate_smiles_endpoint(
    request: ValidateSmilesRequest,
    current_user: User = Depends(get_current_active_user)
):
    """
    Validate a SMILES string and get molecule information

    Args:
        request: Request containing SMILES string
        current_user: Current authenticated user

    Returns:
        Validation result with molecule information
    """
    result = validate_smiles(request.smiles)
    logger.info(f"SMILES validation by {current_user.username}: {request.smiles} -> valid={result.get('valid', False)}")
    return result


@router.post("/new", response_model=ElectrolyteSchema, status_code=status.HTTP_201_CREATED)
def create_electrolyte_new_format(
    electrolyte_data: ElectrolyteCreateNew,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Create a new electrolyte system using new format (concentration-based)

    Args:
        electrolyte_data: Electrolyte creation data (new format)
        db: Database session
        current_user: Current authenticated user

    Returns:
        Electrolyte: Created electrolyte system

    Raises:
        HTTPException: If project not found or validation fails
    """
    # Check if project exists and user has permission
    project = db.query(Project).filter(Project.id == electrolyte_data.project_id).first()
    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Project not found"
        )

    if project.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="You don't have permission to create electrolytes in this project"
        )

    # Convert new format to old format
    try:
        old_format_data = convert_new_to_old_format(electrolyte_data)
    except Exception as e:
        logger.error(f"Error converting electrolyte data: {e}")
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Error converting electrolyte data: {str(e)}"
        )

    # Generate electrolyte name with EL-YYYYMMDD-全局序号 prefix
    from datetime import date, datetime
    from sqlalchemy import func

    today = date.today()
    date_str = today.strftime('%Y%m%d')
    today_start = datetime.combine(today, datetime.min.time())
    today_end = datetime.combine(today, datetime.max.time())

    # Count ALL electrolytes created today (global count, not per-user)
    count = db.query(func.count(ElectrolyteSystem.id)).filter(
        ElectrolyteSystem.created_at >= today_start,
        ElectrolyteSystem.created_at <= today_end
    ).scalar()

    # Next sequential number (starting from 1)
    seq_number = count + 1

    # Format: EL-20251119-0001-description
    description = old_format_data["name"]
    electrolyte_name = f"EL-{date_str}-{seq_number:04d}-{description}"
    old_format_data["name"] = electrolyte_name

    logger.info(f"Generated electrolyte name (global): {electrolyte_name}")

    # Create ElectrolyteCreate object from converted data
    create_data = ElectrolyteCreate(**old_format_data)

    # Calculate hash - convert Pydantic objects to dicts
    hash_key = calculate_system_hash(
        cations=[c.dict() if hasattr(c, 'dict') else c for c in create_data.cations],
        anions=[a.dict() if hasattr(a, 'dict') else a for a in create_data.anions],
        solvents=[s.dict() if hasattr(s, 'dict') else s for s in create_data.solvents] if create_data.solvents else [],
        temperature=create_data.temperature,
        pressure=create_data.pressure
    )

    # Check for duplicate
    existing = db.query(ElectrolyteSystem).filter(
        ElectrolyteSystem.hash_key == hash_key
    ).first()

    if existing:
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail=f"An identical electrolyte system already exists (ID: {existing.id})"
        )

    # Create new electrolyte system
    db_electrolyte = ElectrolyteSystem(
        **create_data.dict(),
        hash_key=hash_key
    )

    db.add(db_electrolyte)
    db.commit()
    db.refresh(db_electrolyte)

    logger.info(f"Created electrolyte system (new format) {db_electrolyte.id} by {current_user.username}")

    return db_electrolyte


@router.post("/", response_model=ElectrolyteSchema, status_code=status.HTTP_201_CREATED)
def create_electrolyte(
    electrolyte_data: ElectrolyteCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Create a new electrolyte system

    Args:
        electrolyte_data: Electrolyte creation data
        db: Database session
        current_user: Current authenticated user

    Returns:
        Electrolyte: Created electrolyte system

    Raises:
        HTTPException: If project not found or duplicate system
    """
    # Check if project exists and user has permission
    project = db.query(Project).filter(Project.id == electrolyte_data.project_id).first()
    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Project not found"
        )

    if project.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not enough permissions"
        )

    # Calculate hash key
    hash_key = calculate_system_hash(
        cations=[c.model_dump() for c in electrolyte_data.cations],
        anions=[a.model_dump() for a in electrolyte_data.anions],
        solvents=[s.model_dump() for s in electrolyte_data.solvents] if electrolyte_data.solvents else [],
        temperature=electrolyte_data.temperature,
        pressure=electrolyte_data.pressure
    )

    # Check for duplicate
    existing = db.query(ElectrolyteSystem).filter(
        ElectrolyteSystem.hash_key == hash_key
    ).first()
    if existing:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Duplicate electrolyte system (ID: {existing.id})"
        )

    # Create electrolyte system
    db_electrolyte = ElectrolyteSystem(
        project_id=electrolyte_data.project_id,
        hash_key=hash_key,
        name=electrolyte_data.name,
        cations=[c.model_dump() for c in electrolyte_data.cations],
        anions=[a.model_dump() for a in electrolyte_data.anions],
        solvents=[s.model_dump() for s in electrolyte_data.solvents] if electrolyte_data.solvents else None,
        temperature=electrolyte_data.temperature,
        pressure=electrolyte_data.pressure,
        density=electrolyte_data.density,
        concentration=electrolyte_data.concentration,
        box_size=electrolyte_data.box_size,
        nsteps_npt=electrolyte_data.nsteps_npt,
        nsteps_nvt=electrolyte_data.nsteps_nvt,
        timestep=electrolyte_data.timestep,
        force_field=electrolyte_data.force_field
    )

    db.add(db_electrolyte)
    db.commit()
    db.refresh(db_electrolyte)

    logger.info(f"Electrolyte system created: {db_electrolyte.name}")
    return db_electrolyte


@router.get("/", response_model=List[ElectrolyteSchema])
def list_electrolytes(
    project_id: int = None,
    skip: int = 0,
    limit: int = 100,
    include_deleted: bool = False,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    List electrolyte systems

    Args:
        project_id: Filter by project ID
        skip: Number of records to skip
        limit: Maximum number of records to return
        include_deleted: 是否包含已删除的配方（仅管理员）
        db: Database session
        current_user: Current authenticated user

    Returns:
        List[Electrolyte]: List of electrolyte systems
    """
    from sqlalchemy import or_

    query = db.query(ElectrolyteSystem)

    # Filter by project if specified
    if project_id:
        project = db.query(Project).filter(Project.id == project_id).first()
        if not project:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Project not found"
            )

        if project.user_id != current_user.id and current_user.role != UserRole.ADMIN:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Not enough permissions"
            )

        query = query.filter(ElectrolyteSystem.project_id == project_id)
    else:
        # Show only user's electrolytes unless admin
        if current_user.role != UserRole.ADMIN:
            user_project_ids = [p.id for p in current_user.projects]
            query = query.filter(ElectrolyteSystem.project_id.in_(user_project_ids))

    # 配方使用硬删除，不需要过滤 is_deleted
    # 已删除的配方会从数据库中移除

    # Order by created_at descending (newest first)
    query = query.order_by(ElectrolyteSystem.created_at.desc())

    electrolytes = query.offset(skip).limit(limit).all()
    return electrolytes


@router.get("/{electrolyte_id}", response_model=ElectrolyteSchema)
def get_electrolyte(
    electrolyte_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Get electrolyte system by ID

    Args:
        electrolyte_id: Electrolyte system ID
        db: Database session
        current_user: Current authenticated user

    Returns:
        Electrolyte: Electrolyte system data

    Raises:
        HTTPException: If not found or no permission
    """
    electrolyte = db.query(ElectrolyteSystem).filter(
        ElectrolyteSystem.id == electrolyte_id
    ).first()

    if not electrolyte:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Electrolyte system not found"
        )

    # Check permission
    if electrolyte.project.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not enough permissions"
        )

    return electrolyte



@router.get("/{electrolyte_id}/editable", response_model=dict)
def get_electrolyte_editable(
    electrolyte_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Get electrolyte system in editable format (new format with concentrations)

    Args:
        electrolyte_id: Electrolyte system ID
        db: Database session
        current_user: Current authenticated user

    Returns:
        dict: Electrolyte data in new format for editing
    """
    electrolyte = db.query(ElectrolyteSystem).filter(
        ElectrolyteSystem.id == electrolyte_id
    ).first()

    if not electrolyte:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Electrolyte system not found"
        )

    # Check permission
    if electrolyte.project.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not enough permissions"
        )

    # Convert to editable format
    editable_data = convert_old_to_new_format(electrolyte)

    logger.info(f"Retrieved editable format for electrolyte {electrolyte_id} by {current_user.username}")

    return editable_data


@router.put("/{electrolyte_id}", response_model=ElectrolyteSchema)
def update_electrolyte(
    electrolyte_id: int,
    electrolyte_data: ElectrolyteCreateNew,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Update an existing electrolyte system using new format

    Args:
        electrolyte_id: Electrolyte system ID
        electrolyte_data: Updated electrolyte data in new format
        db: Database session
        current_user: Current authenticated user

    Returns:
        Electrolyte: Updated electrolyte system
    """
    # Get existing electrolyte
    electrolyte = db.query(ElectrolyteSystem).filter(
        ElectrolyteSystem.id == electrolyte_id
    ).first()

    if not electrolyte:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Electrolyte system not found"
        )

    # Check permission
    if electrolyte.project.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not enough permissions"
        )

    # Convert new format to old format
    try:
        old_format_data = convert_new_to_old_format(electrolyte_data)
    except Exception as e:
        logger.error(f"Error converting electrolyte data: {e}")
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Error converting electrolyte data: {str(e)}"
        )

    # Preserve the EL-YYYYMMDD-序号 prefix from the original name
    # Extract prefix from original name (e.g., "EL-20251119-0001-")
    import re
    original_name = electrolyte.name
    prefix_match = re.match(r'^(EL-\d{8}-\d{4})-', original_name)

    if prefix_match:
        # Keep the original prefix, update only the description part
        prefix = prefix_match.group(1)
        new_description = old_format_data["name"]
        electrolyte_name = f"{prefix}-{new_description}"
        logger.info(f"Updating electrolyte name: {original_name} -> {electrolyte_name}")
    else:
        # Old format name without prefix, keep as is or generate new prefix
        electrolyte_name = old_format_data["name"]
        logger.info(f"Electrolyte has old format name, keeping: {electrolyte_name}")

    old_format_data["name"] = electrolyte_name

    # Create ElectrolyteCreate object from converted data
    create_data = ElectrolyteCreate(**old_format_data)

    # Calculate new hash
    new_hash_key = calculate_system_hash(
        cations=[c.dict() if hasattr(c, 'dict') else c for c in create_data.cations],
        anions=[a.dict() if hasattr(a, 'dict') else a for a in create_data.anions],
        solvents=[s.dict() if hasattr(s, 'dict') else s for s in create_data.solvents] if create_data.solvents else [],
        temperature=create_data.temperature,
        pressure=create_data.pressure
    )

    # Check if new hash conflicts with another system
    if new_hash_key != electrolyte.hash_key:
        existing = db.query(ElectrolyteSystem).filter(
            ElectrolyteSystem.hash_key == new_hash_key,
            ElectrolyteSystem.id != electrolyte_id
        ).first()

        if existing:
            raise HTTPException(
                status_code=status.HTTP_409_CONFLICT,
                detail=f"An identical electrolyte system already exists (ID: {existing.id})"
            )

    # Update electrolyte fields
    electrolyte.project_id = create_data.project_id
    electrolyte.name = create_data.name
    electrolyte.hash_key = new_hash_key
    # Convert Pydantic objects to dicts for JSONB storage
    electrolyte.cations = [c.dict() if hasattr(c, 'dict') else c for c in create_data.cations]
    electrolyte.anions = [a.dict() if hasattr(a, 'dict') else a for a in create_data.anions]
    electrolyte.solvents = [s.dict() if hasattr(s, 'dict') else s for s in create_data.solvents]
    electrolyte.temperature = create_data.temperature
    electrolyte.pressure = create_data.pressure
    electrolyte.box_size = create_data.box_size
    electrolyte.nsteps_npt = create_data.nsteps_npt
    electrolyte.nsteps_nvt = create_data.nsteps_nvt
    electrolyte.timestep = create_data.timestep
    electrolyte.force_field = create_data.force_field

    db.commit()
    db.refresh(electrolyte)

    logger.info(f"Updated electrolyte system {electrolyte_id} by {current_user.username}")

    return electrolyte



@router.delete("/{electrolyte_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_electrolyte(
    electrolyte_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Delete an electrolyte system (软删除 - 保留数据)

    Args:
        electrolyte_id: Electrolyte system ID
        db: Database session
        current_user: Current authenticated user

    Returns:
        None (204 No Content)
    """
    from datetime import datetime

    # Get existing electrolyte
    electrolyte = db.query(ElectrolyteSystem).filter(
        ElectrolyteSystem.id == electrolyte_id
    ).first()

    if not electrolyte:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Electrolyte system not found"
        )

    # Check permission
    if electrolyte.project.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not enough permissions"
        )

    # 检查是否有关联的任务
    has_jobs = db.query(MDJob).filter(
        MDJob.system_id == electrolyte_id,
        MDJob.is_deleted == False
    ).first()

    if has_jobs:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="无法删除配方：该配方有关联的MD任务。请先删除所有关联任务，或者保留配方。"
        )

    # 硬删除电解质配方（配方本身不重要，可以重新创建）
    db.delete(electrolyte)
    db.commit()

    logger.info(f"Deleted electrolyte system {electrolyte_id} by {current_user.username}")

    return None


class BatchUpdateProjectRequest(BaseModel):
    ids: List[int]
    project_id: int


@router.delete("/batch/delete")
def batch_delete_electrolytes(
    ids: List[int],
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    批量删除电解质配方（软删除 - 保留数据）

    只能删除用户自己的配方，或管理员可以删除所有配方
    """
    from datetime import datetime

    deleted_count = 0
    failed_ids = []

    for electrolyte_id in ids:
        electrolyte = db.query(ElectrolyteSystem).filter(ElectrolyteSystem.id == electrolyte_id).first()
        if not electrolyte:
            failed_ids.append(electrolyte_id)
            continue

        # 检查权限 - 通过project关系获取用户ID
        if electrolyte.project.user_id != current_user.id and current_user.role != UserRole.ADMIN:
            failed_ids.append(electrolyte_id)
            continue

        # 检查是否有关联的任务
        has_jobs = db.query(MDJob).filter(
            MDJob.system_id == electrolyte_id,
            MDJob.is_deleted == False
        ).first()

        if has_jobs:
            failed_ids.append(electrolyte_id)
            continue

        # 硬删除
        db.delete(electrolyte)
        deleted_count += 1

    db.commit()

    logger.info(f"Batch deleted {deleted_count} electrolyte systems by {current_user.username}")

    return {
        "deleted_count": deleted_count,
        "failed_ids": failed_ids,
        "message": f"成功删除 {deleted_count} 个配方"
    }


@router.put("/batch/project")
def batch_update_project(
    request: BatchUpdateProjectRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    批量更改电解质配方的项目归属

    只能更改用户自己的配方到自己的项目
    """
    from app.models.project import Project

    # 验证目标项目存在且属于当前用户
    target_project = db.query(Project).filter(Project.id == request.project_id).first()
    if not target_project:
        raise HTTPException(status_code=404, detail="目标项目不存在")

    if target_project.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="无权操作目标项目")

    updated_count = 0
    failed_ids = []

    for electrolyte_id in request.ids:
        electrolyte = db.query(ElectrolyteSystem).filter(ElectrolyteSystem.id == electrolyte_id).first()
        if not electrolyte:
            failed_ids.append(electrolyte_id)
            continue

        # 检查权限
        if electrolyte.project.user_id != current_user.id and current_user.role != UserRole.ADMIN:
            failed_ids.append(electrolyte_id)
            continue

        # 更新项目归属
        electrolyte.project_id = request.project_id
        updated_count += 1

    db.commit()

    logger.info(f"Batch updated {updated_count} electrolyte systems to project {request.project_id} by {current_user.username}")

    return {
        "updated_count": updated_count,
        "failed_ids": failed_ids,
        "message": f"成功移动 {updated_count} 个配方到项目 '{target_project.name}'"
    }
