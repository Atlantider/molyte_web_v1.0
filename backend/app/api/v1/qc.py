"""
Quantum Chemistry (QC) API routes
量子化学计算相关的API接口
"""
import logging
import os
from typing import List, Optional
from datetime import datetime

from fastapi import APIRouter, Depends, HTTPException, Query
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session, selectinload
from sqlalchemy import desc, or_

from app.database import get_db
from app.models.user import User, UserRole
from app.models.qc import QCJob, QCResult, MoleculeQCCache, QCJobStatus as QCJobStatusModel
from app.schemas.qc import (
    QCJobCreate,
    QCJobBatchCreate,
    QCJobUpdate,
    QCJobEdit,
    QCJobRecalculate,
    QCJob as QCJobSchema,
    QCJobWithResults,
    QCJobListResponse,
    QCResult as QCResultSchema,
    MoleculeQCCache as MoleculeQCCacheSchema,
    QCSearchParams,
    QCAccuracyLevel,
    SolventModel,
    GAUSSIAN_SOLVENTS,
    QC_ACCURACY_PRESETS,
    BasisSet,
    Functional,
    DuplicateCheckRequest,
    DuplicateCheckResponse,
    MoleculeCheckResult,
)
from app.dependencies import get_current_active_user

logger = logging.getLogger(__name__)
router = APIRouter()


# ============================================================================
# 重复计算检查 (全局共享)
# ============================================================================

@router.post("/check-duplicates", response_model=DuplicateCheckResponse)
def check_duplicate_calculations(
    request: DuplicateCheckRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    检查分子是否已有相同参数的QC计算结果（全局共享）

    所有QC计算结果是全局共享的，不限于单个用户。
    如果找到相同参数的已完成计算，直接使用已有结果，避免重复计算。

    检查条件（所有参数必须完全匹配）：
    - SMILES
    - 泛函 (functional)
    - 基组 (basis_set)
    - 溶剂模型 (solvent_model)
    - 隐式溶剂名称 (solvent_name)（如果使用隐式溶剂）
    - 电荷 (charge)
    - 自旋多重度 (spin_multiplicity)
    """
    results = []
    existing_count = 0

    for mol in request.molecules:
        result = MoleculeCheckResult(
            smiles=mol.smiles,
            molecule_name=mol.molecule_name,
            has_existing_result=False
        )

        # 构建查询条件 - 查找所有用户的已完成QC计算
        query = db.query(QCJob).filter(
            QCJob.smiles == mol.smiles,
            QCJob.functional == mol.functional,
            QCJob.basis_set == mol.basis_set,
            QCJob.charge == mol.charge,
            QCJob.spin_multiplicity == mol.spin_multiplicity,
            QCJob.status == QCJobStatusModel.COMPLETED
        )

        # 匹配溶剂配置
        # 使用JSONB查询匹配溶剂模型和溶剂名称
        if mol.solvent_model == 'gas':
            # 气相：溶剂模型为gas或者没有solvent_config
            query = query.filter(
                or_(
                    QCJob.config['solvent_config']['model'].astext == 'gas',
                    QCJob.config['solvent_config'].is_(None),
                    ~QCJob.config.has_key('solvent_config')
                )
            )
        else:
            # PCM/SMD：需要匹配溶剂模型和溶剂名称
            query = query.filter(
                QCJob.config['solvent_config']['model'].astext == mol.solvent_model
            )
            if mol.solvent_name:
                query = query.filter(
                    QCJob.config['solvent_config']['solvent_name'].astext == mol.solvent_name
                )

        # 查找最新的已完成计算
        existing_job = query.order_by(desc(QCJob.finished_at)).first()

        if existing_job:
            # 找到已有结果
            result.has_existing_result = True
            result.existing_qc_job_id = existing_job.id
            result.functional = existing_job.functional
            result.basis_set = existing_job.basis_set
            result.solvent_model = existing_job.config.get('solvent_config', {}).get('model', 'gas') if existing_job.config else 'gas'
            result.solvent_name = existing_job.config.get('solvent_config', {}).get('solvent_name') if existing_job.config else None
            result.completed_at = existing_job.finished_at

            # 获取计算结果
            if existing_job.results:
                qc_result = existing_job.results[0]  # 取第一个结果
                result.existing_result_id = qc_result.id
                result.energy_au = qc_result.energy_au
                result.homo_ev = qc_result.homo * 27.2114 if qc_result.homo else None
                result.lumo_ev = qc_result.lumo * 27.2114 if qc_result.lumo else None
                result.homo_lumo_gap_ev = qc_result.homo_lumo_gap

            existing_count += 1
            logger.info(f"Found existing QC result for {mol.smiles[:30]}... (job_id={existing_job.id})")

        results.append(result)

    return DuplicateCheckResponse(
        total_molecules=len(request.molecules),
        existing_count=existing_count,
        new_count=len(request.molecules) - existing_count,
        results=results
    )


# ============================================================================
# QC Job CRUD Operations
# ============================================================================

@router.post("/jobs", response_model=QCJobSchema)
def create_qc_job(
    job_data: QCJobCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    创建QC计算任务

    Args:
        job_data: QC任务创建数据

    Returns:
        创建的QC任务
    """
    # 验证SMILES是否有效
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(job_data.smiles)
        if mol is None:
            raise HTTPException(
                status_code=400,
                detail=f"无效的SMILES: {job_data.smiles}。请检查分子结构是否正确。"
            )
        # 尝试生成3D坐标以验证分子可以被处理
        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result == -1:
            # 尝试使用随机坐标方法
            result = AllChem.EmbedMolecule(mol, useRandomCoords=True, maxAttempts=100, randomSeed=42)
        if result == -1:
            # 某些特殊分子（如PF6-）的UFF力场不支持，需要手动处理
            # 这里只发出警告，不阻止创建任务
            # 实际的3D坐标会在任务执行时通过其他方法生成
            logger.warning(f"无法使用RDKit自动生成3D坐标: {job_data.smiles}，将在任务执行时尝试其他方法")
    except ImportError:
        logger.warning("RDKit not available, skipping SMILES validation")
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=400,
            detail=f"SMILES验证失败: {str(e)}"
        )

    # 提取溶剂模型和溶剂名称
    solvent_model = 'gas'
    solvent_name = None
    if job_data.solvent_config:
        solvent_model = job_data.solvent_config.model or 'gas'
        solvent_name = job_data.solvent_config.solvent_name

    # 检查是否已存在完全相同参数的任务（查重）
    # 检查所有未删除的任务，不仅仅是已完成的
    duplicate_query = db.query(QCJob).filter(
        QCJob.smiles == job_data.smiles,
        QCJob.functional == job_data.functional,
        QCJob.basis_set == job_data.basis_set,
        QCJob.charge == job_data.charge,
        QCJob.spin_multiplicity == job_data.spin_multiplicity,
        QCJob.is_deleted == False
    )

    # 匹配溶剂配置
    if solvent_model == 'gas':
        duplicate_query = duplicate_query.filter(
            or_(
                QCJob.solvent_model == 'gas',
                QCJob.solvent_model.is_(None)
            )
        )
    elif solvent_model == 'custom':
        # 自定义溶剂：需要匹配所有自定义参数
        duplicate_query = duplicate_query.filter(
            QCJob.solvent_model == 'custom'
        )
        # 匹配自定义溶剂的关键参数（eps 是最重要的）
        if job_data.solvent_config:
            eps_val = job_data.solvent_config.eps
            if eps_val is not None:
                duplicate_query = duplicate_query.filter(
                    QCJob.config['solvent_config']['eps'].astext == str(eps_val)
                )
    else:
        duplicate_query = duplicate_query.filter(
            QCJob.solvent_model == solvent_model
        )
        if solvent_name:
            duplicate_query = duplicate_query.filter(
                QCJob.solvent_name == solvent_name
            )

    existing_job = duplicate_query.order_by(desc(QCJob.created_at)).first()

    # 对于自定义溶剂，额外验证所有参数是否完全匹配
    if existing_job and solvent_model == 'custom' and job_data.solvent_config:
        existing_config = existing_job.config.get('solvent_config', {}) if existing_job.config else {}
        # 检查所有关键参数是否匹配
        key_params = ['eps', 'eps_inf', 'hbond_acidity', 'hbond_basicity', 'surface_tension']
        params_match = True
        for key in key_params:
            new_val = getattr(job_data.solvent_config, key, None)
            existing_val = existing_config.get(key)
            if new_val != existing_val:
                params_match = False
                break
        if not params_match:
            existing_job = None  # 参数不完全匹配，不视为重复

    if existing_job:
        # 找到完全相同参数的任务，拒绝创建
        status_text = {
            QCJobStatusModel.CREATED: "已创建",
            QCJobStatusModel.PENDING: "等待中",
            QCJobStatusModel.RUNNING: "运行中",
            QCJobStatusModel.COMPLETED: "已完成",
            QCJobStatusModel.FAILED: "失败",
            QCJobStatusModel.CANCELLED: "已取消"
        }.get(existing_job.status, str(existing_job.status))

        error_msg = (
            f"检测到重复计算！已存在完全相同参数的QC计算任务（ID: {existing_job.id}）。\n"
            f"分子: {existing_job.molecule_name}\n"
            f"泛函: {existing_job.functional}, 基组: {existing_job.basis_set}\n"
            f"溶剂: {solvent_model}"
        )
        if solvent_name:
            error_msg += f" ({solvent_name})"
        error_msg += f"\n任务状态: {status_text}"
        if existing_job.status == QCJobStatusModel.COMPLETED and existing_job.finished_at:
            error_msg += f"\n完成时间: {existing_job.finished_at.strftime('%Y-%m-%d %H:%M:%S')}"
        error_msg += f"\n\n请直接查看已有任务（ID: {existing_job.id}），避免重复计算浪费资源。"

        logger.warning(f"Duplicate QC job detected: {job_data.smiles[:30]}... matches job {existing_job.id} (status: {existing_job.status})")
        raise HTTPException(
            status_code=409,  # 409 Conflict
            detail={
                "message": error_msg,
                "existing_job_id": existing_job.id,
                "existing_job_name": existing_job.molecule_name,
                "existing_job_status": existing_job.status.value if hasattr(existing_job.status, 'value') else str(existing_job.status),
                "completed_at": existing_job.finished_at.isoformat() if existing_job.finished_at else None
            }
        )

    # 构建配置 - 包含精度等级、溶剂配置和Slurm资源配置
    config = job_data.config or {}
    config["accuracy_level"] = job_data.accuracy_level.value if job_data.accuracy_level else "standard"
    config["auto_spin"] = job_data.auto_spin
    if job_data.solvent_config:
        config["solvent_config"] = job_data.solvent_config.model_dump()

    # Slurm 资源配置
    config["slurm_partition"] = job_data.slurm_partition or "cpu"
    config["slurm_cpus"] = job_data.slurm_cpus or 16
    config["slurm_time"] = job_data.slurm_time or 7200

    # 创建QC任务
    db_job = QCJob(
        user_id=current_user.id,
        md_job_id=job_data.md_job_id,
        molecule_name=job_data.molecule_name,
        smiles=job_data.smiles,
        molecule_type=job_data.molecule_type,
        basis_set=job_data.basis_set,
        functional=job_data.functional,
        charge=job_data.charge,
        spin_multiplicity=job_data.spin_multiplicity,
        solvent_model=solvent_model,
        solvent_name=solvent_name,
        accuracy_level=job_data.accuracy_level.value if job_data.accuracy_level else "standard",
        slurm_partition=job_data.slurm_partition or "cpu",
        slurm_cpus=job_data.slurm_cpus or 16,
        slurm_time=job_data.slurm_time or 7200,
        config=config,
        status=QCJobStatusModel.CREATED
    )
    
    db.add(db_job)
    db.commit()
    db.refresh(db_job)
    
    logger.info(f"Created QC job {db_job.id} for molecule {job_data.molecule_name} by user {current_user.username}")
    
    return db_job


@router.post("/jobs/batch", response_model=List[QCJobSchema])
def create_qc_jobs_batch(
    batch_data: QCJobBatchCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    批量创建QC计算任务（带查重）

    Args:
        batch_data: 批量创建数据

    Returns:
        创建的QC任务列表（包括复用的任务）
    """
    created_jobs = []
    reused_count = 0
    skipped_count = 0

    for mol_data in batch_data.molecules:
        # 使用批量参数覆盖单个参数（如果有指定）
        basis_set = mol_data.basis_set or batch_data.basis_set
        functional = mol_data.functional or batch_data.functional
        md_job_id = mol_data.md_job_id or batch_data.md_job_id

        # 提取溶剂配置
        solvent_model = 'gas'
        solvent_name = None
        config = mol_data.config or {}
        if config.get('solvent_config'):
            solvent_model = config['solvent_config'].get('model', 'gas')
            solvent_name = config['solvent_config'].get('solvent_name')

        # ======== 查重逻辑 ========
        duplicate_query = db.query(QCJob).filter(
            QCJob.smiles == mol_data.smiles,
            QCJob.functional == functional,
            QCJob.basis_set == basis_set,
            QCJob.charge == mol_data.charge,
            QCJob.spin_multiplicity == mol_data.spin_multiplicity,
            QCJob.is_deleted == False
        )

        # 匹配溶剂配置
        if solvent_model == 'gas':
            duplicate_query = duplicate_query.filter(
                or_(
                    QCJob.solvent_model == 'gas',
                    QCJob.solvent_model.is_(None)
                )
            )
        elif solvent_model == 'custom':
            # 自定义溶剂：需要匹配所有自定义参数
            duplicate_query = duplicate_query.filter(
                QCJob.solvent_model == 'custom'
            )
            # 匹配自定义溶剂的关键参数（eps 是最重要的）
            custom_solvent_config = config.get('solvent_config', {})
            eps_val = custom_solvent_config.get('eps')
            if eps_val is not None:
                duplicate_query = duplicate_query.filter(
                    QCJob.config['solvent_config']['eps'].astext == str(eps_val)
                )
        else:
            duplicate_query = duplicate_query.filter(
                QCJob.solvent_model == solvent_model
            )
            if solvent_name:
                duplicate_query = duplicate_query.filter(
                    QCJob.solvent_name == solvent_name
                )

        existing_job = duplicate_query.first()

        # 对于自定义溶剂，额外验证所有参数是否完全匹配
        if existing_job and solvent_model == 'custom':
            existing_config = existing_job.config.get('solvent_config', {}) if existing_job.config else {}
            custom_solvent_config = config.get('solvent_config', {})
            # 检查所有关键参数是否匹配
            key_params = ['eps', 'eps_inf', 'hbond_acidity', 'hbond_basicity', 'surface_tension']
            params_match = True
            for key in key_params:
                if custom_solvent_config.get(key) != existing_config.get(key):
                    params_match = False
                    break
            if not params_match:
                existing_job = None  # 参数不完全匹配，不视为重复

        if existing_job:
            if existing_job.status == QCJobStatusModel.COMPLETED:
                # 复用已完成的任务
                logger.info(f"Batch: Reusing completed QC job {existing_job.id} for '{mol_data.molecule_name}'")
                db_job = QCJob(
                    user_id=current_user.id,
                    md_job_id=md_job_id,
                    molecule_name=mol_data.molecule_name,
                    smiles=mol_data.smiles,
                    molecule_type=mol_data.molecule_type,
                    basis_set=basis_set,
                    functional=functional,
                    charge=mol_data.charge,
                    spin_multiplicity=mol_data.spin_multiplicity,
                    solvent_model=solvent_model,
                    solvent_name=solvent_name,
                    config=config,
                    status=QCJobStatusModel.COMPLETED,
                    is_reused=True,
                    reused_from_job_id=existing_job.id
                )
                db.add(db_job)
                created_jobs.append(db_job)
                reused_count += 1
            else:
                # 跳过未完成的重复任务
                logger.info(f"Batch: Skipping duplicate QC job for '{mol_data.molecule_name}' "
                           f"(existing job {existing_job.id} status: {existing_job.status})")
                skipped_count += 1
            continue
        # ======== 查重逻辑结束 ========

        db_job = QCJob(
            user_id=current_user.id,
            md_job_id=md_job_id,
            molecule_name=mol_data.molecule_name,
            smiles=mol_data.smiles,
            molecule_type=mol_data.molecule_type,
            basis_set=basis_set,
            functional=functional,
            charge=mol_data.charge,
            spin_multiplicity=mol_data.spin_multiplicity,
            solvent_model=solvent_model,
            solvent_name=solvent_name,
            config=config,
            status=QCJobStatusModel.CREATED
        )
        db.add(db_job)
        created_jobs.append(db_job)

    db.commit()

    for job in created_jobs:
        db.refresh(job)

    logger.info(f"Batch created {len(created_jobs)} QC jobs (reused: {reused_count}, skipped: {skipped_count}) by user {current_user.username}")

    return created_jobs


@router.get("/jobs", response_model=QCJobListResponse)
def list_qc_jobs(
    status: Optional[str] = Query(None, description="按状态筛选"),
    md_job_id: Optional[int] = Query(None, description="按关联的MD任务筛选"),
    molecule_name: Optional[str] = Query(None, description="按分子名称筛选"),
    smiles: Optional[str] = Query(None, description="按SMILES筛选"),
    functional: Optional[str] = Query(None, description="按泛函筛选"),
    basis_set: Optional[str] = Query(None, description="按基组筛选"),
    include_deleted: bool = Query(False, description="是否包含已删除的任务（仅管理员）"),
    visibility: Optional[str] = Query(None, description="按可见性筛选（PUBLIC/DELAYED/PRIVATE）"),
    skip: int = Query(0, ge=0),
    limit: int = Query(20, ge=1, le=100),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取QC任务列表

    支持多种筛选条件：
    - status: 任务状态
    - md_job_id: 关联的MD任务ID
    - molecule_name: 分子名称（模糊匹配）
    - smiles: SMILES（模糊匹配）
    - functional: 泛函
    - basis_set: 基组
    - visibility: 可见性（PUBLIC/DELAYED/PRIVATE）

    如果指定 visibility=PUBLIC，则搜索所有用户的公开数据
    否则只返回当前用户的数据
    """
    # 如果是搜索公开数据
    if visibility == "PUBLIC":
        query = db.query(QCJob).filter(QCJob.visibility == "PUBLIC")
    else:
        # 只返回当前用户的数据
        query = db.query(QCJob).filter(QCJob.user_id == current_user.id)
        # 如果指定了其他可见性，则筛选
        if visibility:
            query = query.filter(QCJob.visibility == visibility)

    # 排除已删除的任务（管理员可以通过特殊参数查看）
    # 普通用户看不到已删除的数据
    if current_user.role != UserRole.ADMIN:
        query = query.filter(or_(QCJob.is_deleted == False, QCJob.is_deleted.is_(None)))
    elif not include_deleted:
        # 管理员默认也不显示已删除的任务，除非明确要求
        query = query.filter(or_(QCJob.is_deleted == False, QCJob.is_deleted.is_(None)))

    # 其他筛选条件
    if status:
        query = query.filter(QCJob.status == status)
    if md_job_id:
        query = query.filter(QCJob.md_job_id == md_job_id)
    if molecule_name:
        query = query.filter(QCJob.molecule_name.ilike(f"%{molecule_name}%"))
    if smiles:
        query = query.filter(QCJob.smiles.ilike(f"%{smiles}%"))
    if functional:
        query = query.filter(QCJob.functional == functional)
    if basis_set:
        query = query.filter(QCJob.basis_set == basis_set)

    total = query.count()
    jobs = query.options(selectinload(QCJob.results)).order_by(desc(QCJob.created_at)).offset(skip).limit(limit).all()

    # 处理复用任务：如果任务是复用的且没有自己的结果，获取原始任务的结果
    for job in jobs:
        if job.is_reused and job.reused_from_job_id and len(job.results) == 0:
            original_results = db.query(QCResult).filter(
                QCResult.qc_job_id == job.reused_from_job_id
            ).all()
            if original_results:
                job.results = original_results

    return QCJobListResponse(total=total, jobs=jobs)


@router.get("/jobs/{job_id}", response_model=QCJobWithResults)
def get_qc_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取QC任务详情"""
    from datetime import datetime
    from app.models.job import DataVisibility

    job = db.query(QCJob).filter(QCJob.id == job_id).first()

    if not job:
        raise HTTPException(status_code=404, detail="QC job not found")

    # 检查权限（支持公开数据访问）
    is_owner = job.user_id == current_user.id
    is_admin = current_user.role == UserRole.ADMIN
    is_public = job.visibility == "PUBLIC"
    is_delayed_expired = (
        job.visibility == "DELAYED" and
        job.visibility_delay_until and
        job.visibility_delay_until <= datetime.utcnow()
    )

    if not (is_owner or is_admin or is_public or is_delayed_expired):
        raise HTTPException(status_code=403, detail="Permission denied")

    # 如果是复用任务且没有自己的结果，获取原始任务的结果
    if job.is_reused and job.reused_from_job_id and len(job.results) == 0:
        original_results = db.query(QCResult).filter(
            QCResult.qc_job_id == job.reused_from_job_id
        ).all()
        if original_results:
            # 动态添加原始任务的结果到当前任务对象
            # 注意：这不会持久化到数据库，只是为了返回给前端
            job.results = original_results
            logger.info(f"QC job {job_id} is reused from {job.reused_from_job_id}, returning original results")

    return job


@router.get("/jobs/{job_id}/status")
def get_qc_job_status(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取QC任务状态（轻量级轮询接口）"""
    job = db.query(QCJob).filter(QCJob.id == job_id).first()

    if not job:
        raise HTTPException(status_code=404, detail="QC job not found")

    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="Permission denied")

    return {
        "id": job.id,
        "status": job.status.value,
        "progress": job.progress,
        "error_message": job.error_message,
        "slurm_job_id": job.slurm_job_id,
        "updated_at": job.updated_at.isoformat() if job.updated_at else None
    }


@router.put("/jobs/{job_id}")
def update_qc_job(
    job_id: int,
    job_data: QCJobEdit,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """编辑QC任务（仅CREATED状态可编辑）"""
    job = db.query(QCJob).filter(QCJob.id == job_id).first()

    if not job:
        raise HTTPException(status_code=404, detail="QC job not found")

    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="Permission denied")

    # 只有CREATED状态的任务可以编辑
    if job.status != QCJobStatusModel.CREATED:
        raise HTTPException(
            status_code=400,
            detail=f"只有CREATED状态的任务可以编辑，当前状态: {job.status.value}"
        )

    # 更新字段
    update_data = job_data.model_dump(exclude_unset=True)

    for field, value in update_data.items():
        if field == 'solvent_config' and value is not None:
            # 处理溶剂配置
            setattr(job, field, value.model_dump() if hasattr(value, 'model_dump') else value)
        elif field == 'molecule_type' and value is not None:
            setattr(job, field, value)
        elif field == 'accuracy_level' and value is not None:
            setattr(job, field, value)
        else:
            setattr(job, field, value)

    db.commit()
    db.refresh(job)

    logger.info(f"Updated QC job {job_id} by user {current_user.username}")

    return job


@router.delete("/jobs/{job_id}")
def delete_qc_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    删除/取消QC任务

    - 运行中的任务：取消任务但保留记录
    - 未完成/失败的任务：真正删除
    - 已完成的任务：软删除（保留计算结果数据供公开数据库使用）
    """
    from datetime import datetime

    job = db.query(QCJob).filter(QCJob.id == job_id).first()

    if not job:
        raise HTTPException(status_code=404, detail="QC job not found")

    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="Permission denied")

    # 如果任务正在运行，先取消
    if job.status in [QCJobStatusModel.QUEUED, QCJobStatusModel.RUNNING]:
        if job.slurm_job_id:
            try:
                import subprocess
                subprocess.run(["scancel", job.slurm_job_id], check=True)
                logger.info(f"Cancelled Slurm job {job.slurm_job_id}")
            except Exception as e:
                logger.warning(f"Failed to cancel Slurm job: {e}")

        job.status = QCJobStatusModel.CANCELLED
        db.commit()
        return {"message": "QC job cancelled", "id": job_id}

    # 已完成的任务：软删除（保留数据供公开数据库和管理员使用）
    if job.status == QCJobStatusModel.COMPLETED:
        job.is_deleted = True
        job.deleted_at = datetime.now()
        job.deleted_by = current_user.id
        job.delete_reason = "用户主动删除（数据已保留供公开使用）"
        db.commit()
        logger.info(f"Soft deleted completed QC job {job_id} by user {current_user.username}")
        return {"message": "QC job deleted (data preserved)", "id": job_id, "soft_delete": True}

    # 未完成/失败/取消的任务：真正删除
    db.delete(job)
    db.commit()
    logger.info(f"Permanently deleted QC job {job_id} by user {current_user.username}")
    return {"message": "QC job permanently deleted", "id": job_id, "soft_delete": False}


@router.delete("/jobs/batch/delete")
def batch_delete_qc_jobs(
    ids: List[int],
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    批量删除QC任务

    - 运行中的任务：取消
    - 未完成/失败的任务：真正删除
    - 已完成的任务：软删除（保留数据）
    """
    from datetime import datetime

    soft_deleted_count = 0  # 软删除（已完成任务）
    hard_deleted_count = 0  # 真正删除（未完成任务）
    cancelled_count = 0
    failed_ids = []

    for job_id in ids:
        job = db.query(QCJob).filter(QCJob.id == job_id).first()
        if not job:
            failed_ids.append(job_id)
            continue

        # 检查权限
        if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
            failed_ids.append(job_id)
            continue

        # 如果任务正在运行，先取消
        if job.status in [QCJobStatusModel.QUEUED, QCJobStatusModel.RUNNING]:
            if job.slurm_job_id:
                try:
                    import subprocess
                    subprocess.run(["scancel", job.slurm_job_id], check=True)
                except Exception as e:
                    logger.warning(f"Failed to cancel Slurm job: {e}")
            job.status = QCJobStatusModel.CANCELLED
            cancelled_count += 1
        elif job.status == QCJobStatusModel.COMPLETED:
            # 已完成的任务：软删除
            job.is_deleted = True
            job.deleted_at = datetime.now()
            job.deleted_by = current_user.id
            job.delete_reason = "用户批量删除（数据已保留供公开使用）"
            soft_deleted_count += 1
        else:
            # 未完成/失败/取消的任务：真正删除
            db.delete(job)
            hard_deleted_count += 1

    db.commit()

    logger.info(f"Batch delete: soft={soft_deleted_count}, hard={hard_deleted_count}, cancelled={cancelled_count} by {current_user.username}")

    return {
        "deleted_count": soft_deleted_count + hard_deleted_count,
        "soft_deleted_count": soft_deleted_count,
        "hard_deleted_count": hard_deleted_count,
        "cancelled_count": cancelled_count,
        "failed_ids": failed_ids,
        "message": f"成功删除 {soft_deleted_count + hard_deleted_count} 个任务，取消 {cancelled_count} 个任务"
    }


@router.post("/jobs/batch/submit")
def batch_submit_qc_jobs(
    ids: List[int],
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    批量提交QC任务到计算集群

    在混合云架构中，任务由 Polling Worker 自动获取并处理。
    """
    success_count = 0
    failed_count = 0
    errors = []

    for job_id in ids:
        try:
            job = db.query(QCJob).filter(QCJob.id == job_id).first()

            if not job:
                failed_count += 1
                errors.append({"job_id": job_id, "error": "任务不存在"})
                continue

            # 检查权限
            if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
                failed_count += 1
                errors.append({"job_id": job_id, "error": "权限不足"})
                continue

            # 检查状态：只能提交 CREATED 或 FAILED/CANCELLED 的任务
            if job.status not in [QCJobStatusModel.CREATED, QCJobStatusModel.FAILED, QCJobStatusModel.CANCELLED]:
                failed_count += 1
                errors.append({"job_id": job_id, "error": f"任务状态为 {job.status}，无法提交"})
                continue

            # 更新状态为 SUBMITTED，Polling Worker 会拉取
            job.status = QCJobStatusModel.SUBMITTED
            job.config = job.config or {}
            job.config["submitted_at"] = datetime.now().isoformat()
            job.config["submitted_by"] = current_user.username
            job.error_message = None  # 清除之前的错误信息
            db.commit()

            success_count += 1
            logger.info(f"QC job {job_id} status=SUBMITTED, waiting for polling worker")

        except Exception as e:
            failed_count += 1
            errors.append({"job_id": job_id, "error": str(e)})
            logger.error(f"Failed to mark QC job {job_id} for submission: {e}")

    return {
        "success_count": success_count,
        "failed_count": failed_count,
        "errors": errors,
        "message": f"成功提交 {success_count} 个任务，失败 {failed_count} 个"
    }


@router.post("/jobs/batch/cancel")
def batch_cancel_qc_jobs(
    ids: List[int],
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    批量取消QC任务
    """
    import subprocess

    success_count = 0
    failed_count = 0
    errors = []

    for job_id in ids:
        try:
            job = db.query(QCJob).filter(QCJob.id == job_id).first()

            if not job:
                failed_count += 1
                errors.append({"job_id": job_id, "error": "任务不存在"})
                continue

            # 检查权限
            if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
                failed_count += 1
                errors.append({"job_id": job_id, "error": "权限不足"})
                continue

            # 检查状态
            if job.status not in [QCJobStatusModel.CREATED, QCJobStatusModel.PENDING, QCJobStatusModel.QUEUED, QCJobStatusModel.RUNNING]:
                failed_count += 1
                errors.append({"job_id": job_id, "error": f"任务状态为 {job.status}，无法取消"})
                continue

            # 如果有Slurm任务ID，取消Slurm任务
            if job.slurm_job_id:
                try:
                    subprocess.run(["scancel", job.slurm_job_id], check=True)
                    logger.info(f"Cancelled Slurm job {job.slurm_job_id} for QC job {job_id}")
                except Exception as e:
                    logger.warning(f"Failed to cancel Slurm job {job.slurm_job_id}: {e}")

            # 更新任务状态
            job.status = QCJobStatusModel.CANCELLED
            job.error_message = "用户批量取消"
            db.commit()

            success_count += 1
            logger.info(f"QC job {job_id} cancelled by {current_user.username}")

        except Exception as e:
            failed_count += 1
            errors.append({"job_id": job_id, "error": str(e)})
            logger.error(f"Failed to cancel QC job {job_id}: {e}")

    return {
        "success_count": success_count,
        "failed_count": failed_count,
        "errors": errors,
        "message": f"成功取消 {success_count} 个任务，失败 {failed_count} 个"
    }


@router.get("/admin/deleted-jobs")
def list_deleted_qc_jobs(
    skip: int = Query(0, ge=0),
    limit: int = Query(20, ge=1, le=100),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取已删除的QC任务列表（仅管理员）

    这些是用户删除但保留了计算结果的任务，可用于公开数据库
    """
    if current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="只有管理员可以查看已删除的任务")

    query = db.query(QCJob).filter(QCJob.is_deleted == True)
    total = query.count()
    jobs = query.order_by(desc(QCJob.deleted_at)).offset(skip).limit(limit).all()

    # 返回包含删除信息的任务列表
    result = []
    for job in jobs:
        job_dict = {
            "id": job.id,
            "molecule_name": job.molecule_name,
            "smiles": job.smiles,
            "status": job.status.value if job.status else None,
            "functional": job.functional,
            "basis_set": job.basis_set,
            "user_id": job.user_id,
            "deleted_at": job.deleted_at.isoformat() if job.deleted_at else None,
            "deleted_by": job.deleted_by,
            "delete_reason": job.delete_reason,
            "created_at": job.created_at.isoformat() if job.created_at else None,
            "finished_at": job.finished_at.isoformat() if job.finished_at else None,
        }
        result.append(job_dict)

    return {"total": total, "jobs": result}


@router.post("/jobs/{job_id}/restore")
def restore_qc_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """恢复已删除的QC任务（仅管理员）"""
    if current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="只有管理员可以恢复已删除的任务")

    job = db.query(QCJob).filter(QCJob.id == job_id).first()

    if not job:
        raise HTTPException(status_code=404, detail="QC job not found")

    if not job.is_deleted:
        raise HTTPException(status_code=400, detail="任务未被删除，无需恢复")

    job.is_deleted = False
    job.deleted_at = None
    job.deleted_by = None
    job.delete_reason = None
    db.commit()

    logger.info(f"Restored QC job {job_id} by admin {current_user.username}")

    return {"message": "任务已恢复", "id": job_id}


@router.post("/jobs/{job_id}/submit")
def submit_qc_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    提交QC任务到计算集群

    在混合云架构中，任务由 Polling Worker 自动获取并处理。
    这里只需要确保任务状态为 CREATED，Worker 会自动处理。
    """
    job = db.query(QCJob).filter(QCJob.id == job_id).first()

    if not job:
        raise HTTPException(status_code=404, detail="QC job not found")

    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="Permission denied")

    # 只能提交 CREATED 或 FAILED/CANCELLED 的任务
    if job.status not in [QCJobStatusModel.CREATED, QCJobStatusModel.FAILED, QCJobStatusModel.CANCELLED]:
        raise HTTPException(status_code=400, detail=f"Job cannot be submitted in {job.status} status")

    # 更新状态为 SUBMITTED，Polling Worker 会拉取
    job.status = QCJobStatusModel.SUBMITTED
    job.config = job.config or {}
    job.config["submitted_at"] = datetime.now().isoformat()
    job.config["submitted_by"] = current_user.username
    job.error_message = None  # 清除之前的错误信息
    db.commit()

    logger.info(f"QC job {job_id} status=SUBMITTED, waiting for polling worker")

    return {
        "message": "QC任务已提交，等待计算集群处理",
        "job_id": job_id,
        "status": "SUBMITTED"
    }


@router.post("/jobs/{job_id}/recalculate", response_model=QCJobSchema)
def recalculate_qc_job(
    job_id: int,
    recalc_params: QCJobRecalculate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    基于已有QC任务创建新的计算任务（重新计算）

    - 复用原任务的分子信息（SMILES、电荷、自旋多重度）
    - 允许修改计算参数（泛函、基组、溶剂模型）
    - 自动关联到原任务
    - 新任务状态为CREATED，需要手动提交
    """
    # 获取原任务
    original_job = db.query(QCJob).filter(QCJob.id == job_id).first()

    if not original_job:
        raise HTTPException(status_code=404, detail="原任务不存在")

    # 检查权限：只有任务所有者或管理员可以重新计算
    if original_job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="无权限重新计算此任务")

    # 准备新任务的配置
    new_config = original_job.config.copy() if original_job.config else {}
    new_config["recalculated_from"] = job_id
    new_config["recalculated_at"] = datetime.now().isoformat()

    # 使用新参数或原参数
    new_functional = recalc_params.functional or original_job.functional
    new_basis_set = recalc_params.basis_set or original_job.basis_set

    # 溶剂配置
    if recalc_params.solvent_config:
        new_config["solvent_config"] = recalc_params.solvent_config.dict()
    elif "solvent_config" in new_config:
        # 保留原溶剂配置
        pass

    # 创建新任务
    new_job = QCJob(
        user_id=current_user.id,
        md_job_id=original_job.md_job_id,
        molecule_name=f"{original_job.molecule_name}_recalc",
        smiles=original_job.smiles,
        molecule_type=original_job.molecule_type,
        charge=original_job.charge,
        spin_multiplicity=original_job.spin_multiplicity,
        functional=new_functional,
        basis_set=new_basis_set,
        config=new_config,
        status=QCJobStatusModel.CREATED
    )

    db.add(new_job)
    db.commit()
    db.refresh(new_job)

    logger.info(f"Created recalculation job {new_job.id} from original job {job_id} by user {current_user.username}")

    return new_job


# ============================================================================
# QC Results Operations
# ============================================================================

@router.get("/results/{job_id}", response_model=List[QCResultSchema])
def get_qc_results(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取QC任务的计算结果"""
    from datetime import datetime

    job = db.query(QCJob).filter(QCJob.id == job_id).first()

    if not job:
        raise HTTPException(status_code=404, detail="QC job not found")

    # 检查权限（支持公开数据访问）
    is_owner = job.user_id == current_user.id
    is_admin = current_user.role == UserRole.ADMIN
    is_public = job.visibility == "PUBLIC"
    is_delayed_expired = (
        job.visibility == "DELAYED" and
        job.visibility_delay_until and
        job.visibility_delay_until <= datetime.utcnow()
    )

    if not (is_owner or is_admin or is_public or is_delayed_expired):
        raise HTTPException(status_code=403, detail="Permission denied")

    results = db.query(QCResult).filter(QCResult.qc_job_id == job_id).all()

    # 如果是复用任务且没有自己的结果，获取原始任务的结果
    if not results and job.is_reused and job.reused_from_job_id:
        results = db.query(QCResult).filter(
            QCResult.qc_job_id == job.reused_from_job_id
        ).all()
        if results:
            logger.info(f"QC job {job_id} is reused from {job.reused_from_job_id}, returning original results")

    return results


@router.get("/results/by-smiles")
def get_qc_results_by_smiles(
    smiles: str = Query(..., description="SMILES表达式"),
    basis_set: Optional[str] = Query(None, description="基组"),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """根据SMILES查询QC结果"""
    query = db.query(QCResult).filter(QCResult.smiles == smiles)

    if basis_set:
        # 需要join QCJob来过滤basis_set
        query = query.join(QCJob).filter(QCJob.basis_set == basis_set)

    results = query.all()

    return results


# ============================================================================
# Molecule QC Cache Operations
# ============================================================================

@router.get("/cache/{smiles:path}", response_model=MoleculeQCCacheSchema)
def get_molecule_qc_cache(
    smiles: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取分子的QC缓存数据"""
    cache = db.query(MoleculeQCCache).filter(MoleculeQCCache.smiles == smiles).first()

    if not cache:
        raise HTTPException(status_code=404, detail="No QC data found for this molecule")

    return cache


@router.get("/cache")
def search_molecule_qc_cache(
    smiles: Optional[str] = Query(None),
    molecule_name: Optional[str] = Query(None),
    lumo_min: Optional[float] = Query(None),
    lumo_max: Optional[float] = Query(None),
    homo_min: Optional[float] = Query(None),
    homo_max: Optional[float] = Query(None),
    skip: int = Query(0, ge=0),
    limit: int = Query(20, ge=1, le=100),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """搜索分子QC缓存"""
    query = db.query(MoleculeQCCache)

    if smiles:
        query = query.filter(MoleculeQCCache.smiles.ilike(f"%{smiles}%"))
    if molecule_name:
        query = query.filter(MoleculeQCCache.molecule_name.ilike(f"%{molecule_name}%"))
    if lumo_min is not None:
        query = query.filter(MoleculeQCCache.lumo_ev >= lumo_min)
    if lumo_max is not None:
        query = query.filter(MoleculeQCCache.lumo_ev <= lumo_max)
    if homo_min is not None:
        query = query.filter(MoleculeQCCache.homo_ev >= homo_min)
    if homo_max is not None:
        query = query.filter(MoleculeQCCache.homo_ev <= homo_max)

    total = query.count()
    results = query.offset(skip).limit(limit).all()

    return {
        "total": total,
        "data": results
    }


# ============================================================================
# Utility Endpoints
# ============================================================================

@router.get("/basis-sets")
def get_available_basis_sets():
    """获取可用的基组列表"""
    return {
        "basis_sets": [
            {"value": "6-31++g(d,p)", "label": "6-31++G(d,p)", "description": "常用基组，精度适中"},
            {"value": "6-311g(d,p)", "label": "6-311G(d,p)", "description": "三重分裂基组"},
            {"value": "6-311++g(d,p)", "label": "6-311++G(d,p)", "description": "带弥散函数，适合阴离子"},
            {"value": "Def2TZVP", "label": "Def2-TZVP", "description": "高精度基组"},
        ]
    }


@router.get("/functionals")
def get_available_functionals():
    """获取可用的泛函列表"""
    return {
        "functionals": [
            {"value": "B3LYP", "label": "B3LYP", "description": "最常用的杂化泛函"},
            {"value": "M062X", "label": "M06-2X", "description": "适合非共价相互作用"},
            {"value": "wB97XD", "label": "ωB97X-D", "description": "带色散校正"},
            {"value": "PBE0", "label": "PBE0", "description": "无经验参数的杂化泛函"},
        ]
    }


# ============================================================================
# ESP Image Endpoints
# ============================================================================

def get_user_from_token_param(
    token: Optional[str] = Query(None, description="JWT Token for image access"),
    db: Session = Depends(get_db)
) -> Optional[User]:
    """从query参数获取用户（用于图片等资源的访问）"""
    if not token:
        return None
    from app.core.security import decode_access_token
    payload = decode_access_token(token)
    if payload is None:
        return None
    username = payload.get("sub")
    if username is None:
        return None
    return db.query(User).filter(User.username == username).first()


@router.get("/esp-image/{result_id}")
def get_esp_image(
    result_id: int,
    token: Optional[str] = Query(None, description="JWT Token"),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取QC结果的ESP图片"""
    from fastapi.responses import Response
    from datetime import datetime
    import base64

    result = db.query(QCResult).filter(QCResult.id == result_id).first()

    if not result:
        raise HTTPException(status_code=404, detail="QC result not found")

    # 检查权限（支持公开数据访问）
    qc_job = db.query(QCJob).filter(QCJob.id == result.qc_job_id).first()
    if qc_job:
        is_owner = qc_job.user_id == current_user.id
        is_admin = current_user.role == UserRole.ADMIN
        is_public = qc_job.visibility == "PUBLIC"
        is_delayed_expired = (
            qc_job.visibility == "DELAYED" and
            qc_job.visibility_delay_until and
            qc_job.visibility_delay_until <= datetime.utcnow()
        )

        if not (is_owner or is_admin or is_public or is_delayed_expired):
            raise HTTPException(status_code=403, detail="Permission denied")

    # 优先使用数据库中的图片内容（混合云架构）
    if result.esp_image_content:
        try:
            image_data = base64.b64decode(result.esp_image_content)
            return Response(
                content=image_data,
                media_type="image/png",
                headers={"Content-Disposition": f"inline; filename=esp_{result_id}.png"}
            )
        except Exception as e:
            logger.error(f"Failed to decode ESP image content: {e}")

    # 回退到文件路径（本地部署）
    if result.esp_image_path and os.path.exists(result.esp_image_path):
        return FileResponse(
            result.esp_image_path,
            media_type="image/png",
            filename=f"esp_{result_id}.png",
            headers={"Content-Disposition": f"inline; filename=esp_{result_id}.png"}
        )

    raise HTTPException(status_code=404, detail="ESP image not found")


@router.get("/esp-image-download/{result_id}")
def download_esp_image(
    result_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """下载QC结果的ESP图片"""
    result = db.query(QCResult).filter(QCResult.id == result_id).first()

    if not result:
        raise HTTPException(status_code=404, detail="QC result not found")

    # 检查权限
    qc_job = db.query(QCJob).filter(QCJob.id == result.qc_job_id).first()
    if qc_job and qc_job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="Permission denied")

    if not result.esp_image_path or not os.path.exists(result.esp_image_path):
        raise HTTPException(status_code=404, detail="ESP image not found")

    return FileResponse(
        result.esp_image_path,
        media_type="image/png",
        filename=f"esp_{result_id}.png",
        headers={"Content-Disposition": f"attachment; filename=esp_{result_id}.png"}
    )


@router.get("/homo-image/{result_id}")
def get_homo_image(
    result_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取QC结果的HOMO轨道图片"""
    from fastapi.responses import Response
    from datetime import datetime
    import base64

    result = db.query(QCResult).filter(QCResult.id == result_id).first()

    if not result:
        raise HTTPException(status_code=404, detail="QC result not found")

    # 检查权限（支持公开数据访问）
    qc_job = db.query(QCJob).filter(QCJob.id == result.qc_job_id).first()
    if qc_job:
        is_owner = qc_job.user_id == current_user.id
        is_admin = current_user.role == UserRole.ADMIN
        is_public = qc_job.visibility == "PUBLIC"
        is_delayed_expired = (
            qc_job.visibility == "DELAYED" and
            qc_job.visibility_delay_until and
            qc_job.visibility_delay_until <= datetime.utcnow()
        )

        if not (is_owner or is_admin or is_public or is_delayed_expired):
            raise HTTPException(status_code=403, detail="Permission denied")

    # 优先使用数据库中的图片内容（混合云架构）
    if result.homo_image_content:
        try:
            image_data = base64.b64decode(result.homo_image_content)
            return Response(
                content=image_data,
                media_type="image/png",
                headers={"Content-Disposition": f"inline; filename=homo_{result_id}.png"}
            )
        except Exception as e:
            logger.error(f"Failed to decode HOMO image content: {e}")

    # 回退到文件路径（本地部署）
    if result.homo_image_path and os.path.exists(result.homo_image_path):
        return FileResponse(
            result.homo_image_path,
            media_type="image/png",
            filename=f"homo_{result_id}.png",
            headers={"Content-Disposition": f"inline; filename=homo_{result_id}.png"}
        )

    raise HTTPException(status_code=404, detail="HOMO image not found")


@router.get("/lumo-image/{result_id}")
def get_lumo_image(
    result_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取QC结果的LUMO轨道图片"""
    from fastapi.responses import Response
    from datetime import datetime
    import base64

    result = db.query(QCResult).filter(QCResult.id == result_id).first()

    if not result:
        raise HTTPException(status_code=404, detail="QC result not found")

    # 检查权限（支持公开数据访问）
    qc_job = db.query(QCJob).filter(QCJob.id == result.qc_job_id).first()
    if qc_job:
        is_owner = qc_job.user_id == current_user.id
        is_admin = current_user.role == UserRole.ADMIN
        is_public = qc_job.visibility == "PUBLIC"
        is_delayed_expired = (
            qc_job.visibility == "DELAYED" and
            qc_job.visibility_delay_until and
            qc_job.visibility_delay_until <= datetime.utcnow()
        )

        if not (is_owner or is_admin or is_public or is_delayed_expired):
            raise HTTPException(status_code=403, detail="Permission denied")

    # 优先使用数据库中的图片内容（混合云架构）
    if result.lumo_image_content:
        try:
            image_data = base64.b64decode(result.lumo_image_content)
            return Response(
                content=image_data,
                media_type="image/png",
                headers={"Content-Disposition": f"inline; filename=lumo_{result_id}.png"}
            )
        except Exception as e:
            logger.error(f"Failed to decode LUMO image content: {e}")

    # 回退到文件路径（本地部署）
    if result.lumo_image_path and os.path.exists(result.lumo_image_path):
        return FileResponse(
            result.lumo_image_path,
            media_type="image/png",
            filename=f"lumo_{result_id}.png",
            headers={"Content-Disposition": f"inline; filename=lumo_{result_id}.png"}
        )

    raise HTTPException(status_code=404, detail="LUMO image not found")


@router.get("/esp-image-by-smiles/{smiles:path}")
def get_esp_image_by_smiles(
    smiles: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """根据SMILES获取ESP图片"""
    from urllib.parse import unquote
    smiles = unquote(smiles)

    cache = db.query(MoleculeQCCache).filter(
        MoleculeQCCache.smiles == smiles
    ).first()

    if not cache:
        raise HTTPException(status_code=404, detail="No QC cache found for this SMILES")

    if not cache.esp_image_path or not os.path.exists(cache.esp_image_path):
        raise HTTPException(status_code=404, detail="ESP image not found")

    return FileResponse(
        cache.esp_image_path,
        media_type="image/png",
        filename=f"esp_{cache.id}.png",
        headers={"Content-Disposition": f"inline; filename=esp_{cache.id}.png"}
    )


@router.get("/homo-image-by-smiles/{smiles:path}")
def get_homo_image_by_smiles(
    smiles: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """根据SMILES获取HOMO轨道图片"""
    from urllib.parse import unquote
    smiles = unquote(smiles)

    cache = db.query(MoleculeQCCache).filter(
        MoleculeQCCache.smiles == smiles
    ).first()

    if not cache:
        raise HTTPException(status_code=404, detail="No QC cache found for this SMILES")

    if not cache.homo_image_path or not os.path.exists(cache.homo_image_path):
        raise HTTPException(status_code=404, detail="HOMO image not found")

    return FileResponse(
        cache.homo_image_path,
        media_type="image/png",
        filename=f"homo_{cache.id}.png",
        headers={"Content-Disposition": f"inline; filename=homo_{cache.id}.png"}
    )


@router.get("/lumo-image-by-smiles/{smiles:path}")
def get_lumo_image_by_smiles(
    smiles: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """根据SMILES获取LUMO轨道图片"""
    from urllib.parse import unquote
    smiles = unquote(smiles)

    cache = db.query(MoleculeQCCache).filter(
        MoleculeQCCache.smiles == smiles
    ).first()

    if not cache:
        raise HTTPException(status_code=404, detail="No QC cache found for this SMILES")

    if not cache.lumo_image_path or not os.path.exists(cache.lumo_image_path):
        raise HTTPException(status_code=404, detail="LUMO image not found")

    return FileResponse(
        cache.lumo_image_path,
        media_type="image/png",
        filename=f"lumo_{cache.id}.png",
        headers={"Content-Disposition": f"inline; filename=lumo_{cache.id}.png"}
    )


# ============================================================================
# Configuration Endpoints
# ============================================================================

@router.get("/config/accuracy-levels")
def get_accuracy_levels():
    """获取可用的精度等级及其参数"""
    return {
        "levels": [
            {
                "value": QCAccuracyLevel.FAST.value,
                "label": "快速",
                "functional": QC_ACCURACY_PRESETS[QCAccuracyLevel.FAST]["functional"],
                "basis_set": QC_ACCURACY_PRESETS[QCAccuracyLevel.FAST]["basis_set"],
                "description": QC_ACCURACY_PRESETS[QCAccuracyLevel.FAST]["description"],
                "estimated_time": QC_ACCURACY_PRESETS[QCAccuracyLevel.FAST]["estimated_time"],
            },
            {
                "value": QCAccuracyLevel.STANDARD.value,
                "label": "标准",
                "functional": QC_ACCURACY_PRESETS[QCAccuracyLevel.STANDARD]["functional"],
                "basis_set": QC_ACCURACY_PRESETS[QCAccuracyLevel.STANDARD]["basis_set"],
                "description": QC_ACCURACY_PRESETS[QCAccuracyLevel.STANDARD]["description"],
                "estimated_time": QC_ACCURACY_PRESETS[QCAccuracyLevel.STANDARD]["estimated_time"],
            },
            {
                "value": QCAccuracyLevel.ACCURATE.value,
                "label": "精确",
                "functional": QC_ACCURACY_PRESETS[QCAccuracyLevel.ACCURATE]["functional"],
                "basis_set": QC_ACCURACY_PRESETS[QCAccuracyLevel.ACCURATE]["basis_set"],
                "description": QC_ACCURACY_PRESETS[QCAccuracyLevel.ACCURATE]["description"],
                "estimated_time": QC_ACCURACY_PRESETS[QCAccuracyLevel.ACCURATE]["estimated_time"],
            },
            {
                "value": QCAccuracyLevel.CUSTOM.value,
                "label": "自定义",
                "functional": None,
                "basis_set": None,
                "description": "自定义泛函和基组参数",
                "estimated_time": "取决于参数设置",
            },
        ]
    }


@router.get("/config/solvents")
def get_available_solvents():
    """获取可用的溶剂列表"""
    solvents = []
    for name, info in GAUSSIAN_SOLVENTS.items():
        solvents.append({
            "value": name,
            "label": f"{info['description']} ({name})",
            "eps": info["eps"],
            "description": info["description"],
        })
    # 按介电常数排序
    solvents.sort(key=lambda x: x["eps"], reverse=True)
    return {"solvents": solvents}


@router.get("/config/solvent-models")
def get_solvent_models():
    """获取可用的溶剂模型"""
    return {
        "models": [
            {
                "value": SolventModel.GAS.value,
                "label": "气相",
                "description": "无溶剂效应，真空环境计算",
            },
            {
                "value": SolventModel.PCM.value,
                "label": "PCM隐式溶剂",
                "description": "极化连续介质模型，适合大多数溶剂效应计算",
            },
            {
                "value": SolventModel.SMD.value,
                "label": "SMD隐式溶剂",
                "description": "溶剂模型密度，更精确的溶剂化自由能",
            },
            {
                "value": SolventModel.CUSTOM.value,
                "label": "自定义溶剂",
                "description": "自定义溶剂参数（需要提供7个SMD参数）",
            },
        ]
    }


@router.get("/config/basis-sets")
def get_available_basis_sets():
    """获取可用的基组列表"""
    return {
        "basis_sets": [
            {"value": BasisSet.STO3G.value, "label": "STO-3G", "category": "minimal", "description": "最小基组，快速但精度低"},
            {"value": BasisSet.B321G.value, "label": "3-21G", "category": "split-valence", "description": "分裂价层基组"},
            {"value": BasisSet.B631G.value, "label": "6-31G", "category": "split-valence", "description": "常用分裂价层基组"},
            {"value": BasisSet.B631GD.value, "label": "6-31G(d)", "category": "polarized", "description": "带极化函数，推荐用于几何优化"},
            {"value": BasisSet.B631GDP.value, "label": "6-31G(d,p)", "category": "polarized", "description": "带极化函数，适合含氢体系"},
            {"value": BasisSet.B631_PLUSPLUS_GDP.value, "label": "6-31++G(d,p)", "category": "diffuse", "description": "带弥散函数，适合阴离子和弱相互作用"},
            {"value": BasisSet.B6311_GDP.value, "label": "6-311G(d,p)", "category": "triple-zeta", "description": "三重分裂价层"},
            {"value": BasisSet.B6311_PLUSPLUS_GDP.value, "label": "6-311++G(d,p)", "category": "triple-zeta", "description": "高精度计算推荐"},
            {"value": BasisSet.DEF2SVP.value, "label": "Def2-SVP", "category": "def2", "description": "Ahlrichs基组，平衡精度和效率"},
            {"value": BasisSet.DEF2TZVP.value, "label": "Def2-TZVP", "category": "def2", "description": "高精度Ahlrichs基组"},
            {"value": BasisSet.DEF2QZVP.value, "label": "Def2-QZVP", "category": "def2", "description": "极高精度，计算量大"},
            {"value": BasisSet.CCPVDZ.value, "label": "cc-pVDZ", "category": "correlation-consistent", "description": "相关一致基组"},
            {"value": BasisSet.CCPVTZ.value, "label": "cc-pVTZ", "category": "correlation-consistent", "description": "高精度相关一致基组"},
            {"value": BasisSet.AUGCCPVDZ.value, "label": "aug-cc-pVDZ", "category": "correlation-consistent", "description": "带弥散函数的相关一致基组"},
        ]
    }


@router.get("/config/functionals")
def get_available_functionals_v2():
    """获取可用的泛函列表（增强版）"""
    return {
        "functionals": [
            {"value": Functional.HF.value, "label": "HF", "category": "wavefunction", "description": "Hartree-Fock，无电子相关"},
            {"value": Functional.B3LYP.value, "label": "B3LYP", "category": "hybrid", "description": "最常用的杂化泛函，适合大多数体系"},
            {"value": Functional.M062X.value, "label": "M06-2X", "category": "meta-hybrid", "description": "适合非共价相互作用和热化学"},
            {"value": Functional.WB97XD.value, "label": "ωB97X-D", "category": "range-separated", "description": "带色散校正，适合大分子"},
            {"value": Functional.PBE0.value, "label": "PBE0", "category": "hybrid", "description": "无经验参数的杂化泛函"},
            {"value": Functional.CAM_B3LYP.value, "label": "CAM-B3LYP", "category": "range-separated", "description": "长程校正，适合激发态"},
            {"value": Functional.B3PW91.value, "label": "B3PW91", "category": "hybrid", "description": "杂化泛函，适合过渡金属"},
            {"value": Functional.BLYP.value, "label": "BLYP", "category": "gga", "description": "纯GGA泛函，计算快"},
            {"value": Functional.PBE.value, "label": "PBE", "category": "gga", "description": "通用GGA泛函"},
        ]
    }


# ============================================================================
# Common Molecules for QC Calculation
# ============================================================================

@router.get("/config/common-molecules")
def get_common_molecules():
    """获取常用分子列表供QC计算选择"""
    return {
        "categories": [
            {
                "name": "Common Solvents",
                "molecules": [
                    {"name": "Water", "label": "Water", "smiles": "O", "charge": 0},
                    {"name": "Acetonitrile", "label": "Acetonitrile", "smiles": "CC#N", "charge": 0},
                    {"name": "Methanol", "label": "Methanol", "smiles": "CO", "charge": 0},
                    {"name": "Ethanol", "label": "Ethanol", "smiles": "CCO", "charge": 0},
                    {"name": "Acetone", "label": "Acetone", "smiles": "CC(=O)C", "charge": 0},
                    {"name": "DMSO", "label": "DMSO", "smiles": "CS(=O)C", "charge": 0},
                    {"name": "THF", "label": "THF", "smiles": "C1CCOC1", "charge": 0},
                    {"name": "DCM", "label": "DCM", "smiles": "ClCCl", "charge": 0},
                    {"name": "Chloroform", "label": "Chloroform", "smiles": "ClC(Cl)Cl", "charge": 0},
                    {"name": "Benzene", "label": "Benzene", "smiles": "c1ccccc1", "charge": 0},
                    {"name": "Toluene", "label": "Toluene", "smiles": "Cc1ccccc1", "charge": 0},
                    {"name": "DMF", "label": "DMF", "smiles": "CN(C)C=O", "charge": 0},
                ]
            },
            {
                "name": "Cations",
                "molecules": [
                    {"name": "Li", "label": "Li", "smiles": "[Li+]", "charge": 1},
                    {"name": "Na", "label": "Na", "smiles": "[Na+]", "charge": 1},
                    {"name": "K", "label": "K", "smiles": "[K+]", "charge": 1},
                ]
            },
            {
                "name": "Anions",
                "molecules": [
                    {"name": "PF6", "label": "PF6", "smiles": "F[P-](F)(F)(F)(F)F", "charge": -1},
                    {"name": "BF4", "label": "BF4", "smiles": "F[B-](F)(F)F", "charge": -1},
                    {"name": "TFSI", "label": "TFSI", "smiles": "FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F", "charge": -1},
                    {"name": "FSI", "label": "FSI", "smiles": "FS(=O)(=O)[N-]S(=O)(=O)F", "charge": -1},
                    {"name": "DFOB", "label": "DFOB", "smiles": "FB1OC(=O)C(=O)O[B-]1F", "charge": -1},
                    {"name": "ClO4", "label": "ClO4", "smiles": "[O-]Cl(=O)(=O)=O", "charge": -1},
                    {"name": "NO3", "label": "NO3", "smiles": "[O-][N+](=O)[O-]", "charge": -1},
                    {"name": "F", "label": "F", "smiles": "[F-]", "charge": -1},
                    {"name": "Cl", "label": "Cl", "smiles": "[Cl-]", "charge": -1},
                ]
            },
            {
                "name": "Carbonates",
                "molecules": [
                    {"name": "EC", "label": "EC", "smiles": "C1COC(=O)O1", "charge": 0},
                    {"name": "PC", "label": "PC", "smiles": "CC1COC(=O)O1", "charge": 0},
                    {"name": "DMC", "label": "DMC", "smiles": "COC(=O)OC", "charge": 0},
                    {"name": "DEC", "label": "DEC", "smiles": "CCOC(=O)OCC", "charge": 0},
                    {"name": "EMC", "label": "EMC", "smiles": "CCOC(=O)OC", "charge": 0},
                ]
            },
            {
                "name": "Ethers",
                "molecules": [
                    {"name": "DME", "label": "DME", "smiles": "COCCOC", "charge": 0},
                    {"name": "DEGDME", "label": "DEGDME", "smiles": "COCCOCCOC", "charge": 0},
                    {"name": "TEGDME", "label": "TEGDME", "smiles": "COCCOCCOCCOCCOC", "charge": 0},
                    {"name": "DOL", "label": "DOL", "smiles": "C1COCO1", "charge": 0},
                ]
            },
            {
                "name": "Ionic Liquids",
                "molecules": [
                    {"name": "EMIm", "label": "EMIm", "smiles": "CC[n+]1ccn(C)c1", "charge": 1},
                    {"name": "BMIm", "label": "BMIm", "smiles": "CCCC[n+]1ccn(C)c1", "charge": 1},
                    {"name": "Pyr13", "label": "Pyr13", "smiles": "CCC[N+]1(C)CCCC1", "charge": 1},
                    {"name": "TEA", "label": "TEA", "smiles": "CC[N+](CC)(CC)CC", "charge": 1},
                ]
            },
        ]
    }


@router.get("/config/custom-solvent-params")
def get_custom_solvent_params_info():
    """获取自定义溶剂参数说明"""
    return {
        "description": "SMD溶剂模型需要以下7个参数来定义自定义溶剂",
        "parameters": [
            {
                "name": "eps",
                "label": "介电常数 (ε)",
                "description": "静态介电常数，反映溶剂极性",
                "example": 78.3553,
                "unit": "无量纲",
                "range": "1.0 - 200.0"
            },
            {
                "name": "eps_inf",
                "label": "光学介电常数 (n²)",
                "description": "折射率的平方，用于非平衡溶剂化",
                "example": 1.778,
                "unit": "无量纲",
                "range": "1.0 - 5.0"
            },
            {
                "name": "hbond_acidity",
                "label": "氢键酸度 (α)",
                "description": "Abraham氢键酸度参数",
                "example": 0.82,
                "unit": "无量纲",
                "range": "0.0 - 1.0"
            },
            {
                "name": "hbond_basicity",
                "label": "氢键碱度 (β)",
                "description": "Abraham氢键碱度参数",
                "example": 0.35,
                "unit": "无量纲",
                "range": "0.0 - 1.0"
            },
            {
                "name": "surface_tension",
                "label": "表面张力 (γ)",
                "description": "溶剂表面张力",
                "example": 71.99,
                "unit": "cal/mol·Å²",
                "range": "0.0 - 100.0"
            },
            {
                "name": "carbon_aromaticity",
                "label": "芳香碳比例 (φ)",
                "description": "溶剂分子中芳香碳原子的比例",
                "example": 0.0,
                "unit": "无量纲",
                "range": "0.0 - 1.0"
            },
            {
                "name": "halogenicity",
                "label": "卤素比例 (ψ)",
                "description": "溶剂分子中F、Cl、Br原子的比例",
                "example": 0.0,
                "unit": "无量纲",
                "range": "0.0 - 1.0"
            },
        ],
        "example_solvents": {
            "water": {
                "eps": 78.3553,
                "eps_inf": 1.778,
                "hbond_acidity": 0.82,
                "hbond_basicity": 0.35,
                "surface_tension": 71.99,
                "carbon_aromaticity": 0.0,
                "halogenicity": 0.0
            },
            "acetonitrile": {
                "eps": 35.688,
                "eps_inf": 1.806,
                "hbond_acidity": 0.07,
                "hbond_basicity": 0.32,
                "surface_tension": 41.25,
                "carbon_aromaticity": 0.0,
                "halogenicity": 0.0
            }
        }
    }


# ============================================================================
# Spin Multiplicity Calculation
# ============================================================================

@router.post("/calculate-spin")
def calculate_spin_multiplicity(
    smiles: str = Query(..., description="SMILES表达式"),
    charge: int = Query(default=0, description="分子电荷")
):
    """
    根据SMILES和电荷自动计算自旋多重度

    自旋多重度 = 2S + 1，其中S是总自旋量子数
    对于闭壳层分子，S=0，自旋多重度=1
    对于自由基，S=0.5，自旋多重度=2
    """
    try:
        from rdkit import Chem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES")

        # 添加氢原子以获得完整的电子数
        mol = Chem.AddHs(mol)

        # 计算总电子数
        total_electrons = 0
        for atom in mol.GetAtoms():
            total_electrons += atom.GetAtomicNum()

        # 减去电荷（正电荷减少电子，负电荷增加电子）
        total_electrons -= charge

        # 计算未配对电子数
        # 首先检查是否有自由基
        num_radical_electrons = 0
        for atom in mol.GetAtoms():
            num_radical_electrons += atom.GetNumRadicalElectrons()

        # 如果没有显式自由基，根据电子数判断
        if num_radical_electrons == 0:
            # 偶数电子通常是闭壳层（单重态）
            # 奇数电子是双重态
            if total_electrons % 2 == 0:
                spin_multiplicity = 1
            else:
                spin_multiplicity = 2
        else:
            # 有自由基电子
            spin_multiplicity = num_radical_electrons + 1

        return {
            "smiles": smiles,
            "charge": charge,
            "total_electrons": total_electrons,
            "num_radical_electrons": num_radical_electrons,
            "spin_multiplicity": spin_multiplicity,
            "description": f"自旋多重度 = {spin_multiplicity} ({'单重态' if spin_multiplicity == 1 else '双重态' if spin_multiplicity == 2 else '三重态' if spin_multiplicity == 3 else f'{spin_multiplicity}重态'})"
        }
    except ImportError:
        raise HTTPException(status_code=500, detail="RDKit not installed")
    except Exception as e:
        logger.error(f"Error calculating spin multiplicity: {e}")
        raise HTTPException(status_code=400, detail=str(e))
