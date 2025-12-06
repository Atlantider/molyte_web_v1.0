"""
QC 任务复用工具函数

确保只有真正有能量结果的任务才能被复用，
并正确处理复用链条和结果复制。
"""
import logging
from typing import Optional, Tuple, Set
from sqlalchemy.orm import Session, joinedload

from app.models.qc import QCJob, QCResult, QCJobStatus

logger = logging.getLogger(__name__)


def find_root_job_with_result(
    db: Session, 
    job: QCJob, 
    max_depth: int = 50
) -> Optional[QCJob]:
    """
    追溯复用链，找到有实际能量结果的根任务
    
    Args:
        db: 数据库会话
        job: 起始任务
        max_depth: 最大追溯深度，防止无限循环
        
    Returns:
        有能量结果的根任务，如果链条断裂返回 None
    """
    visited: Set[int] = set()
    current = job
    depth = 0
    
    while current and depth < max_depth:
        if current.id in visited:
            logger.warning(f"检测到复用链循环引用: job {current.id}")
            return None
        visited.add(current.id)
        depth += 1
        
        # 检查当前任务是否有能量结果
        result = db.query(QCResult).filter(
            QCResult.qc_job_id == current.id,
            QCResult.energy_au.isnot(None)
        ).first()
        
        if result:
            logger.debug(f"找到根任务 {current.id}，能量={result.energy_au}")
            return current
        
        # 继续向上追溯
        if current.reused_from_job_id:
            current = db.query(QCJob).filter(
                QCJob.id == current.reused_from_job_id
            ).first()
            if current is None:
                logger.warning(f"复用链断裂: reused_from_job_id={job.reused_from_job_id} 不存在")
                return None
        else:
            # 到达链条末端但没有结果
            logger.warning(f"到达复用链末端但无结果: job {current.id}")
            return None
    
    logger.warning(f"复用链过深(>{max_depth}): 起始 job {job.id}")
    return None


def get_energy_from_job(db: Session, job: QCJob) -> Optional[float]:
    """
    获取任务的能量结果，支持复用任务的链式追溯
    
    Args:
        db: 数据库会话
        job: QC 任务
        
    Returns:
        能量值 (Hartree)，如果无法获取返回 None
    """
    # 首先检查任务自身是否有结果
    result = db.query(QCResult).filter(
        QCResult.qc_job_id == job.id,
        QCResult.energy_au.isnot(None)
    ).first()
    
    if result:
        return result.energy_au
    
    # 如果是复用任务，追溯根任务
    if job.is_reused and job.reused_from_job_id:
        root_job = find_root_job_with_result(db, job)
        if root_job:
            root_result = db.query(QCResult).filter(
                QCResult.qc_job_id == root_job.id,
                QCResult.energy_au.isnot(None)
            ).first()
            if root_result:
                return root_result.energy_au
    
    return None


def validate_job_for_reuse(db: Session, job: QCJob) -> Tuple[bool, Optional[QCJob], str]:
    """
    验证任务是否可以被复用
    
    Args:
        db: 数据库会话
        job: 待验证的 QC 任务
        
    Returns:
        (是否可复用, 根任务, 原因说明)
    """
    if job.status != QCJobStatus.COMPLETED:
        return False, None, f"任务状态为 {job.status}，不是 COMPLETED"
    
    if job.is_deleted:
        return False, None, "任务已被删除"
    
    # 检查是否有直接结果
    result = db.query(QCResult).filter(
        QCResult.qc_job_id == job.id,
        QCResult.energy_au.isnot(None)
    ).first()
    
    if result:
        return True, job, "任务有直接能量结果"
    
    # 如果是复用任务，追溯根任务
    if job.is_reused and job.reused_from_job_id:
        root_job = find_root_job_with_result(db, job)
        if root_job:
            return True, root_job, f"复用链有效，根任务={root_job.id}"
        else:
            return False, None, "复用链断裂，无法找到有效结果"
    
    return False, None, "任务没有能量结果且不是复用任务"


def copy_result_for_reused_job(
    db: Session, 
    source_job: QCJob, 
    target_job: QCJob
) -> Optional[QCResult]:
    """
    将源任务的结果复制到目标复用任务
    
    Args:
        db: 数据库会话
        source_job: 源任务（有结果的根任务）
        target_job: 目标任务（复用任务）
        
    Returns:
        新创建的 QCResult，如果失败返回 None
    """
    # 获取源任务的结果
    source_result = db.query(QCResult).filter(
        QCResult.qc_job_id == source_job.id,
        QCResult.energy_au.isnot(None)
    ).first()
    
    if not source_result:
        logger.error(f"源任务 {source_job.id} 没有有效结果")
        return None
    
    # 检查目标任务是否已有结果
    existing = db.query(QCResult).filter(
        QCResult.qc_job_id == target_job.id
    ).first()
    
    if existing:
        logger.info(f"目标任务 {target_job.id} 已有结果，跳过复制")
        return existing

    # 创建新的结果记录
    new_result = QCResult(
        qc_job_id=target_job.id,
        smiles=source_result.smiles,
        energy_au=source_result.energy_au,
        homo=source_result.homo,
        lumo=source_result.lumo,
        homo_lumo_gap=source_result.homo_lumo_gap,
        esp_min_kcal=source_result.esp_min_kcal,
        esp_max_kcal=source_result.esp_max_kcal,
        dipole_moment=source_result.dipole_moment,
        polarizability=source_result.polarizability,
        vip_ev=source_result.vip_ev,
        vea_ev=source_result.vea_ev,
        oxidation_potential_v=source_result.oxidation_potential_v,
        reduction_potential_v=source_result.reduction_potential_v,
        # 图片和文件路径不复制，保持引用源任务的结果
        additional_properties={
            "copied_from_job_id": source_job.id,
            "copied_from_result_id": source_result.id,
            "is_copied": True
        }
    )

    db.add(new_result)
    db.flush()

    logger.info(f"已将任务 {source_job.id} 的结果复制到任务 {target_job.id}")
    return new_result


def fix_reused_job_without_result(db: Session, job: QCJob) -> Tuple[bool, str]:
    """
    修复没有结果的复用任务

    如果能找到有效的根任务，复制结果；
    否则将任务状态重置为待计算。

    Args:
        db: 数据库会话
        job: 需要修复的任务

    Returns:
        (是否成功, 说明)
    """
    if not job.is_reused:
        return False, "不是复用任务"

    # 检查是否已有结果
    existing_result = db.query(QCResult).filter(
        QCResult.qc_job_id == job.id,
        QCResult.energy_au.isnot(None)
    ).first()

    if existing_result:
        return True, "任务已有结果，无需修复"

    # 尝试追溯根任务
    root_job = find_root_job_with_result(db, job)

    if root_job:
        # 找到根任务，复制结果
        new_result = copy_result_for_reused_job(db, root_job, job)
        if new_result:
            # 更新 reused_from_job_id 直接指向根任务
            job.reused_from_job_id = root_job.id
            return True, f"已从根任务 {root_job.id} 复制结果"
        else:
            return False, "复制结果失败"
    else:
        # 复用链断裂，重置任务状态
        job.status = QCJobStatus.SUBMITTED
        job.is_reused = False
        job.reused_from_job_id = None
        job.error_message = "复用链无效，任务需要重新计算"
        return True, "复用链断裂，已重置为待计算状态"


def batch_fix_invalid_reused_jobs(db: Session, dry_run: bool = True) -> dict:
    """
    批量修复无效的复用任务

    Args:
        db: 数据库会话
        dry_run: 是否仅检测不修改

    Returns:
        修复统计信息
    """
    from sqlalchemy import and_

    stats = {
        "total_reused_jobs": 0,
        "valid_jobs": 0,
        "fixed_by_copy": 0,
        "reset_to_submitted": 0,
        "errors": []
    }

    # 查找所有复用任务
    reused_jobs = db.query(QCJob).filter(
        QCJob.is_reused == True,
        QCJob.status == QCJobStatus.COMPLETED,
        QCJob.is_deleted == False
    ).all()

    stats["total_reused_jobs"] = len(reused_jobs)

    for job in reused_jobs:
        # 检查是否有结果
        has_result = db.query(QCResult).filter(
            QCResult.qc_job_id == job.id,
            QCResult.energy_au.isnot(None)
        ).count() > 0

        if has_result:
            stats["valid_jobs"] += 1
            continue

        # 尝试修复
        if not dry_run:
            success, msg = fix_reused_job_without_result(db, job)
            if success:
                if "复制结果" in msg:
                    stats["fixed_by_copy"] += 1
                elif "重置" in msg:
                    stats["reset_to_submitted"] += 1
            else:
                stats["errors"].append({"job_id": job.id, "error": msg})
        else:
            # 仅检测
            root_job = find_root_job_with_result(db, job)
            if root_job:
                stats["fixed_by_copy"] += 1
            else:
                stats["reset_to_submitted"] += 1

    if not dry_run:
        db.commit()

    return stats

