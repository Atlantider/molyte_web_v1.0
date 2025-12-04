"""
热力学循环计算氧化还原电位 API

⚠️ 高风险警告：
- 结果对方法/基组/溶剂模型/构型高度敏感
- 计算量大，经常不收敛
- 数值可能存在数百 mV 的系统性偏差
- 仅供研究参考，不应作为定量预测
"""
from datetime import datetime
from typing import Optional, List
from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from app.database import get_db
from app.dependencies import get_current_active_user
from app.models.user import User
from app.models.job import (
    RedoxPotentialJob, RedoxJobStatus,
    ReorganizationEnergyJob, ReorgEnergyJobStatus
)
from app.schemas.redox import (
    RedoxJobCreate, RedoxJobResponse, RedoxJobListResponse,
    ReorgEnergyJobCreate, ReorgEnergyJobResponse
)

router = APIRouter(prefix="/redox", tags=["redox"])


# ============================================================================
# 热力学循环计算氧化还原电位
# ============================================================================

@router.post("/jobs", response_model=RedoxJobResponse)
def create_redox_job(
    job_data: RedoxJobCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    创建热力学循环计算氧化还原电位任务

    ⚠️ 警告：
    - 结果对方法/基组/溶剂模型高度敏感
    - 计算量大，可能不收敛
    - 建议物种数量不超过 10 个
    """
    # 检查物种数量限制
    species_count = len(job_data.config.species_list)
    if species_count > 20:
        raise HTTPException(
            status_code=400,
            detail=f"物种数量 ({species_count}) 超过限制 (最多 20 个)"
        )

    if species_count > 10:
        # 警告但允许
        pass

    # 创建任务
    job = RedoxPotentialJob(
        md_job_id=job_data.md_job_id,
        user_id=current_user.id,
        status=RedoxJobStatus.CREATED,
        config=job_data.config.model_dump(),
    )

    db.add(job)
    db.commit()
    db.refresh(job)

    return job


@router.get("/jobs", response_model=RedoxJobListResponse)
def list_redox_jobs(
    md_job_id: Optional[int] = Query(None, description="按 MD 任务 ID 筛选"),
    status: Optional[str] = Query(None, description="按状态筛选"),
    skip: int = Query(0, ge=0),
    limit: int = Query(20, ge=1, le=100),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """列出热力学循环任务"""
    query = db.query(RedoxPotentialJob).filter(
        RedoxPotentialJob.user_id == current_user.id
    )

    if md_job_id:
        query = query.filter(RedoxPotentialJob.md_job_id == md_job_id)

    if status:
        try:
            status_enum = RedoxJobStatus(status)
            query = query.filter(RedoxPotentialJob.status == status_enum)
        except ValueError:
            pass

    total = query.count()
    jobs = query.order_by(RedoxPotentialJob.created_at.desc()).offset(skip).limit(limit).all()

    return RedoxJobListResponse(total=total, jobs=jobs)


@router.get("/jobs/{job_id}", response_model=RedoxJobResponse)
def get_redox_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取热力学循环任务详情"""
    job = db.query(RedoxPotentialJob).filter(
        RedoxPotentialJob.id == job_id,
        RedoxPotentialJob.user_id == current_user.id
    ).first()

    if not job:
        raise HTTPException(status_code=404, detail="任务不存在")

    return job


@router.post("/jobs/{job_id}/submit", response_model=RedoxJobResponse)
def submit_redox_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """提交热力学循环任务到计算队列"""
    job = db.query(RedoxPotentialJob).filter(
        RedoxPotentialJob.id == job_id,
        RedoxPotentialJob.user_id == current_user.id
    ).first()

    if not job:
        raise HTTPException(status_code=404, detail="任务不存在")

    if job.status != RedoxJobStatus.CREATED:
        raise HTTPException(
            status_code=400,
            detail=f"任务状态为 {job.status.value}，无法提交"
        )



@router.delete("/jobs/{job_id}")
def delete_redox_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """删除热力学循环任务"""
    job = db.query(RedoxPotentialJob).filter(
        RedoxPotentialJob.id == job_id,
        RedoxPotentialJob.user_id == current_user.id
    ).first()

    if not job:
        raise HTTPException(status_code=404, detail="任务不存在")

    if job.status == RedoxJobStatus.RUNNING:
        raise HTTPException(
            status_code=400,
            detail="任务正在运行中，无法删除"
        )

    db.delete(job)
    db.commit()

    return {"message": "任务已删除", "job_id": job_id}


# ============================================================================
# 重组能计算 (Marcus 理论)
# ============================================================================

@router.post("/reorganization-energy/jobs", response_model=ReorgEnergyJobResponse)
def create_reorg_energy_job(
    job_data: ReorgEnergyJobCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    创建重组能计算任务 (Marcus 理论)

    ⚠️ 极高风险警告：
    - 每个物种至少 2 次优化 + 4 次单点
    - Cluster 体系极易不收敛
    - 构型依赖极强
    - 强烈建议物种数量不超过 5 个
    """
    species_count = len(job_data.config.species_list)
    if species_count > 10:
        raise HTTPException(
            status_code=400,
            detail=f"物种数量 ({species_count}) 超过限制 (最多 10 个)"
        )

    if species_count > 5:
        # 警告但允许（风险由用户承担）
        pass

    job = ReorganizationEnergyJob(
        md_job_id=job_data.md_job_id,
        user_id=current_user.id,
        status=ReorgEnergyJobStatus.CREATED,
        config=job_data.config.model_dump(),
    )

    db.add(job)
    db.commit()
    db.refresh(job)

    return job


@router.get("/reorganization-energy/jobs", response_model=List[ReorgEnergyJobResponse])
def list_reorg_energy_jobs(
    md_job_id: Optional[int] = Query(None, description="按 MD 任务 ID 筛选"),
    status: Optional[str] = Query(None, description="按状态筛选"),
    skip: int = Query(0, ge=0),
    limit: int = Query(20, ge=1, le=100),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """列出重组能计算任务"""
    query = db.query(ReorganizationEnergyJob).filter(
        ReorganizationEnergyJob.user_id == current_user.id
    )

    if md_job_id:
        query = query.filter(ReorganizationEnergyJob.md_job_id == md_job_id)

    if status:
        query = query.filter(ReorganizationEnergyJob.status == status)

    jobs = query.order_by(ReorganizationEnergyJob.created_at.desc()).offset(skip).limit(limit).all()

    return jobs


@router.get("/reorganization-energy/jobs/{job_id}", response_model=ReorgEnergyJobResponse)
def get_reorg_energy_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取重组能任务详情"""
    job = db.query(ReorganizationEnergyJob).filter(
        ReorganizationEnergyJob.id == job_id,
        ReorganizationEnergyJob.user_id == current_user.id
    ).first()

    if not job:
        raise HTTPException(status_code=404, detail="任务不存在")

    return job


@router.post("/reorganization-energy/jobs/{job_id}/submit", response_model=ReorgEnergyJobResponse)
def submit_reorg_energy_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """提交重组能任务到计算队列"""
    job = db.query(ReorganizationEnergyJob).filter(
        ReorganizationEnergyJob.id == job_id,
        ReorganizationEnergyJob.user_id == current_user.id
    ).first()

    if not job:
        raise HTTPException(status_code=404, detail="任务不存在")

    if job.status != ReorgEnergyJobStatus.CREATED:
        raise HTTPException(
            status_code=400,
            detail=f"任务状态为 {job.status.value}，无法提交"
        )

    job.status = ReorgEnergyJobStatus.SUBMITTED
    db.commit()
    db.refresh(job)

    return job


@router.delete("/reorganization-energy/jobs/{job_id}")
def delete_reorg_energy_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """删除重组能任务"""
    job = db.query(ReorganizationEnergyJob).filter(
        ReorganizationEnergyJob.id == job_id,
        ReorganizationEnergyJob.user_id == current_user.id
    ).first()

    if not job:
        raise HTTPException(status_code=404, detail="任务不存在")

    if job.status == ReorgEnergyJobStatus.RUNNING:
        raise HTTPException(
            status_code=400,
            detail="任务正在运行中，无法删除"
        )

    db.delete(job)
    db.commit()

    return {"message": "任务已删除", "job_id": job_id}


# ============================================================================
# 物理常数和参考值
# ============================================================================

@router.get("/constants")
def get_physical_constants():
    """获取物理常数和参考值"""
    from app.schemas.redox import PhysicalConstants

    return {
        "faraday_c_mol": PhysicalConstants.FARADAY_C_MOL,
        "faraday_kcal_mol_v": PhysicalConstants.FARADAY_KCAL_MOL_V,
        "hartree_to_kcal": PhysicalConstants.HARTREE_TO_KCAL,
        "hartree_to_ev": PhysicalConstants.HARTREE_TO_EV,
        "li_absolute_potential_vs_she": PhysicalConstants.LI_ABSOLUTE_POTENTIAL_VS_SHE,
        "she_absolute_potential": PhysicalConstants.SHE_ABSOLUTE_POTENTIAL,
        "notes": {
            "li_reference": "Li+/Li 绝对电位文献有差异 (-2.9 ~ -3.1 V vs SHE)，默认使用 Trasatti 推荐值 -3.04 V",
            "she_reference": "SHE 绝对电位文献有差异 (4.28 ~ 4.44 V vs 真空)，默认使用 4.44 V",
            "warning": "计算结果存在系统性偏差，仅供研究参考"
        }
    }
