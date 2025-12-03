"""
电解液研发 API 端点
用于搜索和展示已完成的电解液计算结果
"""
from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from sqlalchemy import and_, or_
from datetime import datetime
from typing import List, Optional
from app.database import get_db
from app.models.user import User, UserRole
from app.models.electrolyte import ElectrolyteSystem
from app.models.job import MDJob, JobStatus, DataVisibility
from app.models.result import RDFResult, MSDResult, SolvationStructure
from app.dependencies import get_current_active_user, get_optional_current_user
from app.core.logger import logger

router = APIRouter()


def _search_electrolytes_impl(
    cations: Optional[List[str]],
    anions: Optional[List[str]],
    solvents: Optional[List[str]],
    solvent_smiles: Optional[str],
    temp_min: Optional[float],
    temp_max: Optional[float],
    skip: int,
    limit: int,
    db: Session,
    user_id: Optional[int] = None,  # 如果提供，则只搜索该用户的任务
    current_user: Optional[User] = None,  # 当前登录用户（用于可见性过滤）
    public_only: bool = False  # 是否只搜索公开数据
):
    """
    搜索已完成的电解液计算结果

    支持按以下条件搜索：
    - 阴阳离子（多选）
    - 溶剂名称（多选）
    - 溶剂 SMILES（精确搜索）
    - 温度范围

    只返回状态为 COMPLETED 的任务
    可见性规则：
    - 管理员可以看到所有数据
    - 普通用户可以看到：自己的所有数据 + 公开数据 + 已过延期的数据
    - 游客只能看到公开数据
    """
    try:
        # 查询已完成的任务
        query = db.query(MDJob, ElectrolyteSystem, User).join(
            ElectrolyteSystem,
            MDJob.system_id == ElectrolyteSystem.id
        ).join(
            User,
            MDJob.user_id == User.id
        ).filter(
            MDJob.status == JobStatus.COMPLETED
        )

        # 如果指定了用户ID，只搜索该用户的任务（用于"我的数据"）
        if user_id is not None:
            query = query.filter(MDJob.user_id == user_id)
        elif public_only:
            # 公开搜索：应用可见性过滤
            now = datetime.utcnow()
            if current_user and current_user.role == UserRole.ADMIN:
                # 管理员可以看到所有数据
                pass
            elif current_user:
                # 登录用户：自己的数据 + 公开数据 + 已过延期的数据
                query = query.filter(
                    or_(
                        MDJob.user_id == current_user.id,  # 自己的数据
                        MDJob.visibility == DataVisibility.PUBLIC,  # 公开数据
                        and_(  # 已过延期的数据
                            MDJob.visibility == DataVisibility.DELAYED,
                            MDJob.visibility_delay_until <= now
                        )
                    )
                )
            else:
                # 游客：只能看到公开数据和已过延期的数据
                query = query.filter(
                    or_(
                        MDJob.visibility == DataVisibility.PUBLIC,
                        and_(
                            MDJob.visibility == DataVisibility.DELAYED,
                            MDJob.visibility_delay_until <= now
                        )
                    )
                )
        
        # 构建过滤条件
        filters = []
        
        # 阳离子过滤
        if cations and len(cations) > 0:
            cation_filters = []
            for cation in cations:
                # 检查 cations JSON 数组中是否包含该阳离子
                cation_filters.append(
                    ElectrolyteSystem.cations.op('@>')(f'[{{"name": "{cation}"}}]')
                )
            filters.append(or_(*cation_filters))
        
        # 阴离子过滤
        if anions and len(anions) > 0:
            anion_filters = []
            for anion in anions:
                anion_filters.append(
                    ElectrolyteSystem.anions.op('@>')(f'[{{"name": "{anion}"}}]')
                )
            filters.append(or_(*anion_filters))
        
        # 溶剂过滤 - 支持按名称搜索
        if solvents and len(solvents) > 0:
            solvent_filters = []
            for solvent in solvents:
                # 同时支持按 name 或 smiles 搜索
                solvent_filters.append(
                    or_(
                        ElectrolyteSystem.solvents.op('@>')(f'[{{"name": "{solvent}"}}]'),
                        ElectrolyteSystem.solvents.op('@>')(f'[{{"smiles": "{solvent}"}}]')
                    )
                )
            filters.append(or_(*solvent_filters))

        # 溶剂 SMILES 精确搜索（单独的输入框）
        if solvent_smiles and solvent_smiles.strip():
            smiles_clean = solvent_smiles.strip()
            filters.append(
                ElectrolyteSystem.solvents.op('@>')(f'[{{"smiles": "{smiles_clean}"}}]')
            )

        # 温度范围过滤
        if temp_min is not None:
            filters.append(ElectrolyteSystem.temperature >= temp_min)
        if temp_max is not None:
            filters.append(ElectrolyteSystem.temperature <= temp_max)
        
        # 应用所有过滤条件
        if filters:
            query = query.filter(and_(*filters))
        
        # 按创建时间倒序排列
        query = query.order_by(MDJob.created_at.desc())
        
        # 分页
        total = query.count()
        results = query.offset(skip).limit(limit).all()

        # 构建响应数据
        data = []
        for job, system, owner in results:
            # 获取相关的分析结果统计
            rdf_count = db.query(RDFResult).filter(RDFResult.md_job_id == job.id).count()
            msd_count = db.query(MSDResult).filter(MSDResult.md_job_id == job.id).count()
            solvation_count = db.query(SolvationStructure).filter(
                SolvationStructure.md_job_id == job.id
            ).count()

            # 处理匿名公开
            is_anonymous = getattr(job, 'anonymous_public', False)
            is_own_data = current_user and job.user_id == current_user.id
            is_admin = current_user and current_user.role == UserRole.ADMIN

            # 确定是否显示用户信息
            show_user_info = is_own_data or is_admin or not is_anonymous

            # 优先使用 job.config 中保存的 system_snapshot（创建任务时的配方快照）
            # 这样可以避免后续修改配方影响历史任务的显示
            snapshot = job.config.get("system_snapshot") if job.config else None

            # 从快照或当前系统数据中获取配方信息
            if snapshot:
                system_name = snapshot.get("name", system.name)
                cations = snapshot.get("cations", system.cations)
                anions = snapshot.get("anions", system.anions)
                solvents = snapshot.get("solvents", system.solvents)
                temperature = snapshot.get("temperature", system.temperature)
                pressure = snapshot.get("pressure", system.pressure)
            else:
                # 兼容旧任务（没有快照的情况）
                system_name = system.name
                cations = system.cations
                anions = system.anions
                solvents = system.solvents
                temperature = system.temperature
                pressure = system.pressure

            # 获取任务配置中的额外信息
            charge_method = job.config.get("charge_method") if job.config else None
            qc_enabled = job.config.get("qc_enabled", False) if job.config else False
            user_note = job.config.get("user_note") if job.config else None

            item = {
                "job_id": job.id,
                "job_name": job.config.get("job_name") if job.config else None,
                "user_note": user_note,  # 用户备注
                "system_id": system.id,
                "system_name": system_name,
                "cations": cations,
                "anions": anions,
                "solvents": solvents,
                "temperature": temperature,
                "pressure": pressure,
                "density": system.density,
                "concentration": system.concentration,
                # 任务配置信息
                "charge_method": charge_method,  # 电荷计算方式
                "qc_enabled": qc_enabled,  # 是否有QC计算
                "created_at": job.created_at.isoformat() if job.created_at else None,
                "finished_at": job.finished_at.isoformat() if job.finished_at else None,
                # 分析结果统计
                "has_rdf": rdf_count > 0,
                "has_msd": msd_count > 0,
                "has_solvation": solvation_count > 0,
                "rdf_count": rdf_count,
                "msd_count": msd_count,
                "solvation_count": solvation_count,
                # 可见性信息
                "visibility": job.visibility.value if job.visibility else "DELAYED",
                "allow_download": getattr(job, 'allow_download', True),
                "view_count": getattr(job, 'view_count', 0),
                "is_own_data": is_own_data,
            }

            # 用户信息（根据匿名设置）
            if show_user_info:
                item["username"] = owner.username
                item["organization"] = owner.organization
            else:
                item["username"] = "匿名用户"
                item["organization"] = "匿名单位"

            data.append(item)

        return {
            "total": total,
            "skip": skip,
            "limit": limit,
            "data": data
        }
    
    except Exception as e:
        logger.error(f"Failed to search electrolytes: {e}")
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/search")
def search_electrolytes(
    cations: Optional[List[str]] = Query(None, description="阳离子列表（多选）"),
    anions: Optional[List[str]] = Query(None, description="阴离子列表（多选）"),
    solvents: Optional[List[str]] = Query(None, description="溶剂名称列表（多选）"),
    solvent_smiles: Optional[str] = Query(None, description="溶剂 SMILES（精确搜索）"),
    temp_min: Optional[float] = Query(None, description="最低温度 (K)"),
    temp_max: Optional[float] = Query(None, description="最高温度 (K)"),
    skip: int = 0,
    limit: int = 100,
    db: Session = Depends(get_db),
    current_user: Optional[User] = Depends(get_optional_current_user)
):
    """
    搜索已完成的电解液计算结果（公开搜索）

    可见性规则：
    - 管理员可以看到所有数据
    - 登录用户可以看到：自己的所有数据 + 公开数据 + 已过延期的数据
    - 游客只能看到公开数据和已过延期的数据

    支持按以下条件搜索：
    - 阳离子名称（多选）
    - 阴离子名称（多选）
    - 溶剂名称（多选）
    - 溶剂 SMILES（精确搜索，唯一标识符）
    - 温度范围
    """
    return _search_electrolytes_impl(
        cations=cations,
        anions=anions,
        solvents=solvents,
        solvent_smiles=solvent_smiles,
        temp_min=temp_min,
        temp_max=temp_max,
        skip=skip,
        limit=limit,
        db=db,
        user_id=None,
        current_user=current_user,
        public_only=True  # 应用可见性过滤
    )


@router.get("/search/my")
def search_my_electrolytes(
    cations: Optional[List[str]] = Query(None, description="阳离子列表（多选）"),
    anions: Optional[List[str]] = Query(None, description="阴离子列表（多选）"),
    solvents: Optional[List[str]] = Query(None, description="溶剂名称列表（多选）"),
    solvent_smiles: Optional[str] = Query(None, description="溶剂 SMILES（精确搜索）"),
    temp_min: Optional[float] = Query(None, description="最低温度 (K)"),
    temp_max: Optional[float] = Query(None, description="最高温度 (K)"),
    skip: int = 0,
    limit: int = 100,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    搜索当前用户的已完成电解液计算结果（工作台内使用）

    只返回当前登录用户创建的已完成任务（不受可见性限制）
    """
    return _search_electrolytes_impl(
        cations=cations,
        anions=anions,
        solvents=solvents,
        solvent_smiles=solvent_smiles,
        temp_min=temp_min,
        temp_max=temp_max,
        skip=skip,
        limit=limit,
        db=db,
        user_id=current_user.id,  # 只搜索当前用户的任务
        current_user=current_user,
        public_only=False  # 不应用可见性过滤（自己的数据都可见）
    )


@router.post("/job/{job_id}/view")
def record_job_view(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: Optional[User] = Depends(get_optional_current_user)
):
    """记录任务查看次数（用于奖励计算）"""
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="任务不存在")

    # 不记录自己查看自己的数据
    if current_user and job.user_id == current_user.id:
        return {"message": "ok", "view_count": job.view_count or 0}

    # 增加查看次数
    job.view_count = (job.view_count or 0) + 1

    # 给数据所有者发放查看奖励（每次 0.1 核时，每日上限 5 核时）
    # TODO: 实现每日上限检查
    owner = db.query(User).filter(User.id == job.user_id).first()
    if owner and job.visibility == DataVisibility.PUBLIC:
        owner.contribution_points = (owner.contribution_points or 0) + 0.1
        # 暂时不发放核时奖励，只记录贡献积分

    db.commit()

    return {"message": "ok", "view_count": job.view_count}


@router.get("/available-options")
def get_available_search_options(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取数据管理页面的可用搜索选项（从数据库中实际数据提取）

    返回当前用户已完成任务中使用过的：
    - 阳离子列表
    - 阴离子列表
    - 溶剂列表
    """
    # 查询当前用户的已完成任务
    completed_jobs = db.query(ElectrolyteSystem).join(
        MDJob, ElectrolyteSystem.id == MDJob.system_id
    ).filter(
        MDJob.user_id == current_user.id,
        MDJob.status == JobStatus.COMPLETED,
        ElectrolyteSystem.is_deleted == False
    ).all()

    # 提取所有使用过的离子和溶剂
    cations_set = set()
    anions_set = set()
    solvents_set = set()

    for system in completed_jobs:
        # 提取阳离子
        if system.cations:
            for cation in system.cations:
                if isinstance(cation, dict) and 'name' in cation:
                    cations_set.add(cation['name'])

        # 提取阴离子
        if system.anions:
            for anion in system.anions:
                if isinstance(anion, dict) and 'name' in anion:
                    anions_set.add(anion['name'])

        # 提取溶剂
        if system.solvents:
            for solvent in system.solvents:
                if isinstance(solvent, dict) and 'name' in solvent:
                    solvents_set.add(solvent['name'])

    # 转换为排序列表
    cations = sorted(list(cations_set))
    anions = sorted(list(anions_set))
    solvents = sorted(list(solvents_set))

    logger.info(f"Available options for user {current_user.username}: {len(cations)} cations, {len(anions)} anions, {len(solvents)} solvents")

    return {
        "cations": cations,
        "anions": anions,
        "solvents": solvents
    }
