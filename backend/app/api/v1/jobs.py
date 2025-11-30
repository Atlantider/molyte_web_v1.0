"""
Job management API endpoints
"""
from fastapi import APIRouter, Depends, HTTPException, status, BackgroundTasks
from sqlalchemy.orm import Session
from sqlalchemy import func
from typing import List, Dict, Any, Tuple
from datetime import datetime, date
from pathlib import Path
import subprocess
from app.database import get_db
from app.models.user import User, UserRole
from app.models.electrolyte import ElectrolyteSystem
from app.models.job import MDJob, JobStatus
from app.models.qc import QCJob, QCJobStatus
from app.schemas.job import MDJob as MDJobSchema, MDJobCreate, MDJobUpdate, BatchMDJobCreate
from app.dependencies import get_current_active_user
from app.core.logger import logger
from app.core.config import settings
from app.utils.quota import check_user_quota, update_user_usage_stats
from app.workers.molyte_wrapper import MolyteWrapper
from app.workers.molyte_adapter import convert_electrolyte_to_molyte_format
from app.schemas.accuracy_level import (
    AccuracyLevel,
    apply_accuracy_level,
    get_all_accuracy_levels
)

router = APIRouter()


def check_job_permission(job: MDJob, current_user: User, allow_premium: bool = False):
    """
    Check if user has permission to access a job

    Args:
        job: MD job
        current_user: Current user
        allow_premium: Whether to allow PREMIUM users (default: False)

    Raises:
        HTTPException: If user doesn't have permission
    """
    if job.user_id == current_user.id:
        return  # Owner always has permission

    if current_user.role == UserRole.ADMIN:
        return  # Admin always has permission

    if allow_premium and current_user.role == UserRole.PREMIUM:
        return  # Premium user has permission if allowed

    raise HTTPException(
        status_code=status.HTTP_403_FORBIDDEN,
        detail="Not enough permissions"
    )


def _submit_job_to_cluster(job: MDJob, electrolyte: ElectrolyteSystem, db: Session) -> None:
    """
    Internal function to submit a job to Slurm cluster

    Args:
        job: MDJob instance
        electrolyte: ElectrolyteSystem instance
        db: Database session

    Raises:
        Exception: If submission fails
    """
    # Convert to molyte format
    job_data = convert_electrolyte_to_molyte_format(
        job_name=job.config.get('job_name', f'MD-{job.id}'),
        job_config=job.config,
        electrolyte_data={
            "name": electrolyte.name,
            "cations": electrolyte.cations,
            "anions": electrolyte.anions,
            "solvents": electrolyte.solvents,
            "additives": getattr(electrolyte, 'additives', None),  # Optional field
            "box_size": electrolyte.box_size,  # 从电解质模型获取盒子大小
            "temperature": electrolyte.temperature,
            "pressure": electrolyte.pressure,
        }
    )

    logger.info(f"Converted job data for {job_data['name']}")

    # Initialize MolyteWrapper
    wrapper = MolyteWrapper(
        work_base_path=settings.MOLYTE_WORK_BASE_PATH,
        initial_salts_path=settings.MOLYTE_INITIAL_SALTS_PATH,
        ligpargen_path=settings.MOLYTE_LIGPARGEN_PATH,
        packmol_path=settings.MOLYTE_PACKMOL_PATH,
        ltemplify_path=settings.MOLYTE_LTEMPLIFY_PATH,
        moltemplate_path=settings.MOLYTE_MOLTEMPLATE_PATH,
        charge_save_path=settings.MOLYTE_CHARGE_SAVE_PATH
    )

    logger.info(f"Generating LAMMPS input files for {job_data['name']}...")

    # Generate LAMMPS input files
    result = wrapper.generate_lammps_input(
        job_data=job_data,
        generate_atom_mapping=True
    )

    if not result["success"]:
        raise Exception(result.get("error", "Unknown error"))

    logger.info(f"LAMMPS input files generated successfully for {job_data['name']}")

    # Submit to Slurm
    work_dir = result["work_dir"]
    job_script = work_dir / "job.sh"

    if not job_script.exists():
        raise Exception("Job script not found")

    logger.info(f"Submitting job to Slurm: {job_script}")

    # Submit using sbatch
    submit_result = subprocess.run(
        ["sbatch", str(job_script)],
        cwd=str(work_dir),
        capture_output=True,
        text=True
    )

    if submit_result.returncode != 0:
        raise Exception(f"Slurm submission failed: {submit_result.stderr}")

    # Parse Slurm job ID from output
    # Expected output: "Submitted batch job 12345"
    slurm_output = submit_result.stdout.strip()
    slurm_job_id = None

    if "Submitted batch job" in slurm_output:
        slurm_job_id = slurm_output.split()[-1]
        logger.info(f"Slurm job ID: {slurm_job_id}")

    # Update job status and fields
    job.status = JobStatus.QUEUED
    job.slurm_job_id = slurm_job_id  # 更新 slurm_job_id 字段
    job.work_dir = str(work_dir)     # 更新 work_dir 字段

    if job.config is None:
        job.config = {}

    # 同时在 config 中保存（向后兼容）
    job.config["slurm_job_id"] = slurm_job_id
    job.config["work_dir"] = str(work_dir)
    job.config["files"] = result.get("files", {})

    db.commit()

    logger.info(f"Job {job.id} submitted successfully:")
    logger.info(f"  Slurm Job ID: {slurm_job_id}")
    logger.info(f"  Work Directory: {work_dir}")
    logger.info(f"  Status: {job.status}")


def generate_job_name(db: Session, electrolyte_name: str, custom_name: str = None) -> str:
    """
    Generate job name with format: MD-YYYYMMDD-全局序号-配方名-自定义名称（可选）

    Examples:
        - Without custom name: MD-20251119-0001-EL-20251119-0001-Li-PF6-EC-DMC
        - With custom name: MD-20251119-0001-EL-20251119-0001-Li-PF6-EC-DMC-高温测试

    Args:
        db: Database session
        electrolyte_name: Electrolyte system name
        custom_name: Optional custom name suffix

    Returns:
        Generated job name
    """
    import re

    # 清理 electrolyte_name 中的特殊字符（空格、斜杠等）
    # 将空格和斜杠替换为连字符，移除其他特殊字符
    clean_electrolyte_name = re.sub(r'[\s/]+', '-', electrolyte_name)
    clean_electrolyte_name = re.sub(r'[^\w\u4e00-\u9fff-]', '', clean_electrolyte_name)

    today = date.today()
    date_str = today.strftime('%Y%m%d')

    # Count ALL jobs created today (global count, not per-user)
    today_start = datetime.combine(today, datetime.min.time())
    today_end = datetime.combine(today, datetime.max.time())

    count = db.query(func.count(MDJob.id)).filter(
        MDJob.created_at >= today_start,
        MDJob.created_at <= today_end
    ).scalar()

    # Next sequential number (starting from 1)
    seq_number = count + 1

    # Format: MD-20251119-0001-electrolyte_name[-custom_name]
    if custom_name and custom_name.strip():
        # 清理 custom_name 中的特殊字符
        clean_custom_name = re.sub(r'[\s/]+', '-', custom_name.strip())
        clean_custom_name = re.sub(r'[^\w\u4e00-\u9fff-]', '', clean_custom_name)
        job_name = f"MD-{date_str}-{seq_number:04d}-{clean_electrolyte_name}-{clean_custom_name}"
    else:
        job_name = f"MD-{date_str}-{seq_number:04d}-{clean_electrolyte_name}"

    return job_name


@router.get("/accuracy-levels")
def get_accuracy_levels(
    current_user: User = Depends(get_current_active_user)
):
    """
    获取所有精度等级配置

    返回格式:
    {
        "fast": {
            "name": "快速模式",
            "description": "...",
            "charge_method": "ligpargen",
            "nsteps_npt": 100000,
            ...
        },
        "standard": {...},
        "accurate": {...}
    }
    """
    return get_all_accuracy_levels()


def check_daily_job_limit(db: Session, user_id: int, user_role: str) -> Tuple[bool, int, int]:
    """
    Check if user has reached daily job creation limit

    Args:
        db: Database session
        user_id: User ID
        user_role: User role (ADMIN, PREMIUM, USER)

    Returns:
        tuple: (can_create, current_count, limit)
    """
    # Set limits based on user role
    # Convert to uppercase for comparison
    role_upper = str(user_role).upper()

    if role_upper == "ADMIN" or role_upper == "USERROLE.ADMIN":
        limit = 100
    elif role_upper == "PREMIUM" or role_upper == "USERROLE.PREMIUM":
        limit = 100
    else:  # regular user
        limit = 10

    # Count jobs created today by this user
    today = date.today()
    today_start = datetime.combine(today, datetime.min.time())
    today_end = datetime.combine(today, datetime.max.time())

    count = db.query(func.count(MDJob.id)).filter(
        MDJob.user_id == user_id,
        MDJob.created_at >= today_start,
        MDJob.created_at <= today_end
    ).scalar()

    can_create = count < limit

    return can_create, count, limit


@router.post("/", response_model=MDJobSchema, status_code=status.HTTP_201_CREATED)
def create_md_job(
    job_data: MDJobCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Create a new MD job
    
    Args:
        job_data: MD job creation data
        db: Database session
        current_user: Current authenticated user
        
    Returns:
        MDJob: Created MD job
        
    Raises:
        HTTPException: If system not found or no permission
    """
    # Check if electrolyte system exists
    system = db.query(ElectrolyteSystem).filter(
        ElectrolyteSystem.id == job_data.system_id
    ).first()
    
    if not system:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Electrolyte system not found"
        )
    
    # Check permission
    if system.project.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not enough permissions"
        )

    # Check user quota (CPU hours, daily job limit, concurrent job limit)
    quota_check = check_user_quota(current_user, db)
    if not quota_check["allowed"]:
        raise HTTPException(
            status_code=status.HTTP_429_TOO_MANY_REQUESTS,
            detail=quota_check["reason"]
        )

    # 生成任务名称
    # job_data.job_name 作为自定义名称后缀（可选）
    custom_name = job_data.job_name if job_data.job_name and job_data.job_name.strip() else None
    job_name = generate_job_name(db, system.name, custom_name)

    # 构建配置参数
    config = {
        "job_name": job_name,
        "nsteps_npt": job_data.nsteps_npt,
        "nsteps_nvt": job_data.nsteps_nvt,
        "timestep": job_data.timestep,
        "temperature": job_data.temperature,
        "pressure": job_data.pressure,
        "freq_trj_npt": job_data.freq_trj_npt,
        "freq_trj_nvt": job_data.freq_trj_nvt,
        "thermo_freq": job_data.thermo_freq,
        # Slurm 资源配置
        "slurm_partition": job_data.slurm_partition or "cpu",
        "slurm_nodes": job_data.slurm_nodes or 1,
        "slurm_ntasks": job_data.slurm_ntasks or 8,
        "slurm_cpus_per_task": job_data.slurm_cpus_per_task or 8,
        "slurm_time": job_data.slurm_time or 7200,
        # QC计算选项
        "qc_enabled": job_data.qc_options.enabled if job_data.qc_options else False,
        "qc_accuracy_level": job_data.qc_options.accuracy_level if job_data.qc_options and job_data.qc_options.enabled else None,
        # 支持多选泛函、基组、溶剂模型和溶剂
        "qc_functionals": job_data.qc_options.functionals if job_data.qc_options and job_data.qc_options.enabled else None,
        "qc_basis_sets": job_data.qc_options.basis_sets if job_data.qc_options and job_data.qc_options.enabled else None,
        "qc_solvent_models": job_data.qc_options.solvent_models if job_data.qc_options and job_data.qc_options.enabled else None,
        "qc_solvents": job_data.qc_options.solvents if job_data.qc_options and job_data.qc_options.enabled else None,
        # 兼容旧版单选字段
        "qc_functional": job_data.qc_options.functional if job_data.qc_options and job_data.qc_options.enabled else None,
        "qc_basis_set": job_data.qc_options.basis_set if job_data.qc_options and job_data.qc_options.enabled else None,
        "qc_solvent_model": job_data.qc_options.solvent_model if job_data.qc_options and job_data.qc_options.enabled else None,
        "qc_solvent_name": job_data.qc_options.solvent_name if job_data.qc_options and job_data.qc_options.enabled else None,
        "qc_use_recommended_params": job_data.qc_options.use_recommended_params if job_data.qc_options and job_data.qc_options.enabled else None,
    }

    # 应用精度等级配置
    accuracy_level = job_data.accuracy_level or AccuracyLevel.STANDARD
    config = apply_accuracy_level(config, accuracy_level)

    logger.info(f"Applied accuracy level '{accuracy_level.value}' to job {job_name}")
    logger.info(f"QC options: enabled={config.get('qc_enabled')}, qc_options={job_data.qc_options}")
    logger.info(f"  Charge method: {config.get('charge_method')}")
    logger.info(f"  NPT steps: {config.get('nsteps_npt')}, NVT steps: {config.get('nsteps_nvt')}")

    # 根据是否提交到集群设置初始状态
    initial_status = JobStatus.QUEUED if job_data.submit_to_cluster else JobStatus.CREATED

    # Create MD job
    db_job = MDJob(
        system_id=job_data.system_id,
        user_id=current_user.id,
        status=initial_status,
        progress=0.0,
        config=config
    )

    db.add(db_job)
    db.commit()
    db.refresh(db_job)

    # 如果需要提交到集群，立即提交
    if job_data.submit_to_cluster:
        try:
            logger.info(f"Auto-submitting job {db_job.id} to cluster...")

            # 调用提交逻辑（复用 submit_md_job 的代码）
            _submit_job_to_cluster(db_job, system, db)

            logger.info(f"MD job submitted to cluster: ID={db_job.id} by {current_user.username}")
        except Exception as e:
            logger.error(f"Failed to auto-submit job {db_job.id}: {e}")
            # 更新任务状态为失败
            db_job.status = JobStatus.FAILED
            if db_job.config is None:
                db_job.config = {}
            db_job.config["error"] = str(e)
            db.commit()
            db.refresh(db_job)
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail=f"Failed to submit job to cluster: {str(e)}"
            )
    else:
        logger.info(f"MD job created: ID={db_job.id} by {current_user.username}")

    # 如果启用了QC计算，创建QC任务
    if job_data.qc_options and job_data.qc_options.enabled:
        try:
            _create_qc_jobs_for_md(db, db_job, system, current_user, job_data.qc_options)
        except Exception as e:
            logger.error(f"Failed to create QC jobs for MD job {db_job.id}: {e}")
            # QC任务创建失败不影响MD任务

    return db_job


@router.post("/batch")
def batch_create_md_jobs(
    batch_data: BatchMDJobCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    批量创建MD任务

    Args:
        batch_data: 批量创建请求数据（包含system_ids和任务配置）
        db: 数据库会话
        current_user: 当前用户

    Returns:
        批量创建结果
    """
    from fastapi.responses import JSONResponse
    from app.utils.quota import check_user_quota

    # 检查用户配额
    quota_check = check_user_quota(current_user, db)
    if not quota_check["allowed"]:
        return JSONResponse(
            status_code=403,
            content={
                "success": False,
                "message": quota_check["reason"],
                "quota_exceeded": True
            }
        )

    # 检查批量创建是否会超过每日任务限制
    quota_details = quota_check.get("details", {})
    today_jobs = quota_details.get("today_jobs", 0)
    daily_limit = quota_details.get("daily_limit", 999)
    remaining = daily_limit - today_jobs

    if len(batch_data.system_ids) > remaining:
        return JSONResponse(
            status_code=403,
            content={
                "success": False,
                "message": f"批量创建 {len(batch_data.system_ids)} 个任务将超过每日限制（剩余配额：{remaining}）",
                "quota_exceeded": True,
                "requested": len(batch_data.system_ids),
                "remaining": remaining
            }
        )

    results = {
        "success": True,
        "total": len(batch_data.system_ids),
        "success_count": 0,
        "failed_count": 0,
        "success_jobs": [],
        "errors": []
    }

    for system_id in batch_data.system_ids:
        try:
            # 为每个配方创建MD任务
            job_create = MDJobCreate(
                system_id=system_id,
                job_name=batch_data.job_name,
                accuracy_level=batch_data.accuracy_level,
                nsteps_npt=batch_data.nsteps_npt,
                nsteps_nvt=batch_data.nsteps_nvt,
                timestep=batch_data.timestep,
                temperature=batch_data.temperature,
                pressure=batch_data.pressure,
                freq_trj_npt=batch_data.freq_trj_npt,
                freq_trj_nvt=batch_data.freq_trj_nvt,
                thermo_freq=batch_data.thermo_freq,
                submit_to_cluster=batch_data.submit_to_cluster,
                slurm_partition=batch_data.slurm_partition,
                slurm_nodes=batch_data.slurm_nodes,
                slurm_ntasks=batch_data.slurm_ntasks,
                slurm_cpus_per_task=batch_data.slurm_cpus_per_task,
                slurm_time=batch_data.slurm_time,
                qc_options=batch_data.qc_options
            )

            md_job = create_md_job(job_create, db, current_user)

            results["success_count"] += 1
            results["success_jobs"].append({
                "system_id": system_id,
                "job_id": md_job.id,
                "job_name": md_job.config.get('job_name', 'N/A')
            })

        except HTTPException as e:
            results["failed_count"] += 1
            error_detail = e.detail if isinstance(e.detail, str) else str(e.detail)
            results["errors"].append({
                "system_id": system_id,
                "error": error_detail
            })
        except Exception as e:
            results["failed_count"] += 1
            results["errors"].append({
                "system_id": system_id,
                "error": str(e)
            })

    if results["failed_count"] > 0:
        results["success"] = False

    return JSONResponse(content=results)


def _get_recommended_qc_params(mol_type: str, base_options) -> dict:
    """
    根据分子类型获取推荐的QC计算参数

    Args:
        mol_type: 分子类型 (solvent, cation, anion)
        base_options: 用户设置的基础选项

    Returns:
        dict: 包含 basis_set, functional, solvent_model, solvent_name 的字典
    """
    # 基础参数
    params = {
        "basis_set": base_options.basis_set or "6-31++g(d,p)",
        "functional": base_options.functional or "B3LYP",
        "solvent_model": base_options.solvent_model or "pcm",
        "solvent_name": base_options.solvent_name or "water",
        "recommendation_reason": ""
    }

    # 如果不使用推荐参数，直接返回用户设置
    if not getattr(base_options, 'use_recommended_params', True):
        params["recommendation_reason"] = "使用用户自定义参数"
        return params

    if mol_type == "anion":
        # 阴离子需要弥散函数来描述扩散的电子云
        # 推荐使用 6-31++G(d,p) 或 aug-cc-pVDZ
        if "+" not in params["basis_set"]:
            params["basis_set"] = "6-31++g(d,p)"
        params["recommendation_reason"] = "阴离子使用带弥散函数(++)的基组，以更好地描述扩散的电子密度"
        # 阴离子通常需要隐式溶剂模型来稳定
        if params["solvent_model"] == "gas":
            params["solvent_model"] = "pcm"
            params["recommendation_reason"] += "；气相阴离子可能不稳定，建议使用PCM溶剂模型"

    elif mol_type == "cation":
        # 阳离子通常电子更加紧凑，标准基组即可
        # 但仍推荐使用极化函数
        if "d" not in params["basis_set"].lower() and "p" not in params["basis_set"].lower():
            params["basis_set"] = "6-31g(d,p)"
        params["recommendation_reason"] = "阳离子使用带极化函数的基组"

    else:  # solvent
        # 中性溶剂分子使用标准参数
        params["recommendation_reason"] = "中性分子使用标准计算参数"

    return params


def _calculate_spin_multiplicity(smiles: str, charge: int) -> int:
    """
    根据SMILES和电荷计算自旋多重度

    简化规则：
    - 对于常见离子，使用预定义值
    - 对于其他分子，尝试使用RDKit计算
    - 默认返回1（单重态）
    """
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 1

        # 计算总电子数
        total_electrons = sum(atom.GetAtomicNum() for atom in mol.GetAtoms())
        total_electrons -= charge  # 电荷影响电子数

        # 如果电子数为偶数，通常是单重态；如果为奇数，通常是二重态
        if total_electrons % 2 == 0:
            return 1  # 单重态
        else:
            return 2  # 二重态
    except Exception:
        return 1  # 默认单重态


def _create_qc_jobs_for_md(db: Session, md_job: MDJob, system: ElectrolyteSystem,
                           user: User, qc_options):
    """
    为MD任务创建关联的QC任务

    Args:
        db: 数据库会话
        md_job: MD任务
        system: 电解质系统
        user: 用户
        qc_options: QC选项
    """
    from app.tasks.qc_submission import submit_qc_job_task

    # 从电解质配方的JSONB字段中提取分子信息
    # 格式: [{"name": "...", "smiles": "...", ...}, ...]
    molecules_to_calc = []  # [(smiles, name, mol_type, charge), ...]
    seen_smiles = set()

    if qc_options.molecules:
        # 使用用户指定的分子列表
        for smiles in qc_options.molecules:
            if smiles and smiles not in seen_smiles:
                molecules_to_calc.append((smiles, f"custom_{len(molecules_to_calc)}", "custom", 0))
                seen_smiles.add(smiles)
    else:
        # 从电解质配方中提取所有分子
        # 溶剂分子
        if system.solvents:
            for sol in system.solvents:
                smiles = sol.get("smiles")
                name = sol.get("name", "solvent")
                if smiles and smiles not in seen_smiles:
                    molecules_to_calc.append((smiles, name, "solvent", 0))
                    seen_smiles.add(smiles)

        # 阳离子
        if system.cations:
            for cat in system.cations:
                smiles = cat.get("smiles")
                name = cat.get("name", "cation")
                charge = cat.get("charge", 1)  # 尝试从数据中获取电荷，默认+1
                if smiles and smiles not in seen_smiles:
                    molecules_to_calc.append((smiles, name, "cation", charge))
                    seen_smiles.add(smiles)

        # 阴离子
        if system.anions:
            for an in system.anions:
                smiles = an.get("smiles")
                name = an.get("name", "anion")
                charge = an.get("charge", -1)  # 尝试从数据中获取电荷，默认-1
                if smiles and smiles not in seen_smiles:
                    molecules_to_calc.append((smiles, name, "anion", charge))
                    seen_smiles.add(smiles)

    # 获取计算参数列表（支持多选）
    basis_sets = getattr(qc_options, 'basis_sets', None) or [getattr(qc_options, 'basis_set', '6-31++g(d,p)')]
    functionals = getattr(qc_options, 'functionals', None) or [getattr(qc_options, 'functional', 'B3LYP')]
    solvent_models = getattr(qc_options, 'solvent_models', None) or [getattr(qc_options, 'solvent_model', 'pcm')]
    solvents = getattr(qc_options, 'solvents', None) or [getattr(qc_options, 'solvent_name', 'Water')]

    # 构建溶剂组合
    # 如果溶剂模型包含 gas，则只使用 gas（不需要溶剂）
    # 否则，为每个溶剂模型和溶剂的组合创建任务
    solvent_combinations = []
    if 'gas' in solvent_models:
        solvent_combinations.append(('gas', None))

    # 为非气相模型创建溶剂组合
    non_gas_models = [m for m in solvent_models if m != 'gas']
    for model in non_gas_models:
        for solvent in solvents:
            solvent_combinations.append((model, solvent))

    # 如果没有任何组合，使用默认值
    if not solvent_combinations:
        solvent_combinations = [('pcm', 'Water')]

    # 计算总任务数：分子 × 泛函 × 基组 × 溶剂组合
    total_jobs = len(molecules_to_calc) * len(functionals) * len(basis_sets) * len(solvent_combinations)
    logger.info(f"Creating {total_jobs} QC jobs ({len(molecules_to_calc)} molecules × {len(functionals)} functionals × {len(basis_sets)} basis sets × {len(solvent_combinations)} solvent combinations) for MD job {md_job.id}")

    for smiles, mol_name, mol_type, charge in molecules_to_calc:
        # 计算自旋多重度（对所有参数组合都相同）
        spin_multiplicity = _calculate_spin_multiplicity(smiles, charge)

        # 遍历所有参数组合（笛卡尔积）
        for functional in functionals:
            for basis_set in basis_sets:
                for solvent_model, solvent_name in solvent_combinations:
                    # 根据分子类型获取推荐参数
                    params = _get_recommended_qc_params(mol_type, qc_options)

                    # 使用当前循环的参数覆盖推荐参数
                    params["basis_set"] = basis_set
                    params["functional"] = functional
                    params["solvent_model"] = solvent_model
                    params["solvent_name"] = solvent_name

                    # 构建溶剂配置
                    solvent_config = None
                    if params["solvent_model"] != "gas":
                        solvent_config = {
                            "model": params["solvent_model"],
                            "solvent_name": params["solvent_name"]
                        }

                    # 构建任务名称（包含所有参数以区分）
                    # 格式: 分子名_泛函_基组_溶剂模型_溶剂
                    name_parts = [mol_name, functional, basis_set.replace('(', '').replace(')', '').replace('+', 'p')]
                    if solvent_model == 'gas':
                        name_parts.append('gas')
                    else:
                        name_parts.extend([solvent_model, solvent_name])
                    job_mol_name = '_'.join(name_parts)

                    # 创建QC任务（所有参数存储在config中）
                    qc_job = QCJob(
                        user_id=user.id,
                        md_job_id=md_job.id,
                        molecule_name=job_mol_name,
                        smiles=smiles,
                        molecule_type=mol_type,
                        basis_set=params["basis_set"],
                        functional=params["functional"],
                        charge=charge,
                        spin_multiplicity=spin_multiplicity,
                        status=QCJobStatus.CREATED,
                        config={
                            "accuracy_level": getattr(qc_options, 'accuracy_level', 'standard'),
                            "solvent_model": params["solvent_model"],
                            "solvent_name": params["solvent_name"] if params["solvent_model"] != "gas" else None,
                            "solvent_config": solvent_config,
                            "recommendation_reason": params["recommendation_reason"],
                            "auto_params": getattr(qc_options, 'use_recommended_params', True),
                        }
                    )

                    db.add(qc_job)
                    db.commit()
                    db.refresh(qc_job)

                    logger.info(f"Created QC job {qc_job.id} for '{job_mol_name}' (type: {mol_type}, "
                               f"charge: {charge}, spin: {spin_multiplicity}, "
                               f"functional: {functional}, basis: {basis_set}, "
                               f"solvent: {solvent_model}/{solvent_name or 'N/A'})")

                    # 提交QC任务到Celery
                    try:
                        submit_qc_job_task.delay(qc_job.id)
                        logger.info(f"Submitted QC job {qc_job.id} to Celery")
                    except Exception as e:
                        logger.error(f"Failed to submit QC job {qc_job.id} to Celery: {e}")


@router.get("/quota/check")
def check_job_quota(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Check user's daily job creation quota

    Returns:
        dict: {
            "can_create": bool,
            "current_count": int,
            "limit": int,
            "remaining": int
        }
    """
    can_create, current_count, limit = check_daily_job_limit(db, current_user.id, current_user.role)

    return {
        "can_create": can_create,
        "current_count": current_count,
        "limit": limit,
        "remaining": limit - current_count
    }


@router.get("/", response_model=List[MDJobSchema])
def list_md_jobs(
    system_id: int = None,
    status_filter: JobStatus = None,
    skip: int = 0,
    limit: int = 100,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    List MD jobs
    
    Args:
        system_id: Filter by system ID
        status_filter: Filter by job status
        skip: Number of records to skip
        limit: Maximum number of records to return
        db: Database session
        current_user: Current authenticated user
        
    Returns:
        List[MDJob]: List of MD jobs
    """
    query = db.query(MDJob)

    # Filter by user unless admin
    if current_user.role != UserRole.ADMIN:
        query = query.filter(MDJob.user_id == current_user.id)
    
    # Filter by system if specified
    if system_id:
        query = query.filter(MDJob.system_id == system_id)
    
    # Filter by status if specified
    if status_filter:
        query = query.filter(MDJob.status == status_filter)
    
    # Order by creation time (newest first)
    query = query.order_by(MDJob.created_at.desc())

    jobs = query.offset(skip).limit(limit).all()

    # 为管理员添加用户信息
    if current_user.role == UserRole.ADMIN:
        result = []
        for job in jobs:
            job_dict = {
                "id": job.id,
                "system_id": job.system_id,
                "user_id": job.user_id,
                "status": job.status,
                "slurm_job_id": job.slurm_job_id,
                "progress": job.progress,
                "work_dir": job.work_dir,
                "log_file": job.log_file,
                "error_message": job.error_message,
                "config": job.config,
                "created_at": job.created_at,
                "updated_at": job.updated_at,
                "started_at": job.started_at,
                "finished_at": job.finished_at,
            }
            # 添加用户信息
            user = db.query(User).filter(User.id == job.user_id).first()
            if user:
                job_dict["username"] = user.username
                job_dict["user_email"] = user.email
            result.append(job_dict)
        return result

    return jobs


@router.get("/{job_id}", response_model=MDJobSchema)
def get_md_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Get MD job by ID
    
    Args:
        job_id: MD job ID
        db: Database session
        current_user: Current authenticated user
        
    Returns:
        MDJob: MD job data
        
    Raises:
        HTTPException: If not found or no permission
    """
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    
    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="MD job not found"
        )
    
    # Check permission
    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not enough permissions"
        )
    
    return job


@router.get("/{job_id}/qc-jobs")
def get_md_job_qc_jobs(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取MD任务关联的所有QC任务及状态汇总

    Args:
        job_id: MD任务ID
        db: 数据库会话
        current_user: 当前用户

    Returns:
        QC任务列表和状态汇总
    """
    from app.schemas.job import QCJobSummary, QCJobsStatusSummary

    # 获取MD任务
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="MD job not found"
        )

    # 检查权限
    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not enough permissions"
        )

    # 获取关联的QC任务
    qc_jobs = db.query(QCJob).filter(QCJob.md_job_id == job_id).order_by(QCJob.created_at).all()

    # 统计状态
    status_summary = QCJobsStatusSummary(total=len(qc_jobs))
    for qc_job in qc_jobs:
        status_str = qc_job.status.value if hasattr(qc_job.status, 'value') else str(qc_job.status)
        if status_str == "CREATED":
            status_summary.created += 1
        elif status_str == "QUEUED":
            status_summary.queued += 1
        elif status_str == "RUNNING":
            status_summary.running += 1
        elif status_str == "POSTPROCESSING":
            status_summary.postprocessing += 1
        elif status_str == "COMPLETED":
            status_summary.completed += 1
        elif status_str == "FAILED":
            status_summary.failed += 1
        elif status_str == "CANCELLED":
            status_summary.cancelled += 1

    # 转换为Schema
    qc_jobs_data = [
        QCJobSummary(
            id=qc.id,
            molecule_name=qc.molecule_name,
            smiles=qc.smiles,
            molecule_type=qc.molecule_type or "custom",
            status=qc.status.value if hasattr(qc.status, 'value') else str(qc.status),
            progress=qc.progress or 0.0,
            basis_set=qc.basis_set,
            functional=qc.functional or "B3LYP",
            charge=qc.charge if hasattr(qc, 'charge') and qc.charge is not None else 0,
            spin_multiplicity=qc.spin_multiplicity if hasattr(qc, 'spin_multiplicity') and qc.spin_multiplicity is not None else 1,
            solvent_model=qc.solvent_model if hasattr(qc, 'solvent_model') else None,
            solvent_name=qc.solvent_name if hasattr(qc, 'solvent_name') else None,
            accuracy_level=qc.accuracy_level if hasattr(qc, 'accuracy_level') else None,
            is_reused=qc.is_reused if hasattr(qc, 'is_reused') and qc.is_reused else False,
            reused_from_job_id=qc.reused_from_job_id if hasattr(qc, 'reused_from_job_id') else None,
            slurm_job_id=qc.slurm_job_id,
            work_dir=qc.work_dir,
            created_at=qc.created_at,
            started_at=qc.started_at,
            finished_at=qc.finished_at,
            error_message=qc.error_message
        )
        for qc in qc_jobs
    ]

    return {
        "md_job_id": job_id,
        "qc_jobs": qc_jobs_data,
        "status_summary": status_summary,
        "qc_enabled": len(qc_jobs) > 0
    }


@router.put("/{job_id}", response_model=MDJobSchema)
def update_md_job(
    job_id: int,
    job_update: MDJobUpdate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Update MD job

    Args:
        job_id: MD job ID
        job_update: MD job update data
        db: Database session
        current_user: Current authenticated user
        
    Returns:
        MDJob: Updated MD job
        
    Raises:
        HTTPException: If not found or no permission
    """
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    
    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="MD job not found"
        )
    
    # Check permission (only owner or admin can update)
    check_job_permission(job, current_user)
    
    # Update fields
    update_data = job_update.model_dump(exclude_unset=True)
    for field, value in update_data.items():
        setattr(job, field, value)
    
    db.commit()
    db.refresh(job)
    
    logger.info(f"MD job updated: ID={job.id}")
    return job


@router.put("/{job_id}/config", response_model=MDJobSchema)
def update_md_job_config(
    job_id: int,
    config: Dict[str, Any],
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Update MD job configuration

    Args:
        job_id: MD job ID
        config: Configuration data
        db: Database session
        current_user: Current authenticated user

    Returns:
        MDJob: Updated MD job

    Raises:
        HTTPException: If not found, no permission, or job already submitted
    """
    job = db.query(MDJob).filter(MDJob.id == job_id).first()

    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="MD job not found"
        )

    # Check permission
    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not enough permissions"
        )

    # Allow updating config for CREATED, FAILED, CANCELLED, and COMPLETED jobs
    # (for resubmission purposes)
    allowed_statuses = [JobStatus.CREATED, JobStatus.FAILED, JobStatus.CANCELLED, JobStatus.COMPLETED]
    if job.status not in allowed_statuses:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Cannot update config for jobs in {job.status} status. Only CREATED, FAILED, CANCELLED, or COMPLETED jobs can be updated."
        )

    # Update config
    job.config = config
    db.commit()
    db.refresh(job)

    logger.info(f"MD job config updated: ID={job.id}, Status={job.status}")
    return job


@router.post("/{job_id}/submit", response_model=MDJobSchema)
def submit_md_job(
    job_id: int,
    use_celery: bool = True,  # 新增参数：是否使用 Celery 异步提交
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Submit MD job to cluster

    可以选择同步提交（use_celery=False）或异步提交（use_celery=True，默认）

    Args:
        job_id: MD job ID
        use_celery: 是否使用 Celery 异步提交（默认 True）
        db: Database session
        current_user: Current authenticated user

    Returns:
        MDJob: Updated MD job

    Raises:
        HTTPException: If not found, no permission, or job already submitted
    """
    # 检查用户余额（计费系统）
    from app.services.billing import BillingService
    can_submit, reason = BillingService.can_submit_job(db, current_user)
    if not can_submit:
        raise HTTPException(
            status_code=status.HTTP_402_PAYMENT_REQUIRED,
            detail=reason
        )

    job = db.query(MDJob).filter(MDJob.id == job_id).first()

    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="MD job not found"
        )

    # Check permission
    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not enough permissions"
        )

    # Only allow submitting CREATED or CANCELLED jobs
    # CANCELLED jobs can be reconfigured and resubmitted
    if job.status not in [JobStatus.CREATED, JobStatus.CANCELLED]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Cannot submit job with status {job.status}. Only CREATED or CANCELLED jobs can be submitted."
        )

    # Get electrolyte system
    electrolyte = db.query(ElectrolyteSystem).filter(
        ElectrolyteSystem.id == job.system_id
    ).first()

    if not electrolyte:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Electrolyte system not found"
        )

    try:
        # 如果启用了QC计算但还没有创建QC任务，现在创建
        if job.config and job.config.get("qc_enabled"):
            # 检查是否已经有QC任务
            existing_qc_jobs = db.query(QCJob).filter(QCJob.md_job_id == job.id).count()
            if existing_qc_jobs == 0:
                logger.info(f"Creating QC jobs for MD job {job.id} on submit...")
                # 构造QC选项
                from app.schemas.job import MDJobQCOptions
                qc_options = MDJobQCOptions(
                    enabled=True,
                    basis_set=job.config.get("qc_basis_set", "6-31++g(d,p)"),
                    functional=job.config.get("qc_functional", "B3LYP"),
                    molecules=None
                )
                try:
                    _create_qc_jobs_for_md(db, job, electrolyte, current_user, qc_options)
                    logger.info(f"QC jobs created for MD job {job.id}")
                except Exception as qc_error:
                    logger.error(f"Failed to create QC jobs for MD job {job.id}: {qc_error}")
                    # QC任务创建失败不影响MD任务提交

        if use_celery:
            # 使用 Celery 异步提交
            from app.tasks.job_submission import submit_md_job_task

            # 更新状态为 PENDING（等待 Celery 处理）
            job.status = JobStatus.CREATED  # 保持 CREATED，Celery 会更新为 QUEUED
            if job.config is None:
                job.config = {}
            job.config["celery_task_submitted"] = True
            job.config["celery_task_submitted_at"] = datetime.now().isoformat()
            job.config["celery_task_submitted_by"] = current_user.username
            db.commit()

            # 提交到 Celery 队列
            task = submit_md_job_task.delay(job.id)

            # 记录 Celery 任务 ID
            job.config["celery_task_id"] = task.id
            db.commit()
            db.refresh(job)

            logger.info(f"MD job {job.id} submitted to Celery queue: task_id={task.id} by {current_user.username}")
            return job
        else:
            # 同步提交（原有逻辑）
            _submit_job_to_cluster(job, electrolyte, db)
            db.refresh(job)

            # Read from the model field, not config
            slurm_job_id = job.slurm_job_id
            logger.info(f"MD job submitted to cluster (sync): ID={job.id}, Slurm ID={slurm_job_id} by {current_user.username}")
            return job

    except Exception as e:
        logger.error(f"Failed to submit job {job.id}: {e}")

        # Update job status to FAILED
        job.status = JobStatus.FAILED

        if job.config is None:
            job.config = {}

        job.config["error"] = str(e)

        db.commit()
        db.refresh(job)

        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to submit job: {str(e)}"
        )


@router.post("/{job_id}/cancel", response_model=MDJobSchema)
def cancel_md_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Cancel a running or queued MD job

    This endpoint cancels a job that is currently queued or running in Slurm.
    It uses scancel to stop the Slurm job and updates the database status.

    Args:
        job_id: MD job ID
        db: Database session
        current_user: Current authenticated user

    Returns:
        MDJob: Updated MD job

    Raises:
        HTTPException: If not found, no permission, or job cannot be cancelled
    """
    job = db.query(MDJob).filter(MDJob.id == job_id).first()

    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="MD job not found"
        )

    # Check permission
    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not enough permissions"
        )

    # Only allow cancelling QUEUED or RUNNING jobs
    if job.status not in [JobStatus.QUEUED, JobStatus.RUNNING]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Cannot cancel job with status {job.status}. Only QUEUED or RUNNING jobs can be cancelled."
        )

    # Get Slurm job ID
    slurm_job_id = job.slurm_job_id

    if not slurm_job_id:
        # If no Slurm job ID, just update status
        logger.warning(f"Job {job_id} has no Slurm job ID, updating status only")
        job.status = JobStatus.CANCELLED
        job.error_message = "Cancelled by user (no Slurm job ID found)"
        db.commit()
        db.refresh(job)
        return job

    try:
        # Cancel the Slurm job using scancel
        result = subprocess.run(
            ["scancel", str(slurm_job_id)],
            capture_output=True,
            text=True,
            timeout=10
        )

        if result.returncode != 0:
            # scancel might fail if job already finished
            logger.warning(f"scancel failed for job {slurm_job_id}: {result.stderr}")
            # Still update status as user requested cancellation

        # Update job status
        job.status = JobStatus.CANCELLED
        job.error_message = f"Cancelled by user {current_user.username}"

        if job.config is None:
            job.config = {}
        job.config["cancelled_by"] = current_user.username
        job.config["cancelled_at"] = datetime.now().isoformat()

        db.commit()
        db.refresh(job)

        logger.info(f"Job {job_id} (Slurm ID: {slurm_job_id}) cancelled by {current_user.username}")
        return job

    except subprocess.TimeoutExpired:
        raise HTTPException(
            status_code=status.HTTP_504_GATEWAY_TIMEOUT,
            detail="Slurm scancel command timeout"
        )
    except Exception as e:
        logger.error(f"Failed to cancel job {job_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to cancel job: {str(e)}"
        )


@router.post("/{job_id}/sync_status")
def sync_job_status(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Sync job status with Slurm

    Queries Slurm for the current status and updates the database
    """
    # Get job
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # Check permission
    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="Not authorized to access this job")

    # Get Slurm job ID - check both job.slurm_job_id and job.config
    slurm_job_id = job.slurm_job_id
    if not slurm_job_id and job.config:
        slurm_job_id = job.config.get("slurm_job_id")

    if not slurm_job_id:
        raise HTTPException(status_code=400, detail="No Slurm job ID found")

    try:
        # Query Slurm using sacct - get more detailed information
        result = subprocess.run(
            ["sacct", "-j", slurm_job_id, "--format=State,Start,End", "--noheader", "--parsable2"],
            capture_output=True,
            text=True,
            timeout=10
        )

        if result.returncode != 0:
            raise Exception(f"sacct failed: {result.stderr}")

        # Parse status and timestamps
        output_line = result.stdout.strip().split("\n")[0] if result.stdout else ""
        parts = output_line.split("|") if output_line else []

        slurm_status = parts[0] if len(parts) > 0 else ""
        start_time_str = parts[1] if len(parts) > 1 else ""
        end_time_str = parts[2] if len(parts) > 2 else ""

        logger.info(f"Slurm job {slurm_job_id} status: {slurm_status}, start: {start_time_str}, end: {end_time_str}")

        # Map Slurm status to our JobStatus
        status_mapping = {
            "PENDING": JobStatus.QUEUED,
            "RUNNING": JobStatus.RUNNING,
            "COMPLETED": JobStatus.COMPLETED,
            "FAILED": JobStatus.FAILED,
            "CANCELLED": JobStatus.CANCELLED,
            "TIMEOUT": JobStatus.FAILED,
            "NODE_FAIL": JobStatus.FAILED,
            "OUT_OF_MEMORY": JobStatus.FAILED
        }

        new_status = status_mapping.get(slurm_status, job.status)

        logger.info(f"Job {job_id} current status: {job.status}, new status: {new_status}")
        logger.info(f"Job {job_id} current progress: {job.progress}%")

        # Update progress based on status
        old_progress = job.progress
        if new_status == JobStatus.QUEUED:
            job.progress = 0.0
        elif new_status == JobStatus.RUNNING:
            job.progress = 50.0  # Running is halfway
        elif new_status == JobStatus.POSTPROCESSING:
            job.progress = 90.0
        elif new_status == JobStatus.COMPLETED:
            job.progress = 100.0
        elif new_status == JobStatus.FAILED or new_status == JobStatus.CANCELLED:
            # Keep current progress for failed/cancelled jobs
            pass

        logger.info(f"Job {job_id} new progress: {job.progress}%")

        # Update timestamps from Slurm data
        from datetime import datetime
        timestamp_updated = False

        # Parse and update start time
        if start_time_str and start_time_str != "Unknown" and not job.started_at:
            try:
                # Slurm format: 2024-11-27T10:30:45
                job.started_at = datetime.strptime(start_time_str, "%Y-%m-%dT%H:%M:%S")
                timestamp_updated = True
                logger.info(f"Job {job_id} started_at updated from Slurm: {job.started_at}")
            except ValueError as e:
                logger.warning(f"Failed to parse Slurm start time '{start_time_str}': {e}")

        # Parse and update end time
        if end_time_str and end_time_str != "Unknown" and not job.finished_at:
            try:
                job.finished_at = datetime.strptime(end_time_str, "%Y-%m-%dT%H:%M:%S")
                timestamp_updated = True
                logger.info(f"Job {job_id} finished_at updated from Slurm: {job.finished_at}")
            except ValueError as e:
                logger.warning(f"Failed to parse Slurm end time '{end_time_str}': {e}")

        # Update job status if changed
        status_changed = new_status != job.status
        progress_changed = job.progress != old_progress

        if status_changed or progress_changed or timestamp_updated:
            if status_changed:
                old_status = job.status
                job.status = new_status

                # Set timestamps if not already set from Slurm data
                if new_status == JobStatus.RUNNING and not job.started_at:
                    job.started_at = datetime.now()
                    logger.info(f"Job {job_id} started at {job.started_at} (fallback)")

                if new_status in [JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELLED] and not job.finished_at:
                    job.finished_at = datetime.now()
                    logger.info(f"Job {job_id} finished at {job.finished_at} (fallback)")

                logger.info(f"Job {job_id} status updated: {old_status} -> {new_status}")

                # 任务完成后自动结算（计费系统）
                if new_status == JobStatus.COMPLETED and not job.billed:
                    try:
                        from app.services.billing import BillingService
                        success, msg = BillingService.settle_job(db, job)
                        logger.info(f"Job {job_id} billing settled: {msg}")
                    except Exception as billing_error:
                        logger.error(f"Job {job_id} billing failed: {billing_error}")

            if progress_changed:
                logger.info(f"Job {job_id} progress updated: {old_progress}% -> {job.progress}%")

            db.commit()
            db.refresh(job)

        return {
            "job_id": job_id,
            "slurm_job_id": slurm_job_id,
            "slurm_status": slurm_status,
            "job_status": job.status,
            "progress": job.progress,
            "updated": status_changed or progress_changed or timestamp_updated
        }

    except subprocess.TimeoutExpired:
        raise HTTPException(status_code=504, detail="Slurm query timeout")
    except Exception as e:
        logger.error(f"Failed to sync job status: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{job_id}/resubmit", response_model=MDJobSchema)
def resubmit_md_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Resubmit a failed or cancelled MD job to cluster

    This endpoint allows resubmitting jobs that have failed or been cancelled.
    It will regenerate the input files and resubmit to Slurm.

    Args:
        job_id: MD job ID
        db: Database session
        current_user: Current authenticated user

    Returns:
        MDJob: Updated MD job

    Raises:
        HTTPException: If not found, no permission, or job cannot be resubmitted
    """
    job = db.query(MDJob).filter(MDJob.id == job_id).first()

    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="MD job not found"
        )

    # Check permission
    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not enough permissions"
        )

    # Only allow resubmitting FAILED or CANCELLED jobs
    if job.status not in [JobStatus.FAILED, JobStatus.CANCELLED]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Cannot resubmit job with status {job.status}. Only FAILED or CANCELLED jobs can be resubmitted."
        )

    # Get electrolyte system
    electrolyte = db.query(ElectrolyteSystem).filter(
        ElectrolyteSystem.id == job.system_id
    ).first()

    if not electrolyte:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Electrolyte system not found"
        )

    try:
        # Clear previous error message
        if job.config is None:
            job.config = {}
        job.config.pop("error", None)

        # Reset job status to CREATED before resubmission
        job.status = JobStatus.CREATED
        job.error_message = None
        job.progress = 0.0
        db.commit()

        logger.info(f"Resubmitting job {job.id} by {current_user.username}")

        # Use the internal submission function
        _submit_job_to_cluster(job, electrolyte, db)
        db.refresh(job)

        slurm_job_id = job.slurm_job_id
        logger.info(f"MD job resubmitted to cluster: ID={job.id}, Slurm ID={slurm_job_id} by {current_user.username}")
        return job

    except Exception as e:
        logger.error(f"Failed to resubmit job {job.id}: {e}")

        # Update job status back to FAILED
        job.status = JobStatus.FAILED

        if job.config is None:
            job.config = {}

        job.config["error"] = str(e)
        job.error_message = str(e)

        db.commit()
        db.refresh(job)

        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to resubmit job: {str(e)}"
        )


@router.get("/{job_id}/atom_mapping")
def get_atom_mapping(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Get atom mapping for a job

    Returns the atom_mapping.json file content
    """
    # Get job
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # Check permission
    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="Not authorized to access this job")

    # Get work directory from job
    if job.work_dir:
        work_dir = Path(job.work_dir)
    elif job.config and job.config.get("work_dir"):
        work_dir = Path(job.config.get("work_dir"))
    elif job.config and job.config.get("job_name"):
        # Fallback: construct from job_name
        work_dir = settings.MOLYTE_WORK_BASE_PATH / job.config.get("job_name")
    else:
        raise HTTPException(status_code=404, detail="Work directory not found for this job")

    mapping_file = work_dir / "atom_mapping.json"

    if not mapping_file.exists():
        raise HTTPException(status_code=404, detail="Atom mapping file not found")

    try:
        import json
        with open(mapping_file, 'r') as f:
            mapping = json.load(f)

        return mapping
    except Exception as e:
        logger.error(f"Failed to read atom mapping: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{job_id}/available_labels")
def get_available_labels(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Get available atom labels for RDF calculation

    Returns:
        {
            "all": ["Li_Li", "TEP_P00", ...],
            "by_molecule": {"Li": [...], "TEP": [...]},
            "by_element": {"Li": [...], "O": [...]}
        }
    """
    # Get job
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # Check permission
    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="Not authorized to access this job")

    # Get work directory from job
    if job.work_dir:
        work_dir = Path(job.work_dir)
    elif job.config and job.config.get("work_dir"):
        work_dir = Path(job.config.get("work_dir"))
    elif job.config and job.config.get("job_name"):
        # Fallback: construct from job_name
        work_dir = settings.MOLYTE_WORK_BASE_PATH / job.config.get("job_name")
    else:
        raise HTTPException(status_code=404, detail="Work directory not found for this job")

    try:
        from app.workers.rdf_calculator import RDFCalculator

        calculator = RDFCalculator(work_dir)
        labels = calculator.get_available_labels()
        return labels
    except Exception as e:
        logger.error(f"Failed to get available labels: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{job_id}/calculate_rdf")
def calculate_rdf(
    job_id: int,
    center_label: str,
    target_label: str,
    r_max: float = 10.0,
    n_bins: int = 200,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Calculate RDF between center and target atoms

    Args:
        job_id: Job ID
        center_label: Center atom label (e.g., "Li_Li", "Li_*")
        target_label: Target atom label (e.g., "TEP_O01", "*_O*")
        r_max: Maximum distance (Å)
        n_bins: Number of bins

    Returns:
        {
            "r": [0.05, 0.15, ...],
            "g_r": [0.0, 0.1, ...],
            "center_label": "Li_Li",
            "target_label": "TEP_O01",
            "center_atom_count": 50,
            "target_atom_count": 100,
            "frame_count": 100
        }
    """
    # Get job
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # Check permission
    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="Not authorized to access this job")

    # Check job status
    if job.status not in [JobStatus.COMPLETED, JobStatus.POSTPROCESSING]:
        raise HTTPException(
            status_code=400,
            detail=f"Job is not completed yet (status: {job.status})"
        )

    # Get work directory from job
    if job.work_dir:
        work_dir = Path(job.work_dir)
    elif job.config and job.config.get("work_dir"):
        work_dir = Path(job.config.get("work_dir"))
    elif job.config and job.config.get("job_name"):
        # Fallback: construct from job_name
        work_dir = settings.MOLYTE_WORK_BASE_PATH / job.config.get("job_name")
    else:
        raise HTTPException(status_code=404, detail="Work directory not found for this job")

    try:
        from app.workers.rdf_calculator import RDFCalculator

        calculator = RDFCalculator(work_dir)
        result = calculator.calculate_rdf(
            center_label=center_label,
            target_label=target_label,
            r_max=r_max,
            n_bins=n_bins
        )

        logger.info(f"Calculated RDF for job {job_id}: {center_label} -> {target_label}")
        return result

    except Exception as e:
        logger.error(f"Failed to calculate RDF: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{job_id}/msd_results")
def get_msd_results(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Get saved MSD results for a job

    Returns:
        List of MSD results with diffusion coefficients
    """
    from app.tasks.msd_processor import get_msd_results as get_msd_data

    # Get job
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # Check permission
    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="Permission denied")

    # Get MSD results
    results = get_msd_data(db, job_id)

    return results


@router.post("/{job_id}/calculate_msd")
def calculate_msd(
    job_id: int,
    background_tasks: BackgroundTasks,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Calculate MSD for a completed job

    Returns:
        Message indicating MSD calculation started
    """
    from app.tasks.msd_processor import process_msd_data
    from pathlib import Path

    # Get job
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # Check permission
    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="Permission denied")

    # Check job status
    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job must be completed before calculating MSD")

    # Check work directory
    if not job.work_dir or not Path(job.work_dir).exists():
        raise HTTPException(status_code=400, detail="Work directory not found")

    try:
        # 获取任务配置中的温度、盒子体积和离子数量
        temperature = 298.15  # 默认温度
        box_volume = None
        ion_counts = None

        if job.config:
            # 从配置中获取温度
            if 'temperature' in job.config:
                temperature = float(job.config['temperature'])

            # 从配置中获取盒子体积
            if 'box_volume' in job.config:
                box_volume = float(job.config['box_volume'])

            # 从配置中获取离子数量
            if 'ion_counts' in job.config:
                ion_counts = job.config['ion_counts']
            elif 'electrolyte' in job.config:
                # 尝试从电解质配置中提取离子数量
                electrolyte = job.config['electrolyte']
                ion_counts = {}
                if 'cation' in electrolyte and 'cation_count' in electrolyte:
                    ion_counts[electrolyte['cation']] = electrolyte['cation_count']
                if 'anion' in electrolyte and 'anion_count' in electrolyte:
                    ion_counts[electrolyte['anion']] = electrolyte['anion_count']

        # Process MSD data with temperature, box_volume, and ion_counts
        results = process_msd_data(
            db,
            job_id,
            Path(job.work_dir),
            temperature=temperature,
            box_volume=box_volume,
            ion_counts=ion_counts
        )

        return {
            "message": "MSD calculation completed",
            "num_results": len(results),
            "results": [
                {
                    "species": r.species,
                    "diffusion_coefficient": r.diffusion_coefficient,
                    "ionic_conductivity": r.ionic_conductivity,
                    "mobility": r.mobility,
                    "charge": r.charge,
                }
                for r in results
            ]
        }
    except Exception as e:
        logger.error(f"Failed to calculate MSD: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{job_id}/rdf_results")
def get_rdf_results(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Get saved RDF results for a job

    Returns:
        List of RDF results with analysis data
    """
    from app.models.result import RDFResult
    from app.dependencies import check_resource_permission

    # Get job
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # 数据隔离：检查权限（管理员可以访问所有数据，普通用户只能访问自己的数据）
    check_resource_permission(job.user_id, current_user)

    # Get RDF results
    rdf_results = db.query(RDFResult).filter(RDFResult.md_job_id == job_id).all()

    return [
        {
            "id": result.id,
            "center_species": result.center_species,
            "shell_species": result.shell_species,
            "r": result.r_values,
            "g_r": result.g_r_values,
            "coordination_number_values": result.coordination_number_values,  # 配位数数组
            "first_peak_position": result.first_peak_position,
            "first_peak_height": result.first_peak_height,
            "coordination_number": result.coordination_number,  # 最终配位数
            "created_at": result.created_at.isoformat() if result.created_at else None,
        }
        for result in rdf_results
    ]


@router.get("/{job_id}/molecule_templates")
def get_molecule_templates(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Get molecule templates with 3D structures and charges

    Returns:
        List of molecules with PDB structure and charge information
    """
    from app.dependencies import check_resource_permission

    # Get job
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # 数据隔离：检查权限（管理员可以访问所有数据，普通用户只能访问自己的数据）
    check_resource_permission(job.user_id, current_user)

    # Get work directory
    if job.work_dir:
        work_dir = Path(job.work_dir)
    elif job.config and job.config.get("work_dir"):
        work_dir = Path(job.config.get("work_dir"))
    elif job.config and job.config.get("job_name"):
        work_dir = settings.MOLYTE_WORK_BASE_PATH / job.config.get("job_name")
    else:
        raise HTTPException(status_code=404, detail="Work directory not found for this job")

    try:
        import re
        molecules = []
        seen_molecules = set()  # 用于去重

        # 获取电荷计算方法
        charge_method = job.config.get("charge_method", "ligpargen") if job.config else "ligpargen"

        # 首先找到所有 .lt 文件（这些是实际使用的分子）
        lt_files = list(work_dir.glob("*.lt"))

        # 获取任务名称，用于过滤系统级别的 .lt 文件
        job_name = job.config.get("job_name", "") if job.config else ""

        for lt_file in lt_files:
            base_name = lt_file.stem  # 例如 "Li", "FSI", "EC"

            # 跳过系统级别的 .lt 文件（通常是整个任务名称）
            if base_name == job_name or len(base_name) > 50:
                logger.info(f"Skipping system-level .lt file: {base_name}")
                continue

            # 跳过已处理的分子
            if base_name in seen_molecules:
                continue
            seen_molecules.add(base_name)

            # 查找对应的 PDB 文件
            # 根据 molyte_command.py 的逻辑：
            # - 溶剂分子使用 .charmm.pdb（因为 RESP 电荷计算和 Packmol 都使用这个文件）
            # - 离子使用 .pdb（从 initial_salts 目录复制的）
            #
            # 优先级顺序：
            # 1. {base_name}.charmm.pdb (溶剂分子，包含 RESP 电荷)
            # 2. {base_name}.q.pdb (LigParGen 生成的电荷文件)
            # 3. {base_name}.pdb (离子或其他分子)
            pdb_file = None
            pdb_candidates = [
                work_dir / f"{base_name}.charmm.pdb",
                work_dir / f"{base_name}.q.pdb",
                work_dir / f"{base_name}.pdb",
            ]

            for candidate in pdb_candidates:
                if candidate.exists():
                    pdb_file = candidate
                    logger.info(f"Found PDB file for {base_name}: {candidate.name}")
                    break

            if not pdb_file:
                logger.warning(f"PDB file not found for {base_name} in {work_dir}")
                logger.warning(f"Tried: .charmm.pdb, .q.pdb, .pdb")
                continue

            # Read PDB content
            try:
                with open(pdb_file, 'r', encoding='utf-8') as f:
                    pdb_content = f.read()
            except UnicodeDecodeError:
                # Try with latin-1 encoding
                with open(pdb_file, 'r', encoding='latin-1') as f:
                    pdb_content = f.read()

            # Parse atoms from PDB
            atoms = []
            for line in pdb_content.split('\n'):
                # Accept both ATOM and HETATM records
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atom_id = int(line[6:11].strip())
                    atom_name = line[12:16].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    element = line[76:78].strip() if len(line) > 76 else atom_name[0]

                    atoms.append({
                        "id": atom_id,
                        "name": atom_name,
                        "element": element,
                        "x": x,
                        "y": y,
                        "z": z,
                        "charge": None  # Will be filled from .lt file
                    })

            # Read charges from .lt file
            with open(lt_file, 'r') as f:
                lt_content = f.read()

            # Extract charges from "In Charges" section or "Data Atoms" section
            charge_map = {}

            # Method 1: Parse from "In Charges" section
            charges_section = re.search(r'write_once\("In Charges"\)\s*\{(.*?)\}', lt_content, re.DOTALL)
            if charges_section:
                for line in charges_section.group(1).split('\n'):
                    match = re.search(r'set type @atom:(\w+)\s+charge\s+([-\d.]+)', line)
                    if match:
                        atom_type = match.group(1)
                        charge = float(match.group(2))
                        charge_map[atom_type] = charge

            # Method 2: Parse from "Data Atoms" section (has charge in 4th column)
            atoms_section = re.search(r'write\("Data Atoms"\)\s*\{(.*?)\}', lt_content, re.DOTALL)
            if atoms_section:
                atom_lines = []
                for line in atoms_section.group(1).split('\n'):
                    line = line.strip()
                    if line and not line.startswith('#'):
                        # Format: $atom:XXX $mol:YYY @atom:ZZZ CHARGE X Y Z
                        match = re.search(r'\$atom:(\w+)\s+\$mol[:\w]*\s+@atom:(\w+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', line)
                        if match:
                            atom_id_or_name = match.group(1)
                            atom_type = match.group(2)
                            charge = float(match.group(3))
                            x = float(match.group(4))
                            y = float(match.group(5))
                            z = float(match.group(6))
                            atom_lines.append({
                                'id_or_name': atom_id_or_name,
                                'type': atom_type,
                                'charge': charge,
                                'x': x,
                                'y': y,
                                'z': z
                            })

                # 尝试按坐标匹配原子
                for i, lt_atom in enumerate(atom_lines):
                    # 方法1: 按坐标匹配（允许小误差）
                    matched = False
                    for pdb_atom in atoms:
                        if (abs(pdb_atom['x'] - lt_atom['x']) < 0.01 and
                            abs(pdb_atom['y'] - lt_atom['y']) < 0.01 and
                            abs(pdb_atom['z'] - lt_atom['z']) < 0.01):
                            pdb_atom['charge'] = lt_atom['charge']
                            matched = True
                            break

                    # 方法2: 如果坐标不匹配，按顺序匹配（假设 PDB 和 LT 的原子顺序一致）
                    if not matched and i < len(atoms):
                        atoms[i]['charge'] = lt_atom['charge']

            # Apply charge_map to atoms by element type
            for atom in atoms:
                if atom['charge'] is None and atom['element'] in charge_map:
                    atom['charge'] = charge_map[atom['element']]

            # Determine molecule type based on total charge
            total_charge = sum(atom['charge'] for atom in atoms if atom['charge'] is not None)

            if total_charge > 0.5:  # Positive charge -> cation
                mol_type = "cation"
            elif total_charge < -0.5:  # Negative charge -> anion
                mol_type = "anion"
            else:  # Neutral -> solvent
                mol_type = "solvent"

            molecules.append({
                "name": base_name,
                "type": mol_type,
                "pdb_content": pdb_content,
                "atoms": atoms,
                "total_charge": total_charge,
                "charge_method": charge_method
            })

        return {"molecules": molecules}

    except Exception as e:
        logger.error(f"Failed to get molecule templates: {e}")
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{job_id}/structure_info")
async def get_structure_info(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取计算后的结构信息（密度、浓度等）
    从 Results/{job_name}/structure_{job_name}.xlsx 文件中读取
    """
    # 获取任务
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # 检查权限（管理员可以访问所有数据）
    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="Not authorized to access this job")

    try:
        # 获取任务名称
        job_name = job.config.get('job_name', f'MD-{job.id}')

        # 构建工作目录
        work_dir = Path(job.work_dir) if job.work_dir else None
        if not work_dir or not work_dir.exists():
            return {"available": False, "message": "Work directory not found"}

        structure_info = {
            "available": False,
            "sample_name": job_name,
            "box_dimensions": None,
            "density": None,
            "concentration": None,
            "initial_density": None,
            "initial_concentration": None,
            "initial_box_dimensions": None,
        }

        # 方法1: 尝试从 Excel 文件读取
        try:
            import openpyxl
            possible_paths = [
                work_dir / "Results" / job_name / f"structure_{job_name}.xlsx",
                work_dir / f"structure_{job_name}.xlsx",
            ]
            structure_file = None
            for path in possible_paths:
                if path.exists():
                    structure_file = path
                    break

            if structure_file:
                wb = openpyxl.load_workbook(structure_file, data_only=True)
                ws = wb.active
                for row_idx, row in enumerate(ws.iter_rows(min_row=1, max_row=5, values_only=True), 1):
                    if len(row) >= 2:
                        label = str(row[0]) if row[0] else ""
                        value = row[1]
                        if "Sample Name" in label:
                            structure_info["sample_name"] = str(value) if value else None
                        elif "Box Dimensions" in label:
                            structure_info["box_dimensions"] = str(value) if value else None
                        elif "Density" in label:
                            try:
                                structure_info["density"] = float(value) if value else None
                            except (ValueError, TypeError):
                                pass
                        elif "Concentration" in label or "concentration" in label.lower():
                            try:
                                structure_info["concentration"] = float(value) if value else None
                            except (ValueError, TypeError):
                                pass
                if structure_info["density"] is not None:
                    structure_info["available"] = True
        except Exception as e:
            logger.warning(f"Failed to read Excel structure file: {e}")

        # 方法2: 如果 Excel 没有数据，尝试从 LAMMPS log 文件读取
        if structure_info["density"] is None:
            log_file = work_dir / f"{job_name}.log"
            if log_file.exists():
                try:
                    with open(log_file, 'r') as f:
                        lines = f.readlines()

                    # 辅助函数：从一行数据中提取密度和盒子尺寸
                    def parse_thermo_line(parts):
                        density_val = None
                        box_vals = []
                        for val in parts:
                            try:
                                fval = float(val)
                                # 密度通常在 0.5 - 3.0 g/cm³ 之间
                                if 0.5 < fval < 3.0 and density_val is None:
                                    density_val = fval
                                # 盒子尺寸通常在 15-100 Å 之间
                                elif 15 < fval < 100:
                                    if len(box_vals) < 3:
                                        box_vals.append(fval)
                            except ValueError:
                                continue
                        return density_val, box_vals

                    # 判断是否为有效数据行
                    skip_prefixes = ('Loop', 'Performance', 'MPI', 'Section',
                                     'Total', 'Pair', 'Bond', 'Kspace', 'Neigh',
                                     'Comm', 'Output', 'Modify', 'Other', 'Nlocal',
                                     'Nghost', 'Neighs', 'Ave', 'Neighbor', 'Dangerous',
                                     'System', 'PPPM', 'WARNING', 'G vector', 'grid',
                                     'stencil', 'estimated', 'using', '3d grid', '-',
                                     'Step', 'Per', 'run', 'thermo', 'fix', 'dump',
                                     'Memory', 'units', 'atom', 'pair', 'kspace',
                                     'Lattice', 'Created', 'Reading', 'Replicate')

                    def is_data_line(line_str):
                        line_str = line_str.strip()
                        if not line_str:
                            return False
                        if line_str.startswith(skip_prefixes):
                            return False
                        parts = line_str.split()
                        if len(parts) < 8:
                            return False
                        try:
                            int(parts[0])  # 第一列应该是时间步（整数）
                            return True
                        except ValueError:
                            return False

                    # 找到第一个 "Step" 表头行的位置（能量最小化阶段的开始，即真正的初始状态）
                    first_step_header_idx = -1
                    for i, line in enumerate(lines):
                        if line.strip().startswith('Step'):
                            first_step_header_idx = i
                            break  # 只找第一个

                    # 找到第一个 "Step" 表头之后的第一行数据（真正的初始状态）
                    first_data_line = None
                    if first_step_header_idx >= 0:
                        for line in lines[first_step_header_idx + 1:]:
                            if is_data_line(line):
                                first_data_line = line.strip().split()
                                break

                    # 找到最后一行数据（最终状态）
                    last_data_line = None
                    for line in reversed(lines[-300:]):
                        if is_data_line(line):
                            last_data_line = line.strip().split()
                            break

                    # 解析初始状态
                    if first_data_line:
                        init_density, init_box = parse_thermo_line(first_data_line)
                        if init_density:
                            structure_info["initial_density"] = round(init_density, 4)
                        if len(init_box) >= 3:
                            structure_info["initial_box_dimensions"] = f"{init_box[0]:.2f} × {init_box[1]:.2f} × {init_box[2]:.2f}"

                    # 解析最终状态
                    if last_data_line:
                        final_density, final_box = parse_thermo_line(last_data_line)
                        if final_density:
                            structure_info["density"] = round(final_density, 4)
                            structure_info["available"] = True
                        if len(final_box) >= 3:
                            structure_info["box_dimensions"] = f"{final_box[0]:.2f} × {final_box[1]:.2f} × {final_box[2]:.2f}"

                        # 计算浓度（初始和最终）
                        try:
                            electrolyte = db.query(ElectrolyteSystem).filter(
                                ElectrolyteSystem.id == job.system_id
                            ).first()

                            if electrolyte:
                                cation_count = sum(c.get('number', 0) for c in (electrolyte.cations or []))
                                avogadro = 6.022e23

                                # 计算初始浓度
                                if init_box and len(init_box) >= 3 and cation_count > 0:
                                    init_volume_L = (init_box[0] * init_box[1] * init_box[2]) * 1e-27
                                    if init_volume_L > 0:
                                        init_conc = (cation_count / avogadro) / init_volume_L
                                        structure_info["initial_concentration"] = round(init_conc, 4)

                                # 计算最终浓度
                                if final_box and len(final_box) >= 3 and cation_count > 0:
                                    final_volume_L = (final_box[0] * final_box[1] * final_box[2]) * 1e-27
                                    if final_volume_L > 0:
                                        final_conc = (cation_count / avogadro) / final_volume_L
                                        structure_info["concentration"] = round(final_conc, 4)
                        except Exception as e:
                            logger.warning(f"Failed to calculate concentration: {e}")

                    logger.info(f"Parsed LAMMPS log for job {job_id}: initial_density={structure_info['initial_density']}, density={structure_info['density']}")

                except Exception as e:
                    logger.warning(f"Failed to parse LAMMPS log file: {e}")
                    import traceback
                    traceback.print_exc()

        if not structure_info["available"]:
            return {"available": False, "message": "Could not extract structure info from job output"}

        logger.info(f"Read structure info for job {job_id}: {structure_info}")
        return structure_info

    except Exception as e:
        logger.error(f"Failed to read structure info: {e}")
        import traceback
        traceback.print_exc()
        return {"available": False, "message": str(e)}


@router.get("/{job_id}/plots/{plot_name}")
async def get_rdf_plot(
    job_id: int,
    plot_name: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取 RDF 图表文件

    Args:
        job_id: 任务 ID
        plot_name: 图表文件名（如 rdf_combined.png, rdf_categorized.png）

    Returns:
        图表文件
    """
    from fastapi.responses import FileResponse

    # 获取任务
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # 检查权限
    if job.user_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized to access this job")

    # 构建图表文件路径
    if not job.work_dir:
        raise HTTPException(status_code=404, detail="Job work directory not found")

    work_dir = Path(job.work_dir)
    plot_file = work_dir / "plots" / plot_name

    if not plot_file.exists():
        raise HTTPException(status_code=404, detail=f"Plot file not found: {plot_name}")

    return FileResponse(
        path=str(plot_file),
        media_type="image/png",
        filename=plot_name
    )


# ============== Slurm 状态查询 API ==============

@router.get("/{job_id}/slurm_status")
def get_job_slurm_status(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取任务的 Slurm 状态

    返回详细的 Slurm 任务状态信息，包括运行时间、CPU 时间等
    """
    from app.services.slurm import get_job_status, normalize_slurm_state

    # 获取任务
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # 检查权限
    if job.user_id != current_user.id and current_user.role != UserRole.ADMIN:
        raise HTTPException(status_code=403, detail="Not authorized to access this job")

    # 获取 Slurm Job ID
    slurm_job_id = job.slurm_job_id or (job.config.get("slurm_job_id") if job.config else None)

    if not slurm_job_id:
        return {
            "job_id": job_id,
            "slurm_job_id": None,
            "status": "NOT_SUBMITTED",
            "message": "任务尚未提交到 Slurm",
            "job_status": job.status.value if job.status else None,
        }

    # 查询 Slurm 状态
    slurm_status = get_job_status(slurm_job_id)

    if not slurm_status:
        return {
            "job_id": job_id,
            "slurm_job_id": slurm_job_id,
            "status": "UNKNOWN",
            "message": "无法获取 Slurm 状态",
            "job_status": job.status.value if job.status else None,
        }

    return {
        "job_id": job_id,
        "slurm_job_id": slurm_job_id,
        "status": normalize_slurm_state(slurm_status.state),
        "raw_state": slurm_status.state,
        "exit_code": slurm_status.exit_code,
        "start_time": slurm_status.start_time,
        "end_time": slurm_status.end_time,
        "elapsed": slurm_status.elapsed,
        "cpu_time": slurm_status.cpu_time,
        "job_status": job.status.value if job.status else None,
    }


# ============== 溶剂化结构分析 API ==============

@router.get("/{job_id}/solvation")
def get_solvation_structures(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取任务的溶剂化结构分析结果

    返回溶剂化结构列表，包括配位数、组成、结构文件路径等
    """
    from app.models.result import SolvationStructure
    from app.dependencies import check_resource_permission

    # 获取任务
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # 数据隔离：检查权限（管理员可以访问所有数据，普通用户只能访问自己的数据）
    check_resource_permission(job.user_id, current_user)

    # 查询溶剂化结构
    structures = db.query(SolvationStructure).filter(
        SolvationStructure.md_job_id == job_id
    ).order_by(SolvationStructure.coordination_num.desc()).all()

    return [
        {
            "id": s.id,
            "center_ion": s.center_ion,
            "structure_type": s.structure_type,
            "coordination_num": s.coordination_num,
            "composition": s.composition,
            "file_path": s.file_path,
            "snapshot_frame": s.snapshot_frame,
            "description": s.description,
            "created_at": s.created_at.isoformat() if s.created_at else None,
        }
        for s in structures
    ]


@router.post("/{job_id}/trigger_postprocess")
def trigger_postprocess(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    手动触发后处理任务（RDF 分析、图表生成等）

    Args:
        job_id: 任务 ID

    Returns:
        Dict with task status
    """
    from app.dependencies import check_resource_permission
    from app.tasks.postprocess import postprocess_md_job_task

    # 获取任务
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # 数据隔离：检查权限
    check_resource_permission(job.user_id, current_user)

    # 检查任务状态
    if job.status not in [JobStatus.COMPLETED, JobStatus.POSTPROCESSING]:
        raise HTTPException(
            status_code=400,
            detail=f"Job must be COMPLETED to trigger postprocessing. Current status: {job.status}"
        )

    try:
        # 触发 Celery 后处理任务
        task = postprocess_md_job_task.delay(job_id)

        logger.info(f"Postprocessing task triggered for job {job_id}: task_id={task.id}")

        return {
            "success": True,
            "message": "Postprocessing task triggered",
            "task_id": task.id,
            "job_id": job_id
        }
    except Exception as e:
        logger.error(f"Failed to trigger postprocessing for job {job_id}: {e}")
        raise HTTPException(
            status_code=500,
            detail=f"Failed to trigger postprocessing: {str(e)}"
        )


@router.post("/{job_id}/solvation/refresh")
def refresh_solvation_structures(
    job_id: int,
    cutoff: float = 3.0,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    重新计算/刷新溶剂化结构分析

    触发溶剂化结构分析，将结果保存到数据库

    Args:
        job_id: 任务 ID
        cutoff: 溶剂化壳层截断距离 (Å)，默认 3.0
    """
    from app.models.result import SolvationStructure
    from app.services.solvation import analyze_solvation_structures, generate_solvation_statistics
    from app.dependencies import check_resource_permission

    # 获取任务
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # 数据隔离：检查权限（管理员可以访问所有数据，普通用户只能访问自己的数据）
    check_resource_permission(job.user_id, current_user)

    # 检查任务状态
    if job.status not in [JobStatus.COMPLETED, JobStatus.POSTPROCESSING]:
        raise HTTPException(
            status_code=400,
            detail=f"Job must be COMPLETED or POSTPROCESSING to analyze solvation structures. Current status: {job.status}"
        )

    # 检查工作目录
    if not job.work_dir:
        raise HTTPException(status_code=400, detail="Job work directory not found")

    # 获取电解液配置
    electrolyte = db.query(ElectrolyteSystem).filter(
        ElectrolyteSystem.id == job.system_id
    ).first()

    if not electrolyte:
        raise HTTPException(status_code=404, detail="Electrolyte system not found")

    electrolyte_data = {
        "cations": electrolyte.cations,
        "anions": electrolyte.anions,
        "solvents": electrolyte.solvents,
    }

    try:
        # 执行溶剂化分析
        logger.info(f"开始溶剂化结构分析: job_id={job_id}, cutoff={cutoff}")
        results = analyze_solvation_structures(
            work_dir=job.work_dir,
            electrolyte_data=electrolyte_data,
            cutoff=cutoff,
        )

        if not results:
            return {
                "success": False,
                "count": 0,
                "message": "未找到溶剂化结构数据，请检查轨迹文件",
            }

        # 删除旧的溶剂化结构记录
        db.query(SolvationStructure).filter(
            SolvationStructure.md_job_id == job_id
        ).delete()

        # 插入新的溶剂化结构记录
        for result in results:
            structure = SolvationStructure(
                md_job_id=job_id,
                center_ion=result['center_ion'],
                structure_type=result['structure_type'],
                coordination_num=result['coordination_num'],
                composition=result['composition'],
                file_path=result['file_path'],
                snapshot_frame=result['snapshot_frame'],
                description=result['description'],
            )
            db.add(structure)

        db.commit()

        # 生成统计信息
        stats = generate_solvation_statistics(results)

        logger.info(f"溶剂化结构分析完成: job_id={job_id}, count={len(results)}")

        return {
            "success": True,
            "count": len(results),
            "statistics": stats,
            "message": f"成功分析 {len(results)} 个溶剂化结构",
        }

    except Exception as e:
        logger.error(f"溶剂化结构分析失败: job_id={job_id}, error={e}")
        db.rollback()
        raise HTTPException(
            status_code=500,
            detail=f"溶剂化结构分析失败: {str(e)}"
        )


@router.get("/{job_id}/solvation/statistics")
def get_solvation_statistics(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取溶剂化结构统计信息

    返回配位数分布、组成分布等统计数据
    """
    from app.models.result import SolvationStructure
    from collections import Counter, defaultdict
    from app.dependencies import check_resource_permission
    import numpy as np

    # 获取任务
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # 数据隔离：检查权限（管理员可以访问所有数据，普通用户只能访问自己的数据）
    check_resource_permission(job.user_id, current_user)

    # 查询溶剂化结构
    structures = db.query(SolvationStructure).filter(
        SolvationStructure.md_job_id == job_id
    ).all()

    if not structures:
        return {
            "total_count": 0,
            "average_coordination_number": 0,
            "coordination_distribution": {},
            "composition_distribution": {},
            "molecule_counts": {},
            "anion_coordination_distribution": {},
        }

    # 配位数分布
    cn_dist = Counter(s.coordination_num for s in structures)

    # 平均配位数
    avg_cn = np.mean([s.coordination_num for s in structures if s.coordination_num])

    # 组成分布
    composition_counter = defaultdict(int)
    molecule_counter = Counter()

    # 获取电解液配置以识别阴离子
    anion_names = set()
    if job.system:
        # 从 ElectrolyteSystem 获取阴离子名称
        if job.system.anions:
            for a in job.system.anions:
                if isinstance(a, dict) and a.get('name'):
                    anion_names.add(a['name'])
                elif isinstance(a, str):
                    anion_names.add(a)

    # 如果没有从电解液配置获取到阴离子，使用常见阴离子名称
    if not anion_names:
        anion_names = {'FSI', 'TFSI', 'PF6', 'BF4', 'ClO4', 'NO3', 'SO4', 'Cl', 'Br', 'I'}

    # 阴离子配位数分布
    anion_cn_dist = Counter()

    for s in structures:
        if s.composition:
            # 生成组成键
            key = "_".join(f"{k}{v}" for k, v in sorted(s.composition.items()) if v and v > 0)
            if key:
                composition_counter[key] += 1

            # 统计各分子数量
            anion_count = 0
            for mol_name, count in s.composition.items():
                if count and count > 0:
                    molecule_counter[mol_name] += count
                    # 统计阴离子数量
                    if mol_name in anion_names:
                        anion_count += count

            anion_cn_dist[anion_count] += 1

    return {
        "total_count": len(structures),
        "average_coordination_number": round(avg_cn, 2) if avg_cn else 0,
        "coordination_distribution": dict(cn_dist),
        "composition_distribution": dict(composition_counter),
        "molecule_counts": dict(molecule_counter),
        "anion_coordination_distribution": dict(anion_cn_dist),
    }


@router.get("/{job_id}/solvation/structure/{structure_id}")
async def get_solvation_structure_file(
    job_id: int,
    structure_id: int,
    format: str = "file",  # "file" 返回文件下载, "content" 返回 XYZ 内容
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取溶剂化结构 XYZ 文件

    Args:
        format: "file" 返回文件下载, "content" 返回 JSON 包含 XYZ 内容
    """
    from fastapi.responses import FileResponse
    from app.models.result import SolvationStructure
    from app.services.solvation import get_structure_xyz_content
    from app.dependencies import check_resource_permission

    # 获取任务
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # 数据隔离：检查权限（管理员可以访问所有数据，普通用户只能访问自己的数据）
    check_resource_permission(job.user_id, current_user)

    # 获取溶剂化结构
    structure = db.query(SolvationStructure).filter(
        SolvationStructure.id == structure_id,
        SolvationStructure.md_job_id == job_id
    ).first()

    if not structure:
        raise HTTPException(status_code=404, detail="Solvation structure not found")

    if not structure.file_path:
        raise HTTPException(status_code=404, detail="Structure file not available")

    file_path = Path(structure.file_path)
    if not file_path.exists():
        raise HTTPException(status_code=404, detail="Structure file not found on disk")

    if format == "content":
        # 返回 XYZ 内容用于 3D 可视化
        xyz_content = get_structure_xyz_content(str(file_path))
        if not xyz_content:
            raise HTTPException(status_code=500, detail="Failed to read structure file")
        return {
            "id": structure.id,
            "center_ion": structure.center_ion,
            "coordination_num": structure.coordination_num,
            "composition": structure.composition,
            "xyz_content": xyz_content,
            "filename": file_path.name,
        }

    return FileResponse(
        path=str(file_path),
        media_type="chemical/x-xyz",
        filename=file_path.name
    )


@router.get("/{job_id}/solvation/system-structure")
def get_system_structure_endpoint(
    job_id: int,
    frame: int = -1,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取整个体系的结构（用于3D可视化）

    Args:
        frame: 帧索引，-1 表示最后一帧
    """
    from app.services.solvation import get_system_structure, get_frame_count
    from app.dependencies import check_resource_permission

    # 获取任务
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # 数据隔离：检查权限（管理员可以访问所有数据，普通用户只能访问自己的数据）
    check_resource_permission(job.user_id, current_user)

    if not job.work_dir:
        raise HTTPException(status_code=400, detail="Job work directory not found")

    result = get_system_structure(job.work_dir, frame)

    if 'error' in result:
        raise HTTPException(status_code=500, detail=result['error'])

    return result


@router.get("/{job_id}/solvation/frame-count")
def get_frame_count_endpoint(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取轨迹文件的帧数
    """
    from app.services.solvation import get_frame_count
    from app.dependencies import check_resource_permission

    # 获取任务
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # 数据隔离：检查权限（管理员可以访问所有数据，普通用户只能访问自己的数据）
    check_resource_permission(job.user_id, current_user)

    if not job.work_dir:
        raise HTTPException(status_code=400, detail="Job work directory not found")

    frame_count = get_frame_count(job.work_dir)

    return {"frame_count": frame_count}


@router.get("/{job_id}/solvation/export-data")
def export_solvation_data(
    job_id: int,
    format: str = "json",  # "json" 或 "csv"
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    导出溶剂化结构数据（用于用户自行绘图）

    Args:
        format: "json" 或 "csv"
    """
    from fastapi.responses import Response
    from app.models.result import SolvationStructure
    from app.dependencies import check_resource_permission
    import csv
    import io

    # 获取任务
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # 数据隔离：检查权限（管理员可以访问所有数据，普通用户只能访问自己的数据）
    check_resource_permission(job.user_id, current_user)

    # 查询溶剂化结构
    structures = db.query(SolvationStructure).filter(
        SolvationStructure.md_job_id == job_id
    ).order_by(SolvationStructure.id).all()

    if not structures:
        raise HTTPException(status_code=404, detail="No solvation data found")

    # 获取所有分子类型
    all_mol_types = set()
    for s in structures:
        if s.composition:
            all_mol_types.update(s.composition.keys())
    all_mol_types = sorted(all_mol_types)

    if format == "csv":
        # 导出为 CSV
        output = io.StringIO()
        writer = csv.writer(output)

        # 表头
        headers = ["id", "center_ion", "coordination_num"] + all_mol_types + ["description"]
        writer.writerow(headers)

        # 数据行
        for s in structures:
            row = [
                s.id,
                s.center_ion,
                s.coordination_num,
            ]
            for mol in all_mol_types:
                row.append(s.composition.get(mol, 0) if s.composition else 0)
            row.append(s.description or "")
            writer.writerow(row)

        csv_content = output.getvalue()
        return Response(
            content=csv_content,
            media_type="text/csv",
            headers={"Content-Disposition": f"attachment; filename=solvation_job{job_id}.csv"}
        )

    # 默认 JSON 格式
    data = {
        "job_id": job_id,
        "total_count": len(structures),
        "molecule_types": all_mol_types,
        "structures": [
            {
                "id": s.id,
                "center_ion": s.center_ion,
                "coordination_num": s.coordination_num,
                "composition": s.composition,
                "description": s.description,
            }
            for s in structures
        ],
        # 统计数据
        "statistics": {
            "coordination_distribution": dict(Counter(s.coordination_num for s in structures)),
            "average_coordination_number": round(
                sum(s.coordination_num for s in structures) / len(structures), 2
            ) if structures else 0,
        }
    }

    return data
