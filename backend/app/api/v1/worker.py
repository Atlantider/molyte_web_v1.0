"""
Worker API 端点

用于轮询 Worker 与云端后端通信
"""

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
from typing import List, Optional, Dict, Any
from datetime import datetime
import logging

from app.database import get_db
from app.models.job import MDJob, JobStatus
from app.models.qc import QCJob, QCJobStatus
from app.models.user import User, UserRole
from app.dependencies import get_current_user
from pydantic import BaseModel


def is_worker_user(user: User) -> bool:
    """检查用户是否是 Worker 用户（ADMIN 角色）"""
    return user.role == UserRole.ADMIN

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/workers", tags=["workers"])


# ==================== 全局缓存 ====================
# 存储 Worker 上报的分区信息
_partition_cache: Dict[str, Any] = {
    "partitions": [],
    "last_updated": None,
    "worker_name": None,
}


def get_cached_partitions() -> List[dict]:
    """获取缓存的分区信息"""
    return _partition_cache["partitions"]


def get_partition_cache_info() -> Dict[str, Any]:
    """获取分区缓存的完整信息"""
    return _partition_cache.copy()


# ==================== Schemas ====================

class PartitionReport(BaseModel):
    """分区信息上报"""
    name: str
    state: str
    total_nodes: int
    available_nodes: int
    total_cpus: int
    available_cpus: int
    max_time: Optional[str] = None


class WorkerHeartbeat(BaseModel):
    """Worker 心跳"""
    worker_name: str
    status: str
    running_jobs: int
    timestamp: str
    partitions: Optional[List[PartitionReport]] = None  # 可选：同时上报分区信息


class JobStatusUpdate(BaseModel):
    """任务状态更新"""
    status: str
    job_type: str  # MD or QC
    worker_name: str
    slurm_job_id: Optional[str] = None
    work_dir: Optional[str] = None
    error_message: Optional[str] = None
    result_files: Optional[List[str]] = None
    progress: Optional[float] = None  # 任务进度 0-100


class PendingJobResponse(BaseModel):
    """待处理任务响应"""
    id: int
    type: str  # MD or QC
    config: dict
    created_at: datetime


# ==================== API Endpoints ====================

@router.post("/heartbeat")
async def worker_heartbeat(
    heartbeat: WorkerHeartbeat,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Worker 心跳接口

    Worker 定期发送心跳，表明自己在线。
    可以同时上报分区信息。
    """
    # 验证是否是 Worker 用户（ADMIN 角色）
    if not is_worker_user(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only worker users can send heartbeat"
        )

    logger.info(
        f"Received heartbeat from {heartbeat.worker_name}: "
        f"{heartbeat.running_jobs} jobs running"
    )

    # 如果心跳中包含分区信息，更新缓存
    if heartbeat.partitions:
        _partition_cache["partitions"] = [p.dict() for p in heartbeat.partitions]
        _partition_cache["last_updated"] = datetime.now().isoformat()
        _partition_cache["worker_name"] = heartbeat.worker_name
        logger.info(f"Updated partition cache with {len(heartbeat.partitions)} partitions from {heartbeat.worker_name}")

    return {"status": "ok", "timestamp": datetime.now().isoformat()}


@router.post("/partitions")
async def report_partitions(
    partitions: List[PartitionReport],
    current_user: User = Depends(get_current_user),
):
    """
    Worker 上报分区信息

    Worker 定期调用此接口上报校园网集群的分区状态。
    云端会缓存这些信息，供前端获取。
    """
    # 验证是否是 Worker 用户（ADMIN 角色）
    if not is_worker_user(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only worker users can report partitions"
        )

    # 更新分区缓存
    _partition_cache["partitions"] = [p.dict() for p in partitions]
    _partition_cache["last_updated"] = datetime.now().isoformat()
    _partition_cache["worker_name"] = current_user.username

    logger.info(f"Received {len(partitions)} partitions from worker {current_user.username}")

    return {
        "status": "ok",
        "partitions_count": len(partitions),
        "timestamp": datetime.now().isoformat()
    }


@router.get("/partitions")
async def get_worker_partitions(
    current_user: User = Depends(get_current_user),
):
    """
    获取 Worker 上报的分区信息

    返回校园网 Worker 最近上报的分区状态。
    """
    cache_info = get_partition_cache_info()

    return {
        "partitions": cache_info["partitions"],
        "last_updated": cache_info["last_updated"],
        "worker_name": cache_info["worker_name"],
    }


@router.get("/jobs/pending", response_model=List[PendingJobResponse])
async def get_pending_jobs(
    job_type: str = "MD",  # MD or QC
    limit: int = 10,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    获取待处理的任务
    
    Worker 轮询此接口获取新任务
    """
    # 验证是否是 Worker 用户（ADMIN 角色）
    if not is_worker_user(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only worker users can fetch pending jobs"
        )
    
    job_type = job_type.upper()
    
    if job_type == "MD":
        # 获取 SUBMITTED 状态的 MD 任务（用户已提交，等待 Worker 处理）
        from app.models.electrolyte import ElectrolyteSystem

        jobs = db.query(MDJob).filter(
            MDJob.status == JobStatus.SUBMITTED
        ).order_by(MDJob.created_at).limit(limit).all()

        result = []
        for job in jobs:
            # 获取电解质系统数据
            electrolyte = db.query(ElectrolyteSystem).filter(
                ElectrolyteSystem.id == job.system_id
            ).first()

            if electrolyte:
                # 合并 config 和电解质数据
                job_name = f"MD-{job.id}"
                if job.config and job.config.get("job_name"):
                    job_name = job.config.get("job_name")

                job_config = {
                    **(job.config or {}),
                    "name": job_name,
                    "cations": electrolyte.cations,
                    "anions": electrolyte.anions,
                    "solvents": electrolyte.solvents,
                    "additives": getattr(electrolyte, 'additives', None),
                    "box_size": electrolyte.box_size,
                    "temperature": electrolyte.temperature,
                    "pressure": electrolyte.pressure,
                }
            else:
                job_config = job.config or {"name": f"MD-{job.id}"}

            result.append(PendingJobResponse(
                id=job.id,
                type="MD",
                config=job_config,
                created_at=job.created_at
            ))

        return result
    
    elif job_type == "QC":
        # 获取 SUBMITTED 状态的 QC 任务（用户已提交，等待 Worker 处理）
        jobs = db.query(QCJob).filter(
            QCJob.status == QCJobStatus.SUBMITTED
        ).order_by(QCJob.created_at).limit(limit).all()

        return [
            PendingJobResponse(
                id=job.id,
                type="QC",
                config={
                    "molecule_name": job.molecule_name,
                    "smiles": job.smiles,
                    "basis_set": job.basis_set,
                    "functional": job.functional,
                    "charge": job.charge,
                    "spin_multiplicity": job.spin_multiplicity,
                    "solvent_model": job.solvent_model,
                    "solvent_name": job.solvent_name,
                    "slurm_partition": job.slurm_partition or "cpu",
                    "slurm_cpus": job.slurm_cpus or 16,
                    "slurm_time": job.slurm_time or 7200,
                    **(job.config or {}),  # 包含额外的配置信息
                },
                created_at=job.created_at
            )
            for job in jobs
        ]
    
    else:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid job type: {job_type}"
        )


@router.put("/jobs/{job_id}/status")
async def update_job_status(
    job_id: int,
    status_update: JobStatusUpdate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    更新任务状态
    
    Worker 在任务状态变化时调用此接口
    """
    # 验证是否是 Worker 用户（ADMIN 角色）
    if not is_worker_user(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only worker users can update job status"
        )
    
    job_type = status_update.job_type.upper()

    # 状态映射：兼容旧版 Worker 发送的 PROCESSING 状态
    status_mapping = {
        "PROCESSING": "QUEUED",  # 旧版 Worker 发送 PROCESSING，映射到 QUEUED
    }
    mapped_status = status_mapping.get(status_update.status, status_update.status)

    # 验证状态值是否有效
    valid_statuses = [s.value for s in JobStatus] if job_type == "MD" else [s.value for s in QCJobStatus]
    if mapped_status not in valid_statuses:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid status: {mapped_status}. Valid statuses: {valid_statuses}"
        )

    if job_type == "MD":
        job = db.query(MDJob).filter(MDJob.id == job_id).first()
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"MD Job {job_id} not found"
            )

        # 更新状态
        job.status = JobStatus[mapped_status]

        if status_update.slurm_job_id:
            job.slurm_job_id = status_update.slurm_job_id

        if status_update.work_dir:
            job.work_dir = status_update.work_dir

        if status_update.error_message:
            job.error_message = status_update.error_message

        if status_update.progress is not None:
            job.progress = status_update.progress

        if mapped_status == "RUNNING" and not job.started_at:
            job.started_at = datetime.now()

        if mapped_status in ["COMPLETED", "FAILED"]:
            job.finished_at = datetime.now()
            # 完成时设置进度为 100%
            if mapped_status == "COMPLETED":
                job.progress = 100.0

        # 如果有结果文件，可以存储到 config 中
        if status_update.result_files:
            if not job.config:
                job.config = {}
            job.config['result_files'] = status_update.result_files

        db.commit()

        logger.info(
            f"MD Job {job_id} status updated to {mapped_status} "
            f"by worker {status_update.worker_name}"
        )

        return {"status": "ok", "job_id": job_id, "new_status": mapped_status}

    elif job_type == "QC":
        job = db.query(QCJob).filter(QCJob.id == job_id).first()
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"QC Job {job_id} not found"
            )

        # 更新状态
        job.status = QCJobStatus[mapped_status]

        if status_update.slurm_job_id:
            job.slurm_job_id = status_update.slurm_job_id

        if status_update.work_dir:
            job.work_dir = status_update.work_dir

        if status_update.error_message:
            job.error_message = status_update.error_message

        if status_update.progress is not None:
            job.progress = status_update.progress

        if mapped_status == "RUNNING" and not job.started_at:
            job.started_at = datetime.now()

        if mapped_status in ["COMPLETED", "FAILED"]:
            job.finished_at = datetime.now()
            # 完成时设置进度为 100%
            if mapped_status == "COMPLETED":
                job.progress = 100.0

        # 如果有结果文件，可以存储到 config 中
        if status_update.result_files:
            if not job.config:
                job.config = {}
            job.config['result_files'] = status_update.result_files

        db.commit()

        logger.info(
            f"QC Job {job_id} status updated to {mapped_status} "
            f"by worker {status_update.worker_name}"
        )

        return {"status": "ok", "job_id": job_id, "new_status": mapped_status}

    else:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid job type: {job_type}"
        )


@router.get("/jobs/{job_id}/check_cancelled")
async def check_job_cancelled(
    job_id: int,
    job_type: str = "MD",
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    检查任务是否已被用户取消

    Worker 定期调用此接口检查任务是否需要取消
    """
    # 验证是否是 Worker 用户（ADMIN 角色）
    if not is_worker_user(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only worker users can check job cancellation"
        )

    job_type = job_type.upper()

    if job_type == "MD":
        job = db.query(MDJob).filter(MDJob.id == job_id).first()
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"MD Job {job_id} not found"
            )

        # 检查是否已被取消
        cancelled = job.status == JobStatus.CANCELLED

        return {
            "job_id": job_id,
            "type": "MD",
            "cancelled": cancelled,
            "status": job.status.value
        }

    elif job_type == "QC":
        job = db.query(QCJob).filter(QCJob.id == job_id).first()
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"QC Job {job_id} not found"
            )

        # 检查是否已被取消
        cancelled = job.status == QCJobStatus.CANCELLED

        return {
            "job_id": job_id,
            "type": "QC",
            "cancelled": cancelled,
            "status": job.status.value
        }

    else:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid job type: {job_type}"
        )


@router.get("/jobs/{job_id}/input")
async def get_job_input_data(
    job_id: int,
    job_type: str = "MD",
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    获取任务输入数据
    
    Worker 可以通过此接口下载任务所需的输入数据
    """
    # 验证是否是 Worker 用户（ADMIN 角色）
    if not is_worker_user(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only worker users can fetch job input data"
        )
    
    job_type = job_type.upper()
    
    if job_type == "MD":
        job = db.query(MDJob).filter(MDJob.id == job_id).first()
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"MD Job {job_id} not found"
            )

        # 获取电解质系统数据
        from app.models.electrolyte import ElectrolyteSystem
        electrolyte = db.query(ElectrolyteSystem).filter(
            ElectrolyteSystem.id == job.system_id
        ).first()

        if not electrolyte:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Electrolyte system {job.system_id} not found"
            )

        # 合并 config 和电解质数据
        job_data = {
            **(job.config or {}),
            "name": job.config.get("job_name") if job.config else f"MD-{job_id}",
            "cations": electrolyte.cations,
            "anions": electrolyte.anions,
            "solvents": electrolyte.solvents,
            "additives": getattr(electrolyte, 'additives', None),
            "box_size": electrolyte.box_size,
            "temperature": electrolyte.temperature,
            "pressure": electrolyte.pressure,
        }

        return {
            "job_id": job_id,
            "type": "MD",
            "config": job_data,
            "system_id": job.system_id
        }
    
    elif job_type == "QC":
        job = db.query(QCJob).filter(QCJob.id == job_id).first()
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"QC Job {job_id} not found"
            )
        
        return {
            "job_id": job_id,
            "type": "QC",
            "molecule_name": job.molecule_name,
            "smiles": job.smiles,
            "basis_set": job.basis_set,
            "functional": job.functional,
            "charge": job.charge,
            "spin_multiplicity": job.spin_multiplicity,
            "solvent_model": job.solvent_model,
            "solvent_name": job.solvent_name,
        }
    
    else:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid job type: {job_type}"
        )


# ==================== QC 结果上传 ====================

class QCResultUpload(BaseModel):
    """QC计算结果上传"""
    energy_au: Optional[float] = None
    homo: Optional[float] = None
    lumo: Optional[float] = None
    homo_lumo_gap: Optional[float] = None
    dipole_moment: Optional[float] = None
    polarizability: Optional[float] = None
    esp_min_kcal: Optional[float] = None
    esp_max_kcal: Optional[float] = None
    # 文件路径（COS/OSS 上的路径）
    fchk_file_path: Optional[str] = None
    log_file_path: Optional[str] = None
    cube_density_path: Optional[str] = None
    cube_esp_path: Optional[str] = None
    cube_homo_path: Optional[str] = None
    cube_lumo_path: Optional[str] = None
    esp_image_path: Optional[str] = None
    homo_image_path: Optional[str] = None
    lumo_image_path: Optional[str] = None
    # 额外属性
    additional_properties: Optional[Dict[str, Any]] = None


@router.post("/jobs/{job_id}/qc_result")
async def upload_qc_result(
    job_id: int,
    result_data: QCResultUpload,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Worker 上传 QC 计算结果

    Worker 解析 Gaussian 输出后调用此接口保存结果到数据库
    """
    from app.models.qc import QCResult

    # 验证是否是 Worker 用户（ADMIN 角色）
    if not is_worker_user(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only worker users can upload QC results"
        )

    # 获取 QC 任务
    job = db.query(QCJob).filter(QCJob.id == job_id).first()
    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"QC Job {job_id} not found"
        )

    # 检查是否已有结果
    existing_result = db.query(QCResult).filter(QCResult.qc_job_id == job_id).first()

    if existing_result:
        # 更新已有结果
        for field, value in result_data.dict(exclude_unset=True).items():
            if value is not None:
                setattr(existing_result, field, value)
        db.commit()
        logger.info(f"Updated QC result for job {job_id}")
        return {"status": "ok", "message": "Result updated", "result_id": existing_result.id}
    else:
        # 创建新结果
        new_result = QCResult(
            qc_job_id=job_id,
            smiles=job.smiles,
            energy_au=result_data.energy_au,
            homo=result_data.homo,
            lumo=result_data.lumo,
            homo_lumo_gap=result_data.homo_lumo_gap,
            dipole_moment=result_data.dipole_moment,
            polarizability=result_data.polarizability,
            esp_min_kcal=result_data.esp_min_kcal,
            esp_max_kcal=result_data.esp_max_kcal,
            fchk_file_path=result_data.fchk_file_path,
            log_file_path=result_data.log_file_path,
            cube_density_path=result_data.cube_density_path,
            cube_esp_path=result_data.cube_esp_path,
            cube_homo_path=result_data.cube_homo_path,
            cube_lumo_path=result_data.cube_lumo_path,
            esp_image_path=result_data.esp_image_path,
            homo_image_path=result_data.homo_image_path,
            lumo_image_path=result_data.lumo_image_path,
            additional_properties=result_data.additional_properties or {},
        )
        db.add(new_result)
        db.commit()
        db.refresh(new_result)
        logger.info(f"Created QC result for job {job_id}: id={new_result.id}")
        return {"status": "ok", "message": "Result created", "result_id": new_result.id}


# ==================== MD 结果上传（RDF、MSD等） ====================

class RDFResultUpload(BaseModel):
    """RDF 结果上传"""
    center_species: str
    shell_species: str
    r_values: List[float]
    g_r_values: List[float]
    coordination_number_values: Optional[List[float]] = None
    first_peak_position: Optional[float] = None
    first_peak_height: Optional[float] = None
    coordination_number: Optional[float] = None


class MSDResultUpload(BaseModel):
    """MSD 结果上传"""
    species: str
    t_values: List[float]
    msd_x_values: Optional[List[float]] = None
    msd_y_values: Optional[List[float]] = None
    msd_z_values: Optional[List[float]] = None
    msd_total_values: List[float]
    labels: Optional[Dict[str, str]] = None
    diffusion_coefficient: Optional[float] = None
    ionic_conductivity: Optional[float] = None
    mobility: Optional[float] = None
    charge: Optional[int] = None


class SolvationStructureUpload(BaseModel):
    """溶剂化结构上传"""
    center_ion: str
    structure_type: Optional[str] = 'first_shell'
    coordination_num: int
    composition: Dict[str, int]  # {"EC": 3, "DMC": 1, "FSI": 0}
    file_path: Optional[str] = None
    snapshot_frame: Optional[int] = None
    description: Optional[str] = None


class MDResultsUpload(BaseModel):
    """MD 任务结果上传（包含 RDF、MSD、溶剂化结构等）"""
    rdf_results: Optional[List[RDFResultUpload]] = None
    msd_results: Optional[List[MSDResultUpload]] = None
    solvation_structures: Optional[List[SolvationStructureUpload]] = None
    # 结果摘要
    final_density: Optional[float] = None
    final_temperature: Optional[float] = None
    final_pressure: Optional[float] = None
    total_energy: Optional[float] = None
    potential_energy: Optional[float] = None
    kinetic_energy: Optional[float] = None
    total_atoms: Optional[int] = None
    total_molecules: Optional[int] = None


@router.post("/jobs/{job_id}/md_results")
async def upload_md_results(
    job_id: int,
    results_data: MDResultsUpload,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Worker 上传 MD 计算结果（RDF、MSD 等）

    Worker 在 MD 任务完成后，解析结果并调用此接口保存到数据库
    """
    from app.models.result import RDFResult, MSDResult, ResultSummary

    # 验证是否是 Worker 用户（ADMIN 角色）
    if not is_worker_user(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only worker users can upload MD results"
        )

    # 获取 MD 任务
    job = db.query(MDJob).filter(MDJob.id == job_id).first()
    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"MD Job {job_id} not found"
        )

    uploaded_counts = {"rdf": 0, "msd": 0, "summary": False}

    # 1. 保存 RDF 结果
    if results_data.rdf_results:
        # 先删除旧的 RDF 结果
        db.query(RDFResult).filter(RDFResult.md_job_id == job_id).delete()

        for rdf_data in results_data.rdf_results:
            rdf_record = RDFResult(
                md_job_id=job_id,
                center_species=rdf_data.center_species,
                shell_species=rdf_data.shell_species,
                r_values=rdf_data.r_values,
                g_r_values=rdf_data.g_r_values,
                coordination_number_values=rdf_data.coordination_number_values,
                first_peak_position=rdf_data.first_peak_position,
                first_peak_height=rdf_data.first_peak_height,
                coordination_number=rdf_data.coordination_number,
            )
            db.add(rdf_record)
            uploaded_counts["rdf"] += 1

    # 2. 保存 MSD 结果
    if results_data.msd_results:
        # 先删除旧的 MSD 结果
        db.query(MSDResult).filter(MSDResult.md_job_id == job_id).delete()

        for msd_data in results_data.msd_results:
            msd_record = MSDResult(
                md_job_id=job_id,
                species=msd_data.species,
                t_values=msd_data.t_values,
                msd_x_values=msd_data.msd_x_values,
                msd_y_values=msd_data.msd_y_values,
                msd_z_values=msd_data.msd_z_values,
                msd_total_values=msd_data.msd_total_values,
                labels=msd_data.labels,
                diffusion_coefficient=msd_data.diffusion_coefficient,
                ionic_conductivity=msd_data.ionic_conductivity,
                mobility=msd_data.mobility,
                charge=msd_data.charge,
            )
            db.add(msd_record)
            uploaded_counts["msd"] += 1

    # 3. 保存溶剂化结构
    if results_data.solvation_structures:
        from app.models.result import SolvationStructure

        # 先删除旧的溶剂化结构
        db.query(SolvationStructure).filter(SolvationStructure.md_job_id == job_id).delete()

        for solv_data in results_data.solvation_structures:
            solv_record = SolvationStructure(
                md_job_id=job_id,
                center_ion=solv_data.center_ion,
                structure_type=solv_data.structure_type,
                coordination_num=solv_data.coordination_num,
                composition=solv_data.composition,
                file_path=solv_data.file_path,
                snapshot_frame=solv_data.snapshot_frame,
                description=solv_data.description,
            )
            db.add(solv_record)
        uploaded_counts["solvation"] = len(results_data.solvation_structures)

    # 4. 保存结果摘要
    if any([
        results_data.final_density,
        results_data.final_temperature,
        results_data.total_energy,
        results_data.total_atoms
    ]):
        # 检查是否已有摘要
        existing_summary = db.query(ResultSummary).filter(
            ResultSummary.md_job_id == job_id
        ).first()

        if existing_summary:
            # 更新已有摘要
            if results_data.final_density is not None:
                existing_summary.final_density = results_data.final_density
            if results_data.final_temperature is not None:
                existing_summary.final_temperature = results_data.final_temperature
            if results_data.final_pressure is not None:
                existing_summary.final_pressure = results_data.final_pressure
            if results_data.total_energy is not None:
                existing_summary.total_energy = results_data.total_energy
            if results_data.potential_energy is not None:
                existing_summary.potential_energy = results_data.potential_energy
            if results_data.kinetic_energy is not None:
                existing_summary.kinetic_energy = results_data.kinetic_energy
            if results_data.total_atoms is not None:
                existing_summary.total_atoms = results_data.total_atoms
            if results_data.total_molecules is not None:
                existing_summary.total_molecules = results_data.total_molecules
        else:
            # 创建新摘要
            summary = ResultSummary(
                md_job_id=job_id,
                final_density=results_data.final_density,
                final_temperature=results_data.final_temperature,
                final_pressure=results_data.final_pressure,
                total_energy=results_data.total_energy,
                potential_energy=results_data.potential_energy,
                kinetic_energy=results_data.kinetic_energy,
                total_atoms=results_data.total_atoms,
                total_molecules=results_data.total_molecules,
            )
            db.add(summary)
        uploaded_counts["summary"] = True

    # 5. 更新任务配置，记录后处理结果
    if job.config is None:
        job.config = {}

    if "postprocess" not in job.config:
        job.config["postprocess"] = {}

    job.config["postprocess"]["rdf_count"] = uploaded_counts["rdf"]
    job.config["postprocess"]["msd_count"] = uploaded_counts["msd"]
    job.config["postprocess"]["solvation_count"] = uploaded_counts.get("solvation", 0)
    job.config["postprocess"]["completed_at"] = datetime.now().isoformat()
    job.config["postprocess"]["uploaded_by"] = "worker"

    db.commit()

    logger.info(
        f"MD results uploaded for job {job_id}: "
        f"RDF={uploaded_counts['rdf']}, MSD={uploaded_counts['msd']}, "
        f"Solvation={uploaded_counts.get('solvation', 0)}, Summary={uploaded_counts['summary']}"
    )

    return {
        "status": "ok",
        "job_id": job_id,
        "uploaded": uploaded_counts
    }

