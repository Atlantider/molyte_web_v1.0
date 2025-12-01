"""
Slurm 服务模块

提供 Slurm 任务状态查询、队列信息获取和资源推荐功能
"""

import logging
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


@dataclass
class SlurmJobStatus:
    """Slurm 任务状态"""
    job_id: str
    state: str
    exit_code: Optional[str] = None
    start_time: Optional[str] = None
    end_time: Optional[str] = None
    elapsed: Optional[str] = None
    cpu_time: Optional[str] = None


@dataclass
class PartitionInfo:
    """Slurm 分区信息"""
    name: str
    state: str
    total_nodes: int
    available_nodes: int
    total_cpus: int
    available_cpus: int
    max_time: Optional[str] = None


@dataclass
class SlurmSuggestion:
    """Slurm 资源推荐"""
    partition: str
    ntasks: int
    cpus_per_task: int
    reason: str


def get_job_status(slurm_job_id: str) -> Optional[SlurmJobStatus]:
    """
    查询 Slurm 任务状态
    
    Args:
        slurm_job_id: Slurm 任务 ID
        
    Returns:
        SlurmJobStatus 对象，如果查询失败返回 None
    """
    try:
        # 使用 sacct 查询任务状态
        cmd = [
            "sacct",
            "-j", slurm_job_id,
            "--format=JobID,State,ExitCode,Start,End,Elapsed,CPUTimeRAW",
            "--parsable2",
            "--noheader",
        ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=10,
        )
        
        if result.returncode != 0:
            logger.warning(f"sacct failed for job {slurm_job_id}: {result.stderr}")
            return None
        
        # 解析输出（取第一行，即主任务）
        lines = result.stdout.strip().split("\n")
        if not lines or not lines[0]:
            logger.warning(f"No sacct output for job {slurm_job_id}")
            return None
        
        fields = lines[0].split("|")
        if len(fields) < 7:
            logger.warning(f"Invalid sacct output for job {slurm_job_id}: {lines[0]}")
            return None
        
        return SlurmJobStatus(
            job_id=fields[0],
            state=fields[1],
            exit_code=fields[2] if fields[2] else None,
            start_time=fields[3] if fields[3] != "Unknown" else None,
            end_time=fields[4] if fields[4] != "Unknown" else None,
            elapsed=fields[5] if fields[5] else None,
            cpu_time=fields[6] if fields[6] else None,
        )
        
    except subprocess.TimeoutExpired:
        logger.error(f"sacct timeout for job {slurm_job_id}")
        return None
    except Exception as e:
        logger.exception(f"Error querying Slurm status for job {slurm_job_id}: {e}")
        return None


def list_partitions() -> List[PartitionInfo]:
    """
    获取所有 Slurm 分区信息

    在云端部署模式下，从 Worker 上报的缓存中获取分区信息。
    如果 Worker 还没有上报，返回默认的分区列表。

    Returns:
        PartitionInfo 列表
    """
    # 尝试从 Worker 缓存获取分区信息
    try:
        from app.api.v1.worker import get_cached_partitions
        cached = get_cached_partitions()
        if cached:
            return [
                PartitionInfo(
                    name=p["name"],
                    state=p["state"],
                    total_nodes=p.get("total_nodes", 0),
                    available_nodes=p.get("available_nodes", 0),
                    total_cpus=p.get("total_cpus", 0),
                    available_cpus=p.get("available_cpus", 0),
                    max_time=p.get("max_time"),
                )
                for p in cached
            ]
    except Exception as e:
        logger.warning(f"Failed to get cached partitions: {e}")

    # 如果没有缓存，返回默认分区（等待 Worker 上报）
    default_partitions = [
        PartitionInfo(
            name="cpu",
            state="up",
            total_nodes=10,
            available_nodes=10,
            total_cpus=320,
            available_cpus=320,
            max_time="7-00:00:00",
        ),
        PartitionInfo(
            name="gpu",
            state="up",
            total_nodes=4,
            available_nodes=4,
            total_cpus=128,
            available_cpus=128,
            max_time="3-00:00:00",
        ),
        PartitionInfo(
            name="debug",
            state="up",
            total_nodes=2,
            available_nodes=2,
            total_cpus=64,
            available_cpus=64,
            max_time="01:00:00",
        ),
    ]

    return default_partitions


def suggest_partition_and_cpus(
    job_type: str = "md",
    expected_runtime_hours: int = 24,
    system_size: int = 1000,
) -> SlurmSuggestion:
    """
    根据任务类型和预期运行时间推荐分区和 CPU 配置

    Args:
        job_type: 任务类型 ("md", "qc", "postprocess")
        expected_runtime_hours: 预期运行时间（小时）
        system_size: 系统大小（原子数）

    Returns:
        SlurmSuggestion 对象
    """
    partitions = list_partitions()

    # 默认配置
    default_partition = "cpu"
    default_ntasks = 8
    default_cpus_per_task = 1

    if not partitions:
        return SlurmSuggestion(
            partition=default_partition,
            ntasks=default_ntasks,
            cpus_per_task=default_cpus_per_task,
            reason="无法获取分区信息，使用默认配置",
        )

    # 根据任务类型选择配置
    if job_type == "postprocess":
        # 后处理任务通常不需要太多资源
        return SlurmSuggestion(
            partition=partitions[0].name,
            ntasks=4,
            cpus_per_task=1,
            reason="后处理任务使用较少资源",
        )

    if job_type == "qc":
        # QC量子化学任务：Gaussian是单节点多核任务
        # 选择可用CPU最多的分区
        best_partition = None
        best_cpus = 0
        for p in partitions:
            if p.state == "up" and p.available_cpus > best_cpus:
                best_cpus = p.available_cpus
                best_partition = p

        if best_partition is None:
            best_partition = partitions[0]

        # QC任务推荐16核，但不超过可用CPU数
        recommended_cpus = min(16, max(best_partition.available_cpus, 4))

        return SlurmSuggestion(
            partition=best_partition.name,
            ntasks=1,  # QC任务通常使用单节点
            cpus_per_task=recommended_cpus,
            reason=f"QC计算推荐使用{recommended_cpus}核，分区{best_partition.name}有{best_partition.available_cpus}个可用CPU",
        )

    # MD 任务：根据系统大小和可用资源推荐
    best_partition = None
    best_score = -1

    for p in partitions:
        if p.state != "up":
            continue

        # 计算分数：优先选择空闲 CPU 多的分区
        score = p.available_cpus

        if score > best_score:
            best_score = score
            best_partition = p

    if best_partition is None:
        best_partition = partitions[0]

    # 根据系统大小推荐 CPU 数
    if system_size < 500:
        ntasks = 4
    elif system_size < 2000:
        ntasks = 8
    elif system_size < 5000:
        ntasks = 16
    else:
        ntasks = 32

    # 确保不超过可用 CPU
    ntasks = min(ntasks, max(best_partition.available_cpus, 4))

    return SlurmSuggestion(
        partition=best_partition.name,
        ntasks=ntasks,
        cpus_per_task=default_cpus_per_task,
        reason=f"基于系统大小({system_size}原子)和分区{best_partition.name}的可用资源({best_partition.available_cpus} CPUs)推荐",
    )


def normalize_slurm_state(slurm_state: str) -> str:
    """
    将 Slurm 状态标准化为统一格式

    Args:
        slurm_state: 原始 Slurm 状态

    Returns:
        标准化后的状态字符串
    """
    slurm_state = slurm_state.upper()

    # 处理 "CANCELLED BY xxx" 格式
    if slurm_state.startswith("CANCELLED"):
        return "CANCELLED"

    state_mapping = {
        "PENDING": "PENDING",
        "CONFIGURING": "PENDING",
        "RESIZING": "PENDING",
        "RUNNING": "RUNNING",
        "COMPLETING": "RUNNING",
        "COMPLETED": "COMPLETED",
        "FAILED": "FAILED",
        "TIMEOUT": "FAILED",
        "OUT_OF_MEMORY": "FAILED",
        "NODE_FAIL": "FAILED",
        "PREEMPTED": "FAILED",
        "BOOT_FAIL": "CANCELLED",
        "DEADLINE": "CANCELLED",
        "REVOKED": "CANCELLED",
    }

    return state_mapping.get(slurm_state, "UNKNOWN")

