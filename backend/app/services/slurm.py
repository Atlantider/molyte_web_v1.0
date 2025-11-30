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
    
    Returns:
        PartitionInfo 列表
    """
    try:
        # 使用 sinfo 获取分区信息
        cmd = [
            "sinfo",
            "--format=%P|%a|%D|%C|%l",
            "--noheader",
        ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=10,
        )
        
        if result.returncode != 0:
            logger.warning(f"sinfo failed: {result.stderr}")
            return []
        
        partitions = []
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            
            fields = line.split("|")
            if len(fields) < 4:
                continue
            
            name = fields[0].rstrip("*")  # 移除默认分区的 * 标记
            state = fields[1]
            total_nodes = int(fields[2]) if fields[2].isdigit() else 0
            
            # 解析 CPU 信息 (格式: A/I/O/T - Allocated/Idle/Other/Total)
            cpu_info = fields[3].split("/")
            if len(cpu_info) >= 4:
                allocated = int(cpu_info[0]) if cpu_info[0].isdigit() else 0
                idle = int(cpu_info[1]) if cpu_info[1].isdigit() else 0
                total = int(cpu_info[3]) if cpu_info[3].isdigit() else 0
            else:
                allocated = idle = total = 0
            
            max_time = fields[4] if len(fields) > 4 else None
            
            partitions.append(PartitionInfo(
                name=name,
                state=state,
                total_nodes=total_nodes,
                available_nodes=total_nodes,  # 简化处理
                total_cpus=total,
                available_cpus=idle,
                max_time=max_time,
            ))
        
        return partitions
        
    except subprocess.TimeoutExpired:
        logger.error("sinfo timeout")
        return []
    except Exception as e:
        logger.exception(f"Error listing partitions: {e}")
        return []


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

