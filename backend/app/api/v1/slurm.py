"""
Slurm 相关 API 端点

提供 Slurm 分区信息查询和资源推荐功能
"""

from fastapi import APIRouter, Depends, HTTPException, Query
from typing import List, Optional
from pydantic import BaseModel

from app.models.user import User, UserRole
from app.dependencies import get_current_active_user
from app.services.slurm import (
    list_partitions,
    suggest_partition_and_cpus,
    PartitionInfo,
    SlurmSuggestion,
)

router = APIRouter()


class PartitionResponse(BaseModel):
    """分区信息响应"""
    name: str
    state: str
    total_nodes: int
    available_nodes: int
    total_cpus: int
    available_cpus: int
    max_time: Optional[str] = None


class SuggestionResponse(BaseModel):
    """资源推荐响应"""
    partition: str
    ntasks: int
    cpus_per_task: int
    reason: str


@router.get("/partitions", response_model=List[PartitionResponse])
def get_partitions(
    current_user: User = Depends(get_current_active_user)
):
    """
    获取用户可访问的 Slurm 分区信息

    - 管理员：返回所有分区
    - 普通用户：只返回 allowed_partitions 中指定的分区

    返回每个分区的节点数、CPU 数、可用资源等信息
    """
    try:
        partitions = list_partitions()
    except Exception as e:
        # 如果 Slurm 命令失败，返回默认分区列表
        partitions = [
            PartitionInfo(
                name="cpu",
                state="up",
                total_nodes=10,
                available_nodes=8,
                total_cpus=320,
                available_cpus=256,
                max_time="7-00:00:00",
            ),
            PartitionInfo(
                name="gpu",
                state="up",
                total_nodes=4,
                available_nodes=3,
                total_cpus=128,
                available_cpus=96,
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

    # 管理员可以看到所有分区
    if current_user.role == UserRole.ADMIN or current_user.allowed_partitions is None:
        filtered_partitions = partitions
    else:
        # 普通用户只能看到被分配的分区
        allowed = current_user.allowed_partitions if isinstance(current_user.allowed_partitions, list) else []
        filtered_partitions = [p for p in partitions if p.name in allowed]

    return [
        PartitionResponse(
            name=p.name,
            state=p.state,
            total_nodes=p.total_nodes,
            available_nodes=p.available_nodes,
            total_cpus=p.total_cpus,
            available_cpus=p.available_cpus,
            max_time=p.max_time,
        )
        for p in filtered_partitions
    ]


@router.get("/suggestion", response_model=SuggestionResponse)
def get_slurm_suggestion(
    job_type: str = Query("md", description="任务类型: md, qc, postprocess"),
    expected_runtime_hours: int = Query(24, description="预期运行时间（小时）"),
    system_size: int = Query(1000, description="系统大小（原子数）"),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取 Slurm 资源推荐配置
    
    根据任务类型、预期运行时间和系统大小推荐最佳的分区和 CPU 配置
    """
    suggestion = suggest_partition_and_cpus(
        job_type=job_type,
        expected_runtime_hours=expected_runtime_hours,
        system_size=system_size,
    )
    
    return SuggestionResponse(
        partition=suggestion.partition,
        ntasks=suggestion.ntasks,
        cpus_per_task=suggestion.cpus_per_task,
        reason=suggestion.reason,
    )


@router.get("/status")
def get_slurm_cluster_status(
    current_user: User = Depends(get_current_active_user)
):
    """
    获取 Slurm 集群整体状态
    
    返回集群的总体资源使用情况
    """
    partitions = list_partitions()
    
    total_cpus = sum(p.total_cpus for p in partitions)
    available_cpus = sum(p.available_cpus for p in partitions)
    total_nodes = sum(p.total_nodes for p in partitions)
    
    return {
        "cluster_status": "online" if partitions else "offline",
        "partition_count": len(partitions),
        "total_nodes": total_nodes,
        "total_cpus": total_cpus,
        "available_cpus": available_cpus,
        "cpu_utilization": round((total_cpus - available_cpus) / total_cpus * 100, 1) if total_cpus > 0 else 0,
        "partitions": [
            {
                "name": p.name,
                "state": p.state,
                "available_cpus": p.available_cpus,
                "total_cpus": p.total_cpus,
            }
            for p in partitions
        ],
    }

