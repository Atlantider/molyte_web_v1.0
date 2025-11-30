"""
计算精度等级配置
"""
from enum import Enum
from typing import Dict, Any


class AccuracyLevel(str, Enum):
    """计算精度等级"""
    FAST = "fast"           # 快速模式：适合快速测试
    STANDARD = "standard"   # 标准模式：适合一般研究
    ACCURATE = "accurate"   # 精确模式：适合发表论文
    CUSTOM = "custom"       # 自定义模式：用户自己设置所有参数


class ChargeMethod(str, Enum):
    """电荷计算方法"""
    LIGPARGEN = "ligpargen"  # LigParGen CM1A 电荷
    RESP = "resp"            # RESP 电荷（需要 Gaussian + RESP）


# 精度等级配置
ACCURACY_CONFIGS: Dict[AccuracyLevel, Dict[str, Any]] = {
    AccuracyLevel.FAST: {
        "name": "快速模式",
        "description": "适合快速测试和预览结果，计算时间约 1 小时",
        "charge_method": ChargeMethod.LIGPARGEN,
        "nsteps_npt": 100_000,      # 0.1 ns (timestep=1fs)
        "nsteps_nvt": 500_000,      # 0.5 ns
        "timestep": 1.0,            # fs
        "temperature": 298.15,      # K
        "pressure": 1.0,            # atm
        "freq_trj_npt": 20_000,     # 每 20 ps 输出一次轨迹（NPT平衡阶段，输出稀疏）
        "freq_trj_nvt": 5_000,      # 每 5 ps 输出一次轨迹（NVT生产阶段，输出较密）
        "thermo_freq": 1_000,       # 每 1 ps 输出一次热力学数据（热力学数据量小，可以密集）
        "estimated_time_hours": 1,
        "recommended_for": "快速测试、参数调试、结果预览",
        "color": "#52c41a",         # 绿色
        "icon": "🚀"
    },
    AccuracyLevel.STANDARD: {
        "name": "标准模式",
        "description": "适合一般研究使用，平衡精度和计算时间，约 12 小时",
        "charge_method": ChargeMethod.LIGPARGEN,
        "nsteps_npt": 5_000_000,    # 5 ns
        "nsteps_nvt": 10_000_000,   # 10 ns
        "timestep": 1.0,            # fs
        "temperature": 298.15,      # K
        "pressure": 1.0,            # atm
        "freq_trj_npt": 5_000,      # 每 5 ps 输出一次轨迹（NPT平衡阶段）
        "freq_trj_nvt": 1_000,      # 每 1 ps 输出一次轨迹（NVT生产阶段）
        "thermo_freq": 500,         # 每 0.5 ps 输出一次热力学数据
        "estimated_time_hours": 12,
        "recommended_for": "一般研究、课题组内部讨论、初步分析",
        "color": "#1890ff",         # 蓝色
        "icon": "⚖️"
    },
    AccuracyLevel.ACCURATE: {
        "name": "精确模式",
        "description": "适合发表论文的高精度计算，使用 RESP 电荷，约 24-48 小时",
        "charge_method": ChargeMethod.RESP,
        "nsteps_npt": 10_000_000,   # 10 ns
        "nsteps_nvt": 20_000_000,   # 20 ns
        "timestep": 1.0,            # fs
        "temperature": 298.15,      # K
        "pressure": 1.0,            # atm
        "freq_trj_npt": 2_000,      # 每 2 ps 输出一次轨迹（NPT平衡阶段，适度输出）
        "freq_trj_nvt": 500,        # 每 0.5 ps 输出一次轨迹（NVT生产阶段，高频输出）
        "thermo_freq": 100,         # 每 0.1 ps 输出一次热力学数据（高精度监控）
        "estimated_time_hours": 36,
        "recommended_for": "论文发表、高精度分析、重要结果验证",
        "color": "#f5222d",         # 红色
        "icon": "🎯"
    },
    AccuracyLevel.CUSTOM: {
        "name": "自定义模式",
        "description": "完全自定义所有参数，适合有特殊需求的高级用户",
        "charge_method": ChargeMethod.LIGPARGEN,  # 默认使用 LigParGen
        # 使用标准模式的参数作为参考默认值（用户可以修改）
        "nsteps_npt": 5_000_000,    # 参考值：标准模式
        "nsteps_nvt": 10_000_000,   # 参考值：标准模式
        "timestep": 1.0,            # 参考值：标准模式
        "temperature": 298.15,      # 参考值：标准模式
        "pressure": 1.0,            # 参考值：标准模式
        "freq_trj_npt": 5_000,      # 参考值：标准模式
        "freq_trj_nvt": 1_000,      # 参考值：标准模式
        "thermo_freq": 500,         # 参考值：标准模式
        "estimated_time_hours": None,  # 自定义模式无法预估时间
        "recommended_for": "特殊研究需求、参数优化、高级用户",
        "color": "#722ed1",         # 紫色
        "icon": "⚙️"
    }
}


def get_accuracy_config(level: AccuracyLevel) -> Dict[str, Any]:
    """
    获取指定精度等级的配置
    
    Args:
        level: 精度等级
        
    Returns:
        配置字典
    """
    return ACCURACY_CONFIGS[level].copy()


def get_all_accuracy_levels() -> Dict[str, Dict[str, Any]]:
    """
    获取所有精度等级的配置（用于前端显示）
    
    Returns:
        {
            "fast": {...},
            "standard": {...},
            "accurate": {...}
        }
    """
    return {
        level.value: config
        for level, config in ACCURACY_CONFIGS.items()
    }


def apply_accuracy_level(
    job_config: Dict[str, Any],
    accuracy_level: AccuracyLevel
) -> Dict[str, Any]:
    """
    应用精度等级配置到任务配置

    Args:
        job_config: 原始任务配置
        accuracy_level: 精度等级

    Returns:
        更新后的任务配置
    """
    config = get_accuracy_config(accuracy_level)

    # 如果是自定义模式，不自动填充参数（用户必须自己指定）
    if accuracy_level != AccuracyLevel.CUSTOM:
        # 更新模拟参数（如果用户没有手动指定）
        if "nsteps_npt" not in job_config or job_config["nsteps_npt"] is None:
            job_config["nsteps_npt"] = config["nsteps_npt"]

        if "nsteps_nvt" not in job_config or job_config["nsteps_nvt"] is None:
            job_config["nsteps_nvt"] = config["nsteps_nvt"]

        if "timestep" not in job_config or job_config["timestep"] is None:
            job_config["timestep"] = config["timestep"]

        if "temperature" not in job_config or job_config["temperature"] is None:
            job_config["temperature"] = config["temperature"]

        if "pressure" not in job_config or job_config["pressure"] is None:
            job_config["pressure"] = config["pressure"]

        if "freq_trj_npt" not in job_config or job_config["freq_trj_npt"] is None:
            job_config["freq_trj_npt"] = config["freq_trj_npt"]

        if "freq_trj_nvt" not in job_config or job_config["freq_trj_nvt"] is None:
            job_config["freq_trj_nvt"] = config["freq_trj_nvt"]

        if "thermo_freq" not in job_config or job_config["thermo_freq"] is None:
            job_config["thermo_freq"] = config["thermo_freq"]

    # 添加精度等级和电荷方法信息
    job_config["accuracy_level"] = accuracy_level.value
    job_config["charge_method"] = config["charge_method"].value

    return job_config

