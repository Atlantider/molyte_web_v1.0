"""
Molyte 适配器

将数据库中的 ElectrolyteSystem 和 MDJob 转换为 molyte 可以理解的格式
"""

from typing import Dict, Any
import logging

logger = logging.getLogger(__name__)


def convert_electrolyte_to_molyte_format(
    job_name: str,
    job_config: Dict[str, Any],
    electrolyte_data: Dict[str, Any]
) -> Dict[str, Any]:
    """
    将数据库中的 ElectrolyteSystem 转换为 molyte 可以理解的格式
    
    Args:
        job_name: 任务名称（如 "MD-20251119-0001-xxx"）
        job_config: 任务配置（从 MDJob.config 读取）
        electrolyte_data: 电解液数据（从 ElectrolyteSystem 读取）
            {
                "cations": [
                    {"name": "Li", "smiles": "[Li+]", "concentration": 1.0, "number": 50}
                ],
                "anions": [
                    {"name": "PF6", "smiles": "F[P-](F)(F)(F)(F)F", "concentration": 1.0, "number": 50}
                ],
                "solvents": [
                    {"name": "EC", "smiles": "C1COC(=O)O1", "molar_ratio": 5.0, "number": 100},
                    {"name": "DMC", "smiles": "COC(=O)OC", "molar_ratio": 5.0, "number": 100}
                ],
                "temperature": 298.15,
                "pressure": 1.0,
                "box_size": 40.0,
                "nsteps_npt": 5000000,
                "nsteps_nvt": 10000000
            }
    
    Returns:
        {
            "name": "MD-20251119-0001-xxx",
            "box_size": 40.0,
            "temperature_npt": 298.15,
            "temperature_nvt": 298.15,
            "pressure": 1.0,
            "cations": [
                {"name": "Li", "number": 50}
            ],
            "anions": [
                {"name": "PF6", "number": 50}
            ],
            "solvents": [
                {"name": "EC", "smiles": "C1COC(=O)O1", "number": 100},
                {"name": "DMC", "smiles": "COC(=O)OC", "number": 100}
            ],
            "nsteps_npt": 5000000,
            "nsteps_nvt": 10000000,
            "freq_trj_npt": 1000000,
            "freq_trj_nvt": 100000,
            "thermo_freq": 1000000,
            "timestep": 1.0
        }
    """
    # 从 job_config 读取参数（如果有），否则使用 electrolyte 的默认值
    config = job_config or {}
    
    # 基础参数
    result = {
        "name": job_name,
        "box_size": electrolyte_data.get("box_size", 40.0),
        "temperature_npt": config.get("temperature", electrolyte_data.get("temperature", 298.15)),
        "temperature_nvt": config.get("temperature", electrolyte_data.get("temperature", 298.15)),
        "pressure": config.get("pressure", electrolyte_data.get("pressure", 1.0)),
    }
    
    # 阳离子
    result["cations"] = [
        {
            "name": c["name"],
            "number": c.get("number", 50)  # 默认 50 个分子
        }
        for c in electrolyte_data.get("cations", [])
    ]

    # 阴离子
    result["anions"] = [
        {
            "name": a["name"],
            "number": a.get("number", 50)  # 默认 50 个分子
        }
        for a in electrolyte_data.get("anions", [])
    ]
    
    # 溶剂
    result["solvents"] = [
        {
            "name": s["name"],
            "smiles": s["smiles"],
            "number": s.get("number", 100)  # 默认 100 个分子
        }
        for s in electrolyte_data.get("solvents", [])
    ]
    
    # 模拟参数
    result["nsteps_npt"] = config.get("nsteps_npt", electrolyte_data.get("nsteps_npt", 5000000))
    result["nsteps_nvt"] = config.get("nsteps_nvt", electrolyte_data.get("nsteps_nvt", 10000000))
    result["freq_trj_npt"] = config.get("freq_trj_npt", 1000000)
    result["freq_trj_nvt"] = config.get("freq_trj_nvt", 100000)
    result["thermo_freq"] = config.get("thermo_freq", 1000000)
    result["timestep"] = config.get("timestep", 1.0)

    # 电荷计算方法
    result["charge_method"] = config.get("charge_method", "ligpargen")

    # Slurm 资源配置 - 从 job_config 中读取
    result["slurm_partition"] = config.get("slurm_partition", "cpu")
    result["slurm_nodes"] = config.get("slurm_nodes", 1)
    result["slurm_ntasks"] = config.get("slurm_ntasks", 8)
    result["slurm_cpus_per_task"] = config.get("slurm_cpus_per_task", 8)
    result["slurm_time"] = config.get("slurm_time", 7200)

    logger.info(f"Converted electrolyte to molyte format: {job_name}")
    logger.debug(f"  Cations: {len(result['cations'])}, Anions: {len(result['anions'])}, "
                f"Solvents: {len(result['solvents'])}")
    logger.debug(f"  Charge method: {result['charge_method']}")
    logger.debug(f"  Slurm partition: {result['slurm_partition']}, ntasks: {result['slurm_ntasks']}, "
                f"cpus_per_task: {result['slurm_cpus_per_task']}")

    return result


def validate_molyte_data(data: Dict[str, Any]) -> bool:
    """
    验证 molyte 数据格式是否正确
    
    Args:
        data: molyte 格式的数据
    
    Returns:
        True 如果数据格式正确，否则 False
    """
    required_fields = [
        "name", "box_size", "temperature_npt", "temperature_nvt", "pressure",
        "cations", "anions", "solvents",
        "nsteps_npt", "nsteps_nvt"
    ]
    
    for field in required_fields:
        if field not in data:
            logger.error(f"Missing required field: {field}")
            return False
    
    # 验证至少有一个阳离子和一个阴离子
    if not data["cations"]:
        logger.error("No cations specified")
        return False
    
    if not data["anions"]:
        logger.error("No anions specified")
        return False
    
    # 验证至少有一个溶剂
    if not data["solvents"]:
        logger.error("No solvents specified")
        return False
    
    return True

