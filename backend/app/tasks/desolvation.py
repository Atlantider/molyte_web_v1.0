"""
去溶剂化能计算任务

计算公式：ΔE_i = E_cluster - (E_cluster_minus_i + E_i)

两阶段处理：
1. 创建所有需要的 QC 任务（cluster, ligands, cluster_minus）
2. 等待所有 QC 任务完成后，计算去溶剂化能

设计原则：
- 充分利用平台已有资源：
  1. 分子 PDB 从 ResultSummary.molecule_structures 获取
  2. 配体原子范围根据 composition 和原子数确定
  3. 单分子能量可以跨任务复用
"""
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime
from collections import defaultdict

import numpy as np

from app.database import SessionLocal
from app.models.job import PostprocessJob, JobStatus, MDJob
from app.models.qc import QCJob, QCJobStatus, QCResult
from app.models.result import SolvationStructure, DesolvationEnergyResult, ResultSummary

logger = logging.getLogger(__name__)

# Hartree to kcal/mol conversion factor
HARTREE_TO_KCAL = 627.509474

# 原子序数映射（用于计算电子数）
ATOMIC_NUMBERS = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
    'K': 19, 'Ca': 20, 'Br': 35, 'I': 53
}


def calculate_spin_multiplicity_from_xyz(xyz_content: str, charge: int) -> int:
    """
    根据 XYZ 内容和电荷自动计算自旋多重度

    逻辑：
    1. 从 XYZ 内容解析原子，计算总电子数
    2. 电子数 = Σ(原子序数) - 电荷
    3. 偶数电子 → 单重态 (spin=1)
    4. 奇数电子 → 双重态 (spin=2)

    Args:
        xyz_content: XYZ 格式的分子坐标
        charge: 分子电荷

    Returns:
        自旋多重度 (1 或 2)
    """
    try:
        lines = xyz_content.strip().split('\n')
        if len(lines) < 3:
            logger.warning(f"XYZ content too short, using default spin=1")
            return 1

        total_electrons = 0
        for line in lines[2:]:  # 跳过原子数和注释行
            parts = line.split()
            if len(parts) >= 4:
                element = parts[0].strip()
                # 处理元素符号（可能有大小写问题）
                element_cap = element.capitalize()
                atomic_num = ATOMIC_NUMBERS.get(element_cap, 0)
                if atomic_num == 0:
                    logger.warning(f"Unknown element: {element}, assuming 0 electrons")
                total_electrons += atomic_num

        # 减去电荷得到实际电子数
        total_electrons -= charge

        # 计算自旋多重度
        if total_electrons % 2 == 0:
            spin = 1  # 闭壳层，单重态
        else:
            spin = 2  # 开壳层，双重态

        logger.debug(f"Spin calculation: {total_electrons} electrons, charge={charge} -> spin={spin}")
        return spin

    except Exception as e:
        logger.warning(f"Failed to calculate spin multiplicity: {e}, using default spin=1")
        return 1


# 分子原子数映射（用于确定配体原子范围）
MOLECULE_ATOM_COUNTS = {
    'Li': 1,
    'Na': 1,
    'K': 1,
    'FSI': 9,       # F2NO4S2 = 2+1+4+2 = 9
    'TFSI': 15,     # C2F6NO4S2 = 2+6+1+4+2 = 15
    'PF6': 7,       # PF6 = 1+6 = 7
    'BF4': 5,       # BF4 = 1+4 = 5
    'ClO4': 5,      # ClO4 = 1+4 = 5
    'EC': 10,       # C3H4O3 = 3+4+3 = 10
    'DMC': 12,      # C3H6O3 = 3+6+3 = 12
    'EMC': 15,      # C4H8O3 = 4+8+3 = 15
    'DEC': 18,      # C5H10O3 = 5+10+3 = 18
    'DME': 16,      # C4H10O2 = 4+10+2 = 16
    'DOL': 11,      # C3H6O2 = 3+6+2 = 11
    'FEC': 10,      # C3H3FO3 = 3+3+1+3 = 10
    'VC': 8,        # C3H2O3 = 3+2+3 = 8
    'TTE': 20,      # C5H4F8O = 实际原子数
}


def run_desolvation_job(job: PostprocessJob, db: SessionLocal) -> Dict[str, Any]:
    """
    执行去溶剂化能计算任务（两阶段处理）

    阶段 1：创建所有需要的 QC 任务
    阶段 2：等待 QC 任务完成后计算去溶剂化能

    Args:
        job: PostprocessJob 对象
        db: 数据库会话

    Returns:
        结果字典 {"success": bool, "job_id": int, ...}
    """
    try:
        logger.info(f"Starting desolvation job {job.id}")

        # 获取配置
        config = job.config or {}
        solvation_structure_id = config.get("solvation_structure_id")
        method_level = config.get("method_level", "standard")

        if not solvation_structure_id:
            raise ValueError("Missing solvation_structure_id in job config")

        # 加载溶剂化结构
        solvation_structure = db.query(SolvationStructure).filter(
            SolvationStructure.id == solvation_structure_id
        ).first()

        if not solvation_structure:
            raise ValueError(f"Solvation structure {solvation_structure_id} not found")

        if not solvation_structure.xyz_content:
            raise ValueError(f"Solvation structure {solvation_structure_id} has no XYZ content")

        logger.info(f"Loaded solvation structure {solvation_structure_id}: {solvation_structure.center_ion}, CN={solvation_structure.coordination_num}")

        # 检查当前阶段
        phase = config.get("phase", 1)

        if phase == 1:
            # 阶段 1：创建 QC 任务
            return _phase1_create_qc_jobs(job, solvation_structure, method_level, db)
        elif phase == 2:
            # 阶段 2：计算去溶剂化能
            return _phase2_calculate_desolvation(job, solvation_structure, method_level, db)
        else:
            raise ValueError(f"Unknown phase: {phase}")

    except Exception as e:
        logger.error(f"Desolvation job {job.id} failed: {e}", exc_info=True)

        job.status = JobStatus.FAILED
        job.finished_at = datetime.now()
        job.error_message = str(e)

        db.commit()

        return {
            "success": False,
            "job_id": job.id,
            "error": str(e)
        }


def _phase1_create_qc_jobs(
    job: PostprocessJob,
    solvation_structure: SolvationStructure,
    method_level: str,
    db: SessionLocal
) -> Dict[str, Any]:
    """
    阶段 1：创建所有需要的 QC 任务

    支持两种去溶剂化模式：
    - stepwise: 逐级去溶剂（每次去掉一个配体）
      - E_cluster: 完整溶剂化簇
      - E_ligand_i: 每种配体分子
      - E_cluster_minus_i: 移除每个配体后的簇
    - full: 全部去溶剂（直接计算中心离子能量）
      - E_cluster: 完整溶剂化簇
      - E_ligand_i: 每种配体分子
      - E_ion: 中心离子能量

    核心改进：
    - 从 ResultSummary.molecule_structures 获取分子 PDB
    - 单分子能量使用标准 PDB 结构，可复用
    - 支持隐式溶剂模型（gas/pcm/smd/custom）
    """
    config = job.config or {}
    desolvation_mode = config.get("desolvation_mode", "stepwise")

    # 获取溶剂配置
    solvent_config = config.get("solvent_config", {})
    solvent_model = solvent_config.get("model", "gas") if solvent_config else "gas"
    solvent_name = solvent_config.get("solvent_name") if solvent_config else None

    logger.info(f"Phase 1: Creating QC jobs for desolvation job {job.id} (mode={desolvation_mode})")

    # 更新状态为 RUNNING
    job.status = JobStatus.RUNNING
    job.started_at = datetime.now()
    db.commit()

    # 获取 MD job 和 user
    md_job_id = job.md_job_id
    md_job = db.query(MDJob).filter(MDJob.id == md_job_id).first()
    if not md_job:
        raise ValueError(f"MD job {md_job_id} not found")
    user_id = md_job.user_id

    # 获取 ResultSummary 中的 molecule_structures
    result_summary = db.query(ResultSummary).filter(
        ResultSummary.md_job_id == md_job_id
    ).first()

    molecule_structures = []
    if result_summary and result_summary.molecule_structures:
        molecule_structures = result_summary.molecule_structures
        logger.info(f"Loaded {len(molecule_structures)} molecule structures from ResultSummary")
    else:
        logger.warning(f"No molecule_structures found in ResultSummary for MD job {md_job_id}")

    # 解析溶剂化结构，传入 molecule_structures
    cluster_data = parse_solvation_cluster(solvation_structure, molecule_structures)
    logger.info(f"Parsed cluster: center={cluster_data['center_ion']}, ligands={len(cluster_data['ligands'])}")

    # 生成更有意义的集群名称
    # 格式: Li_EC1DMC1EMC1PF62_1004 (中心离子_配位组成_结构ID)
    composition = solvation_structure.composition or {}
    composition_str = "".join([
        f"{name}{count}" for name, count in sorted(composition.items()) if count > 0
    ])
    cluster_base_name = f"{cluster_data['center_ion']}_{composition_str}_{solvation_structure.id}"

    # 获取计算参数
    basis_set, functional = get_qc_params_for_method_level(method_level)

    created_qc_jobs = []

    # 1. 创建 E_cluster QC 任务
    cluster_charge = cluster_data['total_charge']
    cluster_spin = calculate_spin_multiplicity_from_xyz(cluster_data['xyz_content'], cluster_charge)
    logger.info(f"Cluster {cluster_base_name}: charge={cluster_charge}, spin={cluster_spin}")

    cluster_qc_job = create_qc_job_for_structure(
        db=db,
        user_id=user_id,
        md_job_id=md_job_id,
        molecule_name=cluster_base_name,
        xyz_content=cluster_data['xyz_content'],
        charge=cluster_charge,
        basis_set=basis_set,
        functional=functional,
        job_type="cluster",
        spin_multiplicity=cluster_spin,
        solvent_model=solvent_model,
        solvent_name=solvent_name,
        solvent_config=solvent_config
    )
    created_qc_jobs.append(cluster_qc_job.id)
    logger.info(f"Created cluster QC job {cluster_qc_job.id}: {cluster_base_name} (charge={cluster_charge}, spin={cluster_spin}, solvent={solvent_model})")

    # 2. 创建每种配体的 QC 任务（去重，支持跨任务复用）
    ligand_qc_jobs = {}
    reused_qc_jobs = []  # 记录复用的 QC 任务

    for ligand in cluster_data['ligands']:
        ligand_type = ligand['ligand_type']
        ligand_charge = ligand['charge']
        ligand_key = f"{ligand_type}_{ligand_charge}"

        if ligand_key not in ligand_qc_jobs:
            # 尝试从 molecule_structures 获取标准 PDB 结构
            ligand_xyz = get_molecule_xyz_from_structures(
                ligand_type,
                molecule_structures,
                fallback_xyz=ligand['xyz_content']
            )

            # 自动计算配体的自旋多重度
            ligand_spin = calculate_spin_multiplicity_from_xyz(ligand_xyz, ligand_charge)
            logger.info(f"Ligand {ligand_type}: charge={ligand_charge}, spin={ligand_spin}")

            # 使用新的 find_or_create 函数，支持跨任务复用
            qc_job_id, is_reused = find_or_create_ligand_qc_job(
                db=db,
                user_id=user_id,
                md_job_id=md_job_id,
                ligand_type=ligand_type,
                ligand_charge=ligand_charge,
                ligand_xyz=ligand_xyz,
                basis_set=basis_set,
                functional=functional,
                spin_multiplicity=ligand_spin,
                solvent_model=solvent_model,
                solvent_name=solvent_name,
                solvent_config=solvent_config
            )

            ligand_qc_jobs[ligand_key] = qc_job_id
            if is_reused:
                reused_qc_jobs.append(qc_job_id)
                logger.info(f"Reused existing QC job {qc_job_id} for {ligand_type} (charge={ligand_charge}, spin={ligand_spin})")
            else:
                created_qc_jobs.append(qc_job_id)
                logger.info(f"Created new ligand QC job {qc_job_id} for {ligand_type} (charge={ligand_charge}, spin={ligand_spin})")

    # 3. 根据模式创建额外的 QC 任务
    cluster_minus_job_ids = []
    center_ion_job_id = None
    ligand_type_mapping = {}  # 配体类型到 cluster_minus job ID 的映射

    if desolvation_mode == "stepwise":
        # 逐级模式：创建每个 cluster_minus_i 的 QC 任务
        # 重要优化：对于等价配体（同类型），只为第一个创建 QC 任务
        # 因为 QC 优化后结构一致，能量相同
        processed_ligand_types = set()  # 已处理的配体类型
        ligand_type_mapping = {}  # 记录每种配体类型对应的 QC job ID

        for i, ligand in enumerate(cluster_data['ligands']):
            ligand_type = ligand['ligand_type']

            # 检查是否已处理过该类型的配体
            if ligand_type in processed_ligand_types:
                # 等价配体，复用已有的 QC 任务
                existing_job_id = ligand_type_mapping.get(ligand_type)
                if existing_job_id:
                    # 记录映射关系，但不创建新任务
                    logger.info(f"Skipping equivalent ligand {ligand['ligand_label']} (same type as processed, will use job {existing_job_id})")
                continue

            processed_ligand_types.add(ligand_type)

            cluster_minus_xyz = generate_cluster_minus_xyz(cluster_data, ligand)
            cluster_minus_charge = cluster_data['total_charge'] - ligand['charge']
            cluster_minus_spin = calculate_spin_multiplicity_from_xyz(cluster_minus_xyz, cluster_minus_charge)

            # 命名格式: Li_EC1DMC1EMC1PF62_1004_minus_EC1
            cluster_minus_name = f"{cluster_base_name}_minus_{ligand['ligand_label']}"
            logger.info(f"Cluster minus {cluster_minus_name}: charge={cluster_minus_charge}, spin={cluster_minus_spin}")

            cluster_minus_qc_job = create_qc_job_for_structure(
                db=db,
                user_id=user_id,
                md_job_id=md_job_id,
                molecule_name=cluster_minus_name,
                xyz_content=cluster_minus_xyz,
                charge=cluster_minus_charge,
                basis_set=basis_set,
                functional=functional,
                job_type="cluster_minus",
                spin_multiplicity=cluster_minus_spin,
                solvent_model=solvent_model,
                solvent_name=solvent_name,
                solvent_config=solvent_config
            )
            cluster_minus_job_ids.append(cluster_minus_qc_job.id)
            created_qc_jobs.append(cluster_minus_qc_job.id)
            ligand_type_mapping[ligand_type] = cluster_minus_qc_job.id
            logger.info(f"Created cluster_minus QC job {cluster_minus_qc_job.id}: {cluster_minus_name} (charge={cluster_minus_charge}, spin={cluster_minus_spin}, representative for {ligand_type})")

    elif desolvation_mode == "full":
        # 全部去溶剂模式：只需要计算中心离子的能量
        center_ion = cluster_data['center_ion']
        center_ion_charge = cluster_data['center_ion_charge']
        center_ion_xyz = generate_center_ion_xyz(cluster_data)

        # 自动计算中心离子的自旋多重度
        center_ion_spin = calculate_spin_multiplicity_from_xyz(center_ion_xyz, center_ion_charge)
        logger.info(f"Center ion {center_ion}: charge={center_ion_charge}, spin={center_ion_spin}")

        # 尝试复用已有的中心离子 QC 计算
        center_ion_job_id, is_reused = find_or_create_ligand_qc_job(
            db=db,
            user_id=user_id,
            md_job_id=md_job_id,
            ligand_type=center_ion,
            ligand_charge=center_ion_charge,
            ligand_xyz=center_ion_xyz,
            basis_set=basis_set,
            functional=functional,
            spin_multiplicity=center_ion_spin,
            solvent_model=solvent_model,
            solvent_name=solvent_name,
            solvent_config=solvent_config
        )

        if is_reused:
            reused_qc_jobs.append(center_ion_job_id)
            logger.info(f"Reused existing QC job {center_ion_job_id} for center ion {center_ion} (charge={center_ion_charge}, spin={center_ion_spin})")
        else:
            created_qc_jobs.append(center_ion_job_id)
            logger.info(f"Created center ion QC job {center_ion_job_id} for {center_ion} (charge={center_ion_charge}, spin={center_ion_spin})")

    # 保存 QC job IDs 到 config
    # 注意：SQLAlchemy 的 JSON 字段需要创建新字典才能检测到变更
    new_config = dict(job.config or {})
    new_config['qc_job_ids'] = created_qc_jobs
    new_config['reused_qc_job_ids'] = reused_qc_jobs  # 复用的 QC 任务
    new_config['cluster_qc_job_id'] = cluster_qc_job.id
    new_config['ligand_qc_jobs'] = ligand_qc_jobs
    new_config['cluster_minus_job_ids'] = cluster_minus_job_ids  # 逐级模式使用
    new_config['cluster_minus_type_mapping'] = ligand_type_mapping if desolvation_mode == "stepwise" else {}  # 配体类型到 cluster_minus job 的映射
    new_config['center_ion_job_id'] = center_ion_job_id  # 全部去溶剂模式使用
    new_config['desolvation_mode'] = desolvation_mode
    new_config['phase'] = 2  # 下次进入阶段 2
    new_config['cluster_data'] = {
        'center_ion': cluster_data['center_ion'],
        'center_ion_charge': cluster_data.get('center_ion_charge', 1),
        'total_charge': cluster_data['total_charge'],
        'ligands': [
            {
                'ligand_id': lig['ligand_id'],
                'ligand_type': lig['ligand_type'],
                'ligand_label': lig['ligand_label'],
                'charge': lig['charge']
            }
            for lig in cluster_data['ligands']
        ]
    }

    # 重新赋值整个 config，确保 SQLAlchemy 检测到变更
    job.config = new_config

    # 更新状态为 POSTPROCESSING（等待 QC 任务完成）
    job.status = JobStatus.POSTPROCESSING
    job.progress = 10.0

    db.commit()

    # 统计日志
    new_ligand_jobs = len(ligand_qc_jobs) - len(reused_qc_jobs)
    if desolvation_mode == "stepwise":
        logger.info(f"Phase 1 completed (stepwise): Created {len(created_qc_jobs)} new QC jobs, reused {len(reused_qc_jobs)} existing jobs")
        logger.info(f"  - 1 cluster, {new_ligand_jobs} new ligand types, {len(reused_qc_jobs)} reused, {len(cluster_minus_job_ids)} cluster_minus")
    else:
        logger.info(f"Phase 1 completed (full): Created {len(created_qc_jobs)} new QC jobs, reused {len(reused_qc_jobs)} existing jobs")
        logger.info(f"  - 1 cluster, {new_ligand_jobs} new ligand types, {len(reused_qc_jobs)} reused, 1 center ion")

    return {
        "success": True,
        "job_id": job.id,
        "phase": 1,
        "desolvation_mode": desolvation_mode,
        "qc_jobs_created": len(created_qc_jobs)
    }


def _phase2_calculate_desolvation(
    job: PostprocessJob,
    solvation_structure: SolvationStructure,
    method_level: str,
    db: SessionLocal
) -> Dict[str, Any]:
    """
    阶段 2：从 QC 结果计算去溶剂化能

    支持两种模式：
    - stepwise: 逐级去溶剂，计算每个配体的去溶剂化能
    - full: 全部去溶剂，计算总去溶剂化能

    前提：所有 QC 任务已完成
    """
    config = job.config or {}
    desolvation_mode = config.get('desolvation_mode', 'stepwise')

    logger.info(f"Phase 2: Calculating desolvation energies for job {job.id} (mode={desolvation_mode})")

    qc_job_ids = config.get('qc_job_ids', [])
    cluster_qc_job_id = config.get('cluster_qc_job_id')
    ligand_qc_jobs = config.get('ligand_qc_jobs', {})
    cluster_minus_job_ids = config.get('cluster_minus_job_ids', [])  # stepwise 模式使用
    center_ion_job_id = config.get('center_ion_job_id')  # full 模式使用
    cluster_data = config.get('cluster_data', {})

    # 检查所有 QC 任务是否完成
    all_completed = True
    any_failed = False
    for qc_job_id in qc_job_ids:
        qc_job = db.query(QCJob).filter(QCJob.id == qc_job_id).first()
        if not qc_job:
            all_completed = False
            logger.error(f"QC job {qc_job_id} not found")
            break
        if qc_job.status == QCJobStatus.FAILED:
            any_failed = True
            logger.error(f"QC job {qc_job_id} failed: {qc_job.error_message}")
            break
        if qc_job.status != QCJobStatus.COMPLETED:
            all_completed = False
            logger.info(f"QC job {qc_job_id} not completed yet (status: {qc_job.status})")
            break

    if any_failed:
        job.status = JobStatus.FAILED
        job.finished_at = datetime.now()
        job.error_message = "One or more QC jobs failed"
        db.commit()
        return {
            "success": False,
            "job_id": job.id,
            "error": "QC job failed"
        }

    if not all_completed:
        logger.info(f"Not all QC jobs completed yet, waiting...")
        return {
            "success": True,
            "job_id": job.id,
            "phase": 2,
            "status": "waiting_for_qc_jobs"
        }

    # 获取 E_cluster
    cluster_result = db.query(QCResult).filter(
        QCResult.qc_job_id == cluster_qc_job_id
    ).first()

    if not cluster_result or cluster_result.total_energy is None:
        raise ValueError(f"Cluster QC result not found or missing energy")

    e_cluster = cluster_result.total_energy
    logger.info(f"E_cluster = {e_cluster:.6f} A.U.")

    per_ligand_results = []
    per_type_summary = []
    e_ion = None  # full 模式使用
    total_desolvation_energy = None  # full 模式使用

    if desolvation_mode == "stepwise":
        # 逐级模式：计算每个配体的去溶剂化能
        # 注意：等价配体（同类型）共享同一个 cluster_minus 能量
        cluster_minus_type_mapping = config.get('cluster_minus_type_mapping', {})

        for i, ligand_info in enumerate(cluster_data['ligands']):
            ligand_type = ligand_info['ligand_type']
            ligand_label = ligand_info['ligand_label']
            ligand_charge = ligand_info['charge']

            # 获取 E_ligand（同类型配体复用相同的能量）
            ligand_key = f"{ligand_type}_{ligand_charge}"
            ligand_qc_job_id = ligand_qc_jobs.get(ligand_key)

            if not ligand_qc_job_id:
                raise ValueError(f"Ligand QC job not found for {ligand_key}")

            ligand_result = db.query(QCResult).filter(
                QCResult.qc_job_id == ligand_qc_job_id
            ).first()

            if not ligand_result or ligand_result.total_energy is None:
                raise ValueError(f"Ligand QC result not found for {ligand_type}")

            e_ligand = ligand_result.total_energy

            # 获取 E_cluster_minus
            # 优先使用类型映射（等价配体优化）
            cluster_minus_qc_job_id = cluster_minus_type_mapping.get(ligand_type)

            if not cluster_minus_qc_job_id:
                # 兼容旧格式：使用索引
                if i < len(cluster_minus_job_ids):
                    cluster_minus_qc_job_id = cluster_minus_job_ids[i]
                else:
                    cluster_minus_qc_job_id = qc_job_ids[1 + len(ligand_qc_jobs) + i]

            cluster_minus_result = db.query(QCResult).filter(
                QCResult.qc_job_id == cluster_minus_qc_job_id
            ).first()

            if not cluster_minus_result or cluster_minus_result.total_energy is None:
                raise ValueError(f"Cluster_minus QC result not found for ligand {ligand_label}")

            e_cluster_minus = cluster_minus_result.total_energy

            # 计算 ΔE_i = E_cluster - (E_cluster_minus + E_ligand)
            delta_e_au = e_cluster - (e_cluster_minus + e_ligand)
            delta_e_kcal = delta_e_au * HARTREE_TO_KCAL

            per_ligand_results.append({
                'ligand_id': ligand_info['ligand_id'],
                'ligand_type': ligand_type,
                'ligand_label': ligand_label,
                'e_ligand': e_ligand,
                'e_cluster_minus': e_cluster_minus,
                'delta_e': delta_e_kcal
            })

            logger.info(f"Ligand {ligand_label}: E_ligand={e_ligand:.6f}, E_cluster_minus={e_cluster_minus:.6f}, ΔE = {delta_e_kcal:.2f} kcal/mol")

        # 按类型汇总
        per_type_summary = summarize_by_type(per_ligand_results)
        logger.info(f"Type summary: {len(per_type_summary)} types")

    elif desolvation_mode == "full":
        # 全部去溶剂模式：计算总去溶剂化能
        # ΔE_total = E_cluster - (E_ion + Σ E_ligand_i)

        # 获取中心离子能量
        if not center_ion_job_id:
            raise ValueError("Center ion QC job ID not found for full mode")

        center_ion_result = db.query(QCResult).filter(
            QCResult.qc_job_id == center_ion_job_id
        ).first()

        if not center_ion_result or center_ion_result.total_energy is None:
            raise ValueError("Center ion QC result not found or missing energy")

        e_ion = center_ion_result.total_energy
        logger.info(f"E_ion = {e_ion:.6f} A.U.")

        # 计算所有配体能量之和
        total_ligand_energy = 0.0
        ligand_energies = {}  # 缓存每种配体的能量

        for ligand_info in cluster_data['ligands']:
            ligand_type = ligand_info['ligand_type']
            ligand_charge = ligand_info['charge']
            ligand_key = f"{ligand_type}_{ligand_charge}"

            if ligand_key not in ligand_energies:
                ligand_qc_job_id = ligand_qc_jobs.get(ligand_key)
                if not ligand_qc_job_id:
                    raise ValueError(f"Ligand QC job not found for {ligand_key}")

                ligand_result = db.query(QCResult).filter(
                    QCResult.qc_job_id == ligand_qc_job_id
                ).first()

                if not ligand_result or ligand_result.total_energy is None:
                    raise ValueError(f"Ligand QC result not found for {ligand_type}")

                ligand_energies[ligand_key] = ligand_result.total_energy

            total_ligand_energy += ligand_energies[ligand_key]

        # 计算总去溶剂化能
        delta_e_total_au = e_cluster - (e_ion + total_ligand_energy)
        total_desolvation_energy = delta_e_total_au * HARTREE_TO_KCAL

        logger.info(f"Total desolvation energy: {total_desolvation_energy:.2f} kcal/mol")
        logger.info(f"  E_cluster = {e_cluster:.6f} A.U.")
        logger.info(f"  E_ion = {e_ion:.6f} A.U.")
        logger.info(f"  Σ E_ligand = {total_ligand_energy:.6f} A.U.")

    # 获取计算参数
    basis_set, functional = get_qc_params_for_method_level(method_level)

    # 保存结果
    result_data = {
        'postprocess_job_id': job.id,
        'solvation_structure_id': solvation_structure.id,
        'method_level': method_level,
        'basis_set': basis_set,
        'functional': functional,
        'e_cluster': e_cluster,
        'per_ligand_results': per_ligand_results,
        'per_type_summary': per_type_summary
    }

    # full 模式额外保存总去溶剂化能
    if desolvation_mode == "full":
        result_data['per_ligand_results'] = [{
            'ligand_id': 0,
            'ligand_type': 'TOTAL',
            'ligand_label': 'Total Desolvation',
            'e_ligand': e_ion,
            'e_cluster_minus': 0.0,  # full 模式不使用
            'delta_e': total_desolvation_energy
        }]
        result_data['per_type_summary'] = [{
            'ligand_type': 'TOTAL',
            'avg_delta_e': total_desolvation_energy,
            'std_delta_e': 0.0,
            'count': 1,
            'min_delta_e': total_desolvation_energy,
            'max_delta_e': total_desolvation_energy
        }]

    desolvation_result = DesolvationEnergyResult(**result_data)

    db.add(desolvation_result)

    # 更新任务状态
    job.status = JobStatus.COMPLETED
    job.finished_at = datetime.now()
    job.progress = 100.0

    db.commit()

    logger.info(f"Desolvation job {job.id} completed successfully (mode={desolvation_mode})")

    return {
        "success": True,
        "job_id": job.id,
        "result_id": desolvation_result.id,
        "desolvation_mode": desolvation_mode
    }


def get_qc_params_for_method_level(method_level: str) -> tuple:
    """根据 method_level 返回 (basis_set, functional)"""
    if method_level == "fast":
        return ("6-31G(d)", "B3LYP")
    elif method_level == "standard":
        return ("6-31++G(d,p)", "B3LYP")
    elif method_level == "accurate":
        return ("6-311++G(2d,2p)", "wB97XD")
    else:
        return ("6-31++G(d,p)", "B3LYP")  # 默认


def find_existing_molecule_qc_job(
    db: SessionLocal,
    molecule_name: str,
    charge: int,
    spin_multiplicity: int,
    basis_set: str,
    functional: str,
    solvent_model: str = "gas",
    solvent_name: Optional[str] = None
) -> Optional[QCJob]:
    """
    查找已完成的相同分子 QC 任务（用于跨任务复用）

    所有计算参数必须完全匹配才能复用：
    - molecule_name: 分子名称
    - charge: 电荷
    - spin_multiplicity: 自旋多重度
    - basis_set: 基组
    - functional: 泛函
    - solvent_model: 溶剂模型 (gas/pcm/smd)
    - solvent_name: 隐式溶剂名称

    Args:
        db: 数据库会话
        molecule_name: 分子名称（如 EC, DMC, FSI）
        charge: 电荷
        spin_multiplicity: 自旋多重度
        basis_set: 基组
        functional: 泛函
        solvent_model: 溶剂模型
        solvent_name: 隐式溶剂名称

    Returns:
        已完成的 QCJob 对象，如果不存在则返回 None
    """
    # 构建查询条件
    query = db.query(QCJob).filter(
        QCJob.molecule_name == molecule_name,
        QCJob.charge == charge,
        QCJob.spin_multiplicity == spin_multiplicity,
        QCJob.basis_set == basis_set,
        QCJob.functional == functional,
        QCJob.solvent_model == solvent_model,
        QCJob.status == QCJobStatus.COMPLETED
    )

    # 处理 solvent_name 匹配（气相时为 None）
    if solvent_model == "gas":
        # 气相时 solvent_name 应为 None 或空
        query = query.filter(
            (QCJob.solvent_name.is_(None)) | (QCJob.solvent_name == "")
        )
    else:
        # 非气相时必须匹配 solvent_name
        query = query.filter(QCJob.solvent_name == solvent_name)

    existing_job = query.first()

    if existing_job:
        # 验证是否有有效的能量结果
        result = db.query(QCResult).filter(
            QCResult.qc_job_id == existing_job.id,
            QCResult.total_energy.isnot(None)
        ).first()

        if result:
            logger.info(f"Found existing completed QC job {existing_job.id} for {molecule_name} (charge={charge})")
            return existing_job

    return None


def create_qc_job_for_structure(
    db: SessionLocal,
    user_id: int,
    md_job_id: int,
    molecule_name: str,
    xyz_content: str,
    charge: int,
    basis_set: str,
    functional: str,
    job_type: str,
    spin_multiplicity: int = 1,
    solvent_model: str = "gas",
    solvent_name: Optional[str] = None,
    solvent_config: Optional[Dict[str, Any]] = None
) -> QCJob:
    """
    创建 QC 任务

    Args:
        db: 数据库会话
        user_id: 用户 ID
        md_job_id: MD 任务 ID
        molecule_name: 分子名称
        xyz_content: XYZ 格式的分子坐标
        charge: 电荷
        basis_set: 基组
        functional: 泛函
        job_type: 任务类型 (cluster/ligand/cluster_minus/center_ion)
        spin_multiplicity: 自旋多重度
        solvent_model: 溶剂模型 (gas/pcm/smd/custom)
        solvent_name: 溶剂名称
        solvent_config: 完整的溶剂配置（包含自定义参数）
    """
    # 从 XYZ 生成 SMILES（简化版，实际应该用 RDKit）
    smiles = "C"  # 占位符

    # 构建 config
    job_config = {
        "xyz_content": xyz_content,
        "desolvation_job_type": job_type,
        "accuracy_level": "standard",
        "solvent_model": solvent_model,
        "solvent_name": solvent_name
    }

    # 如果有自定义溶剂配置，添加到 config
    if solvent_config:
        job_config["solvent_config"] = solvent_config

    qc_job = QCJob(
        user_id=user_id,
        md_job_id=md_job_id,
        molecule_name=molecule_name,
        smiles=smiles,
        molecule_type="custom",
        basis_set=basis_set,
        functional=functional,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        solvent_model=solvent_model,
        solvent_name=solvent_name,
        status=QCJobStatus.CREATED,
        config=job_config
    )

    db.add(qc_job)
    db.commit()
    db.refresh(qc_job)

    return qc_job


def find_or_create_ligand_qc_job(
    db: SessionLocal,
    user_id: int,
    md_job_id: int,
    ligand_type: str,
    ligand_charge: int,
    ligand_xyz: str,
    basis_set: str,
    functional: str,
    spin_multiplicity: int = 1,
    solvent_model: str = "gas",
    solvent_name: Optional[str] = None,
    solvent_config: Optional[Dict[str, Any]] = None
) -> Tuple[int, bool]:
    """
    查找或创建配体分子的 QC 任务（支持跨任务复用）

    复用条件：所有计算参数必须完全匹配

    Args:
        db: 数据库会话
        user_id: 用户 ID
        md_job_id: MD 任务 ID
        ligand_type: 配体类型（如 EC, DMC）
        ligand_charge: 配体电荷
        ligand_xyz: 配体 XYZ 内容
        basis_set: 基组
        functional: 泛函
        spin_multiplicity: 自旋多重度
        solvent_model: 溶剂模型 (gas/pcm/smd/custom)
        solvent_name: 隐式溶剂名称
        solvent_config: 完整的溶剂配置（包含自定义参数）

    Returns:
        (qc_job_id, is_reused): QC 任务 ID 和是否复用标志
    """
    # 1. 尝试查找已完成的相同计算（所有参数必须匹配）
    existing_job = find_existing_molecule_qc_job(
        db=db,
        molecule_name=ligand_type,
        charge=ligand_charge,
        spin_multiplicity=spin_multiplicity,
        basis_set=basis_set,
        functional=functional,
        solvent_model=solvent_model,
        solvent_name=solvent_name
    )

    if existing_job:
        logger.info(f"Reusing existing QC job {existing_job.id} for {ligand_type} "
                   f"(charge={ligand_charge}, spin={spin_multiplicity}, "
                   f"basis={basis_set}, func={functional}, "
                   f"solvent={solvent_model}/{solvent_name})")
        return existing_job.id, True

    # 2. 创建新的 QC 任务
    new_job = create_qc_job_for_structure(
        db=db,
        user_id=user_id,
        md_job_id=md_job_id,
        molecule_name=ligand_type,
        xyz_content=ligand_xyz,
        charge=ligand_charge,
        basis_set=basis_set,
        functional=functional,
        job_type="ligand",
        spin_multiplicity=spin_multiplicity,
        solvent_model=solvent_model,
        solvent_name=solvent_name,
        solvent_config=solvent_config
    )

    logger.info(f"Created new ligand QC job {new_job.id} for {ligand_type}")
    return new_job.id, False


def get_molecule_xyz_from_structures(
    mol_type: str,
    molecule_structures: List[Dict],
    fallback_xyz: Optional[str] = None
) -> str:
    """
    从 molecule_structures 获取分子的 XYZ 内容

    Args:
        mol_type: 分子类型名称（如 EC, DMC, FSI）
        molecule_structures: 分子结构列表（来自 ResultSummary.molecule_structures）
        fallback_xyz: 备选 XYZ（如果找不到 PDB）

    Returns:
        XYZ 格式的分子结构
    """
    # 查找对应的分子结构
    mol_info = None
    for mol in molecule_structures:
        if mol.get('name') == mol_type:
            mol_info = mol
            break

    if not mol_info:
        logger.warning(f"Molecule {mol_type} not found in molecule_structures, using fallback")
        return fallback_xyz or ""

    # 优先使用 atoms 字段（已解析的原子列表）
    atoms = mol_info.get('atoms', [])
    if atoms:
        return atoms_list_to_xyz(atoms, mol_type)

    # 其次尝试解析 pdb_content
    pdb_content = mol_info.get('pdb_content', '')
    if pdb_content:
        return pdb_to_xyz(pdb_content, mol_type)

    logger.warning(f"No atoms or pdb_content for {mol_type}, using fallback")
    return fallback_xyz or ""


def atoms_list_to_xyz(atoms: List[Dict], comment: str = "") -> str:
    """
    将原子列表转换为 XYZ 格式

    Args:
        atoms: 原子列表 [{"element": "C", "x": 0.0, "y": 0.0, "z": 0.0}, ...]
        comment: 注释行

    Returns:
        XYZ 格式字符串
    """
    lines = [str(len(atoms)), comment]

    for atom in atoms:
        element = atom.get('element', 'X')
        x = atom.get('x', 0.0)
        y = atom.get('y', 0.0)
        z = atom.get('z', 0.0)
        lines.append(f"{element} {x:.6f} {y:.6f} {z:.6f}")

    return '\n'.join(lines)


def pdb_to_xyz(pdb_content: str, comment: str = "") -> str:
    """
    将 PDB 内容转换为 XYZ 格式

    Args:
        pdb_content: PDB 文件内容
        comment: 注释行

    Returns:
        XYZ 格式字符串
    """
    atoms = []

    for line in pdb_content.split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # PDB 格式：
            # ATOM      1  C1  MOL     1       0.000   0.000   0.000  1.00  0.00           C
            # 列 1-6: 记录名
            # 列 13-16: 原子名
            # 列 31-38: x 坐标
            # 列 39-46: y 坐标
            # 列 47-54: z 坐标
            # 列 77-78: 元素符号
            try:
                # 尝试从第 77-78 列获取元素符号
                if len(line) >= 78:
                    element = line[76:78].strip()
                else:
                    # 从原子名推断元素
                    atom_name = line[12:16].strip()
                    element = ''.join(c for c in atom_name if c.isalpha())[:2]
                    # 清理元素符号
                    if len(element) > 1:
                        element = element[0].upper() + element[1].lower()
                    else:
                        element = element.upper()

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                atoms.append((element, x, y, z))
            except (ValueError, IndexError) as e:
                logger.warning(f"Failed to parse PDB line: {line}, error: {e}")
                continue

    if not atoms:
        logger.warning("No atoms parsed from PDB content")
        return ""

    lines = [str(len(atoms)), comment]
    for element, x, y, z in atoms:
        lines.append(f"{element} {x:.6f} {y:.6f} {z:.6f}")

    return '\n'.join(lines)


def parse_solvation_cluster(
    solvation_structure: SolvationStructure,
    molecule_structures: Optional[List[Dict]] = None
) -> Dict[str, Any]:
    """
    解析溶剂化结构，提取中心离子和配体信息

    Args:
        solvation_structure: 溶剂化结构对象
        molecule_structures: 分子结构信息列表（从 ResultSummary 获取）

    Returns:
        {
            'center_ion': str,
            'total_charge': int,
            'ligands': [{ligand_id, ligand_type, ligand_label, atom_indices, xyz_content, charge}],
            'xyz_content': str,
            'all_atoms': [(element, x, y, z), ...]
        }
    """
    xyz_content = solvation_structure.xyz_content
    lines = xyz_content.strip().split('\n')

    # 解析 XYZ 格式
    n_atoms = int(lines[0])
    comment = lines[1]

    atoms = []
    for i in range(2, 2 + n_atoms):
        parts = lines[i].split()
        element = parts[0]
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        atoms.append((element, x, y, z))

    # 第一个原子是中心离子
    center_ion = atoms[0][0]
    center_ion_charge = get_ion_charge(center_ion)

    # 获取 mol_order（优先从数据库字段，其次从 XYZ 注释解析）
    mol_order = getattr(solvation_structure, 'mol_order', None)

    # 根据 mol_order、composition 和 molecule_structures 识别配体
    composition = solvation_structure.composition or {}
    ligands = identify_ligands(
        atoms[1:],  # 跳过中心离子
        composition,
        molecule_structures,
        mol_order=mol_order,
        xyz_comment=comment
    )

    # 计算总电荷
    total_charge = center_ion_charge + sum(lig['charge'] for lig in ligands)

    return {
        'center_ion': center_ion,
        'center_ion_charge': center_ion_charge,  # 新增：中心离子电荷
        'total_charge': total_charge,
        'ligands': ligands,
        'xyz_content': xyz_content,
        'all_atoms': atoms
    }


def get_ion_charge(element: str) -> int:
    """获取离子电荷"""
    charge_map = {
        'Li': 1,
        'Na': 1,
        'K': 1,
        'Mg': 2,
        'Ca': 2,
        'Zn': 2,
    }
    return charge_map.get(element, 0)


def get_molecule_charge(
    molecule_type: str,
    smiles: Optional[str] = None,
    molecule_structures: Optional[List[Dict]] = None
) -> int:
    """
    获取分子电荷

    优先级：
    1. 从 SMILES 计算（最准确）
    2. 从 molecule_structures 中查找
    3. 使用硬编码映射（后备）

    Args:
        molecule_type: 分子类型名称（如 EC, FSI）
        smiles: 分子的 SMILES 字符串（可选）
        molecule_structures: 来自 ResultSummary 的分子结构列表（可选）

    Returns:
        分子电荷
    """
    # 1. 尝试从 SMILES 计算电荷（最准确）
    if smiles:
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
                logger.debug(f"Got charge {charge} for {molecule_type} from SMILES")
                return charge
        except Exception as e:
            logger.debug(f"Failed to calculate charge from SMILES for {molecule_type}: {e}")

    # 2. 尝试从 molecule_structures 获取
    if molecule_structures:
        for mol_struct in molecule_structures:
            if mol_struct.get('name') == molecule_type:
                charge = mol_struct.get('total_charge')
                if charge is not None:
                    logger.debug(f"Got charge {charge} for {molecule_type} from molecule_structures")
                    return charge
                # 如果有 SMILES，递归尝试计算
                mol_smiles = mol_struct.get('smiles')
                if mol_smiles:
                    return get_molecule_charge(molecule_type, smiles=mol_smiles)

    # 3. 后备：使用硬编码映射（常见分子）
    # 注意：这只是后备方案，新分子应该通过 SMILES 或 molecule_structures 获取电荷
    charge_map = {
        # 阴离子 (charge = -1)
        'FSI': -1,
        'TFSI': -1,
        'PF6': -1,
        'BF4': -1,
        'ClO4': -1,
        'DCA': -1,
        'NO3': -1,
        'OTf': -1,   # 三氟甲磺酸根
        # 中性溶剂 (charge = 0)
        'EC': 0,
        'DMC': 0,
        'EMC': 0,
        'DEC': 0,
        'DME': 0,
        'DOL': 0,
        'TTE': 0,
        'FEC': 0,
        'VC': 0,
        'PC': 0,     # 碳酸丙烯酯
        'AN': 0,     # 乙腈
        'THF': 0,    # 四氢呋喃
        'DMSO': 0,   # 二甲基亚砜
    }

    charge = charge_map.get(molecule_type, 0)
    if molecule_type not in charge_map:
        logger.warning(f"Unknown molecule type '{molecule_type}', assuming charge=0. "
                       f"Consider providing SMILES for accurate charge calculation.")
    return charge


def identify_ligands(
    atoms: List[Tuple],
    composition: Dict[str, int],
    molecule_structures: Optional[List[Dict]] = None,
    mol_order: Optional[List[Dict]] = None,
    xyz_comment: Optional[str] = None
) -> List[Dict[str, Any]]:
    """
    识别配体分子

    支持两种模式：
    1. 如果提供了 mol_order，直接使用它（最准确）
    2. 否则尝试从 xyz_comment 解析 mol_order
    3. 最后回退到按 composition 顺序（可能不准确）

    Args:
        atoms: 配体原子列表 [(element, x, y, z), ...]，不包含中心离子
        composition: 分子组成 {"EC": 3, "DMC": 1, "FSI": 1}
        molecule_structures: 分子结构信息 [{name, type, pdb_content, atoms, ...}]
        mol_order: XYZ中分子的实际顺序 [{"mol_name": "EC", "atom_count": 10}, ...]
        xyz_comment: XYZ文件的注释行，可能包含 mol_order 信息

    Returns:
        配体列表，每个配体包含 ligand_id, ligand_type, atom_indices, xyz_content 等
    """
    # 尝试从 xyz_comment 解析 mol_order
    if mol_order is None and xyz_comment and "mol_order:" in xyz_comment:
        try:
            # 格式: "... | mol_order:DMC,12;PF6,7;EC,10"
            mol_order_part = xyz_comment.split("mol_order:")[1].strip()
            mol_order = []
            for item in mol_order_part.split(";"):
                if "," in item:
                    mol_name, atom_count = item.split(",")
                    mol_order.append({"mol_name": mol_name.strip(), "atom_count": int(atom_count.strip())})
            logger.info(f"Parsed mol_order from xyz_comment: {mol_order}")
        except Exception as e:
            logger.warning(f"Failed to parse mol_order from xyz_comment: {e}")
            mol_order = None

    # 如果有 mol_order，直接使用它
    if mol_order:
        return _identify_ligands_with_mol_order(atoms, mol_order, molecule_structures)

    # 否则使用旧的基于 composition 的方法（兼容旧数据）
    logger.warning("Using legacy ligand identification (may be inaccurate for existing data)")
    return _identify_ligands_legacy(atoms, composition, molecule_structures)


def _identify_ligands_with_mol_order(
    atoms: List[Tuple],
    mol_order: List[Dict],
    molecule_structures: Optional[List[Dict]] = None
) -> List[Dict[str, Any]]:
    """
    使用 mol_order 精确识别配体分子

    Args:
        atoms: 配体原子列表（不包含中心离子）
        mol_order: [{"mol_name": "EC", "atom_count": 10}, ...]
        molecule_structures: 分子结构信息（用于获取电荷和 SMILES）

    Returns:
        配体列表
    """
    ligands = []
    ligand_id = 0
    atom_idx = 0
    type_counters = defaultdict(int)

    for mol_info in mol_order:
        mol_name = mol_info['mol_name']
        atom_count = mol_info['atom_count']

        ligand_id += 1
        type_counters[mol_name] += 1

        # 确定该配体的原子范围
        start_idx = atom_idx
        end_idx = min(atom_idx + atom_count, len(atoms))
        ligand_atoms = atoms[start_idx:end_idx]

        if len(ligand_atoms) == 0:
            logger.warning(f"No atoms left for ligand {mol_name}_{type_counters[mol_name]}")
            continue

        if len(ligand_atoms) != atom_count:
            logger.warning(f"Ligand {mol_name}_{type_counters[mol_name]}: expected {atom_count} atoms, got {len(ligand_atoms)}")

        # 生成 XYZ 内容
        xyz_lines = [str(len(ligand_atoms)), f"{mol_name}_{type_counters[mol_name]}"]
        for element, x, y, z in ligand_atoms:
            xyz_lines.append(f"{element} {x:.6f} {y:.6f} {z:.6f}")
        xyz_content = '\n'.join(xyz_lines)

        # 获取电荷：优先从 molecule_structures 获取（支持新分子）
        charge = get_molecule_charge(mol_name, molecule_structures=molecule_structures)

        ligands.append({
            'ligand_id': ligand_id,
            'ligand_type': mol_name,
            'ligand_label': f"{mol_name}_{type_counters[mol_name]}",
            'atom_indices': list(range(start_idx, end_idx)),
            'xyz_content': xyz_content,
            'charge': charge
        })

        atom_idx = end_idx

    logger.info(f"Identified {len(ligands)} ligands using mol_order: {dict(type_counters)}")
    return ligands


def _identify_ligands_legacy(
    atoms: List[Tuple],
    composition: Dict[str, int],
    molecule_structures: Optional[List[Dict]] = None
) -> List[Dict[str, Any]]:
    """
    基于连通性分析的分子识别方法

    使用原子间距离判断化学键，然后通过图的连通分量识别独立分子。
    这比基于顺序的方法更可靠。

    注意：输入的 atoms 应该不包含中心离子（Li/Na/K等）。
    中心离子与配体之间是配位键（非共价键），距离约 2-2.5 Å，
    如果包含中心离子，可能会被错误地识别为与某些配体成键。
    """
    import math

    # 共价半径表 (Angstrom)
    # 注意：碱金属（Li/Na/K）和碱土金属（Mg/Ca）设置为 0，
    # 因为它们在溶剂化壳中是中心离子，不应该与配体形成共价键。
    # 这些原子如果出现在配体列表中，会被识别为独立的单原子分子。
    COVALENT_RADII = {
        'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
        'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Br': 1.20, 'I': 1.39,
        'Li': 0.0, 'Na': 0.0, 'K': 0.0,  # 碱金属：不形成共价键
        'Mg': 0.0, 'Ca': 0.0, 'Zn': 0.0,  # 碱土/过渡金属：配位键不是共价键
        'B': 0.84, 'Si': 1.11, 'Se': 1.20, 'Te': 1.38,
    }
    BOND_TOLERANCE = 0.4  # 键长容差

    def get_distance(atom1, atom2):
        """计算两个原子间的距离"""
        _, x1, y1, z1 = atom1
        _, x2, y2, z2 = atom2
        return math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

    def are_bonded(atom1, atom2):
        """判断两个原子是否成键"""
        elem1, _, _, _ = atom1
        elem2, _, _, _ = atom2
        r1 = COVALENT_RADII.get(elem1, 1.5)
        r2 = COVALENT_RADII.get(elem2, 1.5)
        max_bond_length = r1 + r2 + BOND_TOLERANCE
        return get_distance(atom1, atom2) < max_bond_length

    # 使用 Union-Find 算法识别连通分量
    n = len(atoms)
    parent = list(range(n))

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    # 检查所有原子对
    for i in range(n):
        for j in range(i + 1, n):
            if are_bonded(atoms[i], atoms[j]):
                union(i, j)

    # 按连通分量分组
    from collections import defaultdict
    groups = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(i)

    # 将分组转换为配体列表
    molecules = list(groups.values())
    logger.info(f"Connectivity analysis found {len(molecules)} molecules")

    # 构建分子类型 -> 原子数的映射
    mol_atom_counts = {}
    if molecule_structures:
        for mol in molecule_structures:
            mol_name = mol.get('name', '')
            mol_atoms = mol.get('atoms', [])
            if mol_name and mol_atoms:
                mol_atom_counts[mol_name] = len(mol_atoms)

    for mol_type in composition.keys():
        if mol_type not in mol_atom_counts:
            mol_atom_counts[mol_type] = MOLECULE_ATOM_COUNTS.get(mol_type, 15)

    # 根据原子数匹配分子类型
    def identify_molecule_type(mol_atoms_list, atom_indices):
        """根据原子数和组成识别分子类型"""
        n_atoms = len(atom_indices)

        # 计算分子的元素组成
        elements = {}
        for idx in atom_indices:
            elem = mol_atoms_list[idx][0]
            elements[elem] = elements.get(elem, 0) + 1

        # 尝试匹配已知分子类型
        candidates = []
        for mol_type, expected_atoms in mol_atom_counts.items():
            if expected_atoms == n_atoms:
                candidates.append(mol_type)

        if len(candidates) == 1:
            return candidates[0]
        elif len(candidates) > 1:
            # 多个候选，基于元素组成进一步判断
            # 优先选择在 composition 中有需求的
            for mol_type in candidates:
                if composition.get(mol_type, 0) > 0:
                    return mol_type
            return candidates[0]

        # 没有精确匹配，找最接近的
        best_match = None
        min_diff = float('inf')
        for mol_type, expected_atoms in mol_atom_counts.items():
            diff = abs(expected_atoms - n_atoms)
            if diff < min_diff:
                min_diff = diff
                best_match = mol_type

        return best_match or 'Unknown'

    # 构建配体列表
    ligands = []
    type_counters = defaultdict(int)

    for mol_indices in molecules:
        mol_indices = sorted(mol_indices)  # 按索引排序

        # 识别分子类型
        mol_type = identify_molecule_type(atoms, mol_indices)
        type_counters[mol_type] += 1

        # 提取原子
        ligand_atoms = [atoms[i] for i in mol_indices]

        # 生成 XYZ 内容
        xyz_lines = [str(len(ligand_atoms)), f"{mol_type}_{type_counters[mol_type]}"]
        for element, x, y, z in ligand_atoms:
            xyz_lines.append(f"{element} {x:.6f} {y:.6f} {z:.6f}")
        xyz_content = '\n'.join(xyz_lines)

        # 获取电荷：优先从 molecule_structures 获取（支持新分子）
        charge = get_molecule_charge(mol_type, molecule_structures=molecule_structures)

        ligands.append({
            'ligand_id': len(ligands) + 1,
            'ligand_type': mol_type,
            'ligand_label': f"{mol_type}_{type_counters[mol_type]}",
            'atom_indices': mol_indices,
            'xyz_content': xyz_content,
            'charge': charge
        })

    logger.info(f"Identified {len(ligands)} ligands by connectivity: {dict(type_counters)}")
    return ligands



def generate_cluster_minus_xyz(cluster_data: Dict[str, Any], ligand_to_remove: Dict[str, Any]) -> str:
    """
    生成删除指定配体后的簇 XYZ 内容
    """
    all_atoms = cluster_data['all_atoms']
    remove_indices = set(ligand_to_remove['atom_indices'])

    # 保留的原子：中心离子（索引0）+ 其他配体的原子
    kept_atoms = [all_atoms[0]]  # 中心离子

    for i, atom in enumerate(all_atoms[1:], start=1):
        if (i - 1) not in remove_indices:  # i-1 因为 ligand atom_indices 是从 0 开始的（不包括中心离子）
            kept_atoms.append(atom)

    # 生成 XYZ
    xyz_lines = [
        str(len(kept_atoms)),
        f"Cluster minus {ligand_to_remove['ligand_label']}"
    ]

    for element, x, y, z in kept_atoms:
        xyz_lines.append(f"{element} {x:.6f} {y:.6f} {z:.6f}")

    return '\n'.join(xyz_lines)


def generate_center_ion_xyz(cluster_data: Dict[str, Any]) -> str:
    """
    生成中心离子的 XYZ 内容（用于 full 模式）

    中心离子是 cluster 的第一个原子
    """
    all_atoms = cluster_data['all_atoms']

    if not all_atoms:
        raise ValueError("No atoms in cluster data")

    # 中心离子是第一个原子
    center_ion_atom = all_atoms[0]
    element, x, y, z = center_ion_atom

    xyz_lines = [
        "1",
        f"Center ion: {cluster_data.get('center_ion', 'Unknown')}",
        f"{element} {x:.6f} {y:.6f} {z:.6f}"
    ]

    return '\n'.join(xyz_lines)


def summarize_by_type(per_ligand_results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    按配体类型汇总统计

    Returns:
        [{ligand_type, avg_delta_e, std_delta_e, count, min_delta_e, max_delta_e}]
    """
    type_data = defaultdict(list)

    for result in per_ligand_results:
        type_data[result['ligand_type']].append(result['delta_e'])

    summary = []
    for ligand_type, delta_es in type_data.items():
        delta_es_array = np.array(delta_es)

        summary.append({
            'ligand_type': ligand_type,
            'avg_delta_e': float(np.mean(delta_es_array)),
            'std_delta_e': float(np.std(delta_es_array)),
            'count': len(delta_es),
            'min_delta_e': float(np.min(delta_es_array)),
            'max_delta_e': float(np.max(delta_es_array))
        })

    return summary

