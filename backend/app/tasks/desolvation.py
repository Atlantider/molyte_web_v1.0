"""
去溶剂化能计算任务

计算公式：ΔE_i = E_cluster - (E_cluster_minus_i + E_i)

两阶段处理：
1. 创建所有需要的 QC 任务（cluster, ligands, cluster_minus）
2. 等待所有 QC 任务完成后，计算去溶剂化能
"""
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional
from datetime import datetime
from collections import defaultdict

import numpy as np

from app.database import SessionLocal
from app.models.job import PostprocessJob, JobStatus
from app.models.qc import QCJob, QCJobStatus, QCResult
from app.models.result import SolvationStructure, DesolvationEnergyResult

logger = logging.getLogger(__name__)

# Hartree to kcal/mol conversion factor
HARTREE_TO_KCAL = 627.509474


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

    需要创建的 QC 任务：
    1. E_cluster：完整溶剂化簇
    2. E_ligand_i：每个配体分子（去重）
    3. E_cluster_minus_i：移除每个配体后的簇
    """
    logger.info(f"Phase 1: Creating QC jobs for desolvation job {job.id}")

    # 更新状态为 RUNNING
    job.status = JobStatus.RUNNING
    job.started_at = datetime.now()
    db.commit()

    # 解析溶剂化结构
    cluster_data = parse_solvation_cluster(solvation_structure)
    logger.info(f"Parsed cluster: center={cluster_data['center_ion']}, ligands={len(cluster_data['ligands'])}")

    # 获取计算参数
    basis_set, functional = get_qc_params_for_method_level(method_level)

    # 获取 MD job 和 user
    md_job_id = job.md_job_id
    from app.models.job import MDJob
    md_job = db.query(MDJob).filter(MDJob.id == md_job_id).first()
    if not md_job:
        raise ValueError(f"MD job {md_job_id} not found")
    user_id = md_job.user_id

    created_qc_jobs = []

    # 1. 创建 E_cluster QC 任务
    cluster_qc_job = create_qc_job_for_structure(
        db=db,
        user_id=user_id,
        md_job_id=md_job_id,
        molecule_name=f"Cluster_{solvation_structure.id}",
        xyz_content=cluster_data['xyz_content'],
        charge=cluster_data['total_charge'],
        basis_set=basis_set,
        functional=functional,
        job_type="cluster"
    )
    created_qc_jobs.append(cluster_qc_job.id)
    logger.info(f"Created cluster QC job {cluster_qc_job.id}")

    # 2. 创建每个配体的 QC 任务（去重）
    ligand_qc_jobs = {}
    for ligand in cluster_data['ligands']:
        ligand_key = f"{ligand['ligand_type']}_{ligand['charge']}"

        if ligand_key not in ligand_qc_jobs:
            ligand_qc_job = create_qc_job_for_structure(
                db=db,
                user_id=user_id,
                md_job_id=md_job_id,
                molecule_name=f"{ligand['ligand_type']}",
                xyz_content=ligand['xyz_content'],
                charge=ligand['charge'],
                basis_set=basis_set,
                functional=functional,
                job_type="ligand"
            )
            ligand_qc_jobs[ligand_key] = ligand_qc_job.id
            created_qc_jobs.append(ligand_qc_job.id)
            logger.info(f"Created ligand QC job {ligand_qc_job.id} for {ligand['ligand_type']}")

    # 3. 创建每个 cluster_minus_i 的 QC 任务
    for i, ligand in enumerate(cluster_data['ligands']):
        cluster_minus_xyz = generate_cluster_minus_xyz(cluster_data, ligand)
        cluster_minus_charge = cluster_data['total_charge'] - ligand['charge']

        cluster_minus_qc_job = create_qc_job_for_structure(
            db=db,
            user_id=user_id,
            md_job_id=md_job_id,
            molecule_name=f"Cluster_{solvation_structure.id}_minus_{ligand['ligand_label']}",
            xyz_content=cluster_minus_xyz,
            charge=cluster_minus_charge,
            basis_set=basis_set,
            functional=functional,
            job_type="cluster_minus"
        )
        created_qc_jobs.append(cluster_minus_qc_job.id)
        logger.info(f"Created cluster_minus QC job {cluster_minus_qc_job.id} for ligand {ligand['ligand_label']}")

    # 保存 QC job IDs 到 config
    job.config = job.config or {}
    job.config['qc_job_ids'] = created_qc_jobs
    job.config['cluster_qc_job_id'] = cluster_qc_job.id
    job.config['ligand_qc_jobs'] = ligand_qc_jobs
    job.config['phase'] = 2  # 下次进入阶段 2
    job.config['cluster_data'] = {
        'center_ion': cluster_data['center_ion'],
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

    # 更新状态为 POSTPROCESSING（等待 QC 任务完成）
    job.status = JobStatus.POSTPROCESSING
    job.progress = 10.0

    db.commit()

    logger.info(f"Phase 1 completed: Created {len(created_qc_jobs)} QC jobs")

    return {
        "success": True,
        "job_id": job.id,
        "phase": 1,
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

    前提：所有 QC 任务已完成
    """
    logger.info(f"Phase 2: Calculating desolvation energies for job {job.id}")

    config = job.config or {}
    qc_job_ids = config.get('qc_job_ids', [])
    cluster_qc_job_id = config.get('cluster_qc_job_id')
    ligand_qc_jobs = config.get('ligand_qc_jobs', {})
    cluster_data = config.get('cluster_data', {})

    # 检查所有 QC 任务是否完成
    all_completed = True
    for qc_job_id in qc_job_ids:
        qc_job = db.query(QCJob).filter(QCJob.id == qc_job_id).first()
        if not qc_job or qc_job.status != QCJobStatus.COMPLETED:
            all_completed = False
            logger.info(f"QC job {qc_job_id} not completed yet (status: {qc_job.status if qc_job else 'NOT_FOUND'})")
            break

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

    # 计算每个配体的去溶剂化能
    per_ligand_results = []

    for i, ligand_info in enumerate(cluster_data['ligands']):
        ligand_type = ligand_info['ligand_type']
        ligand_label = ligand_info['ligand_label']
        ligand_charge = ligand_info['charge']

        # 获取 E_ligand
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
        cluster_minus_qc_job_id = qc_job_ids[1 + len(ligand_qc_jobs) + i]  # cluster + ligands + cluster_minus_i
        cluster_minus_result = db.query(QCResult).filter(
            QCResult.qc_job_id == cluster_minus_qc_job_id
        ).first()

        if not cluster_minus_result or cluster_minus_result.total_energy is None:
            raise ValueError(f"Cluster_minus QC result not found for ligand {ligand_label}")

        e_cluster_minus = cluster_minus_result.total_energy

        # 计算 ΔE_i
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

        logger.info(f"Ligand {ligand_label}: ΔE = {delta_e_kcal:.2f} kcal/mol")

    # 按类型汇总
    per_type_summary = summarize_by_type(per_ligand_results)
    logger.info(f"Type summary: {len(per_type_summary)} types")

    # 获取计算参数
    basis_set, functional = get_qc_params_for_method_level(method_level)

    # 保存结果
    desolvation_result = DesolvationEnergyResult(
        postprocess_job_id=job.id,
        solvation_structure_id=solvation_structure.id,
        method_level=method_level,
        basis_set=basis_set,
        functional=functional,
        e_cluster=e_cluster,
        per_ligand_results=per_ligand_results,
        per_type_summary=per_type_summary
    )

    db.add(desolvation_result)

    # 更新任务状态
    job.status = JobStatus.COMPLETED
    job.finished_at = datetime.now()
    job.progress = 100.0

    db.commit()

    logger.info(f"Desolvation job {job.id} completed successfully")

    return {
        "success": True,
        "job_id": job.id,
        "result_id": desolvation_result.id
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


def create_qc_job_for_structure(
    db: SessionLocal,
    user_id: int,
    md_job_id: int,
    molecule_name: str,
    xyz_content: str,
    charge: int,
    basis_set: str,
    functional: str,
    job_type: str
) -> QCJob:
    """创建 QC 任务"""
    # 从 XYZ 生成 SMILES（简化版，实际应该用 RDKit）
    smiles = "C"  # 占位符

    qc_job = QCJob(
        user_id=user_id,
        md_job_id=md_job_id,
        molecule_name=molecule_name,
        smiles=smiles,
        molecule_type="custom",
        basis_set=basis_set,
        functional=functional,
        charge=charge,
        spin_multiplicity=1,  # 默认单重态
        solvent_model="gas",
        status=QCJobStatus.CREATED,
        config={
            "xyz_content": xyz_content,
            "desolvation_job_type": job_type,
            "accuracy_level": "standard"
        }
    )

    db.add(qc_job)
    db.commit()
    db.refresh(qc_job)

    return qc_job


def parse_solvation_cluster(solvation_structure: SolvationStructure) -> Dict[str, Any]:
    """
    解析溶剂化结构，提取中心离子和配体信息
    
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
    center_charge = get_ion_charge(center_ion)
    
    # 根据 composition 识别配体
    composition = solvation_structure.composition or {}
    ligands = identify_ligands(atoms[1:], composition)  # 跳过中心离子
    
    # 计算总电荷
    total_charge = center_charge + sum(lig['charge'] for lig in ligands)
    
    return {
        'center_ion': center_ion,
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


def get_molecule_charge(molecule_type: str) -> int:
    """获取分子电荷"""
    charge_map = {
        'FSI': -1,
        'TFSI': -1,
        'PF6': -1,
        'BF4': -1,
        'ClO4': -1,
        'EC': 0,
        'DMC': 0,
        'EMC': 0,
        'DEC': 0,
        'DME': 0,
        'DOL': 0,
        'TTE': 0,
        'FEC': 0,
    }
    return charge_map.get(molecule_type, 0)


def identify_ligands(atoms: List[Tuple], composition: Dict[str, int]) -> List[Dict[str, Any]]:
    """
    识别配体分子

    简化版本：根据 composition 平均分配原子
    """
    ligands = []
    ligand_id = 0
    atom_idx = 0

    for mol_type, count in composition.items():
        if count == 0:
            continue

        # 估算每个分子的原子数
        atoms_per_mol = estimate_atoms_per_molecule(mol_type)

        for i in range(count):
            ligand_id += 1

            # 提取该配体的原子
            start_idx = atom_idx
            end_idx = min(atom_idx + atoms_per_mol, len(atoms))
            ligand_atoms = atoms[start_idx:end_idx]

            # 生成 XYZ 内容
            xyz_lines = [str(len(ligand_atoms)), f"{mol_type}_{i+1}"]
            for element, x, y, z in ligand_atoms:
                xyz_lines.append(f"{element} {x:.6f} {y:.6f} {z:.6f}")
            xyz_content = '\n'.join(xyz_lines)

            ligands.append({
                'ligand_id': ligand_id,
                'ligand_type': mol_type,
                'ligand_label': f"{mol_type}_{i+1}",
                'atom_indices': list(range(start_idx, end_idx)),
                'xyz_content': xyz_content,
                'charge': get_molecule_charge(mol_type)
            })

            atom_idx = end_idx

    return ligands


def estimate_atoms_per_molecule(mol_type: str) -> int:
    """估算分子的原子数"""
    atom_counts = {
        'FSI': 11,      # F2N-SO2-SO2-NF2
        'TFSI': 15,     # (CF3SO2)2N
        'PF6': 7,
        'BF4': 5,
        'EC': 10,       # C3H4O3
        'DMC': 16,      # C3H6O3
        'EMC': 19,
        'DEC': 22,
        'DME': 20,      # C4H10O2
        'DOL': 13,
        'TTE': 29,      # C6H4F8O2
        'FEC': 11,
    }
    return atom_counts.get(mol_type, 15)  # 默认 15 个原子



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

