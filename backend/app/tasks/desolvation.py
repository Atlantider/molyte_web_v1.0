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

    需要创建的 QC 任务：
    1. E_cluster：完整溶剂化簇
    2. E_ligand_i：每种配体分子（去重，使用 PDB 中的结构）
    3. E_cluster_minus_i：移除每个配体后的簇

    核心改进：
    - 从 ResultSummary.molecule_structures 获取分子 PDB
    - 单分子能量使用标准 PDB 结构，可复用
    """
    logger.info(f"Phase 1: Creating QC jobs for desolvation job {job.id}")

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

    # 获取计算参数
    basis_set, functional = get_qc_params_for_method_level(method_level)

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

    # 2. 创建每种配体的 QC 任务（去重，使用 PDB 结构）
    ligand_qc_jobs = {}

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

            ligand_qc_job = create_qc_job_for_structure(
                db=db,
                user_id=user_id,
                md_job_id=md_job_id,
                molecule_name=f"{ligand_type}",
                xyz_content=ligand_xyz,
                charge=ligand_charge,
                basis_set=basis_set,
                functional=functional,
                job_type="ligand"
            )
            ligand_qc_jobs[ligand_key] = ligand_qc_job.id
            created_qc_jobs.append(ligand_qc_job.id)
            logger.info(f"Created ligand QC job {ligand_qc_job.id} for {ligand_type}")

    # 3. 创建每个 cluster_minus_i 的 QC 任务
    cluster_minus_job_ids = []
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
        cluster_minus_job_ids.append(cluster_minus_qc_job.id)
        created_qc_jobs.append(cluster_minus_qc_job.id)
        logger.info(f"Created cluster_minus QC job {cluster_minus_qc_job.id} for ligand {ligand['ligand_label']}")

    # 保存 QC job IDs 到 config
    job.config = job.config or {}
    job.config['qc_job_ids'] = created_qc_jobs
    job.config['cluster_qc_job_id'] = cluster_qc_job.id
    job.config['ligand_qc_jobs'] = ligand_qc_jobs
    job.config['cluster_minus_job_ids'] = cluster_minus_job_ids  # 新增：按顺序保存
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

    logger.info(f"Phase 1 completed: Created {len(created_qc_jobs)} QC jobs (1 cluster, {len(ligand_qc_jobs)} ligand types, {len(cluster_minus_job_ids)} cluster_minus)")

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
    cluster_minus_job_ids = config.get('cluster_minus_job_ids', [])  # 新增：按配体顺序的 job IDs
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

    # 计算每个配体的去溶剂化能
    per_ligand_results = []

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

        # 获取 E_cluster_minus（使用明确的索引）
        if i < len(cluster_minus_job_ids):
            cluster_minus_qc_job_id = cluster_minus_job_ids[i]
        else:
            # 兼容旧格式
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
    center_charge = get_ion_charge(center_ion)

    # 根据 composition 和 molecule_structures 识别配体
    composition = solvation_structure.composition or {}
    ligands = identify_ligands(atoms[1:], composition, molecule_structures)  # 跳过中心离子

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


def identify_ligands(
    atoms: List[Tuple],
    composition: Dict[str, int],
    molecule_structures: Optional[List[Dict]] = None
) -> List[Dict[str, Any]]:
    """
    识别配体分子

    利用平台已有资源：
    1. 从 molecule_structures 获取每种分子的原子数
    2. 根据 composition 和原子数确定每个配体的原子范围
    3. XYZ 中原子是按分子分组排列的

    Args:
        atoms: 配体原子列表 [(element, x, y, z), ...]，不包含中心离子
        composition: 分子组成 {"EC": 3, "DMC": 1, "FSI": 1}
        molecule_structures: 分子结构信息 [{name, type, pdb_content, atoms, ...}]

    Returns:
        配体列表，每个配体包含 ligand_id, ligand_type, atom_indices, xyz_content 等
    """
    # 构建分子类型 -> 原子数的映射
    mol_atom_counts = {}

    if molecule_structures:
        for mol in molecule_structures:
            mol_name = mol.get('name', '')
            mol_atoms = mol.get('atoms', [])
            if mol_name and mol_atoms:
                mol_atom_counts[mol_name] = len(mol_atoms)
                logger.debug(f"Molecule {mol_name}: {len(mol_atoms)} atoms from molecule_structures")

    # 如果没有从 molecule_structures 获取到，使用预定义的原子数
    for mol_type in composition.keys():
        if mol_type not in mol_atom_counts:
            mol_atom_counts[mol_type] = MOLECULE_ATOM_COUNTS.get(mol_type, 15)
            logger.debug(f"Molecule {mol_type}: {mol_atom_counts[mol_type]} atoms from MOLECULE_ATOM_COUNTS")

    ligands = []
    ligand_id = 0
    atom_idx = 0
    type_counters = defaultdict(int)

    for mol_type, count in composition.items():
        if count == 0:
            continue

        atoms_per_mol = mol_atom_counts.get(mol_type, 15)

        for i in range(count):
            ligand_id += 1
            type_counters[mol_type] += 1

            # 确定该配体的原子范围
            start_idx = atom_idx
            end_idx = min(atom_idx + atoms_per_mol, len(atoms))
            ligand_atoms = atoms[start_idx:end_idx]

            if len(ligand_atoms) == 0:
                logger.warning(f"No atoms left for ligand {mol_type}_{type_counters[mol_type]}")
                continue

            # 生成 XYZ 内容
            xyz_lines = [str(len(ligand_atoms)), f"{mol_type}_{type_counters[mol_type]}"]
            for element, x, y, z in ligand_atoms:
                xyz_lines.append(f"{element} {x:.6f} {y:.6f} {z:.6f}")
            xyz_content = '\n'.join(xyz_lines)

            ligands.append({
                'ligand_id': ligand_id,
                'ligand_type': mol_type,
                'ligand_label': f"{mol_type}_{type_counters[mol_type]}",
                'atom_indices': list(range(start_idx, end_idx)),
                'xyz_content': xyz_content,
                'charge': get_molecule_charge(mol_type)
            })

            atom_idx = end_idx

    logger.info(f"Identified {len(ligands)} ligands: {dict(type_counters)}")
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

