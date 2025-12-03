"""
去溶剂化能计算任务

计算公式：ΔE_i = E_cluster - (E_cluster_minus_i + E_i)
"""
import logging
import tempfile
import subprocess
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime
from collections import defaultdict

import numpy as np

from app.database import SessionLocal
from app.models.job import PostprocessJob, JobStatus
from app.models.result import SolvationStructure, DesolvationEnergyResult

logger = logging.getLogger(__name__)

# Hartree to kcal/mol conversion factor
HARTREE_TO_KCAL = 627.509474


def run_desolvation_job(job: PostprocessJob, db: SessionLocal) -> Dict[str, Any]:
    """
    执行去溶剂化能计算任务
    
    Args:
        job: PostprocessJob 对象
        db: 数据库会话
    
    Returns:
        结果字典 {"success": bool, "job_id": int, ...}
    """
    try:
        logger.info(f"Starting desolvation job {job.id}")
        
        # 1. 更新状态为 RUNNING
        job.status = JobStatus.RUNNING
        job.started_at = datetime.now()
        db.commit()
        
        # 2. 获取配置
        config = job.config or {}
        solvation_structure_id = config.get("solvation_structure_id")
        method_level = config.get("method_level", "fast_xtb")
        
        if not solvation_structure_id:
            raise ValueError("Missing solvation_structure_id in job config")
        
        # 3. 加载溶剂化结构
        solvation_structure = db.query(SolvationStructure).filter(
            SolvationStructure.id == solvation_structure_id
        ).first()
        
        if not solvation_structure:
            raise ValueError(f"Solvation structure {solvation_structure_id} not found")
        
        if not solvation_structure.xyz_content:
            raise ValueError(f"Solvation structure {solvation_structure_id} has no XYZ content")
        
        logger.info(f"Loaded solvation structure {solvation_structure_id}: {solvation_structure.center_ion}, CN={solvation_structure.coordination_num}")
        
        # 4. 解析溶剂化结构
        cluster_data = parse_solvation_cluster(solvation_structure)
        logger.info(f"Parsed cluster: center={cluster_data['center_ion']}, ligands={len(cluster_data['ligands'])}")
        
        # 5. 计算 E_cluster
        logger.info("Calculating cluster energy...")
        e_cluster = calculate_energy_with_xtb(
            cluster_data['xyz_content'],
            charge=cluster_data['total_charge']
        )
        logger.info(f"E_cluster = {e_cluster:.6f} A.U.")
        
        # 6. 计算每个配体的去溶剂化能
        per_ligand_results = []
        for i, ligand in enumerate(cluster_data['ligands']):
            logger.info(f"Processing ligand {i+1}/{len(cluster_data['ligands'])}: {ligand['ligand_label']}")
            
            result = calculate_ligand_desolvation(
                cluster_data, ligand, e_cluster, method_level
            )
            per_ligand_results.append(result)
            
            logger.info(f"  ΔE = {result['delta_e']:.2f} kcal/mol")
        
        # 7. 按类型汇总
        per_type_summary = summarize_by_type(per_ligand_results)
        logger.info(f"Type summary: {len(per_type_summary)} types")
        
        # 8. 保存结果
        desolvation_result = DesolvationEnergyResult(
            postprocess_job_id=job.id,
            solvation_structure_id=solvation_structure_id,
            method_level=method_level,
            basis_set="GFN2-xTB" if method_level == "fast_xtb" else None,
            functional="GFN2-xTB" if method_level == "fast_xtb" else None,
            e_cluster=e_cluster,
            per_ligand_results=per_ligand_results,
            per_type_summary=per_type_summary
        )
        
        db.add(desolvation_result)
        
        # 9. 更新任务状态
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


def calculate_energy_with_xtb(xyz_content: str, charge: int = 0) -> float:
    """
    使用 xTB 计算能量

    Args:
        xyz_content: XYZ 格式的结构
        charge: 电荷

    Returns:
        能量 (A.U.)
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)

        # 写入 XYZ 文件
        xyz_file = tmpdir_path / "structure.xyz"
        xyz_file.write_text(xyz_content)

        # 运行 xTB
        cmd = f"xtb {xyz_file.name} --chrg {charge} --gfn 2"

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                cwd=tmpdir,
                capture_output=True,
                text=True,
                timeout=300  # 5 分钟超时
            )

            if result.returncode != 0:
                raise RuntimeError(f"xTB failed with return code {result.returncode}: {result.stderr}")

            # 解析能量
            for line in result.stdout.split('\n'):
                if 'TOTAL ENERGY' in line:
                    # 格式: "TOTAL ENERGY      -42.123456 Eh"
                    parts = line.split()
                    energy_hartree = float(parts[3])
                    return energy_hartree

            raise ValueError("Failed to extract energy from xTB output")

        except subprocess.TimeoutExpired:
            raise RuntimeError("xTB calculation timed out after 5 minutes")
        except Exception as e:
            logger.error(f"xTB calculation failed: {e}")
            raise


def calculate_ligand_desolvation(
    cluster_data: Dict[str, Any],
    ligand: Dict[str, Any],
    e_cluster: float,
    method_level: str
) -> Dict[str, Any]:
    """
    计算单个配体的去溶剂化能

    Returns:
        {
            ligand_id, ligand_type, ligand_label,
            e_ligand, e_cluster_minus, delta_e
        }
    """
    # 1. 计算单分子能量 E_i
    e_ligand = calculate_energy_with_xtb(
        ligand['xyz_content'],
        charge=ligand['charge']
    )

    # 2. 生成 cluster_minus_i 的 XYZ
    cluster_minus_xyz = generate_cluster_minus_xyz(cluster_data, ligand)
    cluster_minus_charge = cluster_data['total_charge'] - ligand['charge']

    # 3. 计算 E_cluster_minus_i
    e_cluster_minus = calculate_energy_with_xtb(
        cluster_minus_xyz,
        charge=cluster_minus_charge
    )

    # 4. 计算 ΔE_i
    delta_e_au = e_cluster - (e_cluster_minus + e_ligand)
    delta_e_kcal = delta_e_au * HARTREE_TO_KCAL

    return {
        'ligand_id': ligand['ligand_id'],
        'ligand_type': ligand['ligand_type'],
        'ligand_label': ligand['ligand_label'],
        'e_ligand': e_ligand,
        'e_cluster_minus': e_cluster_minus,
        'delta_e': delta_e_kcal
    }


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

