"""
MSD 数据处理任务

增强功能：
- 扩散系数计算
- 离子电导率计算 (Nernst-Einstein)
- 离子迁移率计算
- 迁移数计算
"""
import json
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from collections import Counter
from sqlalchemy.orm import Session

from app.models.result import MSDResult
from app.workers.lammps_msd_reader import (
    LAMMPSMSDReader,
    calculate_diffusion_coefficient,
    calculate_ionic_conductivity,
    calculate_mobility,
    calculate_transference_number,
)


# 默认离子电荷
DEFAULT_ION_CHARGES = {
    'Li': 1, 'Na': 1, 'K': 1, 'Mg': 2, 'Ca': 2, 'Zn': 2, 'Al': 3,
    'FSI': -1, 'TFSI': -1, 'PF6': -1, 'BF4': -1, 'ClO4': -1, 'DCA': -1,
}


def get_ion_charge(species: str) -> int:
    """获取离子电荷"""
    for ion, charge in DEFAULT_ION_CHARGES.items():
        if ion in species:
            return charge
    # 默认返回1（假设是阳离子）
    return 1


def extract_box_volume_and_ion_counts(work_dir: Path) -> Tuple[Optional[float], Optional[Dict[str, int]]]:
    """
    从工作目录中提取盒子体积和离子数量

    Args:
        work_dir: 工作目录路径

    Returns:
        (box_volume, ion_counts) - 盒子体积 (Å³) 和离子数量字典
    """
    box_volume = None
    ion_counts = None

    try:
        # 1. 从 atom_mapping.json 提取离子数量
        atom_mapping_file = work_dir / "atom_mapping.json"
        if atom_mapping_file.exists():
            with open(atom_mapping_file) as f:
                atom_mapping = json.load(f)

            if 'molecules' in atom_mapping:
                molecules = atom_mapping['molecules']
                # 统计各类型分子数量
                mol_names = [mol.get('molecule_name', 'unknown') for mol in molecules]
                counts = Counter(mol_names)

                # 只提取离子的数量
                ion_counts = {}
                for name, count in counts.items():
                    if name in DEFAULT_ION_CHARGES:
                        ion_counts[name] = count

                if ion_counts:
                    print(f"  📊 Extracted ion counts from atom_mapping.json: {ion_counts}")

        # 2. 从 LAMMPS data 文件提取盒子体积
        # 优先使用 NVT 后的 data 文件
        data_files = list(work_dir.glob("*_after_nvt.data"))
        if not data_files:
            data_files = list(work_dir.glob("*_after_npt.data"))
        if not data_files:
            data_files = list(work_dir.glob("*.data"))

        if data_files:
            data_file = data_files[0]
            Lx = Ly = Lz = None

            with open(data_file) as f:
                for line in f:
                    if 'xlo xhi' in line:
                        parts = line.split()
                        xlo, xhi = float(parts[0]), float(parts[1])
                        Lx = xhi - xlo
                    elif 'ylo yhi' in line:
                        parts = line.split()
                        ylo, yhi = float(parts[0]), float(parts[1])
                        Ly = yhi - ylo
                    elif 'zlo zhi' in line:
                        parts = line.split()
                        zlo, zhi = float(parts[0]), float(parts[1])
                        Lz = zhi - zlo
                        break

            if Lx and Ly and Lz:
                box_volume = Lx * Ly * Lz
                print(f"  📦 Extracted box volume from {data_file.name}: {box_volume:.2f} Å³")

    except Exception as e:
        print(f"  ⚠️ Failed to extract box/ion info: {e}")

    return box_volume, ion_counts


def process_msd_data(
    db: Session,
    job_id: int,
    work_dir: Path,
    temperature: float = 298.15,
    box_volume: Optional[float] = None,
    ion_counts: Optional[Dict[str, int]] = None,
) -> List[MSDResult]:
    """
    处理 MSD 数据并保存到数据库

    Args:
        db: 数据库会话
        job_id: MD 任务 ID
        work_dir: 工作目录
        temperature: 温度 (K)，默认 298.15
        box_volume: 模拟盒子体积 (Å³)，用于计算电导率
        ion_counts: 各离子数量，用于计算电导率

    Returns:
        MSD 结果列表
    """
    # 删除旧的 MSD 结果
    db.query(MSDResult).filter(MSDResult.md_job_id == job_id).delete()
    db.commit()

    # 如果没有提供 box_volume 或 ion_counts，尝试自动提取
    if box_volume is None or ion_counts is None:
        extracted_volume, extracted_counts = extract_box_volume_and_ion_counts(work_dir)
        if box_volume is None:
            box_volume = extracted_volume
        if ion_counts is None:
            ion_counts = extracted_counts

    # 读取 MSD 数据
    reader = LAMMPSMSDReader(work_dir)
    msd_data_list = reader.read_all_msd()

    if not msd_data_list:
        print(f"No MSD data found in {work_dir}")
        return []

    # 保存到数据库
    results = []
    for msd_data in msd_data_list:
        species = msd_data['species']

        # 计算扩散系数（使用后半段数据）
        diffusion_coeff = calculate_diffusion_coefficient(
            msd_data['time'],
            msd_data['msd_total']
        )

        # 获取离子电荷
        charge = get_ion_charge(species)

        # 计算离子迁移率
        mobility = calculate_mobility(diffusion_coeff, charge, temperature)

        # 计算离子电导率（如果有盒子体积和离子数量）
        ionic_conductivity = None
        if box_volume and ion_counts and diffusion_coeff:
            ion_count = ion_counts.get(species, 0)
            # 尝试模糊匹配
            if ion_count == 0:
                for key, val in ion_counts.items():
                    if key in species or species in key:
                        ion_count = val
                        break

            if ion_count > 0:
                ionic_conductivity = calculate_ionic_conductivity(
                    diffusion_coeff, ion_count, box_volume, charge, temperature
                )

        # 创建 MSD 结果
        msd_result = MSDResult(
            md_job_id=job_id,
            species=species,
            t_values=msd_data['time'],  # 使用 t_values 列名
            msd_x_values=msd_data['msd_x'],
            msd_y_values=msd_data['msd_y'],
            msd_z_values=msd_data['msd_z'],
            msd_total_values=msd_data['msd_total'],
            labels=msd_data['labels'],
            diffusion_coefficient=diffusion_coeff,
            ionic_conductivity=ionic_conductivity,
            mobility=mobility,
            charge=charge,
        )

        db.add(msd_result)
        results.append(msd_result)

    db.commit()

    print(f"✅ Saved {len(results)} MSD results for job {job_id}")
    for result in results:
        D_str = f"{result.diffusion_coefficient:.2e}" if result.diffusion_coefficient else "N/A"
        sigma_str = f"{result.ionic_conductivity:.2e}" if result.ionic_conductivity else "N/A"
        print(f"  - {result.species}: D = {D_str} cm²/s, σ = {sigma_str} S/cm")

    return results


def get_msd_results(db: Session, job_id: int) -> List[Dict[str, Any]]:
    """
    获取 MSD 结果

    Args:
        db: 数据库会话
        job_id: MD 任务 ID

    Returns:
        MSD 结果列表（字典格式）
    """
    results = db.query(MSDResult).filter(MSDResult.md_job_id == job_id).all()

    return [
        {
            'id': r.id,
            'species': r.species,
            'time': r.t_values,  # 使用 t_values 列名
            'msd_x': r.msd_x_values,
            'msd_y': r.msd_y_values,
            'msd_z': r.msd_z_values,
            'msd_total': r.msd_total_values,
            'labels': r.labels,
            'diffusion_coefficient': r.diffusion_coefficient,
            'ionic_conductivity': r.ionic_conductivity,
            'mobility': r.mobility,
            'charge': r.charge,
            'created_at': r.created_at.isoformat() if r.created_at else None,
        }
        for r in results
    ]


def get_transport_properties_summary(db: Session, job_id: int) -> Dict[str, Any]:
    """
    获取传输性质汇总

    Args:
        db: 数据库会话
        job_id: MD 任务 ID

    Returns:
        传输性质汇总字典
    """
    results = db.query(MSDResult).filter(MSDResult.md_job_id == job_id).all()

    if not results:
        return {}

    summary = {
        'species_data': {},
        'total_conductivity': 0.0,
        'transference_numbers': None,
    }

    cation_D = None
    anion_D = None

    for r in results:
        species_data = {
            'diffusion_coefficient': r.diffusion_coefficient,
            'diffusion_coefficient_unit': 'cm²/s',
            'ionic_conductivity': r.ionic_conductivity,
            'ionic_conductivity_unit': 'S/cm',
            'mobility': r.mobility,
            'mobility_unit': 'cm²/(V·s)',
            'charge': r.charge,
        }
        summary['species_data'][r.species] = species_data

        if r.ionic_conductivity:
            summary['total_conductivity'] += r.ionic_conductivity

        # 记录阳离子和阴离子的扩散系数
        if r.charge and r.charge > 0:
            cation_D = r.diffusion_coefficient
        elif r.charge and r.charge < 0:
            anion_D = r.diffusion_coefficient

    # 计算迁移数
    if cation_D and anion_D:
        t_numbers = calculate_transference_number(cation_D, anion_D)
        if t_numbers:
            summary['transference_numbers'] = {
                't_plus': t_numbers[0],
                't_minus': t_numbers[1],
            }

    summary['total_conductivity_unit'] = 'S/cm'

    return summary

