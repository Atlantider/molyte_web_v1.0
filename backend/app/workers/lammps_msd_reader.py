"""
LAMMPS MSD 数据读取器
读取 out_*_msd.dat 文件并解析 MSD 数据

增强功能：
- 扩散系数计算 (D = MSD / 6t)
- 离子电导率计算 (Nernst-Einstein: σ = nq²D / kT)
- 离子迁移率计算 (μ = qD / kT)
- 迁移数计算 (t+ = D+ / (D+ + D-))
"""
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
import re
import math

# 物理常数
BOLTZMANN_CONSTANT = 1.380649e-23  # J/K
ELEMENTARY_CHARGE = 1.602176634e-19  # C
AVOGADRO_NUMBER = 6.02214076e23  # mol^-1


class LAMMPSMSDReader:
    """LAMMPS MSD 数据读取器"""
    
    def __init__(self, work_dir: Path):
        """
        初始化 MSD 读取器
        
        Args:
            work_dir: 工作目录路径
        """
        self.work_dir = Path(work_dir)
    
    def find_msd_files(self) -> List[Path]:
        """
        查找所有 MSD 文件
        
        Returns:
            MSD 文件路径列表
        """
        msd_files = list(self.work_dir.glob("out_*_msd.dat"))
        return sorted(msd_files)
    
    def read_msd_file(self, msd_file: Path) -> Optional[Dict[str, Any]]:
        """
        读取单个 MSD 文件
        
        Args:
            msd_file: MSD 文件路径
            
        Returns:
            MSD 数据字典，包含：
            - species: 物种名称（如 Li, FSI）
            - time: 时间数组（fs）
            - msd_x: x 方向 MSD 数组
            - msd_y: y 方向 MSD 数组
            - msd_z: z 方向 MSD 数组
            - msd_total: 总 MSD 数组
            - labels: 图例标签字典
        """
        if not msd_file.exists():
            return None
        
        try:
            with open(msd_file, 'r') as f:
                lines = f.readlines()
            
            if len(lines) < 3:
                return None
            
            # 第一行：列标题（跳过）
            # 第二行：图例标签
            legend_line = lines[1].strip().split()
            if len(legend_line) < 5:
                return None
            
            # 提取物种名称（从文件名）
            # 文件名格式: out_Li_msd.dat 或 out_FSI_msd.dat
            match = re.search(r'out_(.+?)_msd\.dat', msd_file.name)
            if not match:
                return None
            species = match.group(1)
            
            # 图例标签
            labels = {
                'time': legend_line[0],  # fs
                'x': legend_line[1],     # species_x
                'y': legend_line[2],     # species_y
                'z': legend_line[3],     # species_z
                'total': legend_line[4], # species_total
            }
            
            # 读取数据（从第三行开始）
            time = []
            msd_x = []
            msd_y = []
            msd_z = []
            msd_total = []
            
            for line in lines[2:]:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split()
                if len(parts) < 5:
                    continue
                
                try:
                    time.append(float(parts[0]))
                    msd_x.append(float(parts[1]))
                    msd_y.append(float(parts[2]))
                    msd_z.append(float(parts[3]))
                    msd_total.append(float(parts[4]))
                except (ValueError, IndexError):
                    continue
            
            if not time:
                return None
            
            return {
                'species': species,
                'time': time,
                'msd_x': msd_x,
                'msd_y': msd_y,
                'msd_z': msd_z,
                'msd_total': msd_total,
                'labels': labels,
            }
        
        except Exception as e:
            print(f"Error reading MSD file {msd_file}: {e}")
            return None
    
    def read_all_msd(self) -> List[Dict[str, Any]]:
        """
        读取所有 MSD 文件
        
        Returns:
            MSD 数据列表
        """
        msd_files = self.find_msd_files()
        results = []
        
        for msd_file in msd_files:
            msd_data = self.read_msd_file(msd_file)
            if msd_data:
                results.append(msd_data)
        
        return results


def calculate_diffusion_coefficient(time: List[float], msd: List[float],
                                    start_idx: int = None, end_idx: int = None) -> Optional[float]:
    """
    计算扩散系数（从 MSD 的线性区域）

    D = MSD / (6 * t)

    Args:
        time: 时间数组（fs）
        msd: MSD 数组（Å²）
        start_idx: 线性拟合起始索引
        end_idx: 线性拟合结束索引

    Returns:
        扩散系数（cm²/s）
    """
    if not time or not msd or len(time) != len(msd):
        return None

    # 默认使用后半段数据进行线性拟合
    if start_idx is None:
        start_idx = len(time) // 2
    if end_idx is None:
        end_idx = len(time)

    if start_idx >= end_idx or start_idx < 0 or end_idx > len(time):
        return None

    # 线性拟合 MSD = slope * t
    time_subset = time[start_idx:end_idx]
    msd_subset = msd[start_idx:end_idx]

    if not time_subset or not msd_subset:
        return None

    # 简单线性拟合：slope = (MSD_end - MSD_start) / (t_end - t_start)
    slope = (msd_subset[-1] - msd_subset[0]) / (time_subset[-1] - time_subset[0])

    # D = slope / 6
    # 单位转换：Å²/fs -> cm²/s
    # 1 Å² = 1e-16 cm²
    # 1 fs = 1e-15 s
    # D (cm²/s) = slope (Å²/fs) * 1e-16 / 1e-15 = slope * 1e-1
    D = slope / 6.0 * 1e-1

    return D


def calculate_ionic_conductivity(
    diffusion_coefficient: float,
    ion_count: int,
    box_volume: float,
    charge: int,
    temperature: float = 298.15
) -> Optional[float]:
    """
    使用 Nernst-Einstein 方程计算离子电导率

    σ = n * q² * D / (k_B * T)

    Args:
        diffusion_coefficient: 扩散系数 (cm²/s)
        ion_count: 离子数量
        box_volume: 模拟盒子体积 (Å³)
        charge: 离子电荷数 (如 Li+ 为 1, FSI- 为 -1)
        temperature: 温度 (K)

    Returns:
        离子电导率 (S/cm)
    """
    if diffusion_coefficient is None or diffusion_coefficient <= 0:
        return None
    if ion_count <= 0 or box_volume <= 0:
        return None

    # 计算数密度 n (ions/cm³)
    # box_volume 单位是 Å³, 1 Å³ = 1e-24 cm³
    n = ion_count / (box_volume * 1e-24)

    # 扩散系数单位转换: cm²/s -> m²/s
    D_m2s = diffusion_coefficient * 1e-4

    # 计算电导率
    # σ = n * q² * D / (k_B * T)
    # 单位: (1/m³) * C² * (m²/s) / (J/K * K) = C²/(J·m·s) = S/m
    q = abs(charge) * ELEMENTARY_CHARGE
    sigma_Sm = n * 1e6 * q * q * D_m2s / (BOLTZMANN_CONSTANT * temperature)

    # 转换为 S/cm
    sigma_Scm = sigma_Sm / 100

    return sigma_Scm


def calculate_molar_conductivity(
    ionic_conductivity: float,
    concentration: float
) -> Optional[float]:
    """
    计算摩尔电导率

    Λ = σ / c

    Args:
        ionic_conductivity: 离子电导率 (S/cm)
        concentration: 浓度 (mol/L)

    Returns:
        摩尔电导率 (S·cm²/mol)
    """
    if ionic_conductivity is None or concentration <= 0:
        return None

    # σ (S/cm) / c (mol/L) = σ (S/cm) / c (mol/1000cm³)
    # = σ * 1000 / c (S·cm²/mol)
    return ionic_conductivity * 1000 / concentration


def calculate_mobility(
    diffusion_coefficient: float,
    charge: int,
    temperature: float = 298.15
) -> Optional[float]:
    """
    计算离子迁移率

    μ = q * D / (k_B * T)

    Args:
        diffusion_coefficient: 扩散系数 (cm²/s)
        charge: 离子电荷数
        temperature: 温度 (K)

    Returns:
        离子迁移率 (cm²/(V·s))
    """
    if diffusion_coefficient is None or diffusion_coefficient <= 0:
        return None

    # D 单位转换: cm²/s -> m²/s
    D_m2s = diffusion_coefficient * 1e-4

    # μ = q * D / (k_B * T)
    # 单位: C * m²/s / (J/K * K) = C·m²/(J·s) = m²/(V·s)
    q = abs(charge) * ELEMENTARY_CHARGE
    mu_m2Vs = q * D_m2s / (BOLTZMANN_CONSTANT * temperature)

    # 转换为 cm²/(V·s)
    mu_cm2Vs = mu_m2Vs * 1e4

    return mu_cm2Vs


def calculate_transference_number(
    D_cation: float,
    D_anion: float
) -> Optional[Tuple[float, float]]:
    """
    计算迁移数

    t+ = D+ / (D+ + D-)
    t- = D- / (D+ + D-)

    Args:
        D_cation: 阳离子扩散系数 (cm²/s)
        D_anion: 阴离子扩散系数 (cm²/s)

    Returns:
        (t+, t-) 阳离子和阴离子迁移数
    """
    if D_cation is None or D_anion is None:
        return None
    if D_cation <= 0 and D_anion <= 0:
        return None

    D_total = abs(D_cation) + abs(D_anion)
    if D_total == 0:
        return None

    t_plus = abs(D_cation) / D_total
    t_minus = abs(D_anion) / D_total

    return (t_plus, t_minus)


def calculate_all_transport_properties(
    msd_data_list: List[Dict[str, Any]],
    temperature: float = 298.15,
    box_volume: float = None,
    ion_counts: Dict[str, int] = None,
    ion_charges: Dict[str, int] = None
) -> Dict[str, Any]:
    """
    计算所有传输性质

    Args:
        msd_data_list: MSD 数据列表
        temperature: 温度 (K)
        box_volume: 模拟盒子体积 (Å³)
        ion_counts: 各离子数量 {'Li': 100, 'FSI': 100}
        ion_charges: 各离子电荷 {'Li': 1, 'FSI': -1}

    Returns:
        包含所有传输性质的字典
    """
    results = {
        'temperature': temperature,
        'species_data': {},
        'total_conductivity': None,
        'transference_numbers': None,
    }

    # 默认离子电荷
    default_charges = {
        'Li': 1, 'Na': 1, 'K': 1, 'Mg': 2, 'Ca': 2, 'Zn': 2, 'Al': 3,
        'FSI': -1, 'TFSI': -1, 'PF6': -1, 'BF4': -1, 'ClO4': -1, 'DCA': -1,
    }

    if ion_charges is None:
        ion_charges = default_charges

    cation_D = None
    anion_D = None
    total_conductivity = 0.0

    for msd_data in msd_data_list:
        species = msd_data.get('species', '')
        time = msd_data.get('time', [])
        msd_total = msd_data.get('msd_total', [])

        # 计算扩散系数
        D = calculate_diffusion_coefficient(time, msd_total)

        species_result = {
            'diffusion_coefficient': D,
            'diffusion_coefficient_unit': 'cm²/s',
        }

        # 获取电荷
        charge = ion_charges.get(species, 1)
        for key, val in ion_charges.items():
            if key in species:
                charge = val
                break

        species_result['charge'] = charge

        # 计算迁移率
        mobility = calculate_mobility(D, charge, temperature)
        species_result['mobility'] = mobility
        species_result['mobility_unit'] = 'cm²/(V·s)'

        # 如果有盒子体积和离子数量，计算电导率
        if box_volume and ion_counts and D:
            ion_count = ion_counts.get(species, 0)
            for key, val in ion_counts.items():
                if key in species:
                    ion_count = val
                    break

            if ion_count > 0:
                conductivity = calculate_ionic_conductivity(
                    D, ion_count, box_volume, charge, temperature
                )
                species_result['ionic_conductivity'] = conductivity
                species_result['ionic_conductivity_unit'] = 'S/cm'

                if conductivity:
                    total_conductivity += conductivity

        # 判断是阳离子还是阴离子
        if charge > 0:
            cation_D = D
        else:
            anion_D = D

        results['species_data'][species] = species_result

    # 计算总电导率
    if total_conductivity > 0:
        results['total_conductivity'] = total_conductivity
        results['total_conductivity_unit'] = 'S/cm'

    # 计算迁移数
    if cation_D and anion_D:
        t_numbers = calculate_transference_number(cation_D, anion_D)
        if t_numbers:
            results['transference_numbers'] = {
                't_plus': t_numbers[0],
                't_minus': t_numbers[1],
            }

    return results

