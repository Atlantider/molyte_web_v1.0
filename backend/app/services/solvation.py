"""
溶剂化结构分析服务

封装 PostPro 中的溶剂化分析逻辑，生成溶剂化结构数据并写入数据库
"""

import os
import sys
import json
import logging
from typing import List, Dict, Any, Optional, Tuple
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

# 配置日志
logger = logging.getLogger(__name__)

# 添加 molyte_v1/postPro 到 Python 路径
POSTPRO_PATH = Path(__file__).parent.parent.parent.parent / "molyte_v1" / "postPro"
if str(POSTPRO_PATH) not in sys.path:
    sys.path.insert(0, str(POSTPRO_PATH))


def apply_pbc(atom_pos: np.ndarray, reference_pos: np.ndarray, box_length: np.ndarray) -> np.ndarray:
    """
    根据周期性边界条件，将原子坐标调整到相对于参考原子最近的位置。
    """
    delta = atom_pos - reference_pos
    for i in range(3):
        if delta[i] > 0.5 * box_length[i]:
            atom_pos[i] -= box_length[i]
        elif delta[i] < -0.5 * box_length[i]:
            atom_pos[i] += box_length[i]
    return atom_pos


def load_atom_mapping(work_dir: str) -> Tuple[Dict[int, str], Dict[str, str], Dict[int, List[int]]]:
    """
    从 atom_mapping.json 加载原子到分子的映射。

    Returns:
        atom_to_mol: 原子ID -> 分子名称
        mol_id_to_name: 分子ID -> 分子名称
        mol_id_to_atoms: 分子ID -> 原子ID列表
    """
    mapping_file = Path(work_dir) / "atom_mapping.json"
    if not mapping_file.exists():
        logger.warning(f"atom_mapping.json 不存在: {mapping_file}")
        return {}, {}, {}

    with open(mapping_file, 'r') as f:
        data = json.load(f)

    atom_to_mol = {}
    mol_id_to_name = {}
    mol_id_to_atoms = {}

    for mol in data.get('molecules', []):
        mol_id = mol['molecule_id']
        mol_name = mol['molecule_name']
        mol_id_to_name[mol_id] = mol_name
        mol_id_to_atoms[mol_id] = []

        for atom in mol.get('atoms', []):
            atom_id = atom['atom_id']
            atom_to_mol[atom_id] = mol_name
            mol_id_to_atoms[mol_id].append(atom_id)

    return atom_to_mol, mol_id_to_name, mol_id_to_atoms


def parse_lammps_dump_frame(file_path: str, frame_idx: int = -1) -> Tuple[Dict[int, Dict], np.ndarray]:
    """
    解析 LAMMPS 轨迹文件的指定帧。

    优化：使用 grep 快速定位帧，避免读取整个大文件

    Args:
        file_path: 轨迹文件路径
        frame_idx: 帧索引，-1 表示最后一帧

    Returns:
        atoms: {atom_id: {'element': str, 'mol': int, 'x': float, 'y': float, 'z': float}}
        box: [lx, ly, lz]
    """
    atoms = {}
    box = np.zeros(3)

    try:
        # 使用 grep 快速找到所有帧的行号
        import subprocess
        result = subprocess.run(
            ['grep', '-n', 'ITEM: TIMESTEP', file_path],
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode != 0 or not result.stdout:
            logger.error("未找到任何帧")
            return atoms, box

        # 解析行号
        frame_lines = []
        for line in result.stdout.strip().split('\n'):
            line_num = int(line.split(':')[0])
            frame_lines.append(line_num)

        if not frame_lines:
            logger.error("未找到任何帧")
            return atoms, box

        # 选择帧
        if frame_idx < 0:
            frame_idx = len(frame_lines) + frame_idx
        if frame_idx < 0 or frame_idx >= len(frame_lines):
            frame_idx = len(frame_lines) - 1

        start_line = frame_lines[frame_idx]
        end_line = frame_lines[frame_idx + 1] if frame_idx + 1 < len(frame_lines) else None

        # 使用 sed 提取指定行范围
        if end_line:
            sed_cmd = f"{start_line},{end_line - 1}p"
        else:
            sed_cmd = f"{start_line},$p"

        result = subprocess.run(
            ['sed', '-n', sed_cmd, file_path],
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode != 0:
            logger.error(f"sed 命令失败: {result.stderr}")
            return atoms, box

        lines = result.stdout.split('\n')

    except Exception as e:
        logger.error(f"使用 grep/sed 提取帧失败: {e}，回退到 Python 实现")
        # 回退到简单的 Python 实现（只适用于小文件）
        with open(file_path, 'r') as f:
            lines = f.readlines()

    # 解析提取的帧数据
    i = 0
    parsing_atoms = False
    id_idx = element_idx = mol_idx = x_idx = y_idx = z_idx = None

    while i < len(lines):
        line = lines[i].strip()

        if "ITEM: BOX BOUNDS" in line:
            # 读取盒子尺寸
            for dim in range(3):
                i += 1
                if i < len(lines):
                    parts = lines[i].strip().split()
                    if len(parts) >= 2:
                        lo, hi = float(parts[0]), float(parts[1])
                        box[dim] = hi - lo

        elif "ITEM: ATOMS" in line:
            # 解析原子数据头
            headers = line.split()[2:]
            id_idx = headers.index('id') if 'id' in headers else None
            element_idx = headers.index('element') if 'element' in headers else None
            mol_idx = headers.index('mol') if 'mol' in headers else None
            x_idx = headers.index('x') if 'x' in headers else None
            y_idx = headers.index('y') if 'y' in headers else None
            z_idx = headers.index('z') if 'z' in headers else None
            parsing_atoms = True

        elif parsing_atoms and line and not line.startswith("ITEM:"):
            # 解析原子行
            parts = line.split()
            try:
                if len(parts) >= max(filter(None, [id_idx, element_idx, mol_idx, x_idx, y_idx, z_idx])) + 1:
                    atom_id = int(parts[id_idx]) if id_idx is not None else i
                    atoms[atom_id] = {
                        'element': parts[element_idx] if element_idx is not None else 'X',
                        'mol': int(parts[mol_idx]) if mol_idx is not None else 0,
                        'x': float(parts[x_idx]) if x_idx is not None else 0.0,
                        'y': float(parts[y_idx]) if y_idx is not None else 0.0,
                        'z': float(parts[z_idx]) if z_idx is not None else 0.0,
                    }
            except (ValueError, IndexError) as e:
                logger.warning(f"解析原子行失败: {line}, error: {e}")

        i += 1

    return atoms, box


def get_molecule_info_from_electrolyte(electrolyte_data: Dict) -> Dict[str, List[str]]:
    """
    从电解液配置中提取分子信息。
    返回 {'cations': ['Li'], 'anions': ['FSI', 'NO3'], 'solvents': ['EC', 'DMC']}
    """
    molecules = {
        'cations': [],
        'anions': [],
        'solvents': []
    }

    if electrolyte_data.get('cations'):
        for c in electrolyte_data['cations']:
            if isinstance(c, dict) and c.get('name'):
                molecules['cations'].append(c['name'])
            elif isinstance(c, str):
                molecules['cations'].append(c)

    if electrolyte_data.get('anions'):
        for a in electrolyte_data['anions']:
            if isinstance(a, dict) and a.get('name'):
                molecules['anions'].append(a['name'])
            elif isinstance(a, str):
                molecules['anions'].append(a)

    if electrolyte_data.get('solvents'):
        for s in electrolyte_data['solvents']:
            if isinstance(s, dict) and s.get('name'):
                molecules['solvents'].append(s['name'])
            elif isinstance(s, str):
                molecules['solvents'].append(s)

    return molecules


def analyze_solvation_structures(
    work_dir: str,
    electrolyte_data: Dict,
    cutoff: float = 3.0,
    center_ion: str = None
) -> List[Dict[str, Any]]:
    """
    分析溶剂化结构。

    Args:
        work_dir: MD 任务工作目录
        electrolyte_data: 电解液配置数据
        cutoff: 溶剂化壳层截断距离 (Å)
        center_ion: 中心离子名称，如果为 None 则自动检测阳离子

    Returns:
        溶剂化结构列表
    """
    work_path = Path(work_dir)
    job_name = work_path.name

    # 查找轨迹文件
    nvt_traj = work_path / f"NVT_{job_name}.lammpstrj"
    after_nvt_traj = work_path / f"{job_name}_after_nvt.lammpstrj"

    # 优先选择 after_nvt 文件（只有最后一帧，速度快）
    if after_nvt_traj.exists():
        traj_file = after_nvt_traj
        frame_idx = 0  # after_nvt 只有一帧
    elif nvt_traj.exists():
        traj_file = nvt_traj
        frame_idx = -1  # 使用最后一帧
    else:
        logger.error(f"未找到轨迹文件: {nvt_traj} 或 {after_nvt_traj}")
        return []

    logger.info(f"分析溶剂化结构: 轨迹={traj_file}")

    # 获取分子信息
    molecules = get_molecule_info_from_electrolyte(electrolyte_data)

    # 确定中心离子
    if center_ion is None:
        if molecules['cations']:
            center_ion = molecules['cations'][0]
        else:
            logger.error("未指定中心离子且电解液配置中没有阳离子")
            return []

    # 创建溶剂化结构输出目录
    solvation_dir = work_path / 'solvation_structures'
    solvation_dir.mkdir(exist_ok=True)

    # 加载原子到分子的映射
    atom_to_mol, mol_id_to_name, mol_id_to_atoms = load_atom_mapping(work_dir)
    if not atom_to_mol:
        logger.error("无法加载 atom_mapping.json")
        return []

    # 解析轨迹指定帧
    atoms, box = parse_lammps_dump_frame(str(traj_file), frame_idx=frame_idx)
    if not atoms:
        logger.error("无法解析轨迹文件")
        return []

    logger.info(f"解析了 {len(atoms)} 个原子，盒子尺寸: {box}")

    # 找到所有中心离子
    center_atom_ids = []
    for atom_id, atom_data in atoms.items():
        mol_name = mol_id_to_name.get(atom_data['mol'], '')
        if mol_name == center_ion:
            center_atom_ids.append(atom_id)

    if not center_atom_ids:
        logger.warning(f"未找到中心离子 {center_ion}")
        return []

    logger.info(f"找到 {len(center_atom_ids)} 个 {center_ion} 原子")

    results = []
    solvation_counts = defaultdict(int)

    # 为每个中心离子分析溶剂化结构
    for i, center_id in enumerate(center_atom_ids):
        center_atom = atoms[center_id]
        center_mol_id = center_atom['mol']
        center_pos = np.array([center_atom['x'], center_atom['y'], center_atom['z']])

        # 查找截断距离内的邻近原子（排除中心离子自身分子）
        neighbor_mols = set()  # 记录邻近分子的 mol_id

        for atom_id, atom_data in atoms.items():
            # 排除中心离子所在分子的所有原子
            if atom_data['mol'] == center_mol_id:
                continue

            # 计算距离（考虑周期性边界条件）
            atom_pos = np.array([atom_data['x'], atom_data['y'], atom_data['z']])
            delta = atom_pos - center_pos

            # 应用最小镜像约定
            for dim in range(3):
                if delta[dim] > box[dim] / 2:
                    delta[dim] -= box[dim]
                elif delta[dim] < -box[dim] / 2:
                    delta[dim] += box[dim]

            dist = np.linalg.norm(delta)

            if dist <= cutoff:
                # 只记录邻近分子的 mol_id，稍后获取完整分子
                neighbor_mols.add(atom_data['mol'])

        # 统计各分子类型的数量（按分子计数，不是按原子）
        mol_counts = Counter()
        for mol_id in neighbor_mols:
            mol_name = mol_id_to_name.get(mol_id, 'UNK')
            if mol_name != center_ion:  # 排除其他中心离子
                mol_counts[mol_name] += 1

        # 构建组成信息
        composition = {}
        for mol_name in molecules['anions'] + molecules['solvents']:
            composition[mol_name] = mol_counts.get(mol_name, 0)

        # 计算配位数（邻近分子数量）
        coord_num = sum(mol_counts.values())

        # 生成溶剂化结构名称
        structure_parts = []
        for mol_name, count in sorted(composition.items()):
            if count > 0:
                structure_parts.append(f"{count}{mol_name}")
        structure_name = "-".join(structure_parts) if structure_parts else "empty"

        # 统计该类型结构的出现次数
        solvation_counts[structure_name] += 1

        # 收集完整分子的所有原子，按分子分组
        mol_groups = {}  # mol_id -> {mol_name, atoms}
        for mol_id in neighbor_mols:
            mol_name = mol_id_to_name.get(mol_id, 'UNK')
            if mol_name == center_ion:  # 排除其他中心离子
                continue

            if mol_id not in mol_groups:
                mol_groups[mol_id] = {'mol_name': mol_name, 'atoms': []}

            # 获取该分子的所有原子 ID
            mol_atom_ids = mol_id_to_atoms.get(mol_id, [])
            for atom_id in mol_atom_ids:
                if atom_id in atoms:
                    atom_data = atoms[atom_id]
                    mol_groups[mol_id]['atoms'].append({
                        'id': atom_id,
                        'element': atom_data['element'],
                        'pos': np.array([atom_data['x'], atom_data['y'], atom_data['z']]),
                    })

        # 记录分子写入顺序（用于后续 cluster minus 计算）
        mol_order = []  # [(mol_name, atom_count), ...]
        xyz_content_lines = []
        total_atoms = 1  # 中心离子

        for mol_id, mol_data in mol_groups.items():
            mol_name = mol_data['mol_name']
            mol_atoms = mol_data['atoms']
            mol_order.append({'mol_name': mol_name, 'atom_count': len(mol_atoms)})
            total_atoms += len(mol_atoms)

            # 为每个原子添加坐标（调整到最小镜像位置）
            for atom_info in mol_atoms:
                pos = atom_info['pos'].copy()
                delta = pos - center_pos
                # 应用最小镜像约定
                for dim in range(3):
                    if delta[dim] > box[dim] / 2:
                        delta[dim] -= box[dim]
                    elif delta[dim] < -box[dim] / 2:
                        delta[dim] += box[dim]
                xyz_content_lines.append(f"{atom_info['element']} {delta[0]:.4f} {delta[1]:.4f} {delta[2]:.4f}")

        # 保存结构文件
        xyz_filename = f"{center_ion}_{i+1}_{structure_name}.xyz"
        xyz_path = solvation_dir / xyz_filename

        # 构建带分子顺序信息的注释行
        # 格式: "Li+ solvation shell (CN=5) | mol_order:DMC,12;PF6,7;EC,10;EMC,15;PF6,7"
        mol_order_str = ";".join([f"{m['mol_name']},{m['atom_count']}" for m in mol_order])
        comment_line = f"{center_ion}+ solvation shell (CN={coord_num}) | mol_order:{mol_order_str}"

        try:
            # 写入 XYZ 文件 - 包含完整分子
            with open(xyz_path, 'w') as f:
                f.write(f"{total_atoms}\n")
                f.write(f"{comment_line}\n")
                # 写入中心离子（作为原点参考）
                f.write(f"{center_atom['element']} 0.0000 0.0000 0.0000\n")
                # 写入所有配体原子
                for line in xyz_content_lines:
                    f.write(f"{line}\n")

            # 同时保存 xyz_content 字符串
            xyz_content = f"{total_atoms}\n{comment_line}\n{center_atom['element']} 0.0000 0.0000 0.0000\n"
            xyz_content += "\n".join(xyz_content_lines)

        except Exception as e:
            logger.warning(f"保存 XYZ 文件失败: {e}")
            xyz_path = None
            xyz_content = None

        results.append({
            'center_ion': center_ion,
            'structure_type': 'first_shell',
            'coordination_num': coord_num,
            'composition': composition,
            'mol_order': mol_order,  # 新增：分子写入顺序
            'file_path': str(xyz_path) if xyz_path else None,
            'xyz_content': xyz_content,  # 新增：直接返回 xyz_content
            'snapshot_frame': len(atoms),  # 使用原子数作为帧标识
            'description': f"{center_ion}+ 第一溶剂化壳层结构 (CN={coord_num}): {structure_name}",
            'ion_index': i + 1,
        })

    # 添加统计信息
    logger.info(f"溶剂化结构统计:")
    for struct_name, count in sorted(solvation_counts.items(), key=lambda x: -x[1]):
        logger.info(f"  {struct_name}: {count} 个")

    return results


def generate_solvation_statistics(structures: List[Dict], electrolyte_data: Dict = None) -> Dict[str, Any]:
    """
    生成溶剂化结构统计信息。
    """
    if not structures:
        return {
            'total_count': 0,
            'coordination_distribution': {},
            'composition_distribution': {},
            'structure_types': {},
            'anion_coordination_distribution': {},
            'molecule_counts': {},
        }

    # 配位数分布
    cn_dist = Counter(s['coordination_num'] for s in structures)

    # 组成分布
    composition_counter = defaultdict(int)
    for s in structures:
        if s.get('composition'):
            key = "_".join(f"{k}{v}" for k, v in sorted(s['composition'].items()) if v > 0)
            composition_counter[key] += 1

    # 结构类型分布
    type_dist = Counter(s.get('structure_type', 'unknown') for s in structures)

    # 平均配位数
    avg_cn = np.mean([s['coordination_num'] for s in structures])

    # 统计各分子类型的总数
    molecule_counts = defaultdict(int)
    for s in structures:
        if s.get('composition'):
            for mol_name, count in s['composition'].items():
                molecule_counts[mol_name] += count

    # 阴离子配位数分布
    anion_cn_dist = Counter()
    anion_names = set()
    if electrolyte_data:
        molecules = get_molecule_info_from_electrolyte(electrolyte_data)
        anion_names = set(molecules.get('anions', []))

    # 如果没有从电解液数据获取阴离子，尝试从 composition 中识别
    if not anion_names:
        common_anions = {'FSI', 'TFSI', 'PF6', 'BF4', 'ClO4', 'NO3', 'SO4', 'Cl', 'Br', 'I'}
        anion_names = common_anions

    # 统计每个结构的阴离子数量
    for s in structures:
        if s.get('composition'):
            anion_count = sum(s['composition'].get(anion, 0) for anion in anion_names)
            anion_cn_dist[anion_count] += 1

    # 确保从0开始连续显示阴离子数量分布
    if anion_cn_dist:
        max_anion_cn = max(anion_cn_dist.keys())
        for i in range(max_anion_cn + 1):
            if i not in anion_cn_dist:
                anion_cn_dist[i] = 0

    return {
        'total_count': len(structures),
        'average_coordination_number': round(avg_cn, 2),
        'coordination_distribution': dict(cn_dist),
        'composition_distribution': dict(composition_counter),
        'structure_types': dict(type_dist),
        'anion_coordination_distribution': dict(anion_cn_dist),
        'molecule_counts': dict(molecule_counts),
    }


def get_frame_count(work_dir: str) -> int:
    """
    获取轨迹文件的帧数。

    优化：使用 grep 命令快速计数，避免在 Python 中读取大文件
    """
    work_path = Path(work_dir)
    job_name = work_path.name

    nvt_traj = work_path / f"NVT_{job_name}.lammpstrj"
    after_nvt_traj = work_path / f"{job_name}_after_nvt.lammpstrj"

    if nvt_traj.exists():
        traj_file = nvt_traj
    elif after_nvt_traj.exists():
        traj_file = after_nvt_traj
    else:
        return 0

    try:
        # 使用 grep -c 快速计数（比 Python 读取大文件快得多）
        import subprocess
        result = subprocess.run(
            ['grep', '-c', 'ITEM: TIMESTEP', str(traj_file)],
            capture_output=True,
            text=True,
            timeout=30
        )
        if result.returncode == 0:
            return int(result.stdout.strip())
        else:
            # grep 失败，回退到 Python 实现
            logger.warning(f"grep failed, falling back to Python implementation")
    except Exception as e:
        logger.warning(f"grep failed with error: {e}, falling back to Python implementation")

    # 回退方案：Python 实现（较慢）
    frame_count = 0
    with open(traj_file, 'r') as f:
        for line in f:
            if "ITEM: TIMESTEP" in line:
                frame_count += 1

    return frame_count


def get_system_structure(work_dir: str, frame_idx: int = -1) -> Dict[str, Any]:
    """
    获取整个体系的结构（用于3D可视化）。

    Args:
        work_dir: 工作目录
        frame_idx: 帧索引，-1 表示最后一帧

    Returns:
        结构数据，包含原子坐标和盒子信息
    """
    work_path = Path(work_dir)
    job_name = work_path.name

    nvt_traj = work_path / f"NVT_{job_name}.lammpstrj"
    after_nvt_traj = work_path / f"{job_name}_after_nvt.lammpstrj"

    # 优化：如果请求最后一帧，优先使用 after_nvt 文件（只有一帧，速度快）
    if frame_idx == -1 and after_nvt_traj.exists():
        traj_file = after_nvt_traj
        frame_idx = 0  # after_nvt 文件只有一帧
    elif nvt_traj.exists():
        traj_file = nvt_traj
    elif after_nvt_traj.exists():
        traj_file = after_nvt_traj
        frame_idx = 0
    else:
        return {'error': 'Trajectory file not found'}

    # 加载原子映射
    _, mol_id_to_name, _ = load_atom_mapping(work_dir)

    # 解析指定帧
    atoms, box = parse_lammps_dump_frame(str(traj_file), frame_idx)

    if not atoms:
        return {'error': 'Failed to parse trajectory'}

    # 获取帧数（只在需要时计算，避免不必要的开销）
    if frame_idx == -1 or frame_idx == 0:
        # 如果使用 after_nvt，帧数就是 1
        if traj_file == after_nvt_traj:
            frame_count = 1
            actual_frame_idx = 0
        else:
            frame_count = get_frame_count(work_dir)
            actual_frame_idx = frame_count - 1 if frame_idx == -1 else frame_idx
    else:
        frame_count = get_frame_count(work_dir)
        actual_frame_idx = frame_idx

    # 转换为 XYZ 格式字符串
    xyz_lines = [str(len(atoms)), f"Frame {actual_frame_idx}"]
    for atom_id in sorted(atoms.keys()):
        atom = atoms[atom_id]
        mol_name = mol_id_to_name.get(atom['mol'], 'UNK')
        xyz_lines.append(f"{atom['element']} {atom['x']:.4f} {atom['y']:.4f} {atom['z']:.4f}")

    xyz_content = "\n".join(xyz_lines)

    return {
        'frame_index': actual_frame_idx,
        'total_frames': frame_count,
        'atom_count': len(atoms),
        'box': box.tolist(),
        'xyz_content': xyz_content,
    }


def get_structure_xyz_content(file_path: str) -> Optional[str]:
    """
    读取溶剂化结构 XYZ 文件内容。
    """
    try:
        with open(file_path, 'r') as f:
            return f.read()
    except Exception as e:
        logger.error(f"读取 XYZ 文件失败: {e}")
        return None

