"""
RDF 计算器

基于 atom_mapping.json 计算径向分布函数（RDF）
"""

import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging

logger = logging.getLogger(__name__)


class RDFCalculator:
    """
    RDF 计算器
    
    基于原子标签计算 RDF，支持通配符匹配
    """
    
    def __init__(self, work_dir: Path):
        """
        Args:
            work_dir: 工作目录（包含 atom_mapping.json 和轨迹文件）
        """
        self.work_dir = Path(work_dir)
        self.atom_mapping = None
        self.label_to_atom_ids = {}
        
        # 加载原子映射
        self._load_atom_mapping()
    
    def _load_atom_mapping(self):
        """加载原子映射文件"""
        mapping_file = self.work_dir / "atom_mapping.json"
        
        if not mapping_file.exists():
            logger.warning(f"Atom mapping file not found: {mapping_file}")
            return
        
        try:
            with open(mapping_file, 'r') as f:
                self.atom_mapping = json.load(f)
            
            # 构建标签到原子 ID 的映射
            self._build_label_index()
            
            logger.info(f"Loaded atom mapping: {len(self.atom_mapping['molecules'])} molecules")
        
        except Exception as e:
            logger.error(f"Failed to load atom mapping: {e}")
    
    def _build_label_index(self):
        """构建标签到原子 ID 的索引"""
        self.label_to_atom_ids = {}
        
        for molecule in self.atom_mapping['molecules']:
            for atom in molecule['atoms']:
                label = atom['label']
                atom_id = atom['atom_id']
                
                if label not in self.label_to_atom_ids:
                    self.label_to_atom_ids[label] = []
                
                self.label_to_atom_ids[label].append(atom_id)
        
        logger.info(f"Built label index: {len(self.label_to_atom_ids)} unique labels")
    
    def get_atom_ids_by_label(
        self,
        label_pattern: str,
        use_wildcard: bool = True
    ) -> List[int]:
        """
        根据标签获取原子 ID 列表
        
        Args:
            label_pattern: 标签模式（如 "Li_Li", "TEP_O*", "*_O_carbonyl"）
            use_wildcard: 是否使用通配符匹配
        
        Returns:
            原子 ID 列表
        """
        if not self.atom_mapping:
            logger.error("Atom mapping not loaded")
            return []
        
        atom_ids = []
        
        if use_wildcard and '*' in label_pattern:
            # 通配符匹配
            for label, ids in self.label_to_atom_ids.items():
                if self._match_pattern(label, label_pattern):
                    atom_ids.extend(ids)
        else:
            # 精确匹配
            atom_ids = self.label_to_atom_ids.get(label_pattern, [])
        
        logger.info(f"Label pattern '{label_pattern}' matched {len(atom_ids)} atoms")
        return atom_ids
    
    def _match_pattern(self, label: str, pattern: str) -> bool:
        """
        简单的通配符匹配
        
        Args:
            label: 标签（如 "TEP_O01"）
            pattern: 模式（如 "TEP_O*", "*_O01"）
        
        Returns:
            是否匹配
        """
        # 将模式转换为正则表达式
        import re
        regex_pattern = pattern.replace('*', '.*')
        regex_pattern = f'^{regex_pattern}$'
        
        return bool(re.match(regex_pattern, label))
    
    def get_available_labels(self) -> Dict[str, List[str]]:
        """
        获取所有可用的原子标签
        
        Returns:
            {
                "all": ["Li_Li", "TEP_P00", "TEP_O01", ...],
                "by_molecule": {
                    "Li": ["Li_Li"],
                    "TEP": ["TEP_P00", "TEP_O01", ...],
                    ...
                },
                "by_element": {
                    "Li": ["Li_Li"],
                    "P": ["TEP_P00"],
                    "O": ["TEP_O01", "TEP_O02", ...],
                    ...
                }
            }
        """
        if not self.atom_mapping:
            return {"all": [], "by_molecule": {}, "by_element": {}}
        
        all_labels = set()
        by_molecule = {}
        by_element = {}
        
        for molecule in self.atom_mapping['molecules']:
            mol_name = molecule['molecule_name']
            
            for atom in molecule['atoms']:
                label = atom['label']
                element = atom['element']
                
                # 添加到所有标签
                all_labels.add(label)
                
                # 按分子分组
                if mol_name not in by_molecule:
                    by_molecule[mol_name] = set()
                by_molecule[mol_name].add(label)
                
                # 按元素分组
                if element not in by_element:
                    by_element[element] = set()
                by_element[element].add(label)
        
        # 转换为列表并排序
        return {
            "all": sorted(list(all_labels)),
            "by_molecule": {k: sorted(list(v)) for k, v in by_molecule.items()},
            "by_element": {k: sorted(list(v)) for k, v in by_element.items()}
        }

    def calculate_rdf(
        self,
        center_label: str,
        target_label: str,
        trajectory_file: Optional[Path] = None,
        r_max: float = 10.0,
        n_bins: int = 200,
        use_wildcard: bool = True
    ) -> Dict:
        """
        计算 RDF

        Args:
            center_label: 中心原子标签（如 "Li_Li", "Li_*"）
            target_label: 目标原子标签（如 "TEP_O01", "*_O*"）
            trajectory_file: 轨迹文件路径（如果为 None，自动查找）
            r_max: 最大距离（Å）
            n_bins: 分箱数量
            use_wildcard: 是否使用通配符

        Returns:
            {
                "r": [0.05, 0.15, 0.25, ...],  # 距离（Å）
                "g_r": [0.0, 0.1, 0.5, ...],   # RDF 值
                "center_label": "Li_Li",
                "target_label": "TEP_O01",
                "center_atom_count": 50,
                "target_atom_count": 100,
                "frame_count": 100
            }
        """
        # 获取原子 ID
        center_ids = self.get_atom_ids_by_label(center_label, use_wildcard)
        target_ids = self.get_atom_ids_by_label(target_label, use_wildcard)

        if not center_ids:
            raise ValueError(f"No atoms found for center label: {center_label}")

        if not target_ids:
            raise ValueError(f"No atoms found for target label: {target_label}")

        logger.info(f"Calculating RDF: {center_label} ({len(center_ids)} atoms) -> "
                   f"{target_label} ({len(target_ids)} atoms)")

        # 查找轨迹文件
        if trajectory_file is None:
            trajectory_file = self._find_trajectory_file()

        if not trajectory_file or not trajectory_file.exists():
            raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")

        # 计算 RDF
        r, g_r, frame_count, avg_box_volume = self._compute_rdf_from_trajectory(
            trajectory_file,
            center_ids,
            target_ids,
            r_max,
            n_bins
        )

        return {
            "r": r.tolist(),
            "g_r": g_r.tolist(),
            "center_label": center_label,
            "target_label": target_label,
            "center_atom_count": len(center_ids),
            "target_atom_count": len(target_ids),
            "frame_count": frame_count,
            "box_volume": avg_box_volume
        }

    def _find_trajectory_file(self) -> Optional[Path]:
        """查找轨迹文件"""
        # 查找 NVT 轨迹文件
        nvt_files = list(self.work_dir.glob("NVT_*.lammpstrj"))

        if nvt_files:
            # 使用第一个找到的文件
            return nvt_files[0]

        # 查找其他轨迹文件
        trj_files = list(self.work_dir.glob("*.lammpstrj"))

        if trj_files:
            return trj_files[0]

        logger.warning("No trajectory file found")
        return None

    def _compute_rdf_from_trajectory(
        self,
        trajectory_file: Path,
        center_ids: List[int],
        target_ids: List[int],
        r_max: float,
        n_bins: int
    ) -> Tuple[np.ndarray, np.ndarray, int, float]:
        """
        从轨迹文件计算 RDF

        Returns:
            (r, g_r, frame_count, avg_box_volume)
        """
        # 初始化
        dr = r_max / n_bins
        r = np.linspace(dr/2, r_max - dr/2, n_bins)
        histogram = np.zeros(n_bins)
        frame_count = 0
        total_box_volume = 0.0

        # 读取轨迹文件
        try:
            frames = self._read_lammpstrj(trajectory_file, center_ids, target_ids)

            for frame in frames:
                center_positions = frame['center_positions']
                target_positions = frame['target_positions']
                box_size = frame['box_size']

                # 计算盒子体积
                box_volume = box_size[0] * box_size[1] * box_size[2]
                total_box_volume += box_volume

                # 计算距离并累加到直方图
                for center_pos in center_positions:
                    distances = self._calculate_distances_pbc(
                        center_pos,
                        target_positions,
                        box_size
                    )

                    # 添加到直方图
                    hist, _ = np.histogram(distances, bins=n_bins, range=(0, r_max))
                    histogram += hist

                frame_count += 1

            if frame_count == 0:
                raise ValueError("No frames found in trajectory")

            # 计算平均盒子体积
            avg_box_volume = total_box_volume / frame_count

            # 归一化
            g_r = self._normalize_rdf(histogram, r, dr, len(center_ids), len(target_ids), frame_count, avg_box_volume)

            logger.info(f"Computed RDF from {frame_count} frames, avg box volume: {avg_box_volume:.2f} Å³")

            return r, g_r, frame_count, avg_box_volume

        except Exception as e:
            logger.error(f"Failed to compute RDF: {e}")
            raise

    def _read_lammpstrj(
        self,
        trajectory_file: Path,
        center_ids: List[int],
        target_ids: List[int]
    ):
        """
        读取 LAMMPS 轨迹文件

        Yields:
            {
                'center_positions': np.ndarray,  # shape: (n_center, 3)
                'target_positions': np.ndarray,  # shape: (n_target, 3)
                'box_size': np.ndarray           # shape: (3,)
            }
        """
        center_id_set = set(center_ids)
        target_id_set = set(target_ids)

        with open(trajectory_file, 'r') as f:
            while True:
                # 读取一帧
                try:
                    # ITEM: TIMESTEP
                    line = f.readline()
                    if not line:
                        break

                    if 'ITEM: TIMESTEP' not in line:
                        continue

                    timestep = int(f.readline().strip())

                    # ITEM: NUMBER OF ATOMS
                    f.readline()
                    n_atoms = int(f.readline().strip())

                    # ITEM: BOX BOUNDS
                    f.readline()
                    xlo_xhi = f.readline().split()
                    ylo_yhi = f.readline().split()
                    zlo_zhi = f.readline().split()

                    box_size = np.array([
                        float(xlo_xhi[1]) - float(xlo_xhi[0]),
                        float(ylo_yhi[1]) - float(ylo_yhi[0]),
                        float(zlo_zhi[1]) - float(zlo_zhi[0])
                    ])

                    # ITEM: ATOMS
                    f.readline()

                    center_positions = []
                    target_positions = []

                    for _ in range(n_atoms):
                        line = f.readline().split()
                        atom_id = int(line[0])
                        x, y, z = float(line[5]), float(line[6]), float(line[7])

                        if atom_id in center_id_set:
                            center_positions.append([x, y, z])

                        if atom_id in target_id_set:
                            target_positions.append([x, y, z])

                    if center_positions and target_positions:
                        yield {
                            'center_positions': np.array(center_positions),
                            'target_positions': np.array(target_positions),
                            'box_size': box_size
                        }

                except (ValueError, IndexError) as e:
                    logger.warning(f"Error reading frame: {e}")
                    break

    def _calculate_distances_pbc(
        self,
        center_pos: np.ndarray,
        target_positions: np.ndarray,
        box_size: np.ndarray
    ) -> np.ndarray:
        """
        计算距离（考虑周期性边界条件）

        Args:
            center_pos: 中心原子位置 (3,)
            target_positions: 目标原子位置 (n, 3)
            box_size: 盒子大小 (3,)

        Returns:
            距离数组 (n,)
        """
        # 计算距离向量
        delta = target_positions - center_pos

        # 应用最小镜像约定
        delta = delta - box_size * np.round(delta / box_size)

        # 计算距离
        distances = np.linalg.norm(delta, axis=1)

        return distances

    def _normalize_rdf(
        self,
        histogram: np.ndarray,
        r: np.ndarray,
        dr: float,
        n_center: int,
        n_target: int,
        n_frames: int,
        box_volume: float
    ) -> np.ndarray:
        """
        归一化 RDF

        Args:
            histogram: 距离直方图
            r: 距离数组
            dr: 分箱宽度
            n_center: 中心原子数量
            n_target: 目标原子数量
            n_frames: 帧数
            box_volume: 盒子体积（Å³）

        Returns:
            归一化的 g(r)
        """
        # 计算每个壳层的体积
        shell_volume = 4.0 / 3.0 * np.pi * ((r + dr/2)**3 - (r - dr/2)**3)

        # 计算目标原子的数密度（atoms/Å³）
        rho = n_target / box_volume

        # 归一化 RDF
        # g(r) = (histogram / n_frames / n_center) / (rho * shell_volume)
        # 其中 histogram 是所有中心原子在所有帧中的距离直方图
        g_r = histogram / (n_frames * n_center * rho * shell_volume)

        return g_r

