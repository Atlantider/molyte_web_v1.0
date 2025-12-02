"""
基于 MDAnalysis 的 RDF 计算器

使用 MDAnalysis 的 InterRDF 功能计算径向分布函数
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis.rdf import InterRDF
    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False

logger = logging.getLogger(__name__)


class MDAnalysisRDFCalculator:
    """使用 MDAnalysis 计算 RDF"""
    
    def __init__(self, work_dir: Path):
        self.work_dir = Path(work_dir)
        self.atom_mapping = None
        self.universe = None
        
        if not HAS_MDANALYSIS:
            raise ImportError("MDAnalysis is required for RDF calculation")
        
        # 加载原子映射
        self._load_atom_mapping()
        
        # 加载轨迹
        self._load_trajectory()
    
    def _load_atom_mapping(self):
        """加载原子映射文件"""
        mapping_file = self.work_dir / "atom_mapping.json"
        
        if not mapping_file.exists():
            raise FileNotFoundError(f"Atom mapping file not found: {mapping_file}")
        
        with open(mapping_file, 'r') as f:
            self.atom_mapping = json.load(f)
        
        logger.info(f"Loaded atom mapping: {len(self.atom_mapping['molecules'])} molecules")
    
    def _load_trajectory(self):
        """加载 LAMMPS 轨迹文件"""
        # 查找轨迹文件
        trj_files = list(self.work_dir.glob("NVT_*.lammpstrj"))
        if not trj_files:
            trj_files = list(self.work_dir.glob("*_after_nvt.lammpstrj"))
        if not trj_files:
            trj_files = list(self.work_dir.glob("*.lammpstrj"))

        if not trj_files:
            raise FileNotFoundError("No LAMMPS trajectory file found")

        trj_file = trj_files[0]

        logger.info(f"Loading trajectory: {trj_file}")

        # 直接从轨迹文件加载，不需要拓扑文件
        # MDAnalysis 会从第一帧推断拓扑
        self.universe = mda.Universe(str(trj_file), format='LAMMPSDUMP')

        logger.info(f"Loaded {len(self.universe.atoms)} atoms, {len(self.universe.trajectory)} frames")
    
    def calculate_rdf(
        self,
        center_label: str,
        target_label: str,
        r_max: float = 10.0,
        n_bins: int = 200,
        exclusion_block: Optional[Tuple[int, int]] = None
    ) -> Dict:
        """
        计算 RDF
        
        Args:
            center_label: 中心原子标签（如 "Li_Li"）
            target_label: 目标原子标签（如 "EC_O01"）
            r_max: 最大距离（Å）
            n_bins: 分箱数量
            exclusion_block: 排除块（用于排除同一分子内的原子）
        
        Returns:
            {
                'center_label': str,
                'target_label': str,
                'r': np.ndarray,
                'g_r': np.ndarray,
                'center_atom_count': int,
                'target_atom_count': int
            }
        """
        # 获取原子 ID
        center_ids = self._get_atom_ids_by_label(center_label)
        target_ids = self._get_atom_ids_by_label(target_label)
        
        if not center_ids or not target_ids:
            raise ValueError(f"No atoms found for labels: {center_label}, {target_label}")
        
        logger.info(f"Calculating RDF: {center_label} ({len(center_ids)} atoms) -> {target_label} ({len(target_ids)} atoms)")

        # 创建原子选择
        # MDAnalysis 的 id 从 1 开始，与 LAMMPS 一致
        center_selection = self.universe.select_atoms(f"index {' '.join(map(lambda x: str(x-1), center_ids))}")
        target_selection = self.universe.select_atoms(f"index {' '.join(map(lambda x: str(x-1), target_ids))}")
        
        # 计算 RDF
        rdf = InterRDF(
            center_selection,
            target_selection,
            nbins=n_bins,
            range=(0.0, r_max),
            exclusion_block=exclusion_block
        )
        
        rdf.run()

        # 计算平均盒子体积
        box_volumes = []
        for ts in self.universe.trajectory:
            box = ts.dimensions  # [lx, ly, lz, alpha, beta, gamma]
            if box is not None:
                volume = box[0] * box[1] * box[2]  # 假设是正交盒子
                box_volumes.append(volume)

        avg_box_volume = np.mean(box_volumes) if box_volumes else None

        return {
            'center_label': center_label,
            'target_label': target_label,
            'r': rdf.results.bins,
            'g_r': rdf.results.rdf,
            'center_atom_count': len(center_ids),
            'target_atom_count': len(target_ids),
            'box_volume': avg_box_volume
        }

    def _get_atom_ids_by_label(self, label: str) -> List[int]:
        """
        根据标签获取原子 ID

        支持以下标签格式:
        - "Li_Li": 精确匹配
        - "EC_O01": 精确匹配原子名称
        - "FSI_F*": 通配符匹配所有以F开头的原子
        - "EC_O(carbonyl)": 化学环境匹配（匹配所有环境为carbonyl的O原子）

        Args:
            label: 原子标签

        Returns:
            原子 ID 列表
        """
        atom_ids = []

        # 解析标签
        parts = label.split('_', 1)  # 只分割一次，保留后面的完整部分
        if len(parts) < 2:
            logger.warning(f"Invalid label format: {label}")
            return atom_ids

        molecule_name = parts[0]
        atom_label = parts[1]

        # 检查是否是化学环境匹配模式: 元素(环境)
        import re
        env_match = re.match(r'^([A-Z][a-z]?)\((\w+)\)$', atom_label)

        # 检查是否是通配符模式: 元素*
        wildcard_match = re.match(r'^([A-Z][a-z]?)\*$', atom_label)

        # 遍历所有分子
        for mol in self.atom_mapping.get('molecules', []):
            mol_name = mol.get('molecule_name') or mol.get('name')
            if mol_name != molecule_name:
                continue

            # 遍历分子中的原子
            for atom in mol.get('atoms', []):
                atom_id = atom.get('atom_id') or atom.get('lammps_id')
                if not atom_id:
                    continue

                atom_name = atom.get('atom_name') or atom.get('label') or atom.get('name', '')
                atom_element = atom.get('element', '')
                atom_env = atom.get('environment', '')
                atom_full_label = atom.get('label', '')

                matched = False

                if env_match:
                    # 化学环境匹配: 元素和环境都必须匹配
                    target_element = env_match.group(1)
                    target_env = env_match.group(2)
                    if atom_element == target_element and atom_env == target_env:
                        matched = True
                elif wildcard_match:
                    # 通配符匹配: 只需元素匹配
                    target_element = wildcard_match.group(1)
                    if atom_element == target_element:
                        matched = True
                elif atom_label == '*':
                    # 全匹配
                    matched = True
                else:
                    # 精确匹配：检查原子名称或完整标签
                    if atom_name == atom_label or atom_full_label == atom_label:
                        matched = True

                if matched:
                    atom_ids.append(atom_id)

        return atom_ids

    def get_available_labels(self) -> List[str]:
        """
        获取所有可用的原子标签

        Returns:
            标签列表，格式如 ["Li_Li", "EC_O01", "FSI_F01", ...]
        """
        labels = []

        for mol in self.atom_mapping['molecules']:
            mol_name = mol['name']
            for atom in mol['atoms']:
                label = f"{mol_name}_{atom['label']}"
                if label not in labels:
                    labels.append(label)

        return labels

