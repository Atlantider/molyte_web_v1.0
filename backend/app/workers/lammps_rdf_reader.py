"""
LAMMPS RDF 文件读取器

从 LAMMPS 输出的 out_rdf.dat 文件中读取 RDF 数据
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
import json

logger = logging.getLogger(__name__)


class LAMMPSRDFReader:
    """读取 LAMMPS 输出的 RDF 文件"""

    def __init__(self, work_dir: Path):
        self.work_dir = work_dir
        self.atom_type_to_label = {}  # 原子类型 -> 标签的映射
        self.atom_id_to_label = {}  # 原子 ID -> 标签的映射
        self.atom_mapping = {}  # 完整的原子映射数据
        self.molecule_templates = {}  # 分子模板信息 {molecule_name: {atoms: [...]}}
        
    def read_rdf_file(self) -> Optional[List[Dict]]:
        """
        读取 LAMMPS 输出的 RDF 文件

        Returns:
            List of RDF results, each containing:
            {
                'center_label': str,
                'target_label': str,
                'r': List[float],
                'g_r': List[float],
                'coordination_number': List[float]  # LAMMPS 计算的累积配位数
            }
        """
        rdf_file = self.work_dir / 'out_rdf.dat'
        in_list_file = self.work_dir / f"{self.work_dir.name}.in.list"

        if not rdf_file.exists():
            logger.warning(f"RDF file not found: {rdf_file}")
            return None

        if not in_list_file.exists():
            logger.warning(f"in.list file not found: {in_list_file}")
            return None

        try:
            # 加载原子映射
            self._load_atom_mapping()

            # 解析 .in.list 文件获取 RDF 对的标签
            rdf_pairs = self._parse_in_list(in_list_file)

            if not rdf_pairs:
                logger.warning("No RDF pairs found in .in.list file")
                return None

            # 读取 RDF 数据
            rdf_data = self._read_rdf_data(rdf_file)

            # 组合标签和数据
            results = []
            for i, (center_label, target_label) in enumerate(rdf_pairs):
                # 每个 RDF 对有两列：g(r) 和 CN
                g_r_col = 1 + i * 2  # 第一列是 r (索引0)，从第二列开始是数据 (索引1)
                cn_col = 1 + i * 2 + 1

                if g_r_col < len(rdf_data):
                    results.append({
                        'center_label': center_label,
                        'target_label': target_label,
                        'r': rdf_data[0],  # 第一列是距离
                        'g_r': rdf_data[g_r_col],
                        'coordination_number': rdf_data[cn_col] if cn_col < len(rdf_data) else None
                    })

            logger.info(f"Read {len(results)} RDF pairs from LAMMPS output")
            return results

        except Exception as e:
            logger.error(f"Failed to read RDF file: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def _parse_in_list(self, in_list_file: Path) -> List[Tuple[str, str]]:
        """
        解析 .in.list 文件，获取 RDF 对的标签

        使用 atom_mapping.json 中的标签信息，智能识别每个原子的来源分子和环境

        Returns:
            List of (center_label, target_label) tuples
        """
        element_list = []
        rdf_pair_indices = []

        # 读取 .in.list 文件
        with open(in_list_file, 'r') as f:
            for line in f:
                if 'element_list' in line:
                    # 提取引号内的内容
                    start = line.find('"') + 1
                    end = line.rfind('"')
                    element_str = line[start:end]
                    element_list = element_str.split()
                    logger.info(f"Element list: {element_list}")

                elif 'rdf_pair' in line:
                    # 提取引号内的内容
                    start = line.find('"') + 1
                    end = line.rfind('"')
                    pair_str = line[start:end]
                    indices = [int(x) for x in pair_str.split()]
                    # 将索引两两配对
                    rdf_pair_indices = [(indices[i], indices[i+1]) for i in range(0, len(indices), 2)]
                    logger.info(f"RDF pair indices: {rdf_pair_indices}")

        if not element_list or not rdf_pair_indices:
            logger.error("Failed to parse element_list or rdf_pair from .in.list file")
            return []

        # 生成标签
        rdf_pairs = []
        for center_idx, target_idx in rdf_pair_indices:
            # element_list 中的索引对应的是模板分子中的原子索引
            # 例如：element_list = ["Li", "S", "F", "O", "N", "C", "O", "O", "O", ...]
            # 索引 1 = Li (第1个原子类型)
            # 索引 2 = S (第2个原子类型，属于第一个 FSI 分子的第一个原子)
            # 索引 7 = O (第7个原子类型，属于第一个 EC 分子的第一个 O 原子)

            center_label = self._get_atom_label_by_template_index(center_idx, element_list)
            target_label = self._get_atom_label_by_template_index(target_idx, element_list)

            if center_label and target_label:
                rdf_pairs.append((center_label, target_label))
                logger.info(f"RDF pair: {center_label} -> {target_label}")
            else:
                logger.warning(f"Failed to generate labels for indices {center_idx}, {target_idx}")

        return rdf_pairs

    def _load_atom_mapping(self):
        """
        加载 atom_mapping.json 文件，并构建原子 ID 到标签的映射
        """
        mapping_file = self.work_dir / 'atom_mapping.json'
        if not mapping_file.exists():
            logger.warning(f"atom_mapping.json not found: {mapping_file}")
            return

        try:
            with open(mapping_file, 'r') as f:
                self.atom_mapping = json.load(f)

            # 构建原子 ID -> 标签的映射
            for molecule in self.atom_mapping.get('molecules', []):
                for atom in molecule.get('atoms', []):
                    atom_id = atom['atom_id']
                    label = atom['label']
                    self.atom_id_to_label[atom_id] = label

            logger.info(f"Loaded atom mapping: {len(self.atom_id_to_label)} atoms")

        except Exception as e:
            logger.error(f"Failed to load atom_mapping.json: {e}")
            import traceback
            traceback.print_exc()

    def _get_atom_label_by_template_index(self, template_idx: int, element_list: List[str]) -> Optional[str]:
        """
        根据 element_list 索引获取原子标签

        element_list 中的每个元素对应一个 LAMMPS atom_type。
        例如：["Li", "S", "F", "O", "N", "C", "O", "O", "O", ...]
        索引 1 对应 atom_type 1 (Li)，索引 2 对应 atom_type 2 (S)，等等。

        Args:
            template_idx: element_list 索引（从 1 开始）
            element_list: 元素列表

        Returns:
            原子标签，例如 "Li_Li", "FSI_S", "EC_type2"
        """
        if not self.atom_mapping:
            logger.warning("Atom mapping not loaded")
            return None

        if template_idx < 1 or template_idx > len(element_list):
            logger.warning(f"Invalid template index: {template_idx}")
            return None

        element = element_list[template_idx - 1]
        atom_type = template_idx  # element_list 索引 = atom_type

        # 从 LAMMPS data 文件中找到第一个该 atom_type 的原子
        data_file = self.work_dir / f"{self.work_dir.name}.data"
        if not data_file.exists():
            logger.warning(f"LAMMPS data file not found: {data_file}")
            return f"{element}_{element}"

        try:
            with open(data_file, 'r') as f:
                content = f.read()

            # 查找 Atoms 部分
            if 'Atoms' not in content:
                logger.warning("No 'Atoms' section found in data file")
                return f"{element}_{element}"

            after_atoms = content.split('Atoms')[1]
            # 按 Bonds 或 Velocities 分割
            if 'Bonds' in after_atoms:
                atoms_section = after_atoms.split('Bonds')[0]
            elif 'Velocities' in after_atoms:
                atoms_section = after_atoms.split('Velocities')[0]
            else:
                atoms_section = after_atoms

            lines = [line.strip() for line in atoms_section.strip().split('\n') if line.strip()]

            # 找到第一个该 atom_type 的原子
            for line in lines:
                parts = line.split()
                if len(parts) >= 4:
                    atom_id = int(parts[0])
                    line_atom_type = int(parts[2])

                    # 检查这个 atom_type 是否匹配
                    if line_atom_type == atom_type:
                        # 从 atom_mapping 中查找这个 atom_id 的标签
                        if atom_id in self.atom_id_to_label:
                            label = self.atom_id_to_label[atom_id]
                            logger.debug(f"element_list[{template_idx}] = {element} -> atom_type {atom_type} -> atom_id {atom_id} -> {label}")
                            return label
                        else:
                            logger.warning(f"atom_id {atom_id} not found in atom_mapping")
                            return f"{element}_{element}"

            logger.warning(f"No atom found for atom_type {atom_type} (element_list[{template_idx}] = {element})")
            return f"{element}_{element}"

        except Exception as e:
            logger.error(f"Failed to parse LAMMPS data file: {e}")
            import traceback
            traceback.print_exc()
            return f"{element}_{element}"
    
    def _read_rdf_data(self, rdf_file: Path) -> List[List[float]]:
        """
        读取 RDF 数据文件
        
        Returns:
            List of columns, each column is a list of float values
        """
        data = []
        
        with open(rdf_file, 'r') as f:
            # 跳过前 4 行注释
            for _ in range(4):
                f.readline()
            
            # 读取数据
            for line in f:
                values = [float(x) for x in line.split()[1:]]  # 跳过第一列（行号）
                if not data:
                    data = [[] for _ in range(len(values))]
                for i, val in enumerate(values):
                    data[i].append(val)
        
        return data

