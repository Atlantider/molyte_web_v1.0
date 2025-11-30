"""
Moltemplate .lt 文件解析器

从 .lt 文件中提取原子信息
"""

import re
from pathlib import Path
from typing import Dict, List, Optional
import logging

logger = logging.getLogger(__name__)


class LTParser:
    """解析 Moltemplate .lt 文件"""

    def parse_lt_file(self, lt_file: Path, molecule_name: Optional[str] = None) -> Dict:
        """
        解析 .lt 文件

        Args:
            lt_file: .lt 文件路径
            molecule_name: 分子名称（可选，如果不提供则从文件名提取）

        Returns:
            {
                "molecule_name": "FSI",
                "atoms": [
                    {
                        "index": 1,
                        "name": "S1",
                        "type": "S",
                        "element": "S",
                        "charge": 1.02,
                        "mass": 32.06
                    },
                    ...
                ],
                "bonds": [(1, 2), (1, 3), ...]
            }
        """
        if not lt_file.exists():
            raise FileNotFoundError(f"LT file not found: {lt_file}")

        # 如果没有提供分子名称，从文件名提取
        if molecule_name is None:
            molecule_name = lt_file.stem

        try:
            with open(lt_file, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()

            # 首先解析质量信息，建立 atom_type -> mass 的映射
            type_to_mass = self._parse_masses(content)

            # 解析原子，使用质量信息推断元素
            atoms = self._parse_atoms(content, type_to_mass)

            # 解析键
            bonds = self._parse_bonds(content)

            # 对于简单离子（只有一个原子），使用元素名称作为原子名称
            if len(atoms) == 1:
                atoms[0]["name"] = atoms[0]["element"]
            else:
                # 为每个原子添加化学环境信息
                atoms = self._add_chemical_environment(atoms, bonds)

            logger.info(f"Parsed {lt_file.name}: {len(atoms)} atoms, {len(bonds)} bonds")

            return {
                "molecule_name": molecule_name,
                "atoms": atoms,
                "bonds": bonds
            }

        except Exception as e:
            logger.error(f"Failed to parse LT file {lt_file}: {e}")
            import traceback
            traceback.print_exc()
            raise

    def _parse_masses(self, content: str) -> Dict[str, float]:
        """
        解析质量信息

        从 write_once("Data Masses") 部分提取原子类型到质量的映射

        Returns:
            {atom_type: mass}
        """
        type_to_mass = {}

        # 查找 write_once("Data Masses") 部分
        masses_section = re.search(
            r'write_once\("Data Masses"\)\s*\{(.*?)\}',
            content,
            re.DOTALL
        )

        if not masses_section:
            logger.warning("No 'Data Masses' section found in LT file")
            return type_to_mass

        masses_text = masses_section.group(1)

        # 解析每一行质量
        # 格式：@atom:S 32.060 或 @atom:type1 12.011
        mass_pattern = r'@atom:(\w+)\s+([\d.]+)'

        for match in re.finditer(mass_pattern, masses_text):
            atom_type = match.group(1)
            mass = float(match.group(2))
            type_to_mass[atom_type] = mass

        logger.info(f"Parsed {len(type_to_mass)} atom types with masses")
        return type_to_mass

    def _parse_atoms(self, content: str, type_to_mass: Dict[str, float]) -> List[Dict]:
        """
        解析原子信息

        从 write("Data Atoms") 部分提取原子信息
        格式：$atom:atomname $mol @atom:type charge x y z

        Args:
            content: .lt 文件内容
            type_to_mass: 原子类型到质量的映射
        """
        atoms = []

        # 查找 write("Data Atoms") 部分
        atoms_section = re.search(
            r'write\("Data Atoms"\)\s*\{(.*?)\}',
            content,
            re.DOTALL
        )

        if not atoms_section:
            logger.warning("No 'Data Atoms' section found in LT file")
            return atoms

        atoms_text = atoms_section.group(1)

        # 解析每一行原子
        # 格式1：$atom:S1 $mol @atom:S 1.020000 -1.492 -0.011 0.142
        # 格式2：$atom:type1 $mol:m1 @atom:type1 0.574893 1.50600 -0.15100 0.08600
        # 注意：$mol 后面可能有 :xxx 也可能没有
        atom_pattern = r'\$atom:(\w+)\s+\$mol(?::\w+)?\s+@atom:(\w+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)'

        for idx, match in enumerate(re.finditer(atom_pattern, atoms_text), start=1):
            atom_name = match.group(1)
            atom_type = match.group(2)
            charge = float(match.group(3))

            # 从质量推断元素
            mass = type_to_mass.get(atom_type, 1.0)
            element = self._get_element_from_mass(mass)

            atoms.append({
                "index": idx,
                "name": atom_name,
                "type": atom_type,
                "element": element,
                "charge": charge,
                "mass": mass
            })

        return atoms

    def _parse_bonds(self, content: str) -> List[tuple]:
        """
        解析键信息

        从 write("Data Bonds") 部分提取键信息
        """
        bonds = []

        # 查找 write("Data Bonds") 部分
        bonds_section = re.search(
            r'write\("Data Bonds"\)\s*\{(.*?)\}',
            content,
            re.DOTALL
        )

        if not bonds_section:
            logger.warning("No 'Data Bonds' section found in LT file")
            return bonds

        bonds_text = bonds_section.group(1)

        # 解析每一行键
        # 格式：$bond:b1 @bond:type $atom:S1 $atom:F3
        bond_pattern = r'\$bond:\w+\s+@bond:\w+\s+\$atom:(\w+)\s+\$atom:(\w+)'

        atom_name_to_idx = {}

        # 首先需要建立原子名称到索引的映射
        atoms_section = re.search(
            r'write\("Data Atoms"\)\s*\{(.*?)\}',
            content,
            re.DOTALL
        )

        if atoms_section:
            atoms_text = atoms_section.group(1)
            atom_pattern = r'\$atom:(\w+)'
            for idx, match in enumerate(re.finditer(atom_pattern, atoms_text), start=1):
                atom_name = match.group(1)
                atom_name_to_idx[atom_name] = idx

        # 解析键
        for match in re.finditer(bond_pattern, bonds_text):
            atom1_name = match.group(1)
            atom2_name = match.group(2)

            if atom1_name in atom_name_to_idx and atom2_name in atom_name_to_idx:
                idx1 = atom_name_to_idx[atom1_name]
                idx2 = atom_name_to_idx[atom2_name]
                bonds.append((idx1, idx2))

        return bonds

    def _add_chemical_environment(self, atoms: List[Dict], bonds: List[tuple]) -> List[Dict]:
        """
        为每个原子添加化学环境信息

        Args:
            atoms: 原子列表
            bonds: 键列表 [(idx1, idx2), ...]

        Returns:
            添加了 environment 字段的原子列表
        """
        # 构建邻接表
        neighbors = {atom['index']: [] for atom in atoms}
        for idx1, idx2 in bonds:
            neighbors[idx1].append(idx2)
            neighbors[idx2].append(idx1)

        # 为每个原子推断化学环境
        for atom in atoms:
            idx = atom['index']
            element = atom['element']
            neighbor_indices = neighbors.get(idx, [])
            neighbor_elements = [atoms[n-1]['element'] for n in neighbor_indices if 1 <= n <= len(atoms)]

            # 根据元素和邻居推断环境
            environment = self._infer_environment(element, neighbor_elements)
            atom['environment'] = environment

        return atoms

    def _infer_environment(self, element: str, neighbor_elements: List[str]) -> str:
        """
        根据元素和邻居推断化学环境

        Args:
            element: 原子元素
            neighbor_elements: 邻居元素列表

        Returns:
            化学环境描述
        """
        if element == 'O':
            # 氧原子的环境判断
            if 'C' in neighbor_elements:
                # 简化判断：如果只连接一个原子，可能是羰基氧
                if len(neighbor_elements) == 1:
                    return 'carbonyl'
                else:
                    return 'ether'
            elif 'S' in neighbor_elements:
                return 'sulfonyl'
            else:
                return element
        elif element == 'S':
            # 硫原子的环境判断
            if 'O' in neighbor_elements and 'N' in neighbor_elements:
                return 'sulfonyl'
            else:
                return element
        elif element == 'N':
            # 氮原子的环境判断
            if 'S' in neighbor_elements:
                s_count = neighbor_elements.count('S')
                if s_count >= 2:
                    return 'imide'
                else:
                    return element
            else:
                return element
        elif element == 'F':
            # 氟原子的环境判断
            if 'S' in neighbor_elements:
                return 'fluorosulfonyl'
            elif 'C' in neighbor_elements:
                return 'fluorocarbon'
            else:
                return element
        elif element == 'C':
            # 碳原子的环境判断
            o_count = neighbor_elements.count('O')
            if o_count >= 2:
                return 'carbonyl'
            else:
                return element
        else:
            return element


    def _extract_element_from_name(self, atom_name: str) -> str:
        """
        从原子名称提取元素符号

        例如：S1 -> S, O01 -> O, C00 -> C, H06 -> H
        """
        # 移除数字，保留字母
        element = re.sub(r'\d+', '', atom_name)

        # 如果是两个字母，首字母大写，第二个字母小写
        if len(element) >= 2:
            element = element[0].upper() + element[1:].lower()
        else:
            element = element.upper()

        return element

    def _get_element_from_mass(self, mass: float) -> str:
        """
        根据原子质量推断元素符号

        Args:
            mass: 原子质量

        Returns:
            元素符号
        """
        # 常见元素的质量范围（允许一定误差）
        mass_to_element = [
            (1.008, "H", 0.1),
            (6.941, "Li", 0.5),
            (12.011, "C", 0.5),
            (14.007, "N", 0.5),
            (15.999, "O", 0.5),
            (18.998, "F", 0.5),
            (22.990, "Na", 0.5),
            (30.974, "P", 0.5),
            (32.06, "S", 0.5),
            (35.45, "Cl", 0.5),
            (39.098, "K", 0.5),
        ]

        for ref_mass, element, tolerance in mass_to_element:
            if abs(mass - ref_mass) < tolerance:
                return element

        # 如果找不到匹配，返回 X
        logger.warning(f"Unknown element for mass {mass}")
        return "X"

