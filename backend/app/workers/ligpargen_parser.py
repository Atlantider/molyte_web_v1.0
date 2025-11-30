"""
LigParGen 文件解析器

解析 LigParGen 生成的 .openmm.xml 文件，提取原子信息和键连接信息
"""

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging

logger = logging.getLogger(__name__)


class LigParGenParser:
    """解析 LigParGen 生成的文件，提取原子信息"""
    
    def parse_openmm_xml(self, xml_file: Path, molecule_name: Optional[str] = None) -> Dict:
        """
        解析 .openmm.xml 文件

        Args:
            xml_file: XML 文件路径
            molecule_name: 分子名称（可选，如果不提供则从文件名提取）

        Returns:
            {
                "molecule_name": "TEP",
                "atoms": [
                    {
                        "index": 1,
                        "name": "P00",
                        "type": "opls_800",
                        "element": "P",
                        "mass": 30.974
                    },
                    ...
                ],
                "bonds": [(1, 2), (1, 3), ...]
            }

        Raises:
            FileNotFoundError: XML 文件不存在
            ValueError: XML 文件格式错误
        """
        if not xml_file.exists():
            raise FileNotFoundError(f"XML file not found: {xml_file}")

        # 如果没有提供分子名称，从文件名提取
        if molecule_name is None:
            # 例如: TEP.openmm.xml -> TEP
            molecule_name = xml_file.stem.replace(".openmm", "")

        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()

            # 解析原子类型
            atom_types = {}
            for atom_type in root.findall(".//AtomTypes/Type"):
                type_name = atom_type.get("name")
                atom_types[type_name] = {
                    "element": atom_type.get("element"),
                    "mass": float(atom_type.get("mass"))
                }

            # 解析原子
            atoms = []
            residue = root.find(".//Residues/Residue")
            if residue is None:
                raise ValueError("No Residue found in XML file")
            
            for idx, atom in enumerate(residue.findall("Atom"), start=1):
                atom_name = atom.get("name")
                atom_type = atom.get("type")
                
                if atom_type not in atom_types:
                    logger.warning(f"Unknown atom type: {atom_type}")
                    continue
                
                type_info = atom_types[atom_type]
                
                atoms.append({
                    "index": idx,
                    "name": atom_name,
                    "type": atom_type,
                    "element": type_info["element"],
                    "mass": type_info["mass"]
                })
            
            # 解析键
            bonds = []
            for bond in residue.findall("Bond"):
                from_idx = int(bond.get("from")) + 1  # XML 从 0 开始，我们从 1 开始
                to_idx = int(bond.get("to")) + 1
                bonds.append((from_idx, to_idx))
            
            logger.info(f"Parsed {xml_file.name}: {len(atoms)} atoms, {len(bonds)} bonds")
            
            return {
                "molecule_name": molecule_name,
                "atoms": atoms,
                "bonds": bonds
            }
            
        except ET.ParseError as e:
            raise ValueError(f"Failed to parse XML file: {e}")
    
    def generate_atom_labels(
        self,
        molecule_name: str,
        atoms: List[Dict]
    ) -> Dict[int, str]:
        """
        生成原子标签
        
        Args:
            molecule_name: 分子名称（如 "TEP", "CO2"）
            atoms: 原子列表
        
        Returns:
            {
                1: "TEP_P00",
                2: "TEP_O01",
                3: "TEP_O02",
                ...
            }
        """
        labels = {}
        for atom in atoms:
            # 标签格式：分子名_原子名
            label = f"{molecule_name}_{atom['name']}"
            labels[atom['index']] = label
        
        return labels
    
    def create_simple_ion_info(
        self,
        ion_name: str,
        element: str,
        mass: float
    ) -> Dict:
        """
        为简单离子（Li, Na, K 等）创建原子信息
        
        Args:
            ion_name: 离子名称（如 "Li", "Na"）
            element: 元素符号
            mass: 原子质量
        
        Returns:
            与 parse_openmm_xml() 相同格式的字典
        """
        return {
            "molecule_name": ion_name,
            "atoms": [
                {
                    "index": 1,
                    "name": ion_name,
                    "type": "simple_ion",
                    "element": element,
                    "mass": mass
                }
            ],
            "bonds": []
        }

