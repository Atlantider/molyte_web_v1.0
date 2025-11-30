"""
原子映射生成器

基于 .lt 文件生成 atom_mapping.json
"""

import json
from pathlib import Path
from typing import Dict, List, Optional
import logging

from .lt_parser import LTParser

logger = logging.getLogger(__name__)


class AtomMappingGenerator:
    """生成原子映射文件"""
    
    def __init__(self, initial_salts_path: Path):
        """
        Args:
            initial_salts_path: 阳离子/阴离子初始结构文件路径
        """
        self.parser = LTParser()
        self.initial_salts_path = initial_salts_path
    
    def generate_atom_mapping(
        self,
        job_data: Dict,
        work_path: Path
    ) -> Dict:
        """
        生成原子映射文件
        
        Args:
            job_data: 任务数据，格式：
                {
                    "name": "MD-20251119-0001-xxx",
                    "cations": [{"name": "Li", "number": 50}],
                    "anions": [{"name": "PF6", "number": 50}],
                    "solvents": [
                        {"name": "TEP", "smiles": "...", "number": 100},
                        {"name": "CO2", "smiles": "...", "number": 10}
                    ]
                }
            work_path: 工作目录路径
        
        Returns:
            {
              "molecules": [
                {
                  "molecule_id": 1,
                  "molecule_name": "Li",
                  "molecule_type": "cation",
                  "atoms": [
                    {
                      "atom_id": 1,
                      "atom_index": 1,
                      "atom_name": "Li",
                      "label": "Li_Li",
                      "element": "Li",
                      "type": "simple_ion"
                    }
                  ]
                },
                ...
              ]
            }
        """
        mapping = {"molecules": []}
        
        atom_id_counter = 1
        molecule_id_counter = 1
        
        # 处理阳离子
        for cation in job_data.get("cations", []):
            name = cation["name"]
            number = cation["number"]
            
            # 获取分子信息
            mol_info = self._get_molecule_info(name, "cation", work_path)
            
            # 为每个分子实例创建映射
            for i in range(number):
                mol_entry = self._create_molecule_entry(
                    molecule_id_counter,
                    name,
                    "cation",
                    mol_info,
                    atom_id_counter
                )
                
                mapping["molecules"].append(mol_entry)
                atom_id_counter += len(mol_info["atoms"])
                molecule_id_counter += 1
        
        # 处理阴离子
        for anion in job_data.get("anions", []):
            name = anion["name"]
            number = anion["number"]
            
            mol_info = self._get_molecule_info(name, "anion", work_path)
            
            for i in range(number):
                mol_entry = self._create_molecule_entry(
                    molecule_id_counter,
                    name,
                    "anion",
                    mol_info,
                    atom_id_counter
                )
                
                mapping["molecules"].append(mol_entry)
                atom_id_counter += len(mol_info["atoms"])
                molecule_id_counter += 1
        
        # 处理溶剂
        for solvent in job_data.get("solvents", []):
            name = solvent["name"]
            number = solvent["number"]
            
            mol_info = self._get_molecule_info(name, "solvent", work_path)
            
            for i in range(number):
                mol_entry = self._create_molecule_entry(
                    molecule_id_counter,
                    name,
                    "solvent",
                    mol_info,
                    atom_id_counter
                )
                
                mapping["molecules"].append(mol_entry)
                atom_id_counter += len(mol_info["atoms"])
                molecule_id_counter += 1
        
        # 保存到文件
        mapping_file = work_path / "atom_mapping.json"
        with open(mapping_file, 'w') as f:
            json.dump(mapping, f, indent=2)
        
        logger.info(f"Generated atom mapping: {len(mapping['molecules'])} molecules, "
                   f"{atom_id_counter - 1} atoms")

        return mapping

    def _get_molecule_info(
        self,
        name: str,
        mol_type: str,
        work_path: Path
    ) -> Dict:
        """
        获取分子信息（从 .lt 文件读取）

        Args:
            name: 分子名称
            mol_type: 分子类型（cation, anion, solvent）
            work_path: 工作目录

        Returns:
            分子信息字典
        """
        # 尝试从工作目录查找 LT 文件（优先）
        lt_file = work_path / f"{name}.lt"
        if lt_file.exists():
            try:
                logger.info(f"Parsing LT file for {name}: {lt_file}")
                return self.parser.parse_lt_file(lt_file, name)
            except Exception as e:
                logger.error(f"Failed to parse LT file {lt_file}: {e}")

        # 尝试从初始盐目录查找 LT 文件（阳离子/阴离子）
        lt_file = self.initial_salts_path / f"{name}.lt"
        if lt_file.exists():
            try:
                logger.info(f"Parsing LT file for {name}: {lt_file}")
                return self.parser.parse_lt_file(lt_file, name)
            except Exception as e:
                logger.error(f"Failed to parse LT file {lt_file}: {e}")

        # 如果找不到 LT 文件，记录警告并创建简单结构
        logger.warning(f"No LT file found for {name}, creating simple structure")
        return {
            "molecule_name": name,
            "atoms": [
                {
                    "index": 1,
                    "name": name,
                    "type": "unknown",
                    "element": name,
                    "charge": 0.0,
                    "mass": 1.0
                }
            ],
            "bonds": []
        }

    def _create_molecule_entry(
        self,
        molecule_id: int,
        molecule_name: str,
        molecule_type: str,
        mol_info: Dict,
        atom_id_start: int
    ) -> Dict:
        """
        创建单个分子的映射条目

        Args:
            molecule_id: 分子 ID
            molecule_name: 分子名称
            molecule_type: 分子类型
            mol_info: 分子信息
            atom_id_start: 起始原子 ID

        Returns:
            分子映射条目
        """
        mol_entry = {
            "molecule_id": molecule_id,
            "molecule_name": molecule_name,
            "molecule_type": molecule_type,
            "atoms": []
        }

        atom_id = atom_id_start
        for atom in mol_info["atoms"]:
            # 生成标签：优先使用环境信息
            environment = atom.get("environment", "")
            element = atom["element"]

            if environment and environment != element:
                # 有化学环境信息，使用 "分子名_元素(环境)" 格式
                label = f"{molecule_name}_{element}({environment})"
            else:
                # 没有环境信息，使用原子名称
                label = f"{molecule_name}_{atom['name']}"

            mol_entry["atoms"].append({
                "atom_id": atom_id,
                "atom_index": atom["index"],
                "atom_name": atom["name"],
                "label": label,
                "element": atom["element"],
                "environment": environment,
                "type": atom.get("type", "unknown"),
                "charge": atom.get("charge", 0.0),
                "mass": atom.get("mass", 1.0)
            })
            atom_id += 1

        return mol_entry

