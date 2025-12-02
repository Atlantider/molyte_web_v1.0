"""
智能 RDF 对生成器 - 专家级版本

根据电解液配方和原子映射，智能生成需要计算的 RDF 对
支持全面的化学环境识别、优先级排序和智能数量控制

作者: Molyte Web Team
版本: 2.0 (Expert Edition)
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Set
from collections import defaultdict

logger = logging.getLogger(__name__)


class SmartRDFGenerator:
    """智能 RDF 对生成器 - 通用版本

    核心思路：
    1. 不预设化学环境，从 atom_mapping.json 自动提取所有配位原子
    2. 基于元素类型和分子类型智能分组
    3. 使用优先级系统：水 > 阴离子主原子 > 溶剂配位原子
    4. 自动处理新分子和新元素，无需修改代码

    优势：
    - ✅ 支持任意新阴离子（Cl, Br, OH, I, ClO4等）
    - ✅ 支持任意新溶剂（DMF, NMP, DMA等）
    - ✅ 自动识别所有配位原子（O, F, N, S, P, Cl, Br, I等）
    - ✅ 无需维护庞大的预设表格
    """

    # 配位元素的优先级（基于电负性和配位能力）
    ELEMENT_PRIORITY = {
        'O': 1,   # 氧：最常见的配位原子
        'N': 2,   # 氮：腈基、酰胺、胺类
        'F': 3,   # 氟：阴离子和溶剂
        'S': 4,   # 硫：砜类、亚砜
        'Cl': 5,  # 氯：氯化物、高氯酸根
        'Br': 6,  # 溴：溴化物
        'P': 7,   # 磷：磷酸酯、PF6
        'I': 8,   # 碘：碘化物
        'B': 9,   # 硼：硼酸酯、BF4
    }

    # 水分子的识别关键词
    WATER_KEYWORDS = ['H2O', 'WAT', 'HOH', 'TIP3', 'SPC', 'TIP4', 'OPC']

    # 每个分子的最大配位原子数（避免生成过多 RDF 对）
    MAX_ATOMS_PER_MOLECULE = 5

    # 总 RDF 对数限制
    MAX_TOTAL_RDF_PAIRS = 30
    
    def __init__(self, work_dir: Path):
        """初始化生成器
        
        Args:
            work_dir: 工作目录，包含 atom_mapping.json
        """
        self.work_dir = Path(work_dir)
        self.atom_mapping = self._load_atom_mapping()
    
    def _load_atom_mapping(self) -> Dict:
        """加载原子映射文件"""
        mapping_file = self.work_dir / "atom_mapping.json"

        if not mapping_file.exists():
            logger.warning(f"Atom mapping file not found: {mapping_file}")
            return {}

        try:
            with open(mapping_file, 'r') as f:
                data = json.load(f)
                logger.info(f"Loaded atom mapping from {mapping_file}")
                return data
        except Exception as e:
            logger.error(f"Failed to load atom mapping: {e}")
            return {}

    def _extract_cations(self, composition: Dict) -> List[str]:
        """从配方中提取阳离子"""
        cations = []
        if 'cations' in composition and composition['cations']:
            for cation_info in composition['cations']:
                if isinstance(cation_info, dict) and 'name' in cation_info:
                    cations.append(cation_info['name'])
                elif isinstance(cation_info, str):
                    cations.append(cation_info)
        return cations

    def _extract_anions(self, composition: Dict) -> List[str]:
        """从配方中提取阴离子"""
        anions = []
        if 'anions' in composition and composition['anions']:
            for anion_info in composition['anions']:
                if isinstance(anion_info, dict) and 'name' in anion_info:
                    anions.append(anion_info['name'])
                elif isinstance(anion_info, str):
                    anions.append(anion_info)
        return anions

    def _extract_solvents(self, composition: Dict) -> List[str]:
        """从配方中提取溶剂"""
        solvents = []
        if 'solvents' in composition and composition['solvents']:
            for solvent_info in composition['solvents']:
                if isinstance(solvent_info, dict) and 'name' in solvent_info:
                    solvents.append(solvent_info['name'])
                elif isinstance(solvent_info, str):
                    solvents.append(solvent_info)
        return solvents

    def _is_water_molecule(self, molecule_name: str) -> bool:
        """判断是否为水分子"""
        return any(keyword in molecule_name.upper() for keyword in self.WATER_KEYWORDS)

    def _get_molecule_info(self, molecule_name: str) -> Dict:
        """获取分子的所有信息

        Returns:
            {
                'name': str,
                'type': str,  # 'cation', 'anion', 'solvent', 'water'
                'atoms': [
                    {
                        'label': str,
                        'element': str,
                        'atom_name': str,
                        ...
                    },
                    ...
                ]
            }
        """
        if not self.atom_mapping or 'molecules' not in self.atom_mapping:
            return None

        for mol in self.atom_mapping['molecules']:
            if mol.get('molecule_name') == molecule_name:
                return mol

        return None

    def _extract_coordination_atoms(self, molecule_name: str) -> Dict[str, List[Dict]]:
        """提取分子中的所有配位原子（通用方法）

        **重要**: 同一分子中相同化学环境的原子会被合并
        - 根据 environment 字段判断化学环境
        - 相同环境的原子合并为一个RDF对（如FSI的6个F都是fluorosulfonyl，合并为一个）
        - 不同环境的原子分开（如EC的O01是carbonyl，O02/O03是ether，不合并）

        Returns:
            {
                'O': [
                    {'label': 'EC_O(carbonyl)', 'environment': 'carbonyl', 'count': 1, ...},
                    {'label': 'EC_O(ether)', 'environment': 'ether', 'count': 2, ...}
                ],
                'F': [{'label': 'FSI_F(fluorosulfonyl)', 'environment': 'fluorosulfonyl', 'count': 6, ...}],
                ...
            }
        """
        mol_info = self._get_molecule_info(molecule_name)
        if not mol_info or 'atoms' not in mol_info:
            logger.warning(f"No atom info found for molecule: {molecule_name}")
            return {}

        # 按 (元素, 环境) 统计原子数量
        # key: (element, environment), value: count
        env_counts = defaultdict(int)
        env_mol_type = {}

        for atom in mol_info['atoms']:
            element = atom.get('element', '')

            # 只关注常见的配位元素
            if element in self.ELEMENT_PRIORITY:
                # 获取化学环境，默认为元素本身
                environment = atom.get('environment', element)
                if not environment:
                    environment = element

                key = (element, environment)
                env_counts[key] += 1
                env_mol_type[key] = mol_info.get('molecule_type', 'unknown')

        # 生成合并后的配位原子信息（按化学环境分组）
        coord_atoms = defaultdict(list)

        for (element, environment), count in env_counts.items():
            # 根据环境生成标签
            if environment and environment != element:
                # 有具体化学环境，使用 "分子名_元素(环境)" 格式
                merged_label = f"{molecule_name}_{element}({environment})"
                display_env = environment
            else:
                # 没有特定环境，使用 "分子名_元素*" 通配符格式
                merged_label = f"{molecule_name}_{element}*"
                display_env = element

            coord_atoms[element].append({
                'label': merged_label,
                'atom_name': f'{element}({display_env})',
                'element': element,
                'environment': environment,
                'molecule': molecule_name,
                'molecule_type': env_mol_type[(element, environment)],
                'count': count  # 记录原子数量，用于日志和描述
            })
            logger.debug(f"Merged {count} {element}({environment}) atoms in {molecule_name} -> {merged_label}")

        return dict(coord_atoms)

    def _calculate_priority(self, atom_info: Dict, is_water: bool = False) -> int:
        """计算原子的优先级

        优先级规则：
        1. 水分子：0（最高）
        2. 阴离子：10 + 元素优先级
        3. 溶剂：20 + 元素优先级
        4. 阳离子：30 + 元素优先级
        """
        if is_water:
            return 0

        element = atom_info.get('element', '')
        mol_type = atom_info.get('molecule_type', 'unknown')

        element_priority = self.ELEMENT_PRIORITY.get(element, 99)

        if mol_type == 'anion':
            return 10 + element_priority
        elif mol_type == 'solvent':
            return 20 + element_priority
        elif mol_type == 'cation':
            return 30 + element_priority
        else:
            return 90 + element_priority

    def _generate_description(self, cation: str, atom_info: Dict, is_water: bool = False) -> str:
        """生成 RDF 对的描述

        格式: {阳离子}-{元素}({环境})_{分子名} x{数量}
        例如: Li-O(carbonyl)_EC x1, Li-O(ether)_EC x2, Na-F(fluorosulfonyl)_FSI x6
        """
        if is_water:
            return f"{cation}-O(water)_H2O"

        element = atom_info.get('element', '')
        molecule = atom_info.get('molecule', '')
        environment = atom_info.get('environment', element)
        count = atom_info.get('count', 1)

        # 使用化学环境信息
        if environment and environment != element:
            env_hint = f"({environment})"
        else:
            env_hint = ""

        # 如果有多个等价原子，显示数量
        count_hint = f" x{count}" if count > 1 else ""

        return f"{cation}-{element}{env_hint}_{molecule}{count_hint}"

    def generate_rdf_pairs(self, composition: Dict) -> List[Tuple[str, str, str]]:
        """生成 RDF 对（主方法 - 通用版本）

        策略：
        1. 从 atom_mapping.json 自动提取所有配位原子
        2. 优先处理水分子（如果存在）
        3. 按优先级排序：水 > 阴离子 > 溶剂
        4. 限制数量：每个分子最多5个原子，总数最多30个

        Args:
            composition: 电解液配方

        Returns:
            [(center_label, shell_label, description), ...]
        """
        logger.info("=" * 80)
        logger.info("开始生成 RDF 对（通用智能模式）")
        logger.info("=" * 80)

        # 1. 提取阴阳离子和溶剂
        cations = self._extract_cations(composition)
        anions = self._extract_anions(composition)
        solvents = self._extract_solvents(composition)

        logger.info(f"阳离子: {cations}")
        logger.info(f"阴离子: {anions}")
        logger.info(f"溶剂: {solvents}")

        if not cations:
            logger.warning("未找到阳离子，无法生成 RDF 对")
            return []

        # 2. 收集所有候选 RDF 对（带优先级）
        rdf_candidates = []  # [(cation_label, shell_label, description, priority), ...]

        for cation in cations:
            # 获取阳离子的标签
            cation_label = self._get_cation_label(cation)
            if not cation_label:
                logger.warning(f"未找到阳离子 {cation} 的标签")
                continue

            # 2.1 检查水分子（最高优先级）
            water_atoms = self._get_water_atoms()
            if water_atoms:
                logger.info(f"✓ 检测到水分子，添加 {len(water_atoms)} 个水-氧 RDF 对")
                for atom_info in water_atoms[:self.MAX_ATOMS_PER_MOLECULE]:
                    desc = self._generate_description(cation, atom_info, is_water=True)
                    priority = 0  # 最高优先级
                    rdf_candidates.append((cation_label, atom_info['label'], desc, priority))

            # 2.2 处理阴离子
            for anion in anions:
                coord_atoms = self._extract_coordination_atoms(anion)

                for element, atoms in coord_atoms.items():
                    # 限制每个分子的原子数
                    for atom_info in atoms[:self.MAX_ATOMS_PER_MOLECULE]:
                        desc = self._generate_description(cation, atom_info)
                        priority = self._calculate_priority(atom_info)
                        rdf_candidates.append((cation_label, atom_info['label'], desc, priority))

                logger.info(f"  阴离子 {anion}: 添加了 {sum(len(atoms) for atoms in coord_atoms.values())} 个配位原子")

            # 2.3 处理溶剂
            for solvent in solvents:
                # 跳过水分子（已在上面处理）
                if self._is_water_molecule(solvent):
                    continue

                coord_atoms = self._extract_coordination_atoms(solvent)

                for element, atoms in coord_atoms.items():
                    # 限制每个分子的原子数
                    for atom_info in atoms[:self.MAX_ATOMS_PER_MOLECULE]:
                        desc = self._generate_description(cation, atom_info)
                        priority = self._calculate_priority(atom_info)
                        rdf_candidates.append((cation_label, atom_info['label'], desc, priority))

                logger.info(f"  溶剂 {solvent}: 添加了 {sum(len(atoms) for atoms in coord_atoms.values())} 个配位原子")

        # 3. 按优先级排序
        rdf_candidates.sort(key=lambda x: x[3])

        # 4. 限制总数
        rdf_candidates = rdf_candidates[:self.MAX_TOTAL_RDF_PAIRS]

        # 5. 返回结果（去掉优先级）
        rdf_pairs = [(center, shell, desc) for center, shell, desc, priority in rdf_candidates]

        logger.info("=" * 80)
        logger.info(f"✓ 生成了 {len(rdf_pairs)} 个 RDF 对:")
        for i, (center, shell, desc) in enumerate(rdf_pairs, 1):
            logger.info(f"  {i:2d}. {desc:30s} : {center} -> {shell}")
        logger.info("=" * 80)

        return rdf_pairs

    def _get_cation_label(self, cation_name: str) -> str:
        """获取阳离子的标签"""
        mol_info = self._get_molecule_info(cation_name)
        if not mol_info or 'atoms' not in mol_info:
            # 尝试简单的命名规则
            return f"{cation_name}_{cation_name}"

        # 返回第一个原子的标签（通常阳离子只有一个原子）
        if mol_info['atoms']:
            return mol_info['atoms'][0].get('label', f"{cation_name}_{cation_name}")

        return f"{cation_name}_{cation_name}"

    def _get_water_atoms(self) -> List[Dict]:
        """获取所有水分子的氧原子"""
        water_atoms = []

        if not self.atom_mapping or 'molecules' not in self.atom_mapping:
            return water_atoms

        for mol in self.atom_mapping['molecules']:
            mol_name = mol.get('molecule_name', '')

            if self._is_water_molecule(mol_name):
                # 提取水分子中的氧原子
                for atom in mol.get('atoms', []):
                    if atom.get('element') == 'O':
                        water_atoms.append({
                            'label': atom.get('label', ''),
                            'atom_name': atom.get('atom_name', ''),
                            'element': 'O',
                            'molecule': mol_name,
                            'molecule_type': 'water'
                        })

        return water_atoms
