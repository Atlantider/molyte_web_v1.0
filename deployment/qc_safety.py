"""
QC 计算安全机制模块

提供多层安全机制确保 QC 计算成功完成：
1. 结构预检查 - 检查原子距离、键长等
2. 收敛辅助关键词 - 根据体系特点选择合适的 SCF 策略
3. 错误检测与自动重试 - 分析错误类型并应用修复策略
"""

import re
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class QCErrorType(Enum):
    """QC 计算错误类型"""
    SCF_NOT_CONVERGED = "scf_not_converged"
    OPT_NOT_CONVERGED = "opt_not_converged"
    LINEAR_DEPENDENCY = "linear_dependency"
    MEMORY_ERROR = "memory_error"
    DISK_ERROR = "disk_error"
    CHARGE_SPIN_ERROR = "charge_spin_error"
    GEOMETRY_ERROR = "geometry_error"
    BASIS_SET_ERROR = "basis_set_error"
    UNKNOWN = "unknown"


@dataclass
class StructureCheckResult:
    """结构检查结果"""
    is_valid: bool
    warnings: List[str]
    errors: List[str]
    suggestions: List[str]
    min_distance: float = 0.0
    max_distance: float = 0.0
    atom_count: int = 0


@dataclass
class QCErrorAnalysis:
    """QC 错误分析结果"""
    error_type: QCErrorType
    error_message: str
    suggestions: List[str]
    can_retry: bool
    retry_keywords: Optional[str] = None
    retry_config: Optional[Dict] = None


# 元素共价半径 (Angstrom)
COVALENT_RADII = {
    'H': 0.31, 'He': 0.28,
    'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
    'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
    'K': 2.03, 'Ca': 1.76, 'Mn': 1.39, 'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22,
    'Br': 1.20, 'I': 1.39,
}

# 元素范德华半径 (Angstrom)
VDW_RADII = {
    'H': 1.20, 'He': 1.40,
    'Li': 1.82, 'Be': 1.53, 'B': 1.92, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'Ne': 1.54,
    'Na': 2.27, 'Mg': 1.73, 'Al': 1.84, 'Si': 2.10, 'P': 1.80, 'S': 1.80, 'Cl': 1.75, 'Ar': 1.88,
    'K': 2.75, 'Ca': 2.31, 'Mn': 2.05, 'Fe': 2.05, 'Co': 2.00, 'Ni': 1.63, 'Cu': 1.40, 'Zn': 1.39,
    'Br': 1.85, 'I': 1.98,
}

# SCF 收敛策略（按优先级排序）
SCF_CONVERGENCE_STRATEGIES = [
    {
        'name': 'default',
        'keywords': '',
        'description': '默认设置'
    },
    {
        'name': 'tight_scf',
        'keywords': 'scf=(maxcycle=200,xqc)',
        'description': '增加 SCF 迭代次数，使用二次收敛'
    },
    {
        'name': 'qc_algo',
        'keywords': 'scf=(maxcycle=300,qc)',
        'description': '使用二次收敛算法'
    },
    {
        'name': 'novaracc',
        'keywords': 'scf=(maxcycle=300,xqc,novaracc)',
        'description': '禁用可变精度，增加稳定性'
    },
    {
        'name': 'vshift',
        'keywords': 'scf=(maxcycle=300,xqc,vshift=500)',
        'description': '添加虚拟轨道移位，帮助收敛'
    },
    {
        'name': 'fermi',
        'keywords': 'scf=(maxcycle=400,xqc,fermi)',
        'description': '使用 Fermi 展宽，适合金属体系'
    },
    {
        'name': 'nodamping',
        'keywords': 'scf=(maxcycle=400,qc,nodamping,novaracc)',
        'description': '禁用阻尼，最后手段'
    }
]

# 几何优化策略
OPT_STRATEGIES = [
    {
        'name': 'default',
        'keywords': 'opt',
        'description': '默认优化'
    },
    {
        'name': 'calcfc',
        'keywords': 'opt=(calcfc,maxcycles=200)',
        'description': '计算初始力常数矩阵'
    },
    {
        'name': 'calcall',
        'keywords': 'opt=(calcall,maxcycles=200)',
        'description': '每步计算力常数矩阵'
    },
    {
        'name': 'gdiis',
        'keywords': 'opt=(gdiis,maxcycles=300)',
        'description': '使用 GDIIS 优化算法'
    },
    {
        'name': 'tight',
        'keywords': 'opt=(tight,calcfc,maxcycles=300)',
        'description': '严格收敛标准'
    }
]


def parse_xyz_coordinates(xyz_content: str) -> List[Tuple[str, float, float, float]]:
    """解析 XYZ 坐标内容"""
    atoms = []
    lines = xyz_content.strip().split('\n')

    for line in lines:
        line = line.strip()
        if not line:
            continue

        parts = line.split()
        if len(parts) >= 4:
            element = parts[0]
            # 移除元素后的数字（如 C1 -> C）
            element = ''.join(c for c in element if not c.isdigit())
            try:
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                atoms.append((element, x, y, z))
            except ValueError:
                continue

    return atoms


def calculate_distance(atom1: Tuple[str, float, float, float],
                       atom2: Tuple[str, float, float, float]) -> float:
    """计算两个原子之间的距离"""
    import math
    dx = atom1[1] - atom2[1]
    dy = atom1[2] - atom2[2]
    dz = atom1[3] - atom2[3]
    return math.sqrt(dx*dx + dy*dy + dz*dz)


def check_structure(xyz_content: str, charge: int = 0,
                    spin_multiplicity: int = 1) -> StructureCheckResult:
    """
    检查分子结构的合理性

    Args:
        xyz_content: XYZ 坐标内容
        charge: 分子电荷
        spin_multiplicity: 自旋多重度

    Returns:
        StructureCheckResult 结构检查结果
    """
    warnings = []
    errors = []
    suggestions = []

    atoms = parse_xyz_coordinates(xyz_content)

    if not atoms:
        return StructureCheckResult(
            is_valid=False,
            warnings=[],
            errors=["无法解析坐标，结构为空"],
            suggestions=["请检查 XYZ 坐标格式"],
            atom_count=0
        )

    atom_count = len(atoms)
    min_distance = float('inf')
    max_distance = 0.0

    # 检查原子间距离
    too_close_pairs = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            dist = calculate_distance(atoms[i], atoms[j])
            min_distance = min(min_distance, dist)
            max_distance = max(max_distance, dist)

            elem1, elem2 = atoms[i][0], atoms[j][0]

            # 检查是否太近（小于共价半径之和的 0.5 倍）
            r1 = COVALENT_RADII.get(elem1, 1.5)
            r2 = COVALENT_RADII.get(elem2, 1.5)
            min_allowed = (r1 + r2) * 0.5

            if dist < min_allowed:
                too_close_pairs.append((i, j, elem1, elem2, dist, min_allowed))

    if too_close_pairs:
        for i, j, e1, e2, dist, min_d in too_close_pairs[:5]:  # 最多显示5对
            errors.append(f"原子 {e1}({i+1}) 和 {e2}({j+1}) 距离过近: {dist:.3f} Å < {min_d:.3f} Å")
        if len(too_close_pairs) > 5:
            errors.append(f"... 还有 {len(too_close_pairs)-5} 对原子距离过近")
        suggestions.append("建议：使用分子力学或半经验方法先进行结构优化")

    # 检查体系大小
    if atom_count > 200:
        warnings.append(f"体系较大 ({atom_count} 原子)，计算可能较慢")
        suggestions.append("建议：考虑使用较小的基组或半经验方法")
    elif atom_count > 100:
        warnings.append(f"体系中等大小 ({atom_count} 原子)")

    # 检查电荷和自旋多重度的一致性
    total_electrons = 0
    for elem, _, _, _ in atoms:
        atomic_number = get_atomic_number(elem)
        total_electrons += atomic_number

    total_electrons -= charge

    # 检查自旋多重度
    unpaired = spin_multiplicity - 1
    if (total_electrons - unpaired) % 2 != 0:
        errors.append(f"电荷 ({charge}) 和自旋多重度 ({spin_multiplicity}) 与电子数 ({total_electrons}) 不匹配")
        correct_spin = (total_electrons % 2) + 1
        suggestions.append(f"建议：尝试自旋多重度 = {correct_spin}")

    is_valid = len(errors) == 0

    return StructureCheckResult(
        is_valid=is_valid,
        warnings=warnings,
        errors=errors,
        suggestions=suggestions,
        min_distance=min_distance if min_distance != float('inf') else 0.0,
        max_distance=max_distance,
        atom_count=atom_count
    )


def get_atomic_number(element: str) -> int:
    """获取元素的原子序数"""
    atomic_numbers = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
        'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
        'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26,
        'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34,
        'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42,
        'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
        'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56,
    }
    return atomic_numbers.get(element.capitalize(), 0)


def analyze_gaussian_error(log_content: str) -> QCErrorAnalysis:
    """
    分析 Gaussian 输出文件中的错误

    Args:
        log_content: Gaussian 日志文件内容

    Returns:
        QCErrorAnalysis 错误分析结果
    """
    error_type = QCErrorType.UNKNOWN
    error_message = ""
    suggestions = []
    can_retry = False
    retry_keywords = None
    retry_config = None

    # SCF 不收敛
    if re.search(r'Convergence failure|SCF has not converged|Convergence criterion not met',
                 log_content, re.IGNORECASE):
        error_type = QCErrorType.SCF_NOT_CONVERGED
        error_message = "SCF 迭代未收敛"
        suggestions = [
            "增加 SCF 迭代次数",
            "使用二次收敛算法 (scf=qc)",
            "添加虚拟轨道移位 (scf=vshift=500)",
            "对于金属体系使用 Fermi 展宽"
        ]
        can_retry = True
        retry_keywords = "scf=(maxcycle=300,xqc,novaracc)"

    # 几何优化不收敛
    elif re.search(r'Optimization stopped|Number of steps exceeded|'
                   r'Convergence failure -- run terminated', log_content, re.IGNORECASE):
        if 'opt' in log_content.lower()[:1000]:
            error_type = QCErrorType.OPT_NOT_CONVERGED
            error_message = "几何优化未收敛"
            suggestions = [
                "增加优化步数 (opt=maxcycles=300)",
                "使用 GDIIS 优化算法",
                "计算初始力常数矩阵 (opt=calcfc)",
                "检查初始结构是否合理"
            ]
            can_retry = True
            retry_keywords = "opt=(calcfc,maxcycles=300)"

    # 线性依赖问题
    elif re.search(r'Linear dependency|basis set is linearly dependent',
                   log_content, re.IGNORECASE):
        error_type = QCErrorType.LINEAR_DEPENDENCY
        error_message = "基组线性依赖"
        suggestions = [
            "使用较小的基组",
            "移除弥散函数",
            "使用 int=nobasistransform 关键词"
        ]
        can_retry = True
        retry_keywords = "int=nobasistransform"

    # 内存不足
    elif re.search(r'Out of memory|galloc|malloc failed|insufficient memory',
                   log_content, re.IGNORECASE):
        error_type = QCErrorType.MEMORY_ERROR
        error_message = "内存不足"
        suggestions = [
            "增加分配的内存",
            "使用较小的基组",
            "使用分子分段方法"
        ]
        can_retry = True
        retry_config = {'mem': '16GB'}

    # 磁盘空间不足
    elif re.search(r'No space left|disk quota|write error', log_content, re.IGNORECASE):
        error_type = QCErrorType.DISK_ERROR
        error_message = "磁盘空间不足"
        suggestions = [
            "清理磁盘空间",
            "检查 scratch 目录配额"
        ]
        can_retry = False

    # 电荷/自旋错误
    elif re.search(r'Charge and multiplicity|Wrong number of electrons|'
                   r'The combination of multiplicity', log_content, re.IGNORECASE):
        error_type = QCErrorType.CHARGE_SPIN_ERROR
        error_message = "电荷或自旋多重度设置错误"
        suggestions = [
            "检查分子电荷设置",
            "检查自旋多重度是否与电子数匹配"
        ]
        can_retry = False

    # 几何结构错误
    elif re.search(r'Atom too close|Bad geometry|Small interatomic distance',
                   log_content, re.IGNORECASE):
        error_type = QCErrorType.GEOMETRY_ERROR
        error_message = "分子几何结构有问题"
        suggestions = [
            "检查初始结构是否有原子重叠",
            "使用分子力学预优化结构"
        ]
        can_retry = False

    # 基组错误
    elif re.search(r'Unknown basis set|basis set not found', log_content, re.IGNORECASE):
        error_type = QCErrorType.BASIS_SET_ERROR
        error_message = "基组不支持"
        suggestions = [
            "检查基组名称是否正确",
            "该元素可能不支持此基组"
        ]
        can_retry = False

    # 尝试提取具体错误信息
    error_match = re.search(r'Error termination.*?(?:via|in)\s+(\w+)', log_content)
    if error_match:
        error_message += f" (in {error_match.group(1)})"

    # 如果没有识别出具体错误，检查是否有 Error termination
    if error_type == QCErrorType.UNKNOWN and 'Error termination' in log_content:
        error_message = "Gaussian 计算异常终止"
        # 提取错误行
        for line in log_content.split('\n'):
            if 'Error' in line and len(line) < 200:
                error_message = line.strip()
                break

    return QCErrorAnalysis(
        error_type=error_type,
        error_message=error_message,
        suggestions=suggestions,
        can_retry=can_retry,
        retry_keywords=retry_keywords,
        retry_config=retry_config
    )


def get_robust_keywords(atom_count: int, has_metal: bool = False,
                        retry_level: int = 0) -> str:
    """
    根据体系特点和重试级别生成稳健的计算关键词

    Args:
        atom_count: 原子数
        has_metal: 是否含有金属
        retry_level: 重试级别 (0-6)

    Returns:
        SCF 关键词字符串
    """
    if retry_level >= len(SCF_CONVERGENCE_STRATEGIES):
        retry_level = len(SCF_CONVERGENCE_STRATEGIES) - 1

    strategy = SCF_CONVERGENCE_STRATEGIES[retry_level]
    keywords = strategy['keywords']

    # 大体系默认增加迭代次数
    if atom_count > 50 and retry_level == 0:
        keywords = "scf=(maxcycle=200)"

    # 金属体系使用更稳健的设置
    if has_metal and retry_level < 4:
        if retry_level == 0:
            keywords = "scf=(maxcycle=200,xqc)"
        elif retry_level == 1:
            keywords = "scf=(maxcycle=300,xqc,vshift=300)"

    return keywords


def generate_retry_gjf(original_gjf_path: Path, retry_level: int,
                       error_analysis: QCErrorAnalysis = None) -> str:
    """
    生成重试用的 GJF 文件内容

    Args:
        original_gjf_path: 原始 GJF 文件路径
        retry_level: 重试级别
        error_analysis: 错误分析结果

    Returns:
        新的 GJF 文件内容
    """
    with open(original_gjf_path, 'r') as f:
        content = f.read()

    lines = content.split('\n')
    new_lines = []

    # 找到关键词行（以 # 开头）
    for i, line in enumerate(lines):
        if line.strip().startswith('#'):
            # 解析现有关键词
            keywords_line = line.strip()

            # 如果有错误分析建议的关键词，使用它
            if error_analysis and error_analysis.retry_keywords:
                # 检查是否已有 scf 或 opt 关键词
                if 'scf=' in error_analysis.retry_keywords:
                    # 移除原有的 scf 设置
                    keywords_line = re.sub(r'scf=\([^)]+\)', '', keywords_line)
                    keywords_line = re.sub(r'scf=\w+', '', keywords_line)
                    keywords_line += f' {error_analysis.retry_keywords}'
                elif 'opt=' in error_analysis.retry_keywords:
                    # 移除原有的 opt 设置并替换
                    keywords_line = re.sub(r'opt=\([^)]+\)', '', keywords_line)
                    keywords_line = re.sub(r'\bopt\b', '', keywords_line)
                    keywords_line += f' {error_analysis.retry_keywords}'
                else:
                    keywords_line += f' {error_analysis.retry_keywords}'
            else:
                # 使用通用的重试策略
                strategy = SCF_CONVERGENCE_STRATEGIES[min(retry_level, len(SCF_CONVERGENCE_STRATEGIES)-1)]
                if strategy['keywords']:
                    keywords_line = re.sub(r'scf=\([^)]+\)', '', keywords_line)
                    keywords_line = re.sub(r'scf=\w+', '', keywords_line)
                    keywords_line += f' {strategy["keywords"]}'

            # 清理多余空格
            keywords_line = ' '.join(keywords_line.split())
            new_lines.append(keywords_line)
        else:
            new_lines.append(line)

    # 更新内存设置（如果需要）
    if error_analysis and error_analysis.retry_config:
        if 'mem' in error_analysis.retry_config:
            for i, line in enumerate(new_lines):
                if line.strip().startswith('%mem='):
                    new_lines[i] = f"%mem={error_analysis.retry_config['mem']}"
                    break

    return '\n'.join(new_lines)


class QCRetryManager:
    """QC 任务重试管理器"""

    MAX_RETRIES = 3

    def __init__(self, work_dir: Path, logger=None):
        self.work_dir = work_dir
        self.logger = logger or logging.getLogger(__name__)
        self.retry_count = 0
        self.retry_history = []

    def should_retry(self, log_content: str) -> Tuple[bool, Optional[str]]:
        """
        判断是否应该重试

        Returns:
            (should_retry, new_gjf_content)
        """
        if self.retry_count >= self.MAX_RETRIES:
            self.logger.warning(f"已达到最大重试次数 ({self.MAX_RETRIES})")
            return False, None

        # 分析错误
        error_analysis = analyze_gaussian_error(log_content)
        self.retry_history.append(error_analysis)

        if not error_analysis.can_retry:
            self.logger.info(f"错误类型 {error_analysis.error_type.value} 不可重试")
            return False, None

        # 找到原始 GJF 文件
        gjf_files = list(self.work_dir.glob('*.gjf'))
        if not gjf_files:
            return False, None

        original_gjf = gjf_files[0]

        # 生成重试 GJF
        self.retry_count += 1
        new_content = generate_retry_gjf(original_gjf, self.retry_count, error_analysis)

        self.logger.info(f"准备第 {self.retry_count} 次重试，策略: {error_analysis.retry_keywords or 'default'}")

        return True, new_content

    def get_retry_summary(self) -> Dict[str, Any]:
        """获取重试摘要"""
        return {
            'total_retries': self.retry_count,
            'max_retries': self.MAX_RETRIES,
            'history': [
                {
                    'error_type': h.error_type.value,
                    'error_message': h.error_message,
                    'suggestions': h.suggestions
                }
                for h in self.retry_history
            ]
        }

