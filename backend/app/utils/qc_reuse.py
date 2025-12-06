"""
QC 任务复用工具函数

确保只有真正有能量结果的任务才能被复用，
并正确处理复用链条和结果复制。

复用规则设计：
1. 坐标指纹匹配（最可靠）：相同几何结构的分子可复用
2. 精确匹配复用：SMILES + charge + spin + functional + basis_set + solvent
3. 名称标准化复用：处理各种命名后缀（#数字、_opt、_neutral 等）
4. 等价分子复用：处理同义名称（Li = Li+ = lithium）
5. 等价基组复用：6-31G(d) = 6-31G*
6. 跨计算类型复用：redox neutral_gas = binding ligand
7. Structure ID 匹配：同一 MD 任务 + 同一 structure_id 的 cluster 可复用

特殊处理：
- Cluster 没有 SMILES，使用 structure_id + md_job_id 匹配
- 阴离子 SMILES 可能不可靠，优先使用 molecule_name 匹配
- 从 cluster 提取的配体使用坐标指纹匹配
"""
import logging
import re
import hashlib
from typing import Optional, Tuple, Set, List, Dict
from sqlalchemy.orm import Session, joinedload

from app.models.qc import QCJob, QCResult, QCJobStatus

logger = logging.getLogger(__name__)


# ============================================================================
# 名称标准化规则
# ============================================================================

# 需要去掉的后缀模式（按顺序从右到左处理）
NAME_SUFFIX_PATTERNS = [
    r'#\d+$',                              # EC#1 -> EC
    r'_opt$', r'_sp$', r'_freq$',          # Li-EC_opt -> Li-EC
    r'_neutral$', r'_charged$', r'_oxidized$', r'_reduced$',  # 状态后缀
    r'_sol$', r'_gas$',                    # 溶剂环境后缀
    r'_at_neutral$', r'_at_ox$', r'_at_red$',  # SP 计算后缀
    r'_neutral_opt$', r'_oxidized_opt$', r'_reduced_opt$',  # 组合后缀
    r'_sp_ox_at_neutral$', r'_sp_neutral_at_ox$',  # Redox 特殊后缀
]

# 等价分子名称映射（标准名称 -> 所有等价名称）
EQUIVALENT_NAMES: Dict[str, List[str]] = {
    # 阳离子
    'Li': ['Li', 'Li+', 'lithium', 'Li_ion', 'li', 'LI'],
    'Na': ['Na', 'Na+', 'sodium', 'Na_ion', 'na', 'NA'],
    'K': ['K', 'K+', 'potassium', 'K_ion', 'k'],
    'Mg': ['Mg', 'Mg2+', 'Mg+2', 'magnesium'],
    'Ca': ['Ca', 'Ca2+', 'Ca+2', 'calcium'],
    'Zn': ['Zn', 'Zn2+', 'Zn+2', 'zinc'],

    # 阴离子
    'PF6': ['PF6', 'PF6-', 'hexafluorophosphate', 'pf6'],
    'TFSI': ['TFSI', 'TFSI-', 'bis(trifluoromethanesulfonyl)imide', 'NTf2', 'NTf2-', 'tfsi'],
    'FSI': ['FSI', 'FSI-', 'bis(fluorosulfonyl)imide', 'fsi'],
    'BF4': ['BF4', 'BF4-', 'tetrafluoroborate', 'bf4'],
    'ClO4': ['ClO4', 'ClO4-', 'perchlorate', 'clo4'],

    # 常见溶剂
    'EC': ['EC', 'ethylene carbonate', 'EthyleneCarbonate', 'ec'],
    'DMC': ['DMC', 'dimethyl carbonate', 'DimethylCarbonate', 'dmc'],
    'EMC': ['EMC', 'ethyl methyl carbonate', 'EthylMethylCarbonate', 'emc'],
    'PC': ['PC', 'propylene carbonate', 'PropyleneCarbonate', 'pc'],
    'DEC': ['DEC', 'diethyl carbonate', 'DiethylCarbonate', 'dec'],
    'VC': ['VC', 'vinylene carbonate', 'VinyleneCarbonate', 'vc'],
    'FEC': ['FEC', 'fluoroethylene carbonate', 'FluoroethyleneCarbonate', 'fec'],
    'DOL': ['DOL', 'dioxolane', '1,3-dioxolane', 'dol'],
    'DME': ['DME', 'dimethoxyethane', '1,2-dimethoxyethane', 'dme'],
}

# 反向映射：任意名称 -> 标准名称
_CANONICAL_NAME_MAP: Dict[str, str] = {}
for canonical, equivalents in EQUIVALENT_NAMES.items():
    for name in equivalents:
        _CANONICAL_NAME_MAP[name.lower()] = canonical

# 等价基组映射
EQUIVALENT_BASIS_SETS: Dict[str, List[str]] = {
    "6-31G(d)": ["6-31G(d)", "6-31G*"],
    "6-31G*": ["6-31G(d)", "6-31G*"],
    "6-31G(d,p)": ["6-31G(d,p)", "6-31G**"],
    "6-31G**": ["6-31G(d,p)", "6-31G**"],
    "6-31+G(d)": ["6-31+G(d)", "6-31+G*"],
    "6-31+G*": ["6-31+G(d)", "6-31+G*"],
    "6-31+G(d,p)": ["6-31+G(d,p)", "6-31+G**"],
    "6-31+G**": ["6-31+G(d,p)", "6-31+G**"],
    "6-311G(d,p)": ["6-311G(d,p)", "6-311G**"],
    "6-311G**": ["6-311G(d,p)", "6-311G**"],
    "6-311+G(d,p)": ["6-311+G(d,p)", "6-311+G**"],
    "6-311+G**": ["6-311+G(d,p)", "6-311+G**"],
    "6-311++G(d,p)": ["6-311++G(d,p)", "6-311++G**"],
    "6-311++G**": ["6-311++G(d,p)", "6-311++G**"],
}

# 不应跨结构复用的任务类型（依赖具体原子坐标）
# 注意：这些任务可以在同一 MD 任务 + 同一 structure_id 下复用
STRUCTURE_DEPENDENT_TASK_TYPES = [
    'cluster',
    'cluster_minus_',
    'intermediate',
    'reorg_cluster_',  # 重组能 cluster
]

# 阴离子列表 - 这些分子的 SMILES 可能不可靠，优先用名称匹配
ANION_NAMES = ['PF6', 'PF6-', 'FSI', 'FSI-', 'TFSI', 'TFSI-', 'BF4', 'BF4-',
               'ClO4', 'ClO4-', 'DCA', 'DCA-', 'OTf', 'OTf-']


# ============================================================================
# 坐标指纹功能
# ============================================================================

def compute_coordinate_fingerprint(xyz_content: str, precision: int = 3) -> Optional[str]:
    """
    计算分子坐标的指纹，用于识别相同几何结构的分子

    算法：
    1. 解析 xyz_content 获取元素和坐标
    2. 按元素排序，再按 x, y, z 坐标排序（消除原子顺序影响）
    3. 将坐标取到指定精度
    4. 计算 hash

    Args:
        xyz_content: XYZ 格式的坐标内容
        precision: 坐标精度（小数位数），默认 3 位

    Returns:
        指纹 hash 字符串，如果解析失败返回 None
    """
    if not xyz_content:
        return None

    try:
        lines = xyz_content.strip().split('\n')

        # 跳过 XYZ 格式的前两行（原子数和注释）
        coord_lines = []
        for i, line in enumerate(lines):
            parts = line.split()
            if len(parts) >= 4:
                # 尝试解析为 元素 x y z 格式
                try:
                    element = parts[0].strip()
                    x = round(float(parts[1]), precision)
                    y = round(float(parts[2]), precision)
                    z = round(float(parts[3]), precision)
                    coord_lines.append((element, x, y, z))
                except (ValueError, IndexError):
                    continue

        if not coord_lines:
            return None

        # 排序：先按元素，再按 x, y, z
        coord_lines.sort(key=lambda c: (c[0], c[1], c[2], c[3]))

        # 构建用于 hash 的字符串
        fingerprint_str = ';'.join(
            f"{e},{x:.{precision}f},{y:.{precision}f},{z:.{precision}f}"
            for e, x, y, z in coord_lines
        )

        # 计算 MD5 hash
        return hashlib.md5(fingerprint_str.encode()).hexdigest()[:16]

    except Exception as e:
        logger.warning(f"计算坐标指纹失败: {e}")
        return None


def is_anion_molecule(molecule_name: str, smiles: Optional[str] = None) -> bool:
    """
    判断是否是阴离子分子

    阴离子的 SMILES 可能不可靠，复用时优先使用名称匹配
    """
    if not molecule_name:
        return False

    # 标准化名称后检查
    normalized = normalize_molecule_name(molecule_name)
    if normalized in ANION_NAMES or f"{normalized}-" in ANION_NAMES:
        return True

    # 检查原始名称
    for anion in ANION_NAMES:
        if anion.replace('-', '') in molecule_name.upper():
            return True

    return False


def normalize_molecule_name(name: str) -> str:
    """
    标准化分子名称，去掉各种后缀和变体

    示例：
    - EC#1 -> EC
    - Li-PF6_neutral_gas -> Li-PF6
    - lithium -> Li
    - 6-31G* -> 6-31G(d) [基组用专门函数]

    Args:
        name: 原始分子名称

    Returns:
        标准化后的名称
    """
    if not name:
        return name

    result = name.strip()

    # 1. 去掉各种后缀（循环处理，直到没有变化）
    changed = True
    while changed:
        changed = False
        for pattern in NAME_SUFFIX_PATTERNS:
            new_result = re.sub(pattern, '', result)
            if new_result != result:
                result = new_result
                changed = True

    # 2. 处理复合名称（如 Li-PF6）：拆分、标准化、重组
    # 支持的分隔符：-, +, _
    if '-' in result or '+' in result:
        # 分隔符处理
        parts = re.split(r'[-+]', result)
        normalized_parts = []
        for part in parts:
            part = part.strip()
            if part:
                # 尝试找到标准名称
                canonical = _CANONICAL_NAME_MAP.get(part.lower())
                normalized_parts.append(canonical if canonical else part)

        # 重新组合（使用 - 作为标准分隔符）
        if len(normalized_parts) > 1:
            # 对于 Li-X 形式，保持阳离子在前
            result = '-'.join(normalized_parts)
        elif normalized_parts:
            result = normalized_parts[0]
    else:
        # 单一分子名称
        canonical = _CANONICAL_NAME_MAP.get(result.lower())
        if canonical:
            result = canonical

    return result


def get_equivalent_basis_sets(basis_set: str) -> List[str]:
    """
    获取等价基组列表

    Args:
        basis_set: 原始基组名称

    Returns:
        所有等价基组名称的列表
    """
    return EQUIVALENT_BASIS_SETS.get(basis_set, [basis_set])


def is_structure_dependent_task(task_type: str) -> bool:
    """
    判断任务是否依赖具体结构（不应跨结构复用）

    Args:
        task_type: 任务类型

    Returns:
        True 如果依赖具体结构
    """
    if not task_type:
        return False

    for pattern in STRUCTURE_DEPENDENT_TASK_TYPES:
        if task_type.startswith(pattern) or task_type == pattern:
            return True

    return False


def extract_base_task_info(task_type: str) -> Tuple[str, Optional[str], Optional[str]]:
    """
    从 task_type 提取基础信息，用于跨计算类型复用

    Args:
        task_type: 任务类型，如 'redox_mol_PF6_neutral_gas', 'ligand_EC', 'dimer_Li-PF6'

    Returns:
        (base_type, molecule_name, state):
        - base_type: 'ion', 'ligand', 'dimer', 'cluster', 'redox_mol', 'redox_dimer'
        - molecule_name: 提取的分子名称
        - state: 状态信息（neutral_gas, charged_sol 等，如果有的话）
    """
    if not task_type:
        return ('unknown', None, None)

    # ion 类型
    if task_type == 'ion':
        return ('ion', 'Li', None)

    # ligand_{name}
    if task_type.startswith('ligand_'):
        mol_name = task_type.replace('ligand_', '')
        return ('ligand', normalize_molecule_name(mol_name), None)

    # dimer_{name}
    if task_type.startswith('dimer_'):
        mol_name = task_type.replace('dimer_', '')
        return ('dimer', normalize_molecule_name(mol_name), None)

    # redox_mol_{name}_{state}
    match = re.match(r'redox_mol_(.+)_(neutral|charged|oxidized|reduced)_(gas|sol)$', task_type)
    if match:
        mol_name = match.group(1)
        state = f"{match.group(2)}_{match.group(3)}"
        return ('redox_mol', normalize_molecule_name(mol_name), state)

    # redox_dimer_{name}_{state}
    match = re.match(r'redox_dimer_(.+)_(neutral|charged|oxidized|reduced)_(gas|sol)$', task_type)
    if match:
        mol_name = match.group(1)
        state = f"{match.group(2)}_{match.group(3)}"
        return ('redox_dimer', normalize_molecule_name(mol_name), state)

    # cluster 类型
    if task_type.startswith('cluster'):
        return ('cluster', task_type, None)

    # 其他未知类型
    return ('other', normalize_molecule_name(task_type), None)


def can_reuse_across_calc_types(
    source_task_type: str,
    target_task_type: str,
    source_charge: int,
    target_charge: int,
    source_solvent_model: str,
    target_solvent_model: str
) -> bool:
    """
    判断两个不同计算类型的任务是否可以复用

    复用规则：
    1. ligand_{name} 可复用于 redox_mol_{name}_neutral_gas（如果都是气相、中性）
    2. dimer_{name} 可复用于 redox_dimer_{name}_neutral_gas
    3. ion 可跨所有类型复用

    Args:
        source_task_type: 源任务类型
        target_task_type: 目标任务类型
        source_charge: 源任务电荷
        target_charge: 目标任务电荷
        source_solvent_model: 源任务溶剂模型
        target_solvent_model: 目标任务溶剂模型

    Returns:
        True 如果可以复用
    """
    # 电荷必须相同
    if source_charge != target_charge:
        return False

    source_base, source_mol, source_state = extract_base_task_info(source_task_type)
    target_base, target_mol, target_state = extract_base_task_info(target_task_type)

    # 分子名称必须匹配
    if source_mol != target_mol:
        return False

    # ion 类型可跨所有类型复用
    if source_base == 'ion' and target_base == 'ion':
        return True

    # ligand <-> redox_mol neutral_gas
    if source_base == 'ligand' and target_base == 'redox_mol':
        # ligand 是气相中性的，redox_mol_neutral_gas 也是
        if target_state == 'neutral_gas' and target_solvent_model in ('gas', None):
            return source_solvent_model in ('gas', None)

    if source_base == 'redox_mol' and target_base == 'ligand':
        if source_state == 'neutral_gas' and source_solvent_model in ('gas', None):
            return target_solvent_model in ('gas', None)

    # dimer <-> redox_dimer neutral_gas
    if source_base == 'dimer' and target_base == 'redox_dimer':
        if target_state == 'neutral_gas' and target_solvent_model in ('gas', None):
            return source_solvent_model in ('gas', None)

    if source_base == 'redox_dimer' and target_base == 'dimer':
        if source_state == 'neutral_gas' and source_solvent_model in ('gas', None):
            return target_solvent_model in ('gas', None)

    # 相同类型可复用
    if source_base == target_base:
        if source_state == target_state:
            return True

    return False


def find_root_job_with_result(
    db: Session, 
    job: QCJob, 
    max_depth: int = 50
) -> Optional[QCJob]:
    """
    追溯复用链，找到有实际能量结果的根任务
    
    Args:
        db: 数据库会话
        job: 起始任务
        max_depth: 最大追溯深度，防止无限循环
        
    Returns:
        有能量结果的根任务，如果链条断裂返回 None
    """
    visited: Set[int] = set()
    current = job
    depth = 0
    
    while current and depth < max_depth:
        if current.id in visited:
            logger.warning(f"检测到复用链循环引用: job {current.id}")
            return None
        visited.add(current.id)
        depth += 1
        
        # 检查当前任务是否有能量结果
        result = db.query(QCResult).filter(
            QCResult.qc_job_id == current.id,
            QCResult.energy_au.isnot(None)
        ).first()
        
        if result:
            logger.debug(f"找到根任务 {current.id}，能量={result.energy_au}")
            return current
        
        # 继续向上追溯
        if current.reused_from_job_id:
            current = db.query(QCJob).filter(
                QCJob.id == current.reused_from_job_id
            ).first()
            if current is None:
                logger.warning(f"复用链断裂: reused_from_job_id={job.reused_from_job_id} 不存在")
                return None
        else:
            # 到达链条末端但没有结果
            logger.warning(f"到达复用链末端但无结果: job {current.id}")
            return None
    
    logger.warning(f"复用链过深(>{max_depth}): 起始 job {job.id}")
    return None


def get_energy_from_job(db: Session, job: QCJob) -> Optional[float]:
    """
    获取任务的能量结果，支持复用任务的链式追溯
    
    Args:
        db: 数据库会话
        job: QC 任务
        
    Returns:
        能量值 (Hartree)，如果无法获取返回 None
    """
    # 首先检查任务自身是否有结果
    result = db.query(QCResult).filter(
        QCResult.qc_job_id == job.id,
        QCResult.energy_au.isnot(None)
    ).first()
    
    if result:
        return result.energy_au
    
    # 如果是复用任务，追溯根任务
    if job.is_reused and job.reused_from_job_id:
        root_job = find_root_job_with_result(db, job)
        if root_job:
            root_result = db.query(QCResult).filter(
                QCResult.qc_job_id == root_job.id,
                QCResult.energy_au.isnot(None)
            ).first()
            if root_result:
                return root_result.energy_au
    
    return None


def validate_job_for_reuse(db: Session, job: QCJob) -> Tuple[bool, Optional[QCJob], str]:
    """
    验证任务是否可以被复用
    
    Args:
        db: 数据库会话
        job: 待验证的 QC 任务
        
    Returns:
        (是否可复用, 根任务, 原因说明)
    """
    if job.status != QCJobStatus.COMPLETED:
        return False, None, f"任务状态为 {job.status}，不是 COMPLETED"
    
    if job.is_deleted:
        return False, None, "任务已被删除"
    
    # 检查是否有直接结果
    result = db.query(QCResult).filter(
        QCResult.qc_job_id == job.id,
        QCResult.energy_au.isnot(None)
    ).first()
    
    if result:
        return True, job, "任务有直接能量结果"
    
    # 如果是复用任务，追溯根任务
    if job.is_reused and job.reused_from_job_id:
        root_job = find_root_job_with_result(db, job)
        if root_job:
            return True, root_job, f"复用链有效，根任务={root_job.id}"
        else:
            return False, None, "复用链断裂，无法找到有效结果"
    
    return False, None, "任务没有能量结果且不是复用任务"


def copy_result_for_reused_job(
    db: Session, 
    source_job: QCJob, 
    target_job: QCJob
) -> Optional[QCResult]:
    """
    将源任务的结果复制到目标复用任务
    
    Args:
        db: 数据库会话
        source_job: 源任务（有结果的根任务）
        target_job: 目标任务（复用任务）
        
    Returns:
        新创建的 QCResult，如果失败返回 None
    """
    # 获取源任务的结果
    source_result = db.query(QCResult).filter(
        QCResult.qc_job_id == source_job.id,
        QCResult.energy_au.isnot(None)
    ).first()
    
    if not source_result:
        logger.error(f"源任务 {source_job.id} 没有有效结果")
        return None
    
    # 检查目标任务是否已有结果
    existing = db.query(QCResult).filter(
        QCResult.qc_job_id == target_job.id
    ).first()
    
    if existing:
        logger.info(f"目标任务 {target_job.id} 已有结果，跳过复制")
        return existing

    # 创建新的结果记录
    new_result = QCResult(
        qc_job_id=target_job.id,
        smiles=source_result.smiles,
        energy_au=source_result.energy_au,
        homo=source_result.homo,
        lumo=source_result.lumo,
        homo_lumo_gap=source_result.homo_lumo_gap,
        esp_min_kcal=source_result.esp_min_kcal,
        esp_max_kcal=source_result.esp_max_kcal,
        dipole_moment=source_result.dipole_moment,
        polarizability=source_result.polarizability,
        vip_ev=source_result.vip_ev,
        vea_ev=source_result.vea_ev,
        oxidation_potential_v=source_result.oxidation_potential_v,
        reduction_potential_v=source_result.reduction_potential_v,
        # 图片和文件路径不复制，保持引用源任务的结果
        additional_properties={
            "copied_from_job_id": source_job.id,
            "copied_from_result_id": source_result.id,
            "is_copied": True
        }
    )

    db.add(new_result)
    db.flush()

    logger.info(f"已将任务 {source_job.id} 的结果复制到任务 {target_job.id}")
    return new_result


def fix_reused_job_without_result(db: Session, job: QCJob) -> Tuple[bool, str]:
    """
    修复没有结果的复用任务

    如果能找到有效的根任务，复制结果；
    否则将任务状态重置为待计算。

    Args:
        db: 数据库会话
        job: 需要修复的任务

    Returns:
        (是否成功, 说明)
    """
    if not job.is_reused:
        return False, "不是复用任务"

    # 检查是否已有结果
    existing_result = db.query(QCResult).filter(
        QCResult.qc_job_id == job.id,
        QCResult.energy_au.isnot(None)
    ).first()

    if existing_result:
        return True, "任务已有结果，无需修复"

    # 尝试追溯根任务
    root_job = find_root_job_with_result(db, job)

    if root_job:
        # 找到根任务，复制结果
        new_result = copy_result_for_reused_job(db, root_job, job)
        if new_result:
            # 更新 reused_from_job_id 直接指向根任务
            job.reused_from_job_id = root_job.id
            return True, f"已从根任务 {root_job.id} 复制结果"
        else:
            return False, "复制结果失败"
    else:
        # 复用链断裂，重置任务状态
        job.status = QCJobStatus.SUBMITTED
        job.is_reused = False
        job.reused_from_job_id = None
        job.error_message = "复用链无效，任务需要重新计算"
        return True, "复用链断裂，已重置为待计算状态"


def batch_fix_invalid_reused_jobs(db: Session, dry_run: bool = True) -> dict:
    """
    批量修复无效的复用任务

    Args:
        db: 数据库会话
        dry_run: 是否仅检测不修改

    Returns:
        修复统计信息
    """
    from sqlalchemy import and_

    stats = {
        "total_reused_jobs": 0,
        "valid_jobs": 0,
        "fixed_by_copy": 0,
        "reset_to_submitted": 0,
        "errors": []
    }

    # 查找所有复用任务
    reused_jobs = db.query(QCJob).filter(
        QCJob.is_reused == True,
        QCJob.status == QCJobStatus.COMPLETED,
        QCJob.is_deleted == False
    ).all()

    stats["total_reused_jobs"] = len(reused_jobs)

    for job in reused_jobs:
        # 检查是否有结果
        has_result = db.query(QCResult).filter(
            QCResult.qc_job_id == job.id,
            QCResult.energy_au.isnot(None)
        ).count() > 0

        if has_result:
            stats["valid_jobs"] += 1
            continue

        # 尝试修复
        if not dry_run:
            success, msg = fix_reused_job_without_result(db, job)
            if success:
                if "复制结果" in msg:
                    stats["fixed_by_copy"] += 1
                elif "重置" in msg:
                    stats["reset_to_submitted"] += 1
            else:
                stats["errors"].append({"job_id": job.id, "error": msg})
        else:
            # 仅检测
            root_job = find_root_job_with_result(db, job)
            if root_job:
                stats["fixed_by_copy"] += 1
            else:
                stats["reset_to_submitted"] += 1

    if not dry_run:
        db.commit()

    return stats


def is_valid_smiles(smiles: Optional[str]) -> bool:
    """
    判断 SMILES 是否有效（非空且是有意义的分子表示）

    无效的 SMILES 包括：
    - None 或空字符串
    - 仅包含占位符（如 "N/A", "unknown", "cluster"）
    - 单个原子离子（如 "[Li+]" 虽然有效但不如 molecule_name 可靠）

    Args:
        smiles: SMILES 字符串

    Returns:
        True 如果是有效的 SMILES
    """
    if not smiles:
        return False

    smiles = smiles.strip()

    # 排除明显无效的占位符
    invalid_patterns = [
        'n/a', 'na', 'unknown', 'cluster', 'none', 'null', ''
    ]
    if smiles.lower() in invalid_patterns:
        return False

    # 排除过短的 SMILES（可能是占位符）
    if len(smiles) < 2:
        return False

    return True


def find_reusable_qc_job(
    db: Session,
    smiles: Optional[str],
    molecule_name: str,
    charge: int,
    spin_multiplicity: int,
    functional: str,
    basis_set: str,
    solvent_model: str = 'gas',
    solvent_name: Optional[str] = None,
    task_type: Optional[str] = None,
    calc_mode: Optional[str] = None
) -> Optional[QCJob]:
    """
    增强的 QC 任务复用查找函数

    复用策略（按优先级）：
    1. 标准化 molecule_name 匹配（最可靠，优先使用）
    2. SMILES 精确匹配（作为补充，仅当有有效 SMILES 时）
    3. 跨计算类型复用（redox_mol_neutral_gas ↔ ligand）

    注意：
    - molecule_name 优先级高于 SMILES，因为 cluster 和某些离子没有有效 SMILES
    - 对于结构依赖的任务（cluster 类型），不进行全局复用
    - 等价基组会自动匹配（如 6-31G(d) = 6-31G*）

    Args:
        db: 数据库会话
        smiles: SMILES 字符串（可选，可能无效）
        molecule_name: 分子名称（主要匹配依据）
        charge: 电荷
        spin_multiplicity: 自旋多重度
        functional: 泛函
        basis_set: 基组
        solvent_model: 溶剂模型
        solvent_name: 溶剂名称
        task_type: 任务类型（用于跨计算类型复用判断）
        calc_mode: 计算模式（opt/sp，用于区分优化和单点）

    Returns:
        可复用的 QCJob，如果没有找到返回 None
    """
    from sqlalchemy import desc, or_
    from sqlalchemy.orm import joinedload

    # 获取等价基组
    equivalent_basis_sets = get_equivalent_basis_sets(basis_set)

    # 标准化分子名称
    normalized_name = normalize_molecule_name(molecule_name) if molecule_name else None

    # 判断是否是结构依赖的任务
    if is_structure_dependent_task(task_type):
        logger.debug(f"任务类型 {task_type} 依赖具体结构，跳过全局复用")
        return None

    # 构建基础过滤条件
    def build_base_query():
        query = db.query(QCJob).options(
            joinedload(QCJob.results)
        ).filter(
            QCJob.charge == charge,
            QCJob.spin_multiplicity == spin_multiplicity,
            QCJob.functional == functional,
            QCJob.basis_set.in_(equivalent_basis_sets),
            QCJob.status == QCJobStatus.COMPLETED,
            QCJob.is_deleted == False
        )

        # 溶剂匹配
        if solvent_model in ('gas', None, ''):
            query = query.filter(
                or_(
                    QCJob.solvent_model == 'gas',
                    QCJob.solvent_model.is_(None)
                )
            )
        else:
            query = query.filter(QCJob.solvent_model == solvent_model)
            if solvent_name:
                query = query.filter(QCJob.solvent_name == solvent_name)

        # 计算模式匹配（opt vs sp）- 仅在明确指定时过滤
        if calc_mode == 'sp':
            # SP 计算：只匹配 SP 类型
            query = query.filter(QCJob.molecule_name.ilike('%_sp_%'))

        return query

    # ========================================
    # 策略 1: 标准化 molecule_name 匹配（优先）
    # ========================================
    # molecule_name 比 SMILES 更可靠，因为：
    # - Cluster 类型没有有效 SMILES
    # - 某些离子（如 PF6-）的 SMILES 可能不标准
    # - molecule_name 经过标准化可以匹配各种变体
    if normalized_name:
        # 构建匹配模式：匹配基础名称开头的所有变体
        name_pattern = f"^{re.escape(normalized_name)}([#_].*)?$"

        query = build_base_query().filter(
            QCJob.molecule_name.op('~')(name_pattern)
        )
        candidates = query.order_by(desc(QCJob.finished_at)).limit(10).all()

        for candidate in candidates:
            # 验证标准化后的名称是否匹配
            candidate_normalized = normalize_molecule_name(candidate.molecule_name)
            if candidate_normalized != normalized_name:
                continue

            is_valid, root_job, reason = validate_job_for_reuse(db, candidate)
            if is_valid:
                result_job = root_job if root_job else candidate
                logger.info(f"[复用-名称标准化匹配] 找到任务 {result_job.id}: "
                           f"'{molecule_name}' -> '{normalized_name}' 匹配 '{candidate.molecule_name}', {reason}")
                return result_job

    # 策略 3: 跨计算类型复用
    if task_type:
        base_type, mol_name, state = extract_base_task_info(task_type)

        # 确定可复用的源类型
        cross_type_patterns = []
        if base_type == 'ligand' or (base_type == 'redox_mol' and state == 'neutral_gas'):
            # ligand 可从 ligand、redox_mol_neutral_gas 复用
            if mol_name:
                cross_type_patterns.append(f"^ligand_{re.escape(mol_name)}.*$")
                cross_type_patterns.append(f"^redox_mol_{re.escape(mol_name)}_neutral_gas.*$")

        if base_type == 'dimer' or (base_type == 'redox_dimer' and state == 'neutral_gas'):
            if mol_name:
                cross_type_patterns.append(f"^dimer_{re.escape(mol_name)}.*$")
                cross_type_patterns.append(f"^redox_dimer_{re.escape(mol_name)}_neutral_gas.*$")

        if base_type == 'ion':
            cross_type_patterns.append('^ion$')

        for pattern in cross_type_patterns:
            query = build_base_query().filter(
                QCJob.molecule_type.op('~')(pattern) |
                QCJob.molecule_name.op('~')(pattern)
            )
            candidates = query.order_by(desc(QCJob.finished_at)).limit(5).all()

            for candidate in candidates:
                # 检查跨类型复用条件
                candidate_task_type = candidate.config.get('task_type') if candidate.config else None
                if candidate_task_type and can_reuse_across_calc_types(
                    candidate_task_type, task_type,
                    candidate.charge, charge,
                    candidate.solvent_model, solvent_model
                ):
                    is_valid, root_job, reason = validate_job_for_reuse(db, candidate)
                    if is_valid:
                        result_job = root_job if root_job else candidate
                        logger.info(f"[复用-跨类型] 找到任务 {result_job.id}: "
                                   f"{candidate_task_type} -> {task_type}, {reason}")
                        return result_job

    logger.debug(f"未找到可复用任务: {molecule_name} (charge={charge}, solvent={solvent_model})")
    return None


def find_reusable_qc_job_enhanced(
    db: Session,
    smiles: Optional[str],
    molecule_name: str,
    charge: int,
    spin_multiplicity: int,
    functional: str,
    basis_set: str,
    solvent_model: str = 'gas',
    solvent_name: Optional[str] = None,
    task_type: Optional[str] = None,
    calc_mode: Optional[str] = None,
    xyz_content: Optional[str] = None,
    structure_id: Optional[int] = None,
    md_job_id: Optional[int] = None
) -> Optional[QCJob]:
    """
    增强版 QC 任务复用查找函数

    新增复用策略：
    1. 坐标指纹匹配（最高优先级）- 相同几何结构的分子
    2. Structure ID 匹配 - 同一 MD 任务 + 同一结构 ID 的 cluster
    3. 标准化名称匹配（阴离子优先使用）
    4. SMILES 匹配（作为补充）
    5. 跨计算类型复用

    Args:
        xyz_content: XYZ 格式的坐标内容（用于指纹匹配）
        structure_id: 溶剂化结构 ID（用于 cluster 匹配）
        md_job_id: MD 任务 ID（用于限定 cluster 复用范围）
    """
    from sqlalchemy import desc, or_
    from sqlalchemy.orm import joinedload

    # 获取等价基组
    equivalent_basis_sets = get_equivalent_basis_sets(basis_set)

    # 标准化分子名称
    normalized_name = normalize_molecule_name(molecule_name) if molecule_name else None

    # 计算坐标指纹
    coord_fingerprint = compute_coordinate_fingerprint(xyz_content) if xyz_content else None

    # 构建基础过滤条件
    def build_base_query():
        query = db.query(QCJob).options(
            joinedload(QCJob.results)
        ).filter(
            QCJob.charge == charge,
            QCJob.spin_multiplicity == spin_multiplicity,
            QCJob.functional == functional,
            QCJob.basis_set.in_(equivalent_basis_sets),
            QCJob.status == QCJobStatus.COMPLETED,
            QCJob.is_deleted == False
        )

        # 溶剂匹配
        if solvent_model in ('gas', None, ''):
            query = query.filter(
                or_(
                    QCJob.solvent_model == 'gas',
                    QCJob.solvent_model.is_(None)
                )
            )
        else:
            query = query.filter(QCJob.solvent_model == solvent_model)
            if solvent_name:
                query = query.filter(QCJob.solvent_name == solvent_name)

        return query

    # ========================================
    # 策略 1: 坐标指纹匹配（最高优先级）
    # ========================================
    if coord_fingerprint:
        # 需要在 config 中存储 coord_fingerprint，这里查找已有的
        # 注意：需要遍历候选任务检查指纹
        candidates = build_base_query().order_by(desc(QCJob.finished_at)).limit(50).all()

        for candidate in candidates:
            config = candidate.config or {}
            candidate_xyz = config.get('xyz_content') or config.get('initial_xyz')
            if candidate_xyz:
                candidate_fp = compute_coordinate_fingerprint(candidate_xyz)
                if candidate_fp == coord_fingerprint:
                    is_valid, root_job, reason = validate_job_for_reuse(db, candidate)
                    if is_valid:
                        result_job = root_job if root_job else candidate
                        logger.info(f"[复用-坐标指纹] 找到任务 {result_job.id}: "
                                   f"指纹={coord_fingerprint}, {reason}")
                        return result_job

    # ========================================
    # 策略 2: Structure ID 匹配（Cluster 类型）
    # ========================================
    if structure_id and md_job_id and is_structure_dependent_task(task_type):
        # 对于 cluster 类型，只在同一 MD 任务 + 同一 structure_id 下复用
        query = build_base_query().filter(
            QCJob.solvation_structure_id == structure_id,
            QCJob.md_job_id == md_job_id
        )
        candidates = query.order_by(desc(QCJob.finished_at)).limit(10).all()

        for candidate in candidates:
            is_valid, root_job, reason = validate_job_for_reuse(db, candidate)
            if is_valid:
                result_job = root_job if root_job else candidate
                logger.info(f"[复用-StructureID] 找到任务 {result_job.id}: "
                           f"structure_id={structure_id}, md_job={md_job_id}, {reason}")
                return result_job

        # 如果是结构依赖的任务且没找到匹配，不再尝试其他策略
        logger.debug(f"Structure-dependent 任务 {task_type} 未找到复用: "
                    f"structure_id={structure_id}")
        return None

    # 对于非结构依赖的任务，继续使用原有策略
    # 如果是阴离子，跳过 SMILES 匹配
    skip_smiles = is_anion_molecule(molecule_name, smiles)

    # ========================================
    # 策略 3: 标准化名称匹配
    # ========================================
    if normalized_name:
        name_pattern = f"^{re.escape(normalized_name)}([#_].*)?$"

        query = build_base_query().filter(
            QCJob.molecule_name.op('~')(name_pattern)
        )
        candidates = query.order_by(desc(QCJob.finished_at)).limit(10).all()

        for candidate in candidates:
            candidate_normalized = normalize_molecule_name(candidate.molecule_name)
            if candidate_normalized != normalized_name:
                continue

            is_valid, root_job, reason = validate_job_for_reuse(db, candidate)
            if is_valid:
                result_job = root_job if root_job else candidate
                logger.info(f"[复用-名称] 找到任务 {result_job.id}: "
                           f"'{molecule_name}' -> '{normalized_name}', {reason}")
                return result_job

    # ========================================
    # 策略 4: SMILES 匹配（非阴离子）
    # ========================================
    if not skip_smiles and is_valid_smiles(smiles):
        query = build_base_query().filter(
            QCJob.smiles == smiles
        )
        candidates = query.order_by(desc(QCJob.finished_at)).limit(10).all()

        for candidate in candidates:
            is_valid, root_job, reason = validate_job_for_reuse(db, candidate)
            if is_valid:
                result_job = root_job if root_job else candidate
                logger.info(f"[复用-SMILES] 找到任务 {result_job.id}: "
                           f"smiles={smiles}, {reason}")
                return result_job

    # ========================================
    # 策略 5: 跨计算类型复用
    # ========================================
    if task_type:
        base_type, mol_name, state = extract_base_task_info(task_type)

        cross_type_patterns = []
        if base_type == 'ligand' or (base_type == 'redox_mol' and state == 'neutral_gas'):
            if mol_name:
                cross_type_patterns.append(f"^ligand_{re.escape(mol_name)}.*$")
                cross_type_patterns.append(f"^redox_mol_{re.escape(mol_name)}_neutral_gas.*$")

        if base_type == 'dimer' or (base_type == 'redox_dimer' and state == 'neutral_gas'):
            if mol_name:
                cross_type_patterns.append(f"^dimer_{re.escape(mol_name)}.*$")
                cross_type_patterns.append(f"^redox_dimer_{re.escape(mol_name)}_neutral_gas.*$")

        if base_type == 'ion':
            cross_type_patterns.append('^ion$')

        for pattern in cross_type_patterns:
            query = build_base_query().filter(
                QCJob.molecule_type.op('~')(pattern) |
                QCJob.molecule_name.op('~')(pattern)
            )
            candidates = query.order_by(desc(QCJob.finished_at)).limit(5).all()

            for candidate in candidates:
                candidate_task_type = candidate.config.get('task_type') if candidate.config else None
                if candidate_task_type and can_reuse_across_calc_types(
                    candidate_task_type, task_type,
                    candidate.charge, charge,
                    candidate.solvent_model, solvent_model
                ):
                    is_valid, root_job, reason = validate_job_for_reuse(db, candidate)
                    if is_valid:
                        result_job = root_job if root_job else candidate
                        logger.info(f"[复用-跨类型] 找到任务 {result_job.id}: "
                                   f"{candidate_task_type} -> {task_type}, {reason}")
                        return result_job

    logger.debug(f"未找到可复用任务: {molecule_name} (charge={charge}, fingerprint={coord_fingerprint})")
    return None
