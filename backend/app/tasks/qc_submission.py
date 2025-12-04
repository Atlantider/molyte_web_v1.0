"""
QC 任务提交 Celery Worker

负责异步提交 QC 任务到 Slurm 集群
"""

import logging
import os
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional

from celery import Task
from sqlalchemy.orm import Session

from app.celery_app import celery_app
from app.database import SessionLocal
from app.models.qc import QCJob, QCJobStatus
from app.core.config import settings
from app.schemas.qc import (
    QCAccuracyLevel,
    SolventModel,
    QC_ACCURACY_PRESETS,
    GAUSSIAN_SOLVENTS
)

logger = logging.getLogger(__name__)


class DatabaseTask(Task):
    """自动管理数据库会话的任务基类"""
    _db: Session = None

    def after_return(self, *args, **kwargs):
        if self._db is not None:
            self._db.close()
            self._db = None

    @property
    def db(self) -> Session:
        if self._db is None:
            self._db = SessionLocal()
        return self._db


def get_predefined_coordinates(smiles: str, molecule_name: str = ""):
    """
    获取预定义的3D坐标（从initial_salts目录的PDB文件读取）

    返回格式: [(atom_symbol, x, y, z), ...]
    """
    from pathlib import Path
    from app.core.config import settings

    # 尝试从分子名称匹配PDB文件
    initial_salts_path = settings.MOLYTE_INITIAL_SALTS_PATH

    # 清理分子名称（去除+/-符号）
    clean_name = molecule_name.replace("+", "").replace("-", "").strip()

    # 可能的PDB文件路径
    possible_paths = [
        initial_salts_path / f"{clean_name}.pdb",
        initial_salts_path / f"{molecule_name}.pdb",
        initial_salts_path / f"{molecule_name.upper()}.pdb",
        initial_salts_path / f"{molecule_name.lower()}.pdb",
    ]

    for pdb_path in possible_paths:
        if pdb_path.exists():
            try:
                coords = parse_pdb_coordinates(pdb_path)
                if coords:
                    logger.info(f"Loaded coordinates for {molecule_name} from {pdb_path}")
                    return coords
            except Exception as e:
                logger.warning(f"Failed to parse PDB file {pdb_path}: {e}")
                continue

    logger.debug(f"No PDB file found for {molecule_name} in {initial_salts_path}")
    return None


def parse_pdb_coordinates(pdb_path):
    """
    从PDB文件解析原子坐标

    返回格式: [(atom_symbol, x, y, z), ...]
    """
    coords = []

    with open(pdb_path, 'r') as f:
        for line in f:
            # 解析ATOM或HETATM行
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # PDB格式：
                # HETATM    1  P           0      -0.089   0.595   0.000                       P
                # 列位置：
                # 13-16: 原子名称
                # 31-38: x坐标
                # 39-46: y坐标
                # 47-54: z坐标
                # 77-78: 元素符号（如果有）

                try:
                    atom_name = line[12:16].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())

                    # 尝试获取元素符号（最后两列）
                    if len(line) > 76:
                        element = line[76:78].strip()
                    else:
                        # 如果没有元素符号，从原子名称推断
                        element = atom_name[0]

                    # 如果元素符号为空，使用原子名称的第一个字符
                    if not element:
                        element = atom_name[0]

                    coords.append((element, x, y, z))
                except (ValueError, IndexError) as e:
                    logger.warning(f"Failed to parse PDB line: {line.strip()}, error: {e}")
                    continue

    return coords if coords else None


def get_charge_and_spin_from_smiles(smiles: str):
    """
    从SMILES获取电荷和自旋多重度

    返回: (charge, spin_multiplicity)

    注意：必须添加氢原子后再计算电子数！
    """
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0, 1

        # 计算总电荷
        total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())

        # 添加氢原子以获得完整的电子数
        mol_with_h = Chem.AddHs(mol)

        # 计算总电子数
        total_electrons = sum(atom.GetAtomicNum() for atom in mol_with_h.GetAtoms())
        total_electrons -= total_charge  # 减去电荷

        # 检查是否有显式自由基
        unpaired_electrons = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())

        if unpaired_electrons > 0:
            # 有显式自由基
            spin_multiplicity = unpaired_electrons + 1
        else:
            # 根据电子数判断
            # 偶数电子 -> 闭壳层 -> 自旋多重度 = 1
            # 奇数电子 -> 开壳层 -> 自旋多重度 = 2
            if total_electrons % 2 == 0:
                spin_multiplicity = 1
            else:
                spin_multiplicity = 2

        logger.debug(f"Charge and spin from SMILES {smiles[:30]}...: charge={total_charge}, electrons={total_electrons}, radicals={unpaired_electrons}, spin={spin_multiplicity}")
        return total_charge, spin_multiplicity

    except Exception as e:
        logger.warning(f"Failed to parse SMILES {smiles}: {e}")
        return 0, 1


def get_accuracy_params(job: QCJob) -> tuple:
    """
    根据精度等级获取泛函和基组参数

    Returns:
        (functional, basis_set)
    """
    config = job.config or {}
    accuracy_level = config.get("accuracy_level", "standard")

    if accuracy_level == QCAccuracyLevel.CUSTOM.value or accuracy_level == "custom":
        # 自定义模式使用任务中指定的参数
        return job.functional or "B3LYP", job.basis_set or "6-31G(d)"

    # 预设模式
    try:
        level = QCAccuracyLevel(accuracy_level)
        preset = QC_ACCURACY_PRESETS.get(level, QC_ACCURACY_PRESETS[QCAccuracyLevel.STANDARD])
        return preset["functional"], preset["basis_set"]
    except (ValueError, KeyError):
        # 默认使用标准精度
        return "B3LYP", "6-31G(d)"


def build_scrf_keyword(job: QCJob) -> tuple:
    """
    构建SCRF（溶剂效应）关键字

    Returns:
        tuple: (SCRF关键字字符串, 自定义溶剂参数块字符串)
               如果是气相则返回 ("", "")
    """
    config = job.config or {}
    solvent_config = config.get("solvent_config", {})

    if not solvent_config:
        return "", ""

    solvent_model = solvent_config.get("model", "gas")

    if solvent_model == SolventModel.GAS.value or solvent_model == "gas":
        return "", ""

    if solvent_model == SolventModel.PCM.value or solvent_model == "pcm":
        solvent_name = solvent_config.get("solvent_name", "Water")
        return f"SCRF=(PCM,Solvent={solvent_name})", ""

    if solvent_model == SolventModel.SMD.value or solvent_model == "smd":
        solvent_name = solvent_config.get("solvent_name", "Water")
        return f"SCRF=(SMD,Solvent={solvent_name})", ""

    if solvent_model == SolventModel.CUSTOM.value or solvent_model == "custom":
        # 自定义溶剂参数 - 使用 SMD 模型的 Generic 溶剂
        eps = solvent_config.get("eps", 78.3553)
        eps_inf = solvent_config.get("eps_inf", 1.778)
        hbond_acidity = solvent_config.get("hbond_acidity", 0.82)
        hbond_basicity = solvent_config.get("hbond_basicity", 0.35)
        surface_tension = solvent_config.get("surface_tension", 71.99)
        carbon_aromaticity = solvent_config.get("carbon_aromaticity", 0.0)
        halogenicity = solvent_config.get("halogenicity", 0.0)

        # SCRF关键字指定使用SMD模型的Generic溶剂并读取参数
        scrf_keyword = "SCRF=(SMD,Solvent=Generic,Read)"

        # 自定义溶剂参数块（放在分子坐标后面）
        # Gaussian格式：eps=xxx epsinf=xxx HBondAcidity=xxx HBondBasicity=xxx SurfaceTension=xxx AromaticCarbon=xxx Halogen=xxx
        custom_params = (
            f"eps={eps}\n"
            f"epsinf={eps_inf}\n"
            f"HBondAcidity={hbond_acidity}\n"
            f"HBondBasicity={hbond_basicity}\n"
            f"SurfaceTension={surface_tension}\n"
            f"AromaticCarbon={carbon_aromaticity}\n"
            f"Halogen={halogenicity}\n"
        )

        return scrf_keyword, custom_params

    return "", ""


def sanitize_filename(name: str) -> str:
    """
    将文件名转换为安全的ASCII格式
    移除或替换中文和特殊字符
    """
    import re
    import unicodedata

    # 将中文转换为拼音（如果有pypinyin）或使用hash
    try:
        from pypinyin import lazy_pinyin
        safe_name = '_'.join(lazy_pinyin(name))
    except ImportError:
        # 没有pypinyin，使用分子名的hash作为后备
        import hashlib
        # 保留字母数字和下划线
        ascii_chars = ''.join(c if c.isalnum() or c in '_-' else '' for c in name)
        if not ascii_chars:
            # 全是中文，使用hash
            safe_name = f"mol_{hashlib.md5(name.encode()).hexdigest()[:8]}"
        else:
            safe_name = ascii_chars

    # 移除可能导致shell问题的字符：括号、空格、引号等
    safe_name = re.sub(r'[()（）\[\]{}\'\"<>|&;`$!#\s]', '_', safe_name)
    # 移除连续下划线
    safe_name = re.sub(r'_+', '_', safe_name)
    # 移除首尾下划线
    safe_name = safe_name.strip('_')

    return safe_name or f"molecule_{id(name) % 10000}"


def generate_gaussian_input(job: QCJob, work_dir: Path, pdb_content: str = None) -> str:
    """
    生成Gaussian输入文件(.gjf)

    Args:
        job: QC任务对象
        work_dir: 工作目录
        pdb_content: PDB文件内容（如果有）

    Returns:
        gjf文件路径
    """
    # 使用安全的文件名
    safe_name = sanitize_filename(job.molecule_name)
    gjf_path = work_dir / f"{safe_name}.gjf"

    # 使用任务中的电荷和自旋多重度，或者从SMILES计算
    charge = job.charge
    spin = job.spin_multiplicity

    config = job.config or {}
    auto_spin = config.get("auto_spin", True)

    if auto_spin and charge == 0 and spin == 1:
        charge, spin = get_charge_and_spin_from_smiles(job.smiles)

    # 根据精度等级获取泛函和基组
    functional, basis_set = get_accuracy_params(job)

    # 构建SCRF关键字 (返回 tuple: (关键字, 自定义参数块))
    scrf_keyword, custom_solvent_params = build_scrf_keyword(job)

    # 构建计算关键字行
    keywords = f"opt freq {functional}/{basis_set}"

    # 添加色散校正（对于DFT方法）
    if functional.upper() not in ["HF"]:
        keywords += " em=gd3bj"

    # 添加溶剂效应
    if scrf_keyword:
        keywords += f" {scrf_keyword}"

    # 生成gjf内容 - 使用安全文件名
    gjf_content = f"""%nprocshared=16
%mem=8GB
%chk={safe_name}.chk
# {keywords}

{job.molecule_name}

{charge} {spin}
"""

    # 首先检查 config 中是否有 xyz_content（来自去溶剂化任务的团簇坐标）
    xyz_content = config.get("xyz_content")

    if xyz_content:
        # 使用 config 中的 XYZ 坐标（团簇计算）
        lines = xyz_content.strip().split('\n')
        # XYZ 格式：第一行是原子数，第二行是注释，后面是坐标
        for line in lines[2:]:  # 跳过原子数和注释行
            parts = line.split()
            if len(parts) >= 4:
                atom_symbol = parts[0]
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                gjf_content += f" {atom_symbol:<2}  {x:>12.6f}  {y:>12.6f}  {z:>12.6f}\n"
        logger.info(f"Used xyz_content from config: {len(lines) - 2} atoms for {job.molecule_name}")
    elif pdb_content:
        # 如果有PDB内容，提取坐标
        for line in pdb_content.split('\n'):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                if len(parts) >= 8:
                    atom_type = parts[-1] if len(parts[-1]) <= 2 else parts[2][0]
                    x, y, z = parts[5], parts[6], parts[7]
                    gjf_content += f" {atom_type:<2}  {float(x):>12.6f}  {float(y):>12.6f}  {float(z):>12.6f}\n"
    else:
        # 没有PDB内容时，需要从SMILES生成3D坐标
        # 首先检查是否有预定义的坐标
        coords = get_predefined_coordinates(job.smiles, job.molecule_name)

        if coords:
            # 使用预定义坐标
            for atom_symbol, x, y, z in coords:
                gjf_content += f" {atom_symbol:<2}  {x:>12.6f}  {y:>12.6f}  {z:>12.6f}\n"
        else:
            # 从SMILES生成3D坐标
            try:
                from rdkit import Chem
                from rdkit.Chem import AllChem

                mol = Chem.MolFromSmiles(job.smiles)
                if mol:
                    mol = Chem.AddHs(mol)

                    # 尝试多种方法生成3D坐标
                    result = AllChem.EmbedMolecule(mol, randomSeed=42)

                    if result == -1:
                        # 第一次失败，尝试使用随机坐标
                        logger.warning(f"First attempt failed, trying with random coords for {job.smiles}")
                        result = AllChem.EmbedMolecule(mol, useRandomCoords=True, maxAttempts=100, randomSeed=42)

                    if result == -1:
                        # 第二次失败，尝试使用ETKDGv3方法
                        logger.warning(f"Second attempt failed, trying with ETKDGv3 for {job.smiles}")
                        params = AllChem.ETKDGv3()
                        params.randomSeed = 42
                        result = AllChem.EmbedMolecule(mol, params)

                    if result == -1:
                        raise ValueError(f"无法生成3D坐标")

                    # 尝试优化几何结构
                    try:
                        AllChem.MMFFOptimizeMolecule(mol)
                    except Exception as opt_error:
                        logger.warning(f"MMFF optimization failed: {opt_error}, using unoptimized coordinates")

                    conf = mol.GetConformer()
                    for i, atom in enumerate(mol.GetAtoms()):
                        pos = conf.GetAtomPosition(i)
                        symbol = atom.GetSymbol()
                        gjf_content += f" {symbol:<2}  {pos.x:>12.6f}  {pos.y:>12.6f}  {pos.z:>12.6f}\n"
            except Exception as e:
                logger.error(f"Failed to generate 3D coordinates: {e}")
                raise ValueError(f"无法为 {job.molecule_name} (SMILES: {job.smiles}) 生成3D坐标。\n"
                               f"建议：1) 检查SMILES是否正确；2) 为特殊分子上传PDB文件；3) 联系管理员添加预定义坐标。\n"
                               f"错误详情: {str(e)}")

    gjf_content += "\n"

    # 如果有自定义溶剂参数，添加参数块
    if custom_solvent_params:
        gjf_content += custom_solvent_params
        gjf_content += "\n"

    gjf_content += "\n"

    with open(gjf_path, 'w') as f:
        f.write(gjf_content)

    logger.info(f"Generated Gaussian input with: functional={functional}, basis_set={basis_set}, scrf={scrf_keyword or 'none'}, custom_solvent={bool(custom_solvent_params)}")

    # 返回文件路径和安全文件名
    return str(gjf_path), safe_name


def generate_slurm_script(job: QCJob, work_dir: Path, safe_name: str) -> str:
    """
    生成Slurm提交脚本

    Args:
        job: QC任务对象
        work_dir: 工作目录
        safe_name: 安全的文件名（无中文和特殊字符）
    """
    script_path = work_dir / "run_qc.sh"

    # 从配置获取资源参数
    config = job.config or {}
    partition = config.get("slurm_partition", "cpu")
    cpus = config.get("slurm_cpus", 16)
    time_limit = config.get("slurm_time", 7200)  # 分钟

    # 使用安全的job name（截断以符合Slurm限制）
    safe_job_name = f"QC_{safe_name}"[:64]

    script_content = f"""#!/bin/bash
#SBATCH --job-name={safe_job_name}
#SBATCH --output=qc_out.log
#SBATCH --error=qc_err.log
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --time={time_limit}

# 进入工作目录
cd $SLURM_SUBMIT_DIR

# 设置Gaussian环境
export g16root=/public/software
export GAUSS_SCRDIR=/public/software/g16/scratch
source /public/software/g16/bsd/g16.profile

# 运行Gaussian - 使用安全文件名
ulimit -s unlimited
g16 < "{safe_name}.gjf" > "{safe_name}_out.log" 2>&1

# 转换checkpoint文件
if [ -f "{safe_name}.chk" ]; then
    formchk "{safe_name}.chk" "{safe_name}.fchk"
fi

echo "QC calculation completed"
echo "Original molecule name: {job.molecule_name}"
"""

    with open(script_path, 'w') as f:
        f.write(script_content)

    os.chmod(script_path, 0o755)

    return str(script_path)


@celery_app.task(
    bind=True,
    base=DatabaseTask,
    name="app.tasks.qc_submission.submit_qc_job_task",
    max_retries=3,
    default_retry_delay=60,
)
def submit_qc_job_task(self, job_id: int) -> Dict[str, Any]:
    """
    异步提交 QC 任务到 Slurm 集群
    """
    db = self.db

    try:
        logger.info(f"[Task {self.request.id}] Starting QC job submission for job_id={job_id}")

        # 1. 从数据库读取任务信息
        job = db.query(QCJob).filter(QCJob.id == job_id).first()
        if not job:
            error_msg = f"QC Job {job_id} not found in database"
            logger.error(f"[Task {self.request.id}] {error_msg}")
            return {"success": False, "error": error_msg}

        logger.info(f"[Task {self.request.id}] QC Job: {job.molecule_name}, SMILES: {job.smiles[:30]}...")

        # 2. 创建工作目录
        qc_base_path = getattr(settings, 'QC_WORK_BASE_PATH', '/public/home/xiaoji/molyte_web/data/qc_jobs')
        work_dir = Path(qc_base_path) / f"qc_job_{job_id}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        work_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"[Task {self.request.id}] Work directory: {work_dir}")

        # 3. 生成Gaussian输入文件
        try:
            gjf_path, safe_name = generate_gaussian_input(job, work_dir)
            logger.info(f"[Task {self.request.id}] Generated Gaussian input: {gjf_path} (safe_name={safe_name})")
        except Exception as e:
            error_msg = f"Failed to generate Gaussian input: {str(e)}"
            logger.error(f"[Task {self.request.id}] {error_msg}")
            job.status = QCJobStatus.FAILED
            job.error_message = error_msg
            db.commit()
            return {"success": False, "error": error_msg}

        # 4. 生成Slurm脚本 - 传入safe_name
        script_path = generate_slurm_script(job, work_dir, safe_name)
        logger.info(f"[Task {self.request.id}] Generated Slurm script: {script_path}")

        # 5. 提交到Slurm
        try:
            result = subprocess.run(
                ["sbatch", script_path],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                check=True
            )

            # 解析slurm job id
            # 输出格式: "Submitted batch job 12345"
            output = result.stdout.strip()
            slurm_job_id = output.split()[-1]

            logger.info(f"[Task {self.request.id}] Submitted to Slurm with job_id={slurm_job_id}")

        except subprocess.CalledProcessError as e:
            error_msg = f"Slurm submission failed: {e.stderr}"
            logger.error(f"[Task {self.request.id}] {error_msg}")
            job.status = QCJobStatus.FAILED
            job.error_message = error_msg
            db.commit()
            return {"success": False, "error": error_msg}

        # 6. 更新数据库
        job.status = QCJobStatus.QUEUED
        job.slurm_job_id = str(slurm_job_id)
        job.work_dir = str(work_dir)
        job.started_at = datetime.now()
        job.config = job.config or {}
        job.config["submitted_at"] = datetime.now().isoformat()
        job.config["gjf_path"] = gjf_path
        job.config["safe_name"] = safe_name  # 保存安全文件名用于后处理
        db.commit()

        logger.info(f"[Task {self.request.id}] QC Job {job_id} submitted successfully")

        return {
            "success": True,
            "job_id": job_id,
            "slurm_job_id": slurm_job_id,
            "work_dir": str(work_dir),
        }

    except Exception as exc:
        logger.exception(f"[Task {self.request.id}] Exception during QC job submission: {exc}")

        try:
            job = db.query(QCJob).filter(QCJob.id == job_id).first()
            if job:
                job.status = QCJobStatus.FAILED
                job.error_message = f"Task exception: {str(exc)}"
                db.commit()
        except Exception as db_exc:
            logger.error(f"[Task {self.request.id}] Failed to update job status: {db_exc}")

        raise self.retry(exc=exc, countdown=60)

