"""
RESP charge calculation utilities

支持功能:
1. 生成 RESP 计算的 Slurm 脚本
2. 提交 RESP 任务到 Slurm
3. 解析 RESP 电荷结果
4. 修改 LAMMPS 文件的电荷
5. 计算核时数
"""
import os
import re
import subprocess
import logging
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from datetime import datetime

logger = logging.getLogger(__name__)

# 电荷保存目录
DEFAULT_CHARGE_SAVE_PATH = Path("/public/home/xiaoji/molyte_web/data/charges")

# RESP2.sh 脚本路径
RESP2_SCRIPT_PATH = Path("/public/home/xiaoji/molyte_v1/RESP/RESP2.sh")


class RESPCalculator:
    """RESP 电荷计算器"""

    def __init__(self, charge_save_path: Path = None, slurm_partition: str = "cpu"):
        """
        Args:
            charge_save_path: 电荷文件保存目录
            slurm_partition: Slurm 分区
        """
        self.charge_save_path = charge_save_path or DEFAULT_CHARGE_SAVE_PATH
        self.charge_save_path.mkdir(parents=True, exist_ok=True)
        self.slurm_partition = slurm_partition

    def check_existing_charges(self, molecule_name: str) -> Optional[Path]:
        """
        检查是否已有计算好的电荷文件

        Args:
            molecule_name: 分子名称

        Returns:
            电荷文件路径，不存在返回 None
        """
        charge_file = self.charge_save_path / f"{molecule_name}.charmm.chg"
        if charge_file.exists():
            logger.info(f"Found existing charge file: {charge_file}")
            return charge_file
        return None

    def generate_resp_slurm_script(self, work_dir: Path, pdb_file: str,
                                   molecule_name: str, charge: int = 0,
                                   spin_multiplicity: int = 1,
                                   solvent: str = "water",
                                   cpus: int = 16,
                                   time_limit_hours: int = 24) -> Path:
        """
        生成 RESP 计算的 Slurm 提交脚本

        Args:
            work_dir: 工作目录
            pdb_file: PDB 文件名
            molecule_name: 分子名称
            charge: 分子电荷
            spin_multiplicity: 自旋多重度
            solvent: 溶剂名称
            cpus: CPU 核心数
            time_limit_hours: 最大运行时间（小时）

        Returns:
            脚本路径
        """
        script_path = work_dir / "resp_job.sh"
        safe_name = molecule_name.replace(" ", "_").replace("/", "_")[:32]

        # 将电荷保存路径写入脚本
        charge_save_file = self.charge_save_path / f"{molecule_name}.charmm.chg"

        script_content = f"""#!/bin/bash
#SBATCH --job-name=RESP_{safe_name}
#SBATCH --output=resp_out.log
#SBATCH --error=resp_err.log
#SBATCH --partition={self.slurm_partition}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --time={time_limit_hours}:00:00

# 进入工作目录
cd $SLURM_SUBMIT_DIR

echo "========================================"
echo "RESP Calculation Start"
echo "Molecule: {molecule_name}"
echo "PDB File: {pdb_file}"
echo "Charge: {charge}, Spin: {spin_multiplicity}"
echo "Solvent: {solvent}"
echo "Start Time: $(date)"
echo "========================================"

# 设置 Gaussian 环境
export g16root=/public/software
export GAUSS_SCRDIR=/public/software/g16/scratch
source /public/software/g16/bsd/g16.profile

# 设置 Multiwfn 环境 - 确保使用正确的版本
# 首先从 PATH 中移除可能存在的旧 Multiwfn 路径
export PATH=$(echo $PATH | tr ':' '\\n' | grep -v -i multiwfn | tr '\\n' ':' | sed 's/:$//')
export Multiwfnpath=/public/software/Multiwfn_3.8_dev_bin_Linux
export PATH=$Multiwfnpath:$PATH
export OMP_STACKSIZE=200M
# 添加共享库路径 (libXm.so.4 for Multiwfn)
export LD_LIBRARY_PATH=/public/software/libs:$Multiwfnpath:/lib64:/usr/lib64:$LD_LIBRARY_PATH
ulimit -s unlimited

# 验证 Multiwfn 路径
echo "Using Multiwfn: $(which Multiwfn)"

# 运行 RESP2 计算
echo "Running RESP2.sh..."
{RESP2_SCRIPT_PATH} "{pdb_file}" {charge} {spin_multiplicity} {solvent}

RESP_EXIT_CODE=$?

if [ $RESP_EXIT_CODE -eq 0 ]; then
    echo "RESP calculation completed successfully"

    # 保存电荷文件到统一目录（供复用）
    CHG_FILE="{molecule_name}.charmm.chg"
    if [ -f "$CHG_FILE" ]; then
        mkdir -p "{self.charge_save_path}"
        cp "$CHG_FILE" "{charge_save_file}"
        echo "Charge file saved to: {charge_save_file}"
    fi

    echo "End Time: $(date)"
    exit 0
else
    echo "RESP calculation failed with exit code: $RESP_EXIT_CODE"
    echo "End Time: $(date)"
    exit 1
fi
"""

        with open(script_path, 'w') as f:
            f.write(script_content)

        os.chmod(script_path, 0o755)
        logger.info(f"Generated RESP Slurm script: {script_path}")

        return script_path

    def submit_resp_job(self, work_dir: Path, script_path: Path) -> Tuple[bool, Optional[str], Optional[str]]:
        """
        提交 RESP 任务到 Slurm

        Args:
            work_dir: 工作目录
            script_path: 脚本路径

        Returns:
            (成功与否, slurm_job_id, 错误信息)
        """
        try:
            result = subprocess.run(
                ['sbatch', str(script_path)],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                timeout=30
            )

            if result.returncode == 0:
                # 解析 job ID: "Submitted batch job 12345"
                match = re.search(r'Submitted batch job (\d+)', result.stdout)
                if match:
                    slurm_job_id = match.group(1)
                    logger.info(f"RESP job submitted successfully: {slurm_job_id}")
                    return True, slurm_job_id, None
                else:
                    logger.error(f"Could not parse slurm job ID from: {result.stdout}")
                    return False, None, "Could not parse slurm job ID"
            else:
                logger.error(f"sbatch failed: {result.stderr}")
                return False, None, result.stderr

        except subprocess.TimeoutExpired:
            logger.error("sbatch command timed out")
            return False, None, "sbatch timeout"
        except Exception as e:
            logger.error(f"Failed to submit RESP job: {e}")
            return False, None, str(e)

    def check_job_status(self, slurm_job_id: str) -> str:
        """
        检查 Slurm 任务状态

        Args:
            slurm_job_id: Slurm 任务 ID

        Returns:
            状态: PENDING, RUNNING, COMPLETED, FAILED, CANCELLED, UNKNOWN
        """
        try:
            result = subprocess.run(
                ['squeue', '-j', slurm_job_id, '-h', '-o', '%T'],
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode == 0 and result.stdout.strip():
                status = result.stdout.strip()
                return status
            else:
                # 任务不在队列中，检查 sacct
                result = subprocess.run(
                    ['sacct', '-j', slurm_job_id, '-n', '-o', 'State', '-X'],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                if result.returncode == 0 and result.stdout.strip():
                    status = result.stdout.strip().split()[0]
                    return status
                return "UNKNOWN"

        except Exception as e:
            logger.error(f"Failed to check job status: {e}")
            return "UNKNOWN"

    def get_job_cpu_hours(self, slurm_job_id: str) -> float:
        """
        获取 Slurm 任务的 CPU 核时数

        Args:
            slurm_job_id: Slurm 任务 ID

        Returns:
            CPU 核时数
        """
        try:
            # sacct -j JOB_ID -o CPUTimeRAW,NCPUS -n -X
            result = subprocess.run(
                ['sacct', '-j', slurm_job_id, '-o', 'CPUTimeRAW,NCPUS,Elapsed', '-n', '-X'],
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode == 0 and result.stdout.strip():
                parts = result.stdout.strip().split()
                if len(parts) >= 2:
                    cpu_time_seconds = int(parts[0])
                    cpu_hours = cpu_time_seconds / 3600.0
                    logger.info(f"Job {slurm_job_id}: {cpu_hours:.2f} CPU hours")
                    return cpu_hours

            return 0.0

        except Exception as e:
            logger.error(f"Failed to get CPU hours: {e}")
            return 0.0


def modify_lmp_charges(lmp_file: Path, charge_file: Path, output_file: Path) -> bool:
    """
    使用 RESP 电荷修改 LAMMPS 文件

    基于原版 /public/home/xiaoji/molyte_v1/command/modify_lmp_charges.py

    Args:
        lmp_file: 原始 LAMMPS 文件 (如 EC.lammps.lmp)
        charge_file: RESP 电荷文件 (如 EC.charmm.chg)
        output_file: 输出文件 (如 EC.lmp)

    Returns:
        成功返回 True，失败返回 False
    """
    try:
        # 读取 lmp 文件
        with open(lmp_file, 'r') as f:
            lmp_contents = f.readlines()

        # 读取电荷文件
        with open(charge_file, 'r') as f:
            chg_contents = f.readlines()

        # 提取电荷 (最后一列)
        charges = []
        for line in chg_contents:
            parts = line.split()
            if len(parts) >= 5:  # 元素 x y z 电荷
                try:
                    charge = float(parts[-1])
                    charges.append(charge)
                except ValueError:
                    continue

        if not charges:
            logger.error(f"No charges found in {charge_file}")
            return False

        # 找到 Atoms 部分
        atom_section_start = None
        for i, line in enumerate(lmp_contents):
            if line.strip() == "Atoms":
                atom_section_start = i + 2  # 跳过 "Atoms" 行和空行
                break

        if atom_section_start is None:
            logger.error("Atoms section not found in LAMMPS file")
            return False

        # 找到 Atoms 部分的结束位置
        atom_section_end = atom_section_start
        for i in range(atom_section_start, len(lmp_contents)):
            if lmp_contents[i].strip() == '':
                atom_section_end = i
                break
        else:
            atom_section_end = len(lmp_contents)

        # 检查原子数和电荷数是否匹配
        num_atoms = atom_section_end - atom_section_start
        if num_atoms != len(charges):
            logger.error(f"Atom count ({num_atoms}) != charge count ({len(charges)})")
            return False

        # 修改电荷 (第4列，索引3)
        modified_lines = []
        for idx, line in enumerate(lmp_contents[atom_section_start:atom_section_end]):
            parts = line.split()
            if len(parts) >= 4:
                # LAMMPS Atoms 格式: atom-ID mol-ID atom-type charge x y z
                parts[3] = f"{charges[idx]:.10f}"
                modified_line = '\t'.join(parts) + '\n'
                modified_lines.append(modified_line)
            else:
                modified_lines.append(line)

        # 替换 Atoms 部分
        lmp_contents[atom_section_start:atom_section_end] = modified_lines

        # 写入输出文件
        with open(output_file, 'w') as f:
            f.writelines(lmp_contents)

        logger.info(f"Modified LAMMPS file with RESP charges: {output_file}")
        return True

    except Exception as e:
        logger.error(f"Failed to modify LAMMPS charges: {e}")
        return False


def copy_charge_file_if_exists(molecule_name: str, work_dir: Path,
                               charge_save_path: Path = None) -> bool:
    """
    如果电荷文件已存在，复制到工作目录

    Args:
        molecule_name: 分子名称
        work_dir: 工作目录
        charge_save_path: 电荷文件保存目录

    Returns:
        是否成功复制
    """
    charge_save_path = charge_save_path or DEFAULT_CHARGE_SAVE_PATH
    source_file = charge_save_path / f"{molecule_name}.charmm.chg"
    target_file = work_dir / f"{molecule_name}.charmm.chg"

    if source_file.exists():
        try:
            shutil.copy(source_file, target_file)
            logger.info(f"Copied existing charge file: {source_file} -> {target_file}")
            return True
        except Exception as e:
            logger.error(f"Failed to copy charge file: {e}")
            return False
    return False

