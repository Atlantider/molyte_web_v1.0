"""
Molyte Wrapper

封装 molyte_command.py 的核心功能，接受 Python dict 作为输入
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional

from .atom_mapping_generator import AtomMappingGenerator

logger = logging.getLogger(__name__)


class MolyteWrapper:
    """
    Molyte 工作流封装类
    
    将原来基于 Excel 的工作流改为基于 Python dict 的工作流
    """
    
    def __init__(
        self,
        work_base_path: Path,
        initial_salts_path: Path,
        ligpargen_path: Path,
        packmol_path: Path,
        ltemplify_path: Path,
        moltemplate_path: Path,
        charge_save_path: Path
    ):
        """
        Args:
            work_base_path: 工作目录基础路径
            initial_salts_path: 初始盐文件路径
            ligpargen_path: LigParGen 可执行文件路径
            packmol_path: Packmol 可执行文件路径
            ltemplify_path: ltemplify.py 路径
            moltemplate_path: moltemplate.sh 路径
            charge_save_path: 电荷文件保存路径
        """
        self.work_base_path = Path(work_base_path)
        self.initial_salts_path = Path(initial_salts_path)
        self.ligpargen_path = Path(ligpargen_path)
        self.packmol_path = Path(packmol_path)
        self.ltemplify_path = Path(ltemplify_path)
        self.moltemplate_path = Path(moltemplate_path)
        self.charge_save_path = Path(charge_save_path)
        
        # 创建必要的目录
        self.work_base_path.mkdir(parents=True, exist_ok=True)
        self.charge_save_path.mkdir(parents=True, exist_ok=True)

        # 创建原子映射生成器
        self.atom_mapping_generator = AtomMappingGenerator(self.initial_salts_path)

    @staticmethod
    def strip_ion_charge(name: str) -> str:
        """去掉离子名称中的 + 和 - 符号"""
        if name:
            return name.replace('+', '').replace('-', '')
        return name

    def run_command(self, command: str) -> bool:
        """
        运行 shell 命令
        
        Args:
            command: 要执行的命令
        
        Returns:
            True 如果成功，False 如果失败
        """
        try:
            logger.info(f"Running command: {command}")
            result = subprocess.run(
                command,
                shell=True,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            output = result.stdout.decode().strip()
            if output:
                logger.info(output)
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Command failed: {e.stderr.decode().strip()}")
            return False
    
    def generate_lammps_input(
        self,
        job_data: Dict[str, Any],
        generate_atom_mapping: bool = True
    ) -> Dict[str, Any]:
        """
        生成 LAMMPS 输入文件
        
        Args:
            job_data: 任务数据（从 molyte_adapter 转换后的格式）
            generate_atom_mapping: 是否生成原子映射文件
        
        Returns:
            {
                "success": True/False,
                "work_dir": Path,
                "files": {
                    "data": "xxx.data",
                    "input": "xxx.in",
                    "job_script": "job.sh",
                    "atom_mapping": "atom_mapping.json"
                },
                "error": "错误信息（如果失败）"
            }
        """
        job_name = job_data["name"]
        logger.info(f"Generating LAMMPS input for {job_name}")

        # 创建工作目录
        work_dir = self.work_base_path / job_name
        if work_dir.exists():
            logger.info(f"Directory {work_dir} already exists. Removing it.")
            subprocess.run(f'rm -rf {work_dir}', shell=True)

        work_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Created work directory: {work_dir}")

        # 只复制需要的阳离子和阴离子文件
        if not self._copy_required_salt_files(job_data, work_dir):
            return {"success": False, "error": "Failed to copy required salt files"}

        # 切换到工作目录
        original_dir = os.getcwd()
        os.chdir(work_dir)
        
        try:
            # 1. 生成 Slurm 作业脚本
            if not self._generate_job_script(job_name, work_dir, job_data):
                return {"success": False, "error": "Failed to generate job script"}
            
            # 2. 处理溶剂（生成电荷文件）
            charge_method = job_data.get("charge_method", "ligpargen")
            if not self._process_solvents(job_data.get("solvents", []), charge_method):
                return {"success": False, "error": "Failed to process solvents"}
            
            # 3. 生成 Packmol 输入文件并运行
            if not self._generate_and_run_packmol(job_data, work_dir):
                return {"success": False, "error": "Failed to run Packmol"}

            # 4. 生成 Moltemplate LT 文件
            if not self._generate_moltemplate_files(job_data, work_dir):
                return {"success": False, "error": "Failed to generate Moltemplate files"}

            # 5. 运行 Moltemplate
            if not self._run_moltemplate(job_name, work_dir):
                return {"success": False, "error": "Failed to run Moltemplate"}

            # 6. 生成 LAMMPS 输入文件
            if not self._generate_lammps_input_file(job_data, work_dir):
                return {"success": False, "error": "Failed to generate LAMMPS input file"}

            # 7. 生成原子映射文件（如果需要）
            atom_mapping_file = None
            if generate_atom_mapping:
                try:
                    self.atom_mapping_generator.generate_atom_mapping(job_data, work_dir)
                    atom_mapping_file = "atom_mapping.json"
                    logger.info("Generated atom mapping file")
                except Exception as e:
                    logger.warning(f"Failed to generate atom mapping: {e}")

            # 返回成功结果
            return {
                "success": True,
                "work_dir": work_dir,
                "files": {
                    "data": f"{job_name}.data",
                    "input": f"{job_name}.in",
                    "job_script": "job.sh",
                    "atom_mapping": atom_mapping_file
                }
            }

        except Exception as e:
            logger.error(f"Error generating LAMMPS input: {e}")
            import traceback
            traceback.print_exc()
            return {"success": False, "error": str(e)}

        finally:
            # 切换回原始目录
            os.chdir(original_dir)

    def _copy_required_salt_files(self, job_data: Dict, work_dir: Path) -> bool:
        """
        只复制需要的阳离子和阴离子文件

        Args:
            job_data: 任务数据
            work_dir: 工作目录

        Returns:
            bool: 是否成功
        """
        import shutil

        # 收集需要的离子名称（原始名称和去掉+/-的名称）
        required_ions = {}  # {display_name: file_name}

        for cation in job_data.get("cations", []):
            if cation.get('name') and cation.get('number', 0) > 0:
                display_name = cation['name']
                # 去掉 + 和 - 符号用于查找文件
                file_name = display_name.replace('+', '').replace('-', '')
                required_ions[display_name] = file_name

        for anion in job_data.get("anions", []):
            if anion.get('name') and anion.get('number', 0) > 0:
                display_name = anion['name']
                # 去掉 + 和 - 符号用于查找文件
                file_name = display_name.replace('+', '').replace('-', '')
                required_ions[display_name] = file_name

        logger.info(f"Required ions: {required_ions}")

        # 复制每个离子的 .pdb 和 .lt 文件
        for display_name, file_name in required_ions.items():
            for ext in ['.pdb', '.lt']:
                # 优先查找带 +/- 的文件，找不到再查找不带的
                src_file = self.initial_salts_path / f"{display_name}{ext}"
                if not src_file.exists():
                    src_file = self.initial_salts_path / f"{file_name}{ext}"

                # 目标文件名统一使用去掉符号的名称（与 packmol/moltemplate 一致）
                dst_file = work_dir / f"{file_name}{ext}"

                if src_file.exists():
                    try:
                        shutil.copy2(src_file, dst_file)
                        logger.info(f"Copied {src_file.name} to {dst_file.name}")
                    except Exception as e:
                        logger.error(f"Failed to copy {src_file}: {e}")
                        return False
                else:
                    logger.warning(f"Salt file not found: {display_name}{ext} or {file_name}{ext}")
                    # 对于 .lt 文件，如果不存在也继续（可能是简单离子）
                    if ext == '.pdb':
                        return False

        # 复制 job.sh 模板（如果存在）
        job_sh_template = self.initial_salts_path / 'job.sh'
        if job_sh_template.exists():
            try:
                shutil.copy2(job_sh_template, work_dir / 'job.sh.template')
                logger.info("Copied job.sh template")
            except Exception as e:
                logger.warning(f"Failed to copy job.sh template: {e}")

        return True

    def _generate_job_script(self, job_name: str, work_dir: Path, job_data: Dict[str, Any]) -> bool:
        """生成 Slurm 作业脚本"""
        # 从 job_data 中获取 Slurm 资源配置，如果没有则使用默认值
        partition = job_data.get("slurm_partition", "cpu")
        nodes = job_data.get("slurm_nodes", 1)
        ntasks = job_data.get("slurm_ntasks", 8)
        cpus_per_task = job_data.get("slurm_cpus_per_task", 8)
        time_minutes = job_data.get("slurm_time", 7200)

        # 计算 MPI 进程数：ntasks * cpus_per_task
        mpi_processes = ntasks * cpus_per_task

        job_script_content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output=out.dat
#SBATCH --error=err.dat
#SBATCH --partition={partition}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ntasks}
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --time={time_minutes}

cd $SLURM_SUBMIT_DIR
export PATH=/public/software/lammps/mpich_install/bin:$PATH
EXEC=/public/software/lammps/mpich_install/bin/lmp_mpi
mpirun -n {mpi_processes} $EXEC <{job_name}.in > {job_name}.log
"""

        try:
            with open(work_dir / 'job.sh', 'w') as f:
                f.write(job_script_content)
            logger.info(f"Generated job.sh for {job_name} with resources: partition={partition}, nodes={nodes}, ntasks={ntasks}, cpus_per_task={cpus_per_task}, mpi_processes={mpi_processes}")
            return True
        except Exception as e:
            logger.error(f"Failed to generate job script: {e}")
            return False

    def _process_solvents(self, solvents: List[Dict], charge_method: str = "ligpargen") -> bool:
        """
        处理溶剂：生成电荷文件

        Args:
            solvents: 溶剂列表
            charge_method: 电荷计算方法 ("ligpargen" 或 "resp")
        """
        import shutil
        from .resp_calculator import modify_lmp_charges, copy_charge_file_if_exists

        for solvent in solvents:
            if not solvent.get('name'):
                continue

            name = solvent['name']
            smiles = solvent.get('smiles', '')

            if not smiles:
                logger.warning(f"No SMILES for solvent {name}, skipping")
                continue

            # 调用 ligpargen 生成初始文件（PDB、LMP等）
            ligpargen_cmd = f"{self.ligpargen_path}/ligpargen -s '{smiles}' -n {name} -r MOL -c 0 -o 0 -cgen CM1A"
            if not self.run_command(ligpargen_cmd):
                logger.error(f"Failed to run ligpargen for {name}")
                return False

            # 根据电荷方法处理
            if charge_method == "resp":
                # 使用 RESP 电荷
                logger.info(f"Using RESP charge method for {name}")

                # 检查是否有已保存的电荷文件
                charge_file = self.charge_save_path / f"{name}.charmm.chg"

                if charge_file.exists():
                    # 复用已有的电荷文件
                    logger.info(f"Reusing existing RESP charges from {charge_file}")
                    lammps_lmp = Path(f"{name}.lammps.lmp")
                    target_lmp = Path(f"{name}.lmp")

                    if lammps_lmp.exists():
                        if modify_lmp_charges(lammps_lmp, charge_file, target_lmp):
                            logger.info(f"Applied RESP charges to {target_lmp}")
                        else:
                            logger.error(f"Failed to apply RESP charges to {name}")
                            return False
                    else:
                        logger.error(f"LigParGen output file not found: {lammps_lmp}")
                        return False
                else:
                    # 没有保存的电荷文件，需要返回特殊标记让 worker 处理
                    # 这种情况应该在 worker 层面提前处理
                    logger.warning(f"No existing RESP charges for {name}, fallback to LigParGen")
                    charge_method = "ligpargen"

            if charge_method == "ligpargen":
                # 使用 LigParGen 生成的电荷
                logger.info(f"Using LigParGen CM1A charges for {name}")

                lammps_lmp = Path(f"{name}.lammps.lmp")
                target_lmp = Path(f"{name}.lmp")

                if lammps_lmp.exists():
                    try:
                        shutil.copy(lammps_lmp, target_lmp)
                        logger.info(f"Copied {lammps_lmp} to {target_lmp}")
                    except Exception as e:
                        logger.error(f"Failed to copy {lammps_lmp}: {e}")
                        return False
                else:
                    logger.error(f"LigParGen output file not found: {lammps_lmp}")
                    return False

        return True

    def get_solvents_needing_resp(self, solvents: List[Dict]) -> List[Dict]:
        """
        检查哪些溶剂需要 RESP 计算

        Args:
            solvents: 溶剂列表

        Returns:
            需要 RESP 计算的溶剂列表
        """
        needing_resp = []
        for solvent in solvents:
            name = solvent.get('name')
            if not name:
                continue

            charge_file = self.charge_save_path / f"{name}.charmm.chg"
            if not charge_file.exists():
                needing_resp.append(solvent)
                logger.info(f"Solvent {name} needs RESP calculation")
            else:
                logger.info(f"Solvent {name} has existing RESP charges")

        return needing_resp

    def _generate_and_run_packmol(self, job_data: Dict, work_dir: Path) -> bool:
        """生成 Packmol 输入文件并运行"""
        job_name = job_data["name"]
        box_size = job_data["box_size"]

        # 生成 Packmol 输入
        packmol_input = f"tolerance 2.0\nfiletype pdb\noutput {job_name}.pdb\n"

        # 添加阳离子
        for cation in job_data.get("cations", []):
            if cation.get('name') and cation.get('number', 0) > 0:
                name = self.strip_ion_charge(cation['name'])
                number = cation['number']
                packmol_input += f"""structure {name}.pdb
   number {number}
   inside box 0. 0. 0 {box_size} {box_size} {box_size}
end structure
"""

        # 添加阴离子
        for anion in job_data.get("anions", []):
            if anion.get('name') and anion.get('number', 0) > 0:
                name = self.strip_ion_charge(anion['name'])
                number = anion['number']
                packmol_input += f"""structure {name}.pdb
   number {number}
   inside box 0. 0. 0 {box_size} {box_size} {box_size}
end structure
"""

        # 添加溶剂
        for solvent in job_data.get("solvents", []):
            if solvent.get('name') and solvent.get('number', 0) > 0:
                name = solvent['name']
                number = solvent['number']
                packmol_input += f"""structure {name}.charmm.pdb
   number {number}
   inside box 0. 0. 0 {box_size} {box_size} {box_size}
end structure
"""

        # 写入文件
        packmol_inp_path = work_dir / f"{job_name}.inp"
        try:
            with open(packmol_inp_path, 'w') as f:
                f.write(packmol_input)
            logger.info(f"Generated Packmol input file for {job_name}")
        except Exception as e:
            logger.error(f"Failed to write Packmol input: {e}")
            return False

        # 运行 Packmol
        if not self.run_command(f"{self.packmol_path} < {packmol_inp_path}"):
            logger.error("Failed to run Packmol")
            return False

        return True

    def _generate_moltemplate_files(self, job_data: Dict, work_dir: Path) -> bool:
        """生成 Moltemplate LT 文件"""
        job_name = job_data["name"]
        box_size = job_data["box_size"]

        # 为溶剂生成 LT 文件
        for solvent in job_data.get("solvents", []):
            if solvent.get('name'):
                name = solvent['name']
                lmp_file = work_dir / f"{name}.lmp"
                lt_file = work_dir / f"{name}.lt"

                if lmp_file.exists():
                    cmd = f"{self.ltemplify_path} -name {name} {lmp_file} > {lt_file}"
                    if not self.run_command(cmd):
                        logger.warning(f"Failed to generate LT file for {name}")

        # 生成主 LT 文件
        moltemplate_content = ""

        # 导入所有分子的 LT 文件
        for cation in job_data.get("cations", []):
            if cation.get('name') and cation.get('number', 0) > 0:
                name = self.strip_ion_charge(cation["name"])
                moltemplate_content += f'import "{name}.lt"\n'

        for anion in job_data.get("anions", []):
            if anion.get('name') and anion.get('number', 0) > 0:
                name = self.strip_ion_charge(anion["name"])
                moltemplate_content += f'import "{name}.lt"\n'

        for solvent in job_data.get("solvents", []):
            if solvent.get('name') and solvent.get('number', 0) > 0:
                moltemplate_content += f'import "{solvent["name"]}.lt"\n'

        # 实例化分子
        for cation in job_data.get("cations", []):
            if cation.get('name') and cation.get('number', 0) > 0:
                name = self.strip_ion_charge(cation['name'])
                number = cation['number']
                moltemplate_content += f'{name}s = new {name}[{number}]\n'

        for anion in job_data.get("anions", []):
            if anion.get('name') and anion.get('number', 0) > 0:
                name = self.strip_ion_charge(anion['name'])
                number = anion['number']
                moltemplate_content += f'{name}s = new {name}[{number}]\n'

        for solvent in job_data.get("solvents", []):
            if solvent.get('name') and solvent.get('number', 0) > 0:
                name = solvent['name']
                number = solvent['number']
                moltemplate_content += f'{name}s = new {name}[{number}]\n'

        # 添加边界条件
        moltemplate_content += f"""
write_once("Data Boundary") {{
    0 {box_size} xlo xhi
    0 {box_size} ylo yhi
    0 {box_size} zlo zhi
}}
"""

        # 写入文件
        moltemplate_lt_path = work_dir / f"{job_name}.lt"
        try:
            with open(moltemplate_lt_path, 'w') as f:
                f.write(moltemplate_content)
            logger.info(f"Generated Moltemplate LT file for {job_name}")
            return True
        except Exception as e:
            logger.error(f"Failed to write Moltemplate LT file: {e}")
            return False

    def _run_moltemplate(self, job_name: str, work_dir: Path) -> bool:
        """运行 Moltemplate"""
        pdb_file = work_dir / f"{job_name}.pdb"
        lt_file = work_dir / f"{job_name}.lt"

        if not pdb_file.exists():
            logger.error(f"PDB file not found: {pdb_file}")
            return False

        if not lt_file.exists():
            logger.error(f"LT file not found: {lt_file}")
            return False

        cmd = f"{self.moltemplate_path} -pdb {pdb_file} {lt_file}"
        if not self.run_command(cmd):
            logger.error("Failed to run Moltemplate")
            return False

        # 检查生成的 data 文件
        data_file = work_dir / f"{job_name}.data"
        if not data_file.exists():
            logger.error(f"LAMMPS data file not generated: {data_file}")
            return False

        return True

    def _generate_lammps_input_file(self, job_data: Dict, work_dir: Path) -> bool:
        """生成 LAMMPS 输入文件（.in）"""
        job_name = job_data["name"]

        # 解析 LAMMPS data 文件获取元素列表
        element_list = self._parse_lammps_data_for_elements(work_dir / f"{job_name}.data")

        # 生成 .in.list 文件
        if not self._generate_in_list_file(job_data, work_dir, element_list):
            return False

        # 生成主输入文件
        lmp_in_content = self._build_lammps_input_content(job_data, element_list)

        # 写入文件
        try:
            with open(work_dir / f"{job_name}.in", 'w') as f:
                f.write(lmp_in_content)
            logger.info(f"Generated LAMMPS .in file for {job_name}")
            return True
        except Exception as e:
            logger.error(f"Failed to write LAMMPS input file: {e}")
            return False

    def _parse_lammps_data_for_elements(self, data_file: Path) -> List[str]:
        """解析 LAMMPS data 文件获取元素列表"""
        element_list = []

        try:
            with open(data_file, 'r', encoding='utf-8') as f:
                lines = f.readlines()

            # 查找 Masses 部分
            in_masses = False
            for i, line in enumerate(lines):
                if line.strip().lower() == 'masses':
                    in_masses = True
                    continue

                if in_masses:
                    # 跳过空行
                    if not line.strip():
                        continue

                    # 遇到 Atoms 部分则停止
                    if 'atoms' in line.lower():
                        break

                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            mass = round(float(parts[1]))
                            element = self._mass_to_element(mass)
                            element_list.append(element)
                        except ValueError:
                            continue

            logger.info(f"Parsed {len(element_list)} elements from data file")
            return element_list

        except Exception as e:
            logger.error(f"Failed to parse LAMMPS data file: {e}")
            return []

    def _mass_to_element(self, mass: int) -> str:
        """根据原子质量返回元素符号"""
        mass_to_element_map = {
            1: 'H', 4: 'He', 7: 'Li', 9: 'Be', 11: 'B', 12: 'C',
            14: 'N', 16: 'O', 19: 'F', 20: 'Ne', 23: 'Na', 24: 'Mg',
            27: 'Al', 28: 'Si', 31: 'P', 32: 'S', 35: 'Cl', 39: 'K',
            40: 'Ca', 45: 'Sc', 48: 'Ti', 51: 'V', 52: 'Cr', 55: 'Mn',
            56: 'Fe', 59: 'Co', 64: 'Cu', 65: 'Zn', 70: 'Ga', 73: 'Ge',
            75: 'As', 79: 'Se'
        }
        return mass_to_element_map.get(mass, 'A')

    def _generate_in_list_file(
        self,
        job_data: Dict,
        work_dir: Path,
        element_list: List[str]
    ) -> bool:
        """生成 .in.list 文件"""
        job_name = job_data["name"]

        # 构建 RDF 对列表（使用原子类型索引，不是元素符号）
        # LAMMPS 的 compute rdf 命令需要原子类型索引（从1开始）
        rdf_pair_indices = []

        # 定义阳离子和配位原子
        cation_elements = ['Li', 'Na', 'K', 'Mg', 'Ca', 'Zn', 'Al']
        coord_elements = ['O', 'N', 'P', 'B', 'F', 'S', 'Cl']

        # 遍历元素列表，生成 RDF 对的索引
        for cation_idx, cation in enumerate(element_list, start=1):
            if cation in cation_elements:
                # 对每个阳离子，找到所有配位原子
                for coord_idx, coord in enumerate(element_list, start=1):
                    if coord in coord_elements:
                        # 添加索引对（例如：1 4 表示类型1和类型4之间的RDF）
                        rdf_pair_indices.append(f"{cation_idx} {coord_idx}")

        # 生成内容
        in_list_content = f'variable element_list index "{" ".join(element_list)}"\n'
        in_list_content += f'variable rdf_pair string "{" ".join(rdf_pair_indices)}"\n'

        # 写入文件
        try:
            with open(work_dir / f"{job_name}.in.list", 'w') as f:
                f.write(in_list_content)
            logger.info(f"Generated .in.list file for {job_name}")
            logger.info(f"  Element list: {element_list}")
            logger.info(f"  RDF pairs (indices): {rdf_pair_indices}")
            return True
        except Exception as e:
            logger.error(f"Failed to write .in.list file: {e}")
            return False

    def _build_lammps_input_content(self, job_data: Dict, element_list: List[str]) -> str:
        """构建 LAMMPS 输入文件内容"""
        job_name = job_data["name"]

        # 变量定义部分
        content = "# ----------------- Variable Section -----------------\n"

        # 定义阳离子变量
        for i, cation in enumerate(job_data.get("cations", []), start=1):
            if cation.get('name'):
                content += f'variable Cation{i} index {cation["name"]}\n'

        # 定义阴离子变量
        for i, anion in enumerate(job_data.get("anions", []), start=1):
            if anion.get('name'):
                content += f'variable Anion{i} index {anion["name"]}\n'

        # 定义溶剂变量
        for i, solvent in enumerate(job_data.get("solvents", []), start=1):
            if solvent.get('name'):
                content += f'variable Solvent{i} index {solvent["name"]}\n'

        # 其他变量
        content += f"""
variable infile string {job_name}
variable outname string {job_name}

include {job_name}.in.init
read_data {job_name}.data
include {job_name}.in.settings
include {job_name}.in.list
include {job_name}.in.list_salt

variable Nsteps_NPT equal {job_data.get("nsteps_npt", 5000000)}
variable Nsteps_NVT equal {job_data.get("nsteps_nvt", 10000000)}

variable Freq_trj_npt index {job_data.get("freq_trj_npt", 1000000)}
variable Freq_trj_nvt index {job_data.get("freq_trj_nvt", 100000)}

variable Temp_NPT equal {job_data.get("temperature_npt", 298.15)}
variable Temp_NVT equal {job_data.get("temperature_nvt", 298.15)}

variable thermo_freq equal {job_data.get("thermo_freq", 1000000)}
variable timestep equal {job_data.get("timestep", 1.0)}
"""

        # 定义 group
        for i, cation in enumerate(job_data.get("cations", []), start=1):
            if cation.get('name') and cation.get('number', 0) > 0:
                content += f'group Cation{i} union {cation["name"]}\n'

        for i, anion in enumerate(job_data.get("anions", []), start=1):
            if anion.get('name') and anion.get('number', 0) > 0:
                content += f'group Anion{i} union {anion["name"]}\n'

        # 模拟设置
        content += """
thermo_style custom step cpu cpuremain temp density lx ly lz etotal ke pe evdwl ecoul elong ebond eangle edihed eimp
thermo ${thermo_freq}
timestep ${timestep}

minimize 1.0e-4 1.0e-6 5000 10000

write_data ${infile}_after_minimize.data nocoeff
write_dump all custom ${infile}_after_minimize.lammpstrj id element mol type x y z q modify element ${element_list} sort id

reset_timestep 0

dump trj_npt all custom ${Freq_trj_npt} NPT_${outname}.lammpstrj id element mol type x y z q
dump_modify trj_npt flush yes element ${element_list} sort id

fix fxnpt all npt temp ${Temp_NPT} ${Temp_NPT} $(100.0*dt) iso 0.0 0.0 $(1000.0*dt)
run ${Nsteps_NPT}
unfix fxnpt

undump trj_npt

write_data ${infile}_after_npt.data
write_restart ${infile}_restart_after_npt.data
write_dump all custom ${infile}_after_npt.lammpstrj id element mol type x y z q modify element ${element_list} sort id

reset_timestep 0
"""

        # MSD 计算（阴离子）
        for i, anion in enumerate(job_data.get("anions", []), start=1):
            if anion.get('name') and anion.get('number', 0) > 0:
                name = self.strip_ion_charge(anion['name'])
                content += f'compute An{i} {name} msd com yes\n'
                content += f'fix An{i}msd {name} ave/time ${{Freq_trj_nvt}} 1 ${{Freq_trj_nvt}} c_An{i}[1] c_An{i}[2] c_An{i}[3] c_An{i}[4] file out_{name}_msd.dat title1 "t msd msd msd msd_{name}" title2 "fs {name}_x {name}_y {name}_z {name}_total"\n'

        # MSD 计算（阳离子）
        for i, cation in enumerate(job_data.get("cations", []), start=1):
            if cation.get('name') and cation.get('number', 0) > 0:
                name = self.strip_ion_charge(cation['name'])
                content += f'compute Ca{i} {name} msd com yes\n'
                content += f'fix Ca{i}msd {name} ave/time ${{Freq_trj_nvt}} 1 ${{Freq_trj_nvt}} c_Ca{i}[1] c_Ca{i}[2] c_Ca{i}[3] c_Ca{i}[4] file out_{name}_msd.dat title1 "t msd msd msd msd_{name}" title2 "fs {name}_x {name}_y {name}_z {name}_total"\n'

        # RDF 计算和 NVT 模拟
        content += """
compute rdfc1 all rdf 100 ${rdf_pair}
fix rdff1 all ave/time $(v_Nsteps_NVT/1000) 1000 ${Nsteps_NVT} c_rdfc1[*] file out_rdf.dat mode vector

dump trj_nvt all custom ${Freq_trj_nvt} NVT_${outname}.lammpstrj id element mol type x y z q
dump_modify trj_nvt flush yes element ${element_list} sort id
dump utrj_nvt all custom ${Freq_trj_nvt} NVT_${outname}_un.lammpstrj id element mol type xu yu zu ix iy iz q
dump_modify utrj_nvt flush yes element ${element_list} sort id

fix fxnvt all nvt temp ${Temp_NPT} ${Temp_NPT} $(100.0*dt)
run ${Nsteps_NVT}
unfix fxnvt

undump trj_nvt
undump utrj_nvt

write_data ${infile}_after_nvt.data
write_restart ${infile}_restart_after_nvt.data
write_dump all custom ${infile}_after_nvt.lammpstrj id element mol type x y z q modify element ${element_list} sort id
"""

        return content

    def submit_to_slurm(self, work_dir: Path) -> Dict[str, Any]:
        """
        提交任务到 Slurm 队列

        Args:
            work_dir: 工作目录，包含 job.sh 脚本

        Returns:
            Dict with 'success', 'slurm_job_id' or 'error'
        """
        import subprocess
        import re

        job_script = work_dir / "job.sh"

        if not job_script.exists():
            return {
                'success': False,
                'error': f'Job script not found: {job_script}'
            }

        try:
            # 切换到工作目录并提交任务
            result = subprocess.run(
                ['sbatch', str(job_script)],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                timeout=60
            )

            if result.returncode != 0:
                return {
                    'success': False,
                    'error': f'sbatch failed: {result.stderr}'
                }

            # 解析 Slurm job ID，格式: "Submitted batch job 12345"
            output = result.stdout.strip()
            match = re.search(r'Submitted batch job (\d+)', output)

            if match:
                slurm_job_id = match.group(1)
                logger.info(f"Successfully submitted to Slurm: job_id={slurm_job_id}")
                return {
                    'success': True,
                    'slurm_job_id': slurm_job_id
                }
            else:
                return {
                    'success': False,
                    'error': f'Could not parse Slurm job ID from: {output}'
                }

        except subprocess.TimeoutExpired:
            return {
                'success': False,
                'error': 'sbatch command timed out'
            }
        except Exception as e:
            logger.error(f"Error submitting to Slurm: {e}")
            return {
                'success': False,
                'error': str(e)
            }

