#!/usr/bin/env python3
"""
混合云轮询 Worker

功能：
1. 定期轮询阿里云 API，获取待处理任务
2. 下载任务输入数据
3. 生成 LAMMPS/Gaussian 输入文件
4. 提交到 Slurm 集群
5. 监控任务状态
6. 上传结果到阿里云 OSS
7. 更新任务状态
"""

import os
import sys
import time
import yaml
import logging
import requests
import subprocess
import json
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent / "backend"))

try:
    from qcloud_cos import CosConfig, CosS3Client
except ImportError:
    print("请安装腾讯云 COS SDK: pip install cos-python-sdk-v5")
    sys.exit(1)


class PollingWorker:
    """轮询 Worker 主类"""
    
    def __init__(self, config_path: str = "deployment/polling_worker_config.yaml"):
        """初始化 Worker"""
        # 加载配置
        self.config = self._load_config(config_path)
        
        # 设置日志
        self._setup_logging()
        
        # 初始化 OSS 客户端
        self._init_oss_client()
        
        # 初始化 API 客户端
        self._init_api_client()
        
        # 当前运行的任务
        self.running_jobs: Dict[int, Dict] = {}
        
        self.logger.info(f"Worker '{self.config['worker']['name']}' 已启动")
    
    def _load_config(self, config_path: str) -> Dict:
        """加载配置文件"""
        config_file = Path(config_path)
        if not config_file.exists():
            raise FileNotFoundError(f"配置文件不存在: {config_path}")
        
        with open(config_file, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    
    def _setup_logging(self):
        """设置日志"""
        log_file = self.config['worker']['log_file']
        log_level = getattr(logging, self.config['worker']['log_level'])
        
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s | %(levelname)s | %(name)s | %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger('PollingWorker')
    
    def _init_oss_client(self):
        """初始化 COS 客户端（腾讯云对象存储）"""
        # 支持两种配置：阿里云 OSS 或腾讯云 COS
        if 'cos' in self.config:
            # 腾讯云 COS
            cos_config = self.config['cos']
            config = CosConfig(
                Region=cos_config['region'],
                SecretId=cos_config['secret_id'],
                SecretKey=cos_config['secret_key'],
                Scheme='https'
            )
            self.cos_client = CosS3Client(config)
            self.cos_bucket = cos_config['bucket']
            self.storage_type = 'cos'
            self.logger.info(f"腾讯云 COS 客户端已初始化 (Bucket: {self.cos_bucket})")
        elif 'oss' in self.config:
            # 阿里云 OSS（向后兼容）
            import oss2
            oss_config = self.config['oss']
            auth = oss2.Auth(
                oss_config['access_key_id'],
                oss_config['access_key_secret']
            )
            self.oss_bucket = oss2.Bucket(
                auth,
                oss_config['endpoint'],
                oss_config['bucket_name']
            )
            self.storage_type = 'oss'
            self.logger.info("阿里云 OSS 客户端已初始化")
        else:
            raise ValueError("配置文件中必须包含 'cos' 或 'oss' 配置")
    
    def _init_api_client(self):
        """初始化 API 客户端"""
        self.api_base_url = self.config['api']['base_url']
        self.api_headers = {
            'Authorization': f'Bearer {self.config["api"]["worker_token"]}',
            'Content-Type': 'application/json'
        }
        self.logger.info(f"API 客户端已初始化: {self.api_base_url}")
    
    def run(self):
        """主循环"""
        self.logger.info("开始轮询...")
        
        last_heartbeat = time.time()
        
        while True:
            try:
                # 发送心跳
                if time.time() - last_heartbeat > self.config['worker']['heartbeat_interval']:
                    self._send_heartbeat()
                    last_heartbeat = time.time()
                
                # 检查运行中的任务
                self._check_running_jobs()
                
                # 获取新任务
                if len(self.running_jobs) < self.config['worker']['max_concurrent_jobs']:
                    self._fetch_and_process_new_jobs()
                
                # 等待下一次轮询
                time.sleep(self.config['api']['poll_interval'])
                
            except KeyboardInterrupt:
                self.logger.info("收到中断信号，正在退出...")
                break
            except Exception as e:
                self.logger.error(f"主循环错误: {e}", exc_info=True)
                time.sleep(10)
    
    def _get_slurm_partitions(self) -> List[Dict]:
        """获取 Slurm 分区信息"""
        partitions = []
        try:
            # 使用 sinfo 获取分区信息
            cmd = ["sinfo", "--format=%P|%a|%D|%C|%l", "--noheader"]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)

            if result.returncode != 0:
                self.logger.warning(f"sinfo 命令失败: {result.stderr}")
                return partitions

            for line in result.stdout.strip().split("\n"):
                if not line:
                    continue

                fields = line.split("|")
                if len(fields) < 4:
                    continue

                name = fields[0].rstrip("*")  # 移除默认分区的 * 标记
                state = fields[1]
                total_nodes = int(fields[2]) if fields[2].isdigit() else 0

                # 解析 CPU 信息 (格式: A/I/O/T - Allocated/Idle/Other/Total)
                cpu_info = fields[3].split("/")
                if len(cpu_info) >= 4:
                    allocated = int(cpu_info[0]) if cpu_info[0].isdigit() else 0
                    idle = int(cpu_info[1]) if cpu_info[1].isdigit() else 0
                    total = int(cpu_info[3]) if cpu_info[3].isdigit() else 0
                else:
                    allocated = idle = total = 0

                max_time = fields[4] if len(fields) > 4 else None

                partitions.append({
                    "name": name,
                    "state": state,
                    "total_nodes": total_nodes,
                    "available_nodes": total_nodes,  # 简化处理
                    "total_cpus": total,
                    "available_cpus": idle,
                    "max_time": max_time,
                })

            self.logger.debug(f"获取到 {len(partitions)} 个分区信息")

        except subprocess.TimeoutExpired:
            self.logger.error("sinfo 命令超时")
        except FileNotFoundError:
            self.logger.warning("sinfo 命令未找到（可能不在 Slurm 环境中）")
        except Exception as e:
            self.logger.error(f"获取分区信息失败: {e}")

        return partitions

    def _send_heartbeat(self):
        """发送心跳到云端（同时上报分区信息）"""
        try:
            # 获取当前分区信息
            partitions = self._get_slurm_partitions()

            data = {
                'worker_name': self.config['worker']['name'],
                'status': 'running',
                'running_jobs': len(self.running_jobs),
                'timestamp': datetime.now().isoformat(),
                'partitions': partitions  # 上报分区信息
            }

            response = requests.post(
                f"{self.api_base_url}/workers/heartbeat",
                headers=self.api_headers,
                json=data,
                timeout=self.config['api']['timeout']
            )

            if response.status_code == 200:
                self.logger.debug(f"心跳已发送: {len(self.running_jobs)} 个任务运行中, {len(partitions)} 个分区")
            else:
                self.logger.warning(f"发送心跳失败: {response.status_code} - {response.text}")

        except Exception as e:
            self.logger.warning(f"发送心跳失败: {e}")
    
    def _fetch_and_process_new_jobs(self):
        """获取并处理新任务"""
        try:
            # 获取待处理的 MD 任务
            md_jobs = self._fetch_pending_jobs('md')
            for job in md_jobs:
                if len(self.running_jobs) >= self.config['worker']['max_concurrent_jobs']:
                    break
                self._process_md_job(job)
            
            # 获取待处理的 QC 任务
            qc_jobs = self._fetch_pending_jobs('qc')
            for job in qc_jobs:
                if len(self.running_jobs) >= self.config['worker']['max_concurrent_jobs']:
                    break
                self._process_qc_job(job)
                
        except Exception as e:
            self.logger.error(f"获取新任务失败: {e}", exc_info=True)
    
    def _fetch_pending_jobs(self, job_type: str) -> List[Dict]:
        """从阿里云获取待处理任务"""
        try:
            endpoint = f"{self.api_base_url}/workers/jobs/pending"
            params = {
                'job_type': job_type.upper(),
                'limit': 10
            }
            
            response = requests.get(
                endpoint,
                headers=self.api_headers,
                params=params,
                timeout=self.config['api']['timeout']
            )
            
            if response.status_code == 200:
                jobs = response.json()
                self.logger.info(f"获取到 {len(jobs)} 个待处理的 {job_type.upper()} 任务")
                return jobs
            else:
                self.logger.warning(f"获取任务失败: {response.status_code}")
                return []
                
        except Exception as e:
            self.logger.error(f"获取 {job_type} 任务失败: {e}")
            return []
    
    def _process_md_job(self, job: Dict):
        """处理 MD 任务"""
        job_id = job['id']
        self.logger.info(f"开始处理 MD 任务 {job_id}")

        try:
            # 1. 立即更新任务状态为 QUEUED（表示 Worker 已接收，正在准备）
            self._update_job_status(job_id, 'QUEUED', 'md')
            
            # 2. 导入 MolyteWrapper
            from app.workers.molyte_wrapper import MolyteWrapper
            from app.core.config import settings
            
            # 3. 初始化 MolyteWrapper
            wrapper = MolyteWrapper(
                work_base_path=Path(self.config['local']['work_base_path']),
                initial_salts_path=Path(self.config['local']['initial_salts_path']),
                ligpargen_path=Path(self.config['local']['ligpargen_path']),
                packmol_path=Path(self.config['local']['packmol_path']),
                ltemplify_path=Path(self.config['local']['ltemplify_path']),
                moltemplate_path=Path(self.config['local']['moltemplate_path']),
                charge_save_path=Path(self.config['local']['charge_save_path']),
            )
            
            # 4. 生成 LAMMPS 输入文件
            job_data = job['config']
            result = wrapper.generate_lammps_input(job_data)
            
            if not result['success']:
                raise Exception(result.get('error', 'Unknown error'))
            
            work_dir = result['work_dir']
            
            # 5. 提交到 Slurm
            slurm_result = wrapper.submit_to_slurm(work_dir)
            
            if not slurm_result['success']:
                raise Exception(slurm_result.get('error', 'Failed to submit to Slurm'))
            
            slurm_job_id = slurm_result['slurm_job_id']
            
            # 6. 更新任务状态为 RUNNING
            self._update_job_status(
                job_id, 'RUNNING', 'md',
                slurm_job_id=slurm_job_id,
                work_dir=str(work_dir)
            )
            
            # 7. 添加到运行中的任务列表
            self.running_jobs[job_id] = {
                'type': 'md',
                'slurm_job_id': slurm_job_id,
                'work_dir': work_dir,
                'start_time': time.time()
            }
            
            self.logger.info(f"MD 任务 {job_id} 已提交到 Slurm (Job ID: {slurm_job_id})")
            
        except Exception as e:
            self.logger.error(f"处理 MD 任务 {job_id} 失败: {e}", exc_info=True)
            self._update_job_status(job_id, 'FAILED', 'md', error_message=str(e))
    
    def _process_qc_job(self, job: Dict):
        """处理 QC 任务"""
        job_id = job['id']
        self.logger.info(f"开始处理 QC 任务 {job_id}")

        try:
            # 1. 立即更新任务状态为 QUEUED（表示 Worker 已接收，正在准备）
            self._update_job_status(job_id, 'QUEUED', 'qc')

            # 2. 获取任务配置
            config = job.get('config', {})
            molecule_name = config.get('molecule_name', f'QC_{job_id}')
            smiles = config.get('smiles', '')
            basis_set = config.get('basis_set', '6-31++g(d,p)')
            functional = config.get('functional', 'B3LYP')
            charge = config.get('charge', 0)
            spin_multiplicity = config.get('spin_multiplicity', 1)
            solvent_model = config.get('solvent_model', 'gas')
            solvent_name = config.get('solvent_name', '')
            slurm_partition = config.get('slurm_partition', 'cpu')
            slurm_cpus = config.get('slurm_cpus', 16)
            slurm_time = config.get('slurm_time', 7200)

            # 3. 创建工作目录
            work_base = Path(self.config['local']['work_base_path'])
            work_dir = work_base / "qc_work" / f"QC-{job_id}-{molecule_name}"
            work_dir.mkdir(parents=True, exist_ok=True)

            # 4. 生成安全的文件名
            safe_name = self._sanitize_filename(molecule_name)

            # 5. 生成 Gaussian 输入文件
            gjf_path = work_dir / f"{safe_name}.gjf"
            self._generate_gaussian_input(
                gjf_path, molecule_name, smiles,
                charge, spin_multiplicity,
                functional, basis_set,
                solvent_model, solvent_name
            )
            self.logger.info(f"生成 Gaussian 输入文件: {gjf_path}")

            # 6. 生成 Slurm 作业脚本
            job_script = work_dir / "job.sh"
            self._generate_qc_job_script(
                job_script, safe_name, slurm_partition, slurm_cpus, slurm_time
            )
            self.logger.info(f"生成 Slurm 作业脚本: {job_script}")

            # 7. 提交到 Slurm
            slurm_result = self._submit_to_slurm(work_dir)

            if not slurm_result['success']:
                raise Exception(slurm_result.get('error', 'Failed to submit to Slurm'))

            slurm_job_id = slurm_result['slurm_job_id']

            # 8. 更新任务状态为 RUNNING
            self._update_job_status(
                job_id, 'RUNNING', 'qc',
                slurm_job_id=slurm_job_id,
                work_dir=str(work_dir)
            )

            # 9. 添加到运行中的任务列表
            self.running_jobs[job_id] = {
                'type': 'qc',
                'slurm_job_id': slurm_job_id,
                'work_dir': work_dir,
                'start_time': time.time()
            }

            self.logger.info(f"QC 任务 {job_id} 已提交到 Slurm (Job ID: {slurm_job_id})")

        except Exception as e:
            self.logger.error(f"处理 QC 任务 {job_id} 失败: {e}", exc_info=True)
            self._update_job_status(job_id, 'FAILED', 'qc', error_message=str(e))

    def _sanitize_filename(self, name: str) -> str:
        """清理文件名，去除不安全字符"""
        import re
        # 替换不安全字符
        safe = re.sub(r'[^\w\-.]', '_', name)
        # 去除连续的下划线
        safe = re.sub(r'_+', '_', safe)
        # 去除首尾的下划线
        safe = safe.strip('_')
        return safe or 'molecule'

    def _generate_gaussian_input(self, gjf_path: Path, molecule_name: str, smiles: str,
                                   charge: int, spin_multiplicity: int,
                                   functional: str, basis_set: str,
                                   solvent_model: str, solvent_name: str):
        """生成 Gaussian 输入文件"""
        # 构建计算关键字
        keywords = f"opt freq {functional}/{basis_set}"

        # 添加色散校正（对于DFT方法）
        if functional.upper() not in ["HF"]:
            keywords += " em=gd3bj"

        # 添加溶剂效应
        if solvent_model and solvent_model.lower() != 'gas':
            if solvent_model.lower() == 'pcm':
                keywords += f" scrf=(pcm,solvent={solvent_name or 'water'})"
            elif solvent_model.lower() == 'smd':
                keywords += f" scrf=(smd,solvent={solvent_name or 'water'})"

        safe_name = self._sanitize_filename(molecule_name)

        # 生成 gjf 内容
        gjf_content = f"""%nprocshared=16
%mem=8GB
%chk={safe_name}.chk
# {keywords}

{molecule_name}

{charge} {spin_multiplicity}
"""

        # 尝试从 SMILES 生成 3D 坐标
        coords = self._get_3d_coordinates(smiles, molecule_name)
        if coords:
            for atom, x, y, z in coords:
                gjf_content += f" {atom:<2}  {x:>12.8f}  {y:>12.8f}  {z:>12.8f}\n"
        else:
            # 如果无法生成坐标，使用 SMILES 作为注释
            gjf_content += f"! SMILES: {smiles}\n"
            gjf_content += "! 请手动添加分子坐标\n"

        gjf_content += "\n"  # 空行结尾

        with open(gjf_path, 'w') as f:
            f.write(gjf_content)

    def _get_3d_coordinates(self, smiles: str, molecule_name: str):
        """从 SMILES 或 PDB 文件获取 3D 坐标"""
        # 首先尝试从 initial_salts 目录加载 PDB 文件
        initial_salts_path = Path(self.config['local']['initial_salts_path'])
        clean_name = molecule_name.replace("+", "").replace("-", "").strip()

        possible_paths = [
            initial_salts_path / f"{clean_name}.pdb",
            initial_salts_path / f"{molecule_name}.pdb",
        ]

        for pdb_path in possible_paths:
            if pdb_path.exists():
                coords = self._parse_pdb_coordinates(pdb_path)
                if coords:
                    self.logger.info(f"从 PDB 文件加载坐标: {pdb_path}")
                    return coords

        # 尝试使用 RDKit 从 SMILES 生成坐标
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem

            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(mol)

                coords = []
                conf = mol.GetConformer()
                for i, atom in enumerate(mol.GetAtoms()):
                    pos = conf.GetAtomPosition(i)
                    coords.append((atom.GetSymbol(), pos.x, pos.y, pos.z))

                self.logger.info(f"从 SMILES 生成 3D 坐标")
                return coords
        except ImportError:
            self.logger.warning("RDKit 未安装，无法从 SMILES 生成坐标")
        except Exception as e:
            self.logger.warning(f"从 SMILES 生成坐标失败: {e}")

        return None

    def _parse_pdb_coordinates(self, pdb_path: Path):
        """解析 PDB 文件坐标"""
        coords = []
        try:
            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        # PDB 格式: ATOM/HETATM, atom_num, atom_name, ..., x, y, z
                        atom_name = line[12:16].strip()
                        # 提取元素符号（去掉数字）
                        element = ''.join(c for c in atom_name if c.isalpha())
                        if not element:
                            element = atom_name[0]
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append((element, x, y, z))
            return coords if coords else None
        except Exception as e:
            self.logger.warning(f"解析 PDB 文件失败: {e}")
            return None

    def _generate_qc_job_script(self, script_path: Path, safe_name: str,
                                  partition: str, cpus: int, time_limit: int):
        """生成 QC 任务的 Slurm 作业脚本"""
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

# 设置 Gaussian 环境
export g16root=/public/software
export GAUSS_SCRDIR=/public/software/g16/scratch
source /public/software/g16/bsd/g16.profile

# 运行 Gaussian
ulimit -s unlimited
g16 < "{safe_name}.gjf" > "{safe_name}_out.log" 2>&1

# 转换 checkpoint 文件
if [ -f "{safe_name}.chk" ]; then
    formchk "{safe_name}.chk" "{safe_name}.fchk"
fi

echo "QC calculation completed"
"""

        with open(script_path, 'w') as f:
            f.write(script_content)

        import os
        os.chmod(script_path, 0o755)

    def _submit_to_slurm(self, work_dir: Path) -> Dict[str, Any]:
        """提交任务到 Slurm"""
        import re

        job_script = work_dir / "job.sh"

        if not job_script.exists():
            return {'success': False, 'error': f'Job script not found: {job_script}'}

        try:
            result = subprocess.run(
                ['sbatch', str(job_script)],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                timeout=60
            )

            if result.returncode != 0:
                return {'success': False, 'error': f'sbatch failed: {result.stderr}'}

            # 解析 Slurm job ID
            output = result.stdout.strip()
            match = re.search(r'Submitted batch job (\d+)', output)

            if match:
                slurm_job_id = match.group(1)
                return {'success': True, 'slurm_job_id': slurm_job_id}
            else:
                return {'success': False, 'error': f'Could not parse Slurm job ID from: {output}'}

        except subprocess.TimeoutExpired:
            return {'success': False, 'error': 'sbatch command timed out'}
        except Exception as e:
            return {'success': False, 'error': str(e)}
    
    def _check_running_jobs(self):
        """检查运行中的任务状态"""
        completed_jobs = []

        for job_id, job_info in self.running_jobs.items():
            try:
                slurm_job_id = job_info['slurm_job_id']

                # 检查任务是否被用户取消
                if self._check_if_cancelled(job_id, job_info['type']):
                    self.logger.info(f"任务 {job_id} 已被用户取消，执行 scancel")
                    self._cancel_slurm_job(slurm_job_id)
                    completed_jobs.append(job_id)
                    continue

                status = self._check_slurm_status(slurm_job_id)

                if status == 'QUEUED':
                    # Slurm 任务还在排队，更新数据库状态为 QUEUED
                    self._update_job_status(job_id, 'QUEUED', job_info['type'], slurm_job_id=slurm_job_id)

                elif status == 'RUNNING':
                    # Slurm 任务正在运行，更新数据库状态为 RUNNING
                    self._update_job_status(job_id, 'RUNNING', job_info['type'], slurm_job_id=slurm_job_id)

                elif status == 'COMPLETED':
                    self.logger.info(f"任务 {job_id} (Slurm: {slurm_job_id}) 已完成")
                    self._handle_job_completion(job_id, job_info)
                    completed_jobs.append(job_id)

                elif status == 'FAILED':
                    self.logger.error(f"任务 {job_id} (Slurm: {slurm_job_id}) 失败")
                    self._update_job_status(job_id, 'FAILED', job_info['type'])
                    completed_jobs.append(job_id)

                elif status == 'CANCELLED':
                    self.logger.info(f"任务 {job_id} (Slurm: {slurm_job_id}) 已取消")
                    completed_jobs.append(job_id)

            except Exception as e:
                self.logger.error(f"检查任务 {job_id} 状态失败: {e}")

        # 移除已完成的任务
        for job_id in completed_jobs:
            del self.running_jobs[job_id]

    def _check_if_cancelled(self, job_id: int, job_type: str) -> bool:
        """检查任务是否被用户取消或删除（通过 API 查询）"""
        try:
            endpoint = f"{self.api_base_url}/workers/jobs/{job_id}/check_cancelled"
            params = {'job_type': job_type.upper()}

            response = requests.get(
                endpoint,
                headers=self.api_headers,
                params=params,
                timeout=self.config['api']['timeout']
            )

            if response.status_code == 200:
                data = response.json()
                return data.get('cancelled', False)

            # 如果返回 404，说明任务已被删除，应该取消 Slurm 任务
            if response.status_code == 404:
                self.logger.warning(f"任务 {job_id} 在数据库中不存在（可能已被删除），将取消 Slurm 任务")
                return True

            return False

        except Exception as e:
            self.logger.warning(f"检查任务 {job_id} 取消状态失败: {e}")
            return False

    def _cancel_slurm_job(self, slurm_job_id: str):
        """取消 Slurm 任务"""
        try:
            cmd = f"scancel {slurm_job_id}"
            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=10
            )

            if result.returncode == 0:
                self.logger.info(f"成功取消 Slurm 任务 {slurm_job_id}")
            else:
                self.logger.warning(f"取消 Slurm 任务失败: {result.stderr}")

        except Exception as e:
            self.logger.error(f"取消 Slurm 任务 {slurm_job_id} 失败: {e}")

    def _check_slurm_status(self, slurm_job_id: str) -> str:
        """检查 Slurm 任务状态"""
        try:
            cmd = f"squeue -j {slurm_job_id} -h -o %T"
            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=10
            )

            if result.returncode == 0 and result.stdout.strip():
                status = result.stdout.strip()
                # Slurm 状态: PENDING, RUNNING, COMPLETED, FAILED, CANCELLED 等
                if status == 'PENDING':
                    return 'QUEUED'  # Slurm 排队中
                elif status == 'RUNNING':
                    return 'RUNNING'  # Slurm 正在运行
                elif status == 'COMPLETING':
                    return 'RUNNING'  # 正在完成，视为运行中
                elif status == 'COMPLETED':
                    return 'COMPLETED'
                elif status == 'CANCELLED':
                    return 'CANCELLED'
                else:
                    return 'FAILED'
            else:
                # 任务不在队列中，检查是否已完成
                cmd = f"sacct -j {slurm_job_id} -n -o State"
                result = subprocess.run(
                    cmd, shell=True, capture_output=True, text=True, timeout=10
                )
                if result.returncode == 0 and result.stdout.strip():
                    status = result.stdout.strip().split()[0]
                    if 'COMPLETED' in status:
                        return 'COMPLETED'
                    elif 'CANCELLED' in status:
                        return 'CANCELLED'
                    else:
                        return 'FAILED'
                return 'UNKNOWN'

        except Exception as e:
            self.logger.error(f"检查 Slurm 状态失败: {e}")
            return 'UNKNOWN'

    def _handle_job_completion(self, job_id: int, job_info: Dict):
        """处理任务完成"""
        try:
            work_dir = Path(job_info['work_dir'])
            job_type = job_info['type']

            self.logger.info(f"开始处理任务 {job_id} ({job_type}) 的结果")

            # 1. 上传结果文件到 OSS/COS
            uploaded_files = self._upload_results_to_oss(job_id, work_dir)

            # 2. 针对不同任务类型执行后处理
            if job_type == 'qc':
                # QC 任务：解析 Gaussian 输出并上传结果
                qc_result = self._parse_gaussian_output(work_dir)
                if qc_result:
                    # 添加文件路径到结果
                    qc_result['log_file_path'] = next(
                        (f for f in uploaded_files if f.endswith('.log')), None
                    )
                    qc_result['fchk_file_path'] = next(
                        (f for f in uploaded_files if f.endswith('.fchk')), None
                    )
                    # 上传 QC 结果到数据库
                    self._upload_qc_result(job_id, qc_result)
                else:
                    self.logger.warning(f"QC 任务 {job_id} 未能解析 Gaussian 输出")

            elif job_type == 'md':
                # MD 任务：解析 RDF、MSD 等结果并上传
                md_results = self._parse_md_results(work_dir)
                if md_results:
                    self._upload_md_results(job_id, md_results)
                else:
                    self.logger.warning(f"MD 任务 {job_id} 未能解析结果数据")

            # 3. 更新任务状态为 COMPLETED
            self._update_job_status(
                job_id, 'COMPLETED', job_type,
                result_files=uploaded_files,
                progress=100.0
            )

            self.logger.info(f"任务 {job_id} 处理完成，上传了 {len(uploaded_files)} 个文件")

        except Exception as e:
            self.logger.error(f"处理任务 {job_id} 完成失败: {e}", exc_info=True)
            self._update_job_status(job_id, 'FAILED', job_info['type'], error_message=str(e))

    def _parse_gaussian_output(self, work_dir: Path) -> Optional[Dict[str, Any]]:
        """解析 Gaussian 输出文件，提取能量、HOMO、LUMO 等"""
        import re

        try:
            # 查找 .log 文件
            log_files = list(work_dir.glob("*.log"))
            if not log_files:
                self.logger.warning(f"未找到 Gaussian 日志文件: {work_dir}")
                return None

            log_file = log_files[0]
            self.logger.info(f"解析 Gaussian 输出: {log_file.name}")

            content = log_file.read_text(errors='ignore')

            result = {}

            # 解析 SCF 能量 (最后一个)
            energy_matches = re.findall(r'SCF Done:.*?=\s*([-\d.]+)', content)
            if energy_matches:
                result['energy_au'] = float(energy_matches[-1])
                self.logger.info(f"SCF 能量: {result['energy_au']} Hartree")

            # 解析 HOMO 和 LUMO
            # 查找 Alpha 轨道能量
            orbital_section = re.search(
                r'Population analysis.*?Alpha\s+occ\.\s+eigenvalues.*?([\s\S]*?)(?:Alpha virt\.|Beta)',
                content
            )
            if orbital_section:
                occ_text = orbital_section.group(1)
                # 提取所有占据轨道能量
                occ_energies = re.findall(r'[-\d.]+', occ_text)
                if occ_energies:
                    result['homo'] = float(occ_energies[-1])
                    self.logger.info(f"HOMO: {result['homo']} Hartree")

            # 查找虚轨道能量（LUMO）
            virt_section = re.search(
                r'Alpha virt\.\s+eigenvalues.*?([-\d.\s]+)',
                content
            )
            if virt_section:
                virt_text = virt_section.group(1)
                virt_energies = re.findall(r'[-\d.]+', virt_text)
                if virt_energies:
                    result['lumo'] = float(virt_energies[0])
                    self.logger.info(f"LUMO: {result['lumo']} Hartree")

            # 计算 HOMO-LUMO gap (转换为 eV)
            if 'homo' in result and 'lumo' in result:
                gap_hartree = result['lumo'] - result['homo']
                result['homo_lumo_gap'] = gap_hartree * 27.2114  # Hartree to eV
                self.logger.info(f"HOMO-LUMO Gap: {result['homo_lumo_gap']:.3f} eV")

            # 解析偶极矩
            dipole_match = re.search(r'Dipole moment.*?Tot=\s*([\d.]+)', content, re.DOTALL)
            if dipole_match:
                result['dipole_moment'] = float(dipole_match.group(1))
                self.logger.info(f"偶极矩: {result['dipole_moment']} Debye")

            # 检查计算是否正常结束
            if 'Normal termination' not in content:
                self.logger.warning("Gaussian 计算未正常结束")
                # 仍然返回已解析的结果

            return result if result else None

        except Exception as e:
            self.logger.error(f"解析 Gaussian 输出失败: {e}", exc_info=True)
            return None

    def _upload_qc_result(self, job_id: int, result: Dict[str, Any]):
        """上传 QC 计算结果到后端 API"""
        try:
            endpoint = f"{self.api_base_url}/workers/jobs/{job_id}/qc_result"

            response = requests.post(
                endpoint,
                headers=self.api_headers,
                json=result,
                timeout=self.config['api']['timeout']
            )

            if response.status_code == 200:
                data = response.json()
                self.logger.info(f"QC 结果上传成功: job_id={job_id}, result_id={data.get('result_id')}")
            else:
                self.logger.error(f"QC 结果上传失败: {response.status_code} - {response.text}")

        except Exception as e:
            self.logger.error(f"上传 QC 结果失败: {e}", exc_info=True)

    def _parse_md_results(self, work_dir: Path) -> Optional[Dict[str, Any]]:
        """解析 MD 任务的结果数据（RDF、MSD 等）"""
        import re

        results = {
            'rdf_results': [],
            'msd_results': [],
        }

        try:
            # 1. 解析 RDF 数据（从 out_rdf.dat 文件）
            rdf_file = work_dir / "out_rdf.dat"
            if rdf_file.exists():
                rdf_data = self._parse_lammps_rdf(rdf_file)
                if rdf_data:
                    results['rdf_results'] = rdf_data
                    self.logger.info(f"解析到 {len(rdf_data)} 个 RDF 数据对")

            # 2. 解析 MSD 数据（从 msd_*.dat 文件）
            msd_files = list(work_dir.glob("msd_*.dat"))
            for msd_file in msd_files:
                msd_data = self._parse_lammps_msd(msd_file)
                if msd_data:
                    results['msd_results'].append(msd_data)

            if results['msd_results']:
                self.logger.info(f"解析到 {len(results['msd_results'])} 个 MSD 数据")

            # 3. 解析日志文件获取能量、密度等
            log_file = work_dir / "log.lammps"
            if log_file.exists():
                summary = self._parse_lammps_log(log_file)
                results.update(summary)

            return results if (results['rdf_results'] or results['msd_results']) else None

        except Exception as e:
            self.logger.error(f"解析 MD 结果失败: {e}", exc_info=True)
            return None

    def _parse_lammps_rdf(self, rdf_file: Path) -> List[Dict[str, Any]]:
        """解析 LAMMPS RDF 输出文件"""
        import re

        rdf_results = []

        try:
            content = rdf_file.read_text()
            lines = content.strip().split('\n')

            # 解析头部获取原子对信息
            # 格式: # Row c_rdf[1] c_rdf[2] c_rdf[3] c_rdf[4] ...
            header_line = None
            for line in lines:
                if line.startswith('# Row'):
                    header_line = line
                    break

            if not header_line:
                self.logger.warning("RDF 文件缺少头部信息")
                return []

            # 提取列信息
            # c_rdf[1], c_rdf[2] 是第一对的 r 和 g(r)
            # c_rdf[3], c_rdf[4] 是第二对的 r 和 g(r)，以此类推
            columns = header_line.split()[2:]  # 跳过 "# Row"
            num_pairs = len(columns) // 2

            # 读取数据部分（最后一个时间步的数据）
            data_blocks = []
            current_block = []

            for line in lines:
                if line.startswith('#'):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    # 检查是否是新的时间步（第一列是行号，从 1 开始）
                    try:
                        row_num = int(parts[0])
                        if row_num == 1 and current_block:
                            data_blocks.append(current_block)
                            current_block = []
                        current_block.append([float(x) for x in parts[1:]])
                    except ValueError:
                        continue

            if current_block:
                data_blocks.append(current_block)

            if not data_blocks:
                return []

            # 使用最后一个数据块
            last_block = data_blocks[-1]

            # 尝试从配置文件或工作目录名获取原子对标签
            pair_labels = self._get_rdf_pair_labels(rdf_file.parent)

            # 为每个原子对提取 RDF 数据
            for i in range(num_pairs):
                r_col = i * 2
                g_col = i * 2 + 1

                if g_col >= len(last_block[0]):
                    break

                r_values = [row[r_col] for row in last_block if len(row) > g_col]
                g_r_values = [row[g_col] for row in last_block if len(row) > g_col]

                # 计算配位数（对 g(r) 积分）
                coord_numbers = self._calculate_coordination_number(r_values, g_r_values)

                # 查找第一峰
                first_peak_pos, first_peak_height = self._find_first_peak(r_values, g_r_values)

                # 获取原子对标签
                if pair_labels and i < len(pair_labels):
                    center, shell = pair_labels[i]
                else:
                    center = f"Type{i*2+1}"
                    shell = f"Type{i*2+2}"

                rdf_results.append({
                    'center_species': center,
                    'shell_species': shell,
                    'r_values': r_values,
                    'g_r_values': g_r_values,
                    'coordination_number_values': coord_numbers,
                    'first_peak_position': first_peak_pos,
                    'first_peak_height': first_peak_height,
                    'coordination_number': coord_numbers[-1] if coord_numbers else None,
                })

            return rdf_results

        except Exception as e:
            self.logger.error(f"解析 RDF 文件失败: {e}", exc_info=True)
            return []

    def _get_rdf_pair_labels(self, work_dir: Path) -> List[tuple]:
        """从工作目录获取 RDF 原子对标签"""
        # 尝试从 rdf_pairs.json 文件读取
        pairs_file = work_dir / "rdf_pairs.json"
        if pairs_file.exists():
            try:
                import json
                with open(pairs_file) as f:
                    data = json.load(f)
                return [(p['center'], p['target']) for p in data]
            except:
                pass
        return []

    def _calculate_coordination_number(self, r: List[float], g_r: List[float]) -> List[float]:
        """计算配位数（g(r) 的积分）"""
        import math

        coord_numbers = []
        integral = 0.0

        for i in range(1, len(r)):
            dr = r[i] - r[i-1]
            # 使用梯形法则
            avg_g = (g_r[i] + g_r[i-1]) / 2
            avg_r = (r[i] + r[i-1]) / 2
            # N(r) = 4π ∫ r² g(r) ρ dr，这里假设 ρ=1
            integral += 4 * math.pi * avg_r * avg_r * avg_g * dr
            coord_numbers.append(integral)

        return coord_numbers

    def _find_first_peak(self, r: List[float], g_r: List[float]) -> tuple:
        """查找 g(r) 的第一个峰"""
        if len(r) < 3:
            return None, None

        # 找第一个峰（g(r) > 1 的第一个局部最大值）
        for i in range(1, len(g_r) - 1):
            if g_r[i] > g_r[i-1] and g_r[i] > g_r[i+1] and g_r[i] > 1.0:
                return r[i], g_r[i]

        return None, None

    def _parse_lammps_msd(self, msd_file: Path) -> Optional[Dict[str, Any]]:
        """解析 LAMMPS MSD 输出文件"""
        try:
            # 从文件名获取物种名称
            # 格式: msd_Li.dat, msd_PF6.dat 等
            species = msd_file.stem.replace('msd_', '')

            content = msd_file.read_text()
            lines = content.strip().split('\n')

            t_values = []
            msd_x = []
            msd_y = []
            msd_z = []
            msd_total = []

            for line in lines:
                if line.startswith('#'):
                    continue
                parts = line.strip().split()
                if len(parts) >= 5:
                    try:
                        t_values.append(float(parts[0]))
                        msd_x.append(float(parts[1]))
                        msd_y.append(float(parts[2]))
                        msd_z.append(float(parts[3]))
                        msd_total.append(float(parts[4]))
                    except ValueError:
                        continue

            if not t_values:
                return None

            # 计算扩散系数（从 MSD 斜率，D = MSD / (6t)）
            diffusion_coeff = None
            if len(t_values) > 10:
                # 使用后半部分数据线性拟合
                mid = len(t_values) // 2
                t_fit = t_values[mid:]
                msd_fit = msd_total[mid:]

                if len(t_fit) > 2:
                    # 简单线性回归
                    n = len(t_fit)
                    sum_t = sum(t_fit)
                    sum_msd = sum(msd_fit)
                    sum_t_msd = sum(t * m for t, m in zip(t_fit, msd_fit))
                    sum_t2 = sum(t * t for t in t_fit)

                    slope = (n * sum_t_msd - sum_t * sum_msd) / (n * sum_t2 - sum_t * sum_t)
                    # D = slope / 6 (3D), 单位转换 Å²/fs -> cm²/s
                    diffusion_coeff = slope / 6.0 * 1e-4  # Å²/fs 转 cm²/s

            return {
                'species': species,
                't_values': t_values,
                'msd_x_values': msd_x,
                'msd_y_values': msd_y,
                'msd_z_values': msd_z,
                'msd_total_values': msd_total,
                'diffusion_coefficient': diffusion_coeff,
                'labels': {
                    'time': 'fs',
                    'x': f'{species}_x',
                    'y': f'{species}_y',
                    'z': f'{species}_z',
                    'total': f'{species}_total'
                }
            }

        except Exception as e:
            self.logger.error(f"解析 MSD 文件 {msd_file} 失败: {e}", exc_info=True)
            return None

    def _parse_lammps_log(self, log_file: Path) -> Dict[str, Any]:
        """从 LAMMPS 日志文件解析能量、温度等"""
        import re

        result = {}

        try:
            content = log_file.read_text()

            # 查找最后的 thermo 输出
            # 通常格式: Step Temp Press E_total KinEng PotEng Density ...
            thermo_pattern = r'Step\s+Temp\s+.*?\n([\s\S]*?)Loop time'
            matches = re.findall(thermo_pattern, content)

            if matches:
                last_thermo = matches[-1].strip().split('\n')
                if last_thermo:
                    last_line = last_thermo[-1].split()
                    if len(last_line) >= 6:
                        try:
                            result['final_temperature'] = float(last_line[1])
                            result['final_pressure'] = float(last_line[2])
                            result['total_energy'] = float(last_line[3])
                            result['kinetic_energy'] = float(last_line[4])
                            result['potential_energy'] = float(last_line[5])
                            if len(last_line) >= 7:
                                result['final_density'] = float(last_line[6])
                        except (ValueError, IndexError):
                            pass

        except Exception as e:
            self.logger.error(f"解析 LAMMPS 日志失败: {e}")

        return result

    def _upload_md_results(self, job_id: int, results: Dict[str, Any]):
        """上传 MD 结果到后端 API"""
        try:
            endpoint = f"{self.api_base_url}/workers/jobs/{job_id}/md_results"

            response = requests.post(
                endpoint,
                headers=self.api_headers,
                json=results,
                timeout=self.config['api']['timeout']
            )

            if response.status_code == 200:
                data = response.json()
                self.logger.info(
                    f"MD 结果上传成功: job_id={job_id}, "
                    f"RDF={data.get('uploaded', {}).get('rdf', 0)}, "
                    f"MSD={data.get('uploaded', {}).get('msd', 0)}"
                )
            else:
                self.logger.error(f"MD 结果上传失败: {response.status_code} - {response.text}")

        except Exception as e:
            self.logger.error(f"上传 MD 结果失败: {e}", exc_info=True)

    def _extract_last_frame(self, work_dir: Path, job_id: int) -> Optional[Path]:
        """从 LAMMPS 轨迹文件中提取最后一帧"""
        try:
            # 查找 dump 文件
            dump_files = list(work_dir.glob("*.dump"))
            if not dump_files:
                self.logger.warning(f"未找到轨迹文件: {work_dir}")
                return None

            dump_file = dump_files[0]
            self.logger.info(f"提取最后一帧: {dump_file.name}")

            # 读取 dump 文件，找到最后一帧
            last_frame_lines = []
            current_frame = []
            in_frame = False

            with open(dump_file, 'r') as f:
                for line in f:
                    if line.startswith('ITEM: TIMESTEP'):
                        # 新的一帧开始
                        if current_frame:
                            last_frame_lines = current_frame
                        current_frame = [line]
                        in_frame = True
                    elif in_frame:
                        current_frame.append(line)

                # 最后一帧
                if current_frame:
                    last_frame_lines = current_frame

            if not last_frame_lines:
                self.logger.warning("未找到有效帧")
                return None

            # 保存为 PDB 格式（简化版）
            output_file = work_dir / f"final_frame_{job_id}.pdb"

            # 解析 LAMMPS dump 并转换为 PDB
            atoms = []
            atom_section = False
            for line in last_frame_lines:
                if line.startswith('ITEM: ATOMS'):
                    atom_section = True
                    continue
                if atom_section and not line.startswith('ITEM:'):
                    parts = line.split()
                    if len(parts) >= 5:  # id type x y z
                        atoms.append(parts)

            # 写入 PDB 文件
            with open(output_file, 'w') as f:
                f.write("REMARK   Final frame extracted from LAMMPS trajectory\n")
                f.write(f"REMARK   Job ID: {job_id}\n")
                for i, atom in enumerate(atoms, 1):
                    # PDB 格式: ATOM serial name resName chainID resSeq x y z occupancy tempFactor element
                    atom_id = atom[0] if len(atom) > 0 else str(i)
                    atom_type = atom[1] if len(atom) > 1 else "C"
                    x = float(atom[2]) if len(atom) > 2 else 0.0
                    y = float(atom[3]) if len(atom) > 3 else 0.0
                    z = float(atom[4]) if len(atom) > 4 else 0.0

                    f.write(f"ATOM  {i:5d}  {atom_type:3s} MOL A   1    "
                           f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {atom_type:2s}\n")
                f.write("END\n")

            self.logger.info(f"最后一帧已保存: {output_file.name}")
            return output_file

        except Exception as e:
            self.logger.error(f"提取最后一帧失败: {e}", exc_info=True)
            return None

    def _upload_results_to_oss(self, job_id: int, work_dir: Path) -> List[str]:
        """上传结果文件到对象存储（COS 或 OSS）"""
        uploaded_files = []

        try:
            # 1. 处理轨迹文件：提取最后一帧
            if self.config['upload'].get('trajectory_handling', {}).get('enabled', True):
                if self.config['upload']['trajectory_handling'].get('extract_last_frame', True):
                    last_frame_file = self._extract_last_frame(work_dir, job_id)
                    if last_frame_file:
                        # 将最后一帧添加到上传列表
                        self.logger.info(f"将上传最后一帧: {last_frame_file.name}")

            # 2. 获取需要上传的文件
            essential_patterns = self.config['upload']['essential_files']
            optional_patterns = self.config['upload'].get('optional_large_files', [])
            excluded_patterns = self.config['upload'].get('excluded_files', [])

            max_size = self.config['upload']['max_file_size'] * 1024 * 1024  # MB to bytes
            compress_threshold = self.config['upload'].get('compress_threshold', 50) * 1024 * 1024

            # 获取结果文件前缀
            if self.storage_type == 'cos':
                result_prefix = self.config['cos']['result_prefix']
            else:
                result_prefix = self.config['oss']['result_prefix']

            # 3. 上传必要文件
            all_patterns = essential_patterns.copy()
            if self.config['upload']['strategy'].get('upload_optional_on_demand', False):
                # 如果配置了按需上传，这里不上传可选文件
                pass
            else:
                all_patterns.extend(optional_patterns)

            for pattern in all_patterns:
                for file_path in work_dir.glob(pattern):
                    if not file_path.is_file():
                        continue

                    # 检查是否在排除列表中
                    excluded = False
                    for exclude_pattern in excluded_patterns:
                        if file_path.match(exclude_pattern):
                            self.logger.info(f"跳过排除的文件: {file_path.name}")
                            excluded = True
                            break
                    if excluded:
                        continue

                    # 检查文件大小
                    file_size = file_path.stat().st_size
                    if file_size > max_size:
                        self.logger.warning(f"文件 {file_path.name} ({file_size/1024/1024:.1f}MB) 超过大小限制，跳过")
                        continue

                    # 构建对象存储 Key
                    object_key = f"{result_prefix}{job_id}/{file_path.name}"

                    self.logger.info(f"上传文件: {file_path.name} ({file_size/1024/1024:.1f}MB)")

                    # 根据存储类型上传
                    if self.storage_type == 'cos':
                        # 腾讯云 COS
                        with open(file_path, 'rb') as f:
                            self.cos_client.put_object(
                                Bucket=self.cos_bucket,
                                Body=f,
                                Key=object_key
                            )
                    else:
                        # 阿里云 OSS
                        self.oss_bucket.put_object_from_file(object_key, str(file_path))

                    uploaded_files.append(object_key)
                    self.logger.info(f"✅ 上传成功: {file_path.name}")

            self.logger.info(f"共上传 {len(uploaded_files)} 个文件")
            return uploaded_files

        except Exception as e:
            self.logger.error(f"上传结果文件失败: {e}", exc_info=True)
            raise

    def _update_job_status(
        self,
        job_id: int,
        status: str,
        job_type: str,
        slurm_job_id: Optional[str] = None,
        work_dir: Optional[str] = None,
        error_message: Optional[str] = None,
        result_files: Optional[List[str]] = None,
        progress: Optional[float] = None
    ):
        """更新任务状态到阿里云"""
        try:
            endpoint = f"{self.api_base_url}/workers/jobs/{job_id}/status"

            # 确保状态值有效
            valid_statuses = ["CREATED", "SUBMITTED", "QUEUED", "RUNNING", "POSTPROCESSING", "COMPLETED", "FAILED", "CANCELLED"]
            if status not in valid_statuses:
                self.logger.error(f"无效的状态值: {status}，有效值: {valid_statuses}")
                return

            data = {
                'status': status,
                'job_type': job_type.upper(),
                'worker_name': self.config['worker']['name']
            }

            if slurm_job_id:
                data['slurm_job_id'] = slurm_job_id
            if work_dir:
                data['work_dir'] = work_dir
            if error_message:
                # 截断错误消息，防止过长
                data['error_message'] = str(error_message)[:500]
            if result_files:
                data['result_files'] = result_files
            if progress is not None:
                data['progress'] = progress

            response = requests.put(
                endpoint,
                headers=self.api_headers,
                json=data,
                timeout=self.config['api']['timeout']
            )

            if response.status_code == 200:
                self.logger.info(f"任务 {job_id} 状态已更新为 {status}")
            else:
                self.logger.warning(f"更新任务状态失败: {response.status_code} - {response.text}")

        except Exception as e:
            self.logger.error(f"更新任务 {job_id} 状态失败: {e}", exc_info=True)


def main():
    """主函数"""
    import argparse

    parser = argparse.ArgumentParser(description='混合云轮询 Worker')
    parser.add_argument(
        '--config',
        default='deployment/polling_worker_config.yaml',
        help='配置文件路径'
    )
    args = parser.parse_args()

    # 创建并运行 Worker
    worker = PollingWorker(config_path=args.config)
    worker.run()


if __name__ == '__main__':
    main()

