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
                    # 获取失败原因
                    error_msg = self._get_job_failure_reason(slurm_job_id, job_info.get('work_dir'))
                    self._update_job_status(job_id, 'FAILED', job_info['type'], error_message=error_msg)
                    completed_jobs.append(job_id)

                elif status == 'CANCELLED':
                    self.logger.info(f"任务 {job_id} (Slurm: {slurm_job_id}) 已取消")
                    self._update_job_status(job_id, 'CANCELLED', job_info['type'], error_message="任务被取消")
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

    def _get_job_failure_reason(self, slurm_job_id: str, work_dir: str = None) -> str:
        """获取任务失败的原因"""
        reasons = []

        try:
            # 1. 从 sacct 获取失败原因
            cmd = f"sacct -j {slurm_job_id} -n -o State,ExitCode,Reason --parsable2"
            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0 and result.stdout.strip():
                lines = result.stdout.strip().split('\n')
                for line in lines:
                    parts = line.split('|')
                    if len(parts) >= 3:
                        state, exit_code, reason = parts[0], parts[1], parts[2]
                        if 'FAILED' in state or 'TIMEOUT' in state or 'OUT_OF_MEMORY' in state:
                            reasons.append(f"Slurm状态: {state}, 退出码: {exit_code}")
                            if reason and reason != 'None':
                                reasons.append(f"原因: {reason}")
                            break

            # 2. 检查 slurm 输出文件
            if work_dir:
                work_path = Path(work_dir)
                slurm_out = work_path / f"slurm-{slurm_job_id}.out"
                if slurm_out.exists():
                    # 读取最后 50 行
                    content = slurm_out.read_text()
                    last_lines = content.strip().split('\n')[-50:]

                    # 查找错误信息
                    error_keywords = ['ERROR', 'Error', 'error', 'FATAL', 'Fatal',
                                     'Segmentation fault', 'OOM', 'Out of memory',
                                     'CANCELLED', 'TIME LIMIT', 'srun: error']
                    for line in last_lines:
                        if any(kw in line for kw in error_keywords):
                            reasons.append(f"日志: {line.strip()[:200]}")
                            break

                # 3. 检查 LAMMPS 日志（如果是 MD 任务）
                log_lammps = work_path / "log.lammps"
                if log_lammps.exists():
                    content = log_lammps.read_text()
                    last_lines = content.strip().split('\n')[-30:]
                    for line in last_lines:
                        if 'ERROR' in line or 'error' in line.lower():
                            reasons.append(f"LAMMPS: {line.strip()[:200]}")
                            break

                # 4. 检查 Gaussian 日志（如果是 QC 任务）
                for log_file in work_path.glob("*.log"):
                    content = log_file.read_text()
                    if 'Error termination' in content or 'Convergence failure' in content:
                        # 获取错误行
                        lines = content.split('\n')
                        for i, line in enumerate(lines):
                            if 'Error termination' in line or 'Convergence failure' in line:
                                reasons.append(f"Gaussian: {line.strip()[:200]}")
                                break
                        break

        except Exception as e:
            self.logger.error(f"获取失败原因时出错: {e}")
            reasons.append(f"获取详细原因失败: {str(e)}")

        if reasons:
            return "; ".join(reasons)
        else:
            return f"Slurm 任务 {slurm_job_id} 失败，未能获取详细原因"

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

                    # 生成 QC 可视化图片（ESP、HOMO、LUMO）
                    fchk_files = list(work_dir.glob("*.fchk"))
                    if fchk_files:
                        fchk_file = fchk_files[0]
                        try:
                            vis_result = self._generate_qc_visualizations(work_dir, fchk_file)
                            if vis_result:
                                qc_result.update(vis_result)
                                self.logger.info(f"QC 可视化生成成功: ESP={vis_result.get('esp_image_path')}, HOMO={vis_result.get('homo_image_path')}, LUMO={vis_result.get('lumo_image_path')}")

                                # 上传新生成的PNG图片文件
                                png_files = list(work_dir.glob("*.png"))
                                if png_files:
                                    additional_uploaded = self._upload_additional_files(job_id, png_files)
                                    uploaded_files.extend(additional_uploaded)
                                    self.logger.info(f"额外上传了 {len(additional_uploaded)} 个图片文件")
                        except Exception as e:
                            self.logger.warning(f"QC 可视化生成失败: {e}")

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
            # 查找真正的 Gaussian 输出文件
            # 优先级：1. *_out.log (Gaussian输出) 2. *.out (备选) 3. 其他.log文件
            gaussian_log_file = None

            # 首先查找 *_out.log 格式的文件（这是真正的Gaussian输出）
            out_log_files = list(work_dir.glob("*_out.log"))
            if out_log_files:
                # 排除 qc_out.log（这是Slurm的输出）
                gaussian_files = [f for f in out_log_files if f.name != "qc_out.log"]
                if gaussian_files:
                    gaussian_log_file = gaussian_files[0]

            # 如果没找到，尝试查找 .out 文件
            if not gaussian_log_file:
                out_files = list(work_dir.glob("*.out"))
                if out_files:
                    gaussian_log_file = out_files[0]

            # 最后尝试其他 .log 文件（排除已知的非Gaussian文件）
            if not gaussian_log_file:
                log_files = list(work_dir.glob("*.log"))
                excluded_names = {"qc_out.log", "qc_err.log", "slurm.log"}
                gaussian_files = [f for f in log_files if f.name not in excluded_names]
                if gaussian_files:
                    gaussian_log_file = gaussian_files[0]

            if not gaussian_log_file:
                self.logger.warning(f"未找到 Gaussian 输出文件: {work_dir}")
                return None

            self.logger.info(f"解析 Gaussian 输出: {gaussian_log_file.name}")
            content = gaussian_log_file.read_text(errors='ignore')

            result = {}

            # 使用后端成熟的解析逻辑
            # 正则表达式模式
            energy_pattern = re.compile(r'SCF Done:.*?=\s*([-\d.]+)')
            alpha_occ_pattern = re.compile(r'Alpha\s+occ\.\s+eigenvalues\s+--\s+(.*)')
            alpha_virt_pattern = re.compile(r'Alpha\s+virt\.\s+eigenvalues\s+--\s+(.*)')

            last_energy = None
            last_homo = None
            last_lumo = None

            # 逐行处理，确保配对匹配
            lines = content.splitlines()
            for i in range(len(lines) - 1):
                # 匹配SCF能量
                match_energy = energy_pattern.search(lines[i])
                if match_energy:
                    last_energy = float(match_energy.group(1))

                # 匹配HOMO和LUMO（确保是连续的行）
                match_occ = alpha_occ_pattern.search(lines[i])
                match_virt = alpha_virt_pattern.search(lines[i + 1]) if i + 1 < len(lines) else None

                if match_occ and match_virt:
                    occ_values = match_occ.group(1).split()
                    virt_values = match_virt.group(1).split()

                    if occ_values:
                        last_homo = float(occ_values[-1])
                    if virt_values:
                        last_lumo = float(virt_values[0])

            # 设置结果
            if last_energy is not None:
                result['energy_au'] = last_energy
                self.logger.info(f"SCF 能量: {result['energy_au']} Hartree")

            if last_homo is not None:
                result['homo'] = last_homo
                self.logger.info(f"HOMO: {result['homo']} Hartree")

            if last_lumo is not None:
                result['lumo'] = last_lumo
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

    def _generate_qc_visualizations(self, work_dir: Path, fchk_file: Path) -> Optional[Dict[str, Any]]:
        """
        生成 QC 可视化图片（ESP、HOMO、LUMO）
        使用后端的 qc_postprocess 模块中的函数
        """
        import sys
        import base64

        try:
            # 添加后端模块路径
            backend_path = Path(__file__).parent.parent / "backend"
            if str(backend_path) not in sys.path:
                sys.path.insert(0, str(backend_path))

            from app.tasks.qc_postprocess import (
                generate_esp_visualization,
                generate_orbital_visualization,
                extract_esp_values
            )

            result = {}

            # 获取分子名称
            molecule_name = fchk_file.stem

            # 1. 生成 ESP 图片
            try:
                esp_image_path = generate_esp_visualization(work_dir, molecule_name, str(fchk_file))
                if esp_image_path:
                    result['esp_image_path'] = esp_image_path
                    self.logger.info(f"ESP 图片生成: {esp_image_path}")

                    # 读取图片并转换为base64（用于混合云架构）
                    try:
                        with open(esp_image_path, 'rb') as f:
                            esp_image_data = f.read()
                            result['esp_image_content'] = base64.b64encode(esp_image_data).decode('utf-8')
                            self.logger.info(f"ESP 图片已编码为base64，大小: {len(esp_image_data)} bytes")
                    except Exception as e:
                        self.logger.warning(f"读取ESP图片失败: {e}")

                # 提取 ESP 极值
                surfanalysis_file = work_dir / "surfanalysis.txt"
                esp_min, esp_max = extract_esp_values(str(surfanalysis_file))
                if esp_min is not None:
                    result['esp_min_kcal'] = esp_min
                if esp_max is not None:
                    result['esp_max_kcal'] = esp_max
            except Exception as e:
                self.logger.warning(f"ESP 可视化生成失败: {e}")

            # 2. 生成 HOMO 轨道图片
            try:
                homo_image_path = generate_orbital_visualization(work_dir, str(fchk_file), "HOMO")
                if homo_image_path:
                    result['homo_image_path'] = homo_image_path
                    self.logger.info(f"HOMO 图片生成: {homo_image_path}")

                    # 读取图片并转换为base64（用于混合云架构）
                    try:
                        with open(homo_image_path, 'rb') as f:
                            homo_image_data = f.read()
                            result['homo_image_content'] = base64.b64encode(homo_image_data).decode('utf-8')
                            self.logger.info(f"HOMO 图片已编码为base64，大小: {len(homo_image_data)} bytes")
                    except Exception as e:
                        self.logger.warning(f"读取HOMO图片失败: {e}")
            except Exception as e:
                self.logger.warning(f"HOMO 可视化生成失败: {e}")

            # 3. 生成 LUMO 轨道图片
            try:
                lumo_image_path = generate_orbital_visualization(work_dir, str(fchk_file), "LUMO")
                if lumo_image_path:
                    result['lumo_image_path'] = lumo_image_path
                    self.logger.info(f"LUMO 图片生成: {lumo_image_path}")

                    # 读取图片并转换为base64（用于混合云架构）
                    try:
                        with open(lumo_image_path, 'rb') as f:
                            lumo_image_data = f.read()
                            result['lumo_image_content'] = base64.b64encode(lumo_image_data).decode('utf-8')
                            self.logger.info(f"LUMO 图片已编码为base64，大小: {len(lumo_image_data)} bytes")
                    except Exception as e:
                        self.logger.warning(f"读取LUMO图片失败: {e}")
            except Exception as e:
                self.logger.warning(f"LUMO 可视化生成失败: {e}")

            return result if result else None

        except ImportError as e:
            self.logger.warning(f"无法导入 qc_postprocess 模块: {e}")
            return None
        except Exception as e:
            self.logger.error(f"生成 QC 可视化失败: {e}", exc_info=True)
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
        """
        解析 MD 任务的结果数据（RDF、MSD 等）

        使用后端的 LAMMPS 读取器来正确解析数据，确保：
        1. RDF 标签正确映射到分子/原子名称
        2. MSD 文件正确解析（out_*_msd.dat 格式）
        """
        import sys

        results = {
            'rdf_results': [],
            'msd_results': [],
        }

        try:
            # 添加后端模块路径
            backend_path = Path(__file__).parent.parent / "backend"
            if str(backend_path) not in sys.path:
                sys.path.insert(0, str(backend_path))

            # 1. 解析 RDF 数据（使用后端的 LAMMPSRDFReader）
            try:
                from app.workers.lammps_rdf_reader import LAMMPSRDFReader

                reader = LAMMPSRDFReader(work_dir)
                rdf_data = reader.read_rdf_file()

                if rdf_data:
                    for data in rdf_data:
                        # 转换为 API 上传格式
                        rdf_item = {
                            'center_species': data['center_label'],
                            'shell_species': data['target_label'],
                            'r_values': data['r'],
                            'g_r_values': data['g_r'],
                            'coordination_number_values': data.get('coordination_number'),
                        }
                        # 计算第一峰信息
                        first_peak_pos, first_peak_height = self._find_first_peak(data['r'], data['g_r'])
                        rdf_item['first_peak_position'] = first_peak_pos
                        rdf_item['first_peak_height'] = first_peak_height
                        if data.get('coordination_number'):
                            rdf_item['coordination_number'] = data['coordination_number'][-1]

                        results['rdf_results'].append(rdf_item)

                    self.logger.info(f"解析到 {len(results['rdf_results'])} 个 RDF 数据对")
            except ImportError as e:
                self.logger.warning(f"无法导入 LAMMPSRDFReader，使用简化解析: {e}")
                # 降级到简化解析
                rdf_file = work_dir / "out_rdf.dat"
                if rdf_file.exists():
                    rdf_data = self._parse_lammps_rdf_simple(rdf_file, work_dir)
                    if rdf_data:
                        results['rdf_results'] = rdf_data
                        self.logger.info(f"解析到 {len(rdf_data)} 个 RDF 数据对（简化模式）")

            # 2. 解析 MSD 数据（使用后端的完整计算逻辑）
            try:
                from app.workers.lammps_msd_reader import (
                    LAMMPSMSDReader,
                    calculate_diffusion_coefficient,
                    calculate_ionic_conductivity,
                    calculate_mobility,
                )
                from app.tasks.msd_processor import extract_box_volume_and_ion_counts, get_ion_charge

                msd_reader = LAMMPSMSDReader(work_dir)
                msd_files = msd_reader.find_msd_files()

                # 提取盒子体积和离子数量（用于电导率计算）
                box_volume, ion_counts = extract_box_volume_and_ion_counts(work_dir)
                self.logger.info(f"Box volume: {box_volume}, Ion counts: {ion_counts}")

                temperature = 298.15  # 默认室温

                for msd_file in msd_files:
                    msd_data = msd_reader.read_msd_file(msd_file)
                    if msd_data:
                        species = msd_data['species']
                        time_arr = msd_data['time']
                        msd_total = msd_data['msd_total']

                        # 计算扩散系数
                        diffusion_coeff = calculate_diffusion_coefficient(time_arr, msd_total)

                        # 获取离子电荷
                        charge = get_ion_charge(species)

                        # 计算迁移率
                        mobility = calculate_mobility(diffusion_coeff, charge, temperature)

                        # 计算电导率（需要盒子体积和离子数量）
                        ionic_conductivity = None
                        if box_volume and ion_counts and diffusion_coeff:
                            ion_count = ion_counts.get(species, 0)
                            # 模糊匹配
                            if ion_count == 0:
                                for key, val in ion_counts.items():
                                    if key in species or species in key:
                                        ion_count = val
                                        break

                            if ion_count > 0:
                                ionic_conductivity = calculate_ionic_conductivity(
                                    diffusion_coeff, ion_count, box_volume, charge, temperature
                                )

                        # 转换为 API 上传格式
                        msd_item = {
                            'species': species,
                            't_values': time_arr,
                            'msd_x_values': msd_data.get('msd_x'),
                            'msd_y_values': msd_data.get('msd_y'),
                            'msd_z_values': msd_data.get('msd_z'),
                            'msd_total_values': msd_total,
                            'diffusion_coefficient': diffusion_coeff,
                            'mobility': mobility,
                            'ionic_conductivity': ionic_conductivity,
                            'charge': charge,
                            'labels': msd_data.get('labels'),
                        }
                        results['msd_results'].append(msd_item)
                        self.logger.info(f"MSD {species}: D={diffusion_coeff}, μ={mobility}, σ={ionic_conductivity}")

                if results['msd_results']:
                    self.logger.info(f"解析到 {len(results['msd_results'])} 个 MSD 数据（完整计算）")
            except ImportError as e:
                self.logger.warning(f"无法导入后端 MSD 模块，使用简化解析: {e}")
                # 降级到简化解析 - 正确的文件名模式
                msd_files = list(work_dir.glob("out_*_msd.dat"))
                for msd_file in msd_files:
                    msd_data = self._parse_lammps_msd_simple(msd_file)
                    if msd_data:
                        results['msd_results'].append(msd_data)

                if results['msd_results']:
                    self.logger.info(f"解析到 {len(results['msd_results'])} 个 MSD 数据（简化模式）")

            # 3. 解析日志文件获取能量、密度等
            log_file = work_dir / "log.lammps"
            if log_file.exists():
                summary = self._parse_lammps_log(log_file)
                results.update(summary)

            # 4. 分析溶剂化结构（使用后端的 solvation 服务）
            electrolyte_data = self._load_electrolyte_config(work_dir)

            try:
                from app.services.solvation import analyze_solvation_structures

                if electrolyte_data:
                    self.logger.info(f"开始溶剂化结构分析: {work_dir}")
                    solvation_results = analyze_solvation_structures(
                        work_dir=str(work_dir),
                        electrolyte_data=electrolyte_data,
                        cutoff=3.0,  # 默认截断距离
                    )

                    if solvation_results:
                        # 读取 XYZ 文件内容
                        for solv in solvation_results:
                            if solv.get('file_path'):
                                try:
                                    with open(solv['file_path'], 'r') as f:
                                        solv['xyz_content'] = f.read()
                                except Exception as e:
                                    self.logger.warning(f"读取溶剂化 XYZ 失败: {e}")

                        results['solvation_structures'] = solvation_results
                        self.logger.info(f"溶剂化结构分析完成: {len(solvation_results)} 个结构")
                else:
                    self.logger.warning("未找到电解液配置，跳过溶剂化分析")
            except ImportError as e:
                self.logger.warning(f"无法导入溶剂化分析模块: {e}")
            except Exception as e:
                self.logger.warning(f"溶剂化分析失败: {e}")

            # 5. 提取系统结构（使用后端的 get_system_structure）
            try:
                from app.services.solvation import get_system_structure

                system_result = get_system_structure(str(work_dir), frame_idx=-1)
                if system_result and 'xyz_content' in system_result:
                    results['system_xyz_content'] = system_result['xyz_content']
                    # 从系统结构获取盒子尺寸（如果日志解析没有）
                    if 'box' in system_result and system_result['box']:
                        box = system_result['box']
                        if not results.get('box_x') and len(box) >= 3:
                            results['box_x'] = box[0]
                            results['box_y'] = box[1]
                            results['box_z'] = box[2]
                    self.logger.info(f"系统结构提取完成: {system_result.get('atom_count', 0)} 原子")
            except ImportError as e:
                self.logger.warning(f"无法导入 get_system_structure: {e}")
            except Exception as e:
                self.logger.warning(f"提取系统结构失败: {e}")

            # 6. 提取分子结构（PDB 和电荷信息）
            try:
                molecule_structures = self._extract_molecule_structures(work_dir)
                if molecule_structures:
                    results['molecule_structures'] = molecule_structures
                    self.logger.info(f"分子结构提取完成: {len(molecule_structures)} 个分子")
            except Exception as e:
                self.logger.warning(f"提取分子结构失败: {e}")

            # 7. 计算浓度（如果有盒子尺寸和电解液配置）
            if electrolyte_data and results.get('box_x'):
                try:
                    conc_data = self._calculate_concentration(results, electrolyte_data)
                    results.update(conc_data)
                except Exception as e:
                    self.logger.warning(f"计算浓度失败: {e}")

            return results if (results['rdf_results'] or results['msd_results'] or results.get('solvation_structures') or results.get('system_xyz_content') or results.get('molecule_structures')) else None

        except Exception as e:
            self.logger.error(f"解析 MD 结果失败: {e}", exc_info=True)
            return None

    def _load_electrolyte_config(self, work_dir: Path) -> Optional[Dict[str, Any]]:
        """加载电解液配置"""
        import json

        # 尝试从 electrolyte.json 加载
        electrolyte_file = work_dir / "electrolyte.json"
        if electrolyte_file.exists():
            try:
                with open(electrolyte_file) as f:
                    return json.load(f)
            except Exception as e:
                self.logger.warning(f"读取 electrolyte.json 失败: {e}")

        # 尝试从 job_config.json 加载
        config_file = work_dir / "job_config.json"
        if config_file.exists():
            try:
                with open(config_file) as f:
                    config = json.load(f)
                    if 'electrolyte' in config:
                        return config['electrolyte']
            except Exception as e:
                self.logger.warning(f"读取 job_config.json 失败: {e}")

        # 尝试从 atom_mapping.json 构建基本配置
        atom_mapping_file = work_dir / "atom_mapping.json"
        if atom_mapping_file.exists():
            try:
                with open(atom_mapping_file) as f:
                    atom_mapping = json.load(f)

                # 从 atom_mapping 提取分子信息
                molecules = atom_mapping.get('molecules', [])
                if molecules:
                    # 统计各类型分子
                    from collections import Counter
                    mol_names = [m.get('molecule_name', '') for m in molecules]
                    mol_counts = Counter(mol_names)

                    # 构建基本电解液配置
                    cations = ['Li', 'Na', 'K', 'Mg', 'Ca', 'Zn']
                    anions = ['FSI', 'TFSI', 'PF6', 'BF4', 'ClO4', 'DCA', 'Cl', 'Br', 'I']

                    electrolyte_config = {
                        'cations': [],
                        'anions': [],
                        'solvents': [],
                    }

                    for mol_name, count in mol_counts.items():
                        if any(c in mol_name for c in cations):
                            electrolyte_config['cations'].append({
                                'name': mol_name, 'count': count
                            })
                        elif any(a in mol_name for a in anions):
                            electrolyte_config['anions'].append({
                                'name': mol_name, 'count': count
                            })
                        else:
                            electrolyte_config['solvents'].append({
                                'name': mol_name, 'count': count
                            })

                    if electrolyte_config['cations']:
                        return electrolyte_config
            except Exception as e:
                self.logger.warning(f"从 atom_mapping 构建配置失败: {e}")

        return None

    def _calculate_concentration(self, results: Dict, electrolyte_data: Dict) -> Dict[str, float]:
        """根据盒子尺寸计算浓度"""
        conc_result = {}
        AVOGADRO = 6.022e23

        try:
            # 计算阳离子数量
            cation_count = 0
            for cat in electrolyte_data.get('cations', []):
                cation_count += cat.get('count', cat.get('number', 0))

            if cation_count == 0:
                return conc_result

            # 计算最终浓度
            if results.get('box_x') and results.get('box_y') and results.get('box_z'):
                volume_L = (results['box_x'] * results['box_y'] * results['box_z']) * 1e-27
                if volume_L > 0:
                    conc_result['concentration'] = round((cation_count / AVOGADRO) / volume_L, 4)

            # 计算初始浓度
            if results.get('initial_box_x') and results.get('initial_box_y') and results.get('initial_box_z'):
                init_volume_L = (results['initial_box_x'] * results['initial_box_y'] * results['initial_box_z']) * 1e-27
                if init_volume_L > 0:
                    conc_result['initial_concentration'] = round((cation_count / AVOGADRO) / init_volume_L, 4)

        except Exception as e:
            self.logger.warning(f"浓度计算失败: {e}")

        return conc_result

    def _extract_molecule_structures(self, work_dir: Path) -> List[Dict[str, Any]]:
        """
        从工作目录提取分子结构信息（PDB 内容和电荷）
        这是后端 get_molecule_templates API 的逻辑移植到 Worker
        """
        import json
        import re

        molecules = []
        seen_molecules = set()

        try:
            # 获取任务名称（用于过滤系统级别的 .lt 文件）
            job_name = work_dir.name

            # 首先找到所有 .lt 文件（这些是实际使用的分子）
            lt_files = list(work_dir.glob("*.lt"))

            for lt_file in lt_files:
                base_name = lt_file.stem

                # 跳过系统级别的 .lt 文件
                if base_name == job_name or len(base_name) > 50:
                    continue

                if base_name in seen_molecules:
                    continue
                seen_molecules.add(base_name)

                # 查找对应的 PDB 文件
                pdb_file = None
                pdb_candidates = [
                    work_dir / f"{base_name}.charmm.pdb",
                    work_dir / f"{base_name}.q.pdb",
                    work_dir / f"{base_name}.pdb",
                ]

                for candidate in pdb_candidates:
                    if candidate.exists():
                        pdb_file = candidate
                        break

                if not pdb_file:
                    continue

                # 读取 PDB 内容
                try:
                    pdb_content = pdb_file.read_text(encoding='utf-8')
                except UnicodeDecodeError:
                    pdb_content = pdb_file.read_text(encoding='latin-1')

                # 解析原子
                atoms = []
                for line in pdb_content.split('\n'):
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        try:
                            atom_id = int(line[6:11].strip())
                            atom_name = line[12:16].strip()
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            element = line[76:78].strip() if len(line) > 76 else atom_name[0]

                            atoms.append({
                                "id": atom_id,
                                "name": atom_name,
                                "element": element,
                                "x": x,
                                "y": y,
                                "z": z,
                                "charge": None
                            })
                        except (ValueError, IndexError):
                            continue

                # 从 .lt 文件读取电荷
                lt_content = lt_file.read_text()
                charge_map = {}

                # 方法1: 从 "In Charges" 部分解析
                charges_section = re.search(r'write_once\("In Charges"\)\s*\{(.*?)\}', lt_content, re.DOTALL)
                if charges_section:
                    for line in charges_section.group(1).split('\n'):
                        match = re.search(r'set type @atom:(\w+)\s+charge\s+([-\d.]+)', line)
                        if match:
                            charge_map[match.group(1)] = float(match.group(2))

                # 方法2: 从 "Data Atoms" 部分解析
                atoms_section = re.search(r'write\("Data Atoms"\)\s*\{(.*?)\}', lt_content, re.DOTALL)
                if atoms_section:
                    atom_lines = []
                    for line in atoms_section.group(1).split('\n'):
                        line = line.strip()
                        if line and not line.startswith('#'):
                            match = re.search(
                                r'\$atom:(\w+)\s+\$mol[:\w]*\s+@atom:(\w+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)',
                                line
                            )
                            if match:
                                atom_lines.append({
                                    'charge': float(match.group(3)),
                                    'x': float(match.group(4)),
                                    'y': float(match.group(5)),
                                    'z': float(match.group(6))
                                })

                    # 按坐标或顺序匹配电荷
                    for i, lt_atom in enumerate(atom_lines):
                        matched = False
                        for pdb_atom in atoms:
                            if (abs(pdb_atom['x'] - lt_atom['x']) < 0.01 and
                                abs(pdb_atom['y'] - lt_atom['y']) < 0.01 and
                                abs(pdb_atom['z'] - lt_atom['z']) < 0.01):
                                pdb_atom['charge'] = lt_atom['charge']
                                matched = True
                                break

                        if not matched and i < len(atoms):
                            atoms[i]['charge'] = lt_atom['charge']

                # 应用 charge_map
                for atom in atoms:
                    if atom['charge'] is None and atom['element'] in charge_map:
                        atom['charge'] = charge_map[atom['element']]

                # 确定分子类型
                total_charge = sum(a['charge'] for a in atoms if a['charge'] is not None)

                if total_charge > 0.5:
                    mol_type = "cation"
                elif total_charge < -0.5:
                    mol_type = "anion"
                else:
                    mol_type = "solvent"

                molecules.append({
                    "name": base_name,
                    "type": mol_type,
                    "pdb_content": pdb_content,
                    "atoms": atoms,
                    "total_charge": total_charge,
                    "charge_method": "resp"
                })

            return molecules

        except Exception as e:
            self.logger.error(f"提取分子结构失败: {e}", exc_info=True)
            return []

    def _parse_lammps_rdf_simple(self, rdf_file: Path, work_dir: Path) -> List[Dict[str, Any]]:
        """
        简化版 RDF 解析（当无法导入后端模块时使用）
        尝试从 atom_mapping.json 和 .in.list 获取标签
        """
        import json
        import re

        rdf_results = []

        try:
            # 尝试加载 atom_mapping.json 获取标签信息
            atom_mapping_file = work_dir / "atom_mapping.json"
            type_to_label = {}

            if atom_mapping_file.exists():
                with open(atom_mapping_file) as f:
                    atom_mapping = json.load(f)

                # 从 atom_types 提取标签
                if 'atom_types' in atom_mapping:
                    for at in atom_mapping['atom_types']:
                        type_id = at.get('type_id')
                        label = at.get('label') or at.get('element', f'Type{type_id}')
                        mol_name = at.get('molecule_name', '')
                        if mol_name:
                            type_to_label[type_id] = f"{mol_name}_{label}"
                        else:
                            type_to_label[type_id] = label

            # 尝试从 .in.list 获取 RDF 对
            in_list_file = work_dir / f"{work_dir.name}.in.list"
            rdf_pairs = []

            if in_list_file.exists():
                content = in_list_file.read_text()
                # 解析 compute rdf 命令
                # 格式: compute rdf_compute all rdf 100 1 2 1 3 2 3 ...
                for line in content.split('\n'):
                    if 'compute' in line and 'rdf' in line:
                        parts = line.split()
                        # 找到数字对
                        numbers = []
                        for p in parts:
                            try:
                                numbers.append(int(p))
                            except ValueError:
                                continue
                        # 第一个数字是 bins，后面是原子类型对
                        if len(numbers) > 1:
                            for i in range(1, len(numbers) - 1, 2):
                                type1, type2 = numbers[i], numbers[i + 1]
                                label1 = type_to_label.get(type1, f"Type{type1}")
                                label2 = type_to_label.get(type2, f"Type{type2}")
                                rdf_pairs.append((label1, label2))
                        break

            # 读取 RDF 数据
            content = rdf_file.read_text()
            lines = content.strip().split('\n')

            # 找到数据块
            data_blocks = []
            current_block = []

            for line in lines:
                if line.startswith('#'):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
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
            num_pairs = len(rdf_pairs) if rdf_pairs else (len(last_block[0]) // 2)

            for i in range(num_pairs):
                r_col = i * 2
                g_col = i * 2 + 1

                if g_col >= len(last_block[0]):
                    break

                r_values = [row[r_col] for row in last_block if len(row) > g_col]
                g_r_values = [row[g_col] for row in last_block if len(row) > g_col]

                # 获取标签
                if i < len(rdf_pairs):
                    center, shell = rdf_pairs[i]
                else:
                    center = f"Type{i*2+1}"
                    shell = f"Type{i*2+2}"

                # 计算配位数和第一峰
                coord_numbers = self._calculate_coordination_number(r_values, g_r_values)
                first_peak_pos, first_peak_height = self._find_first_peak(r_values, g_r_values)

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
            self.logger.error(f"简化解析 RDF 文件失败: {e}", exc_info=True)
            return []

    def _parse_lammps_msd_simple(self, msd_file: Path) -> Optional[Dict[str, Any]]:
        """
        简化版 MSD 解析（当无法导入后端模块时使用）
        文件名格式: out_Li_msd.dat, out_FSI_msd.dat 等
        """
        import re

        try:
            # 从文件名提取物种名称
            match = re.search(r'out_(.+?)_msd\.dat', msd_file.name)
            if not match:
                return None
            species = match.group(1)

            content = msd_file.read_text()
            lines = content.strip().split('\n')

            if len(lines) < 3:
                return None

            # 第二行是图例
            legend_line = lines[1].strip().split()
            labels = {
                'time': legend_line[0] if len(legend_line) > 0 else 'fs',
                'x': legend_line[1] if len(legend_line) > 1 else f'{species}_x',
                'y': legend_line[2] if len(legend_line) > 2 else f'{species}_y',
                'z': legend_line[3] if len(legend_line) > 3 else f'{species}_z',
                'total': legend_line[4] if len(legend_line) > 4 else f'{species}_total',
            }

            # 读取数据
            t_values = []
            msd_x = []
            msd_y = []
            msd_z = []
            msd_total = []

            for line in lines[2:]:
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

            # 计算扩散系数
            diffusion_coeff = None
            if len(t_values) > 10:
                mid = len(t_values) // 2
                t_fit = t_values[mid:]
                msd_fit = msd_total[mid:]

                if len(t_fit) > 2:
                    n = len(t_fit)
                    sum_t = sum(t_fit)
                    sum_msd = sum(msd_fit)
                    sum_t_msd = sum(t * m for t, m in zip(t_fit, msd_fit))
                    sum_t2 = sum(t * t for t in t_fit)

                    denom = n * sum_t2 - sum_t * sum_t
                    if abs(denom) > 1e-10:
                        slope = (n * sum_t_msd - sum_t * sum_msd) / denom
                        diffusion_coeff = slope / 6.0 * 1e-4  # Å²/fs -> cm²/s

            # 计算离子电荷、迁移率、电导率
            ion_charge = self._get_ion_charge(species)
            mobility = None
            ionic_conductivity = None

            if diffusion_coeff and ion_charge:
                # 物理常数
                BOLTZMANN = 1.380649e-23  # J/K
                ELEMENTARY_CHARGE = 1.602176634e-19  # C
                temperature = 298.15  # K，默认室温

                # 迁移率 μ = qD / kT (cm²/V·s)
                # D 单位是 cm²/s，需要转换
                mobility = abs(ion_charge) * ELEMENTARY_CHARGE * diffusion_coeff / (BOLTZMANN * temperature)

            return {
                'species': species,
                't_values': t_values,
                'msd_x_values': msd_x,
                'msd_y_values': msd_y,
                'msd_z_values': msd_z,
                'msd_total_values': msd_total,
                'diffusion_coefficient': diffusion_coeff,
                'mobility': mobility,
                'charge': ion_charge,
                'labels': labels,
            }

        except Exception as e:
            self.logger.error(f"简化解析 MSD 文件 {msd_file} 失败: {e}", exc_info=True)
            return None

    def _get_ion_charge(self, species: str) -> int:
        """获取离子电荷"""
        # 常见离子电荷表
        ION_CHARGES = {
            'Li': 1, 'Na': 1, 'K': 1, 'Mg': 2, 'Ca': 2, 'Zn': 2, 'Al': 3,
            'FSI': -1, 'TFSI': -1, 'PF6': -1, 'BF4': -1, 'ClO4': -1, 'DCA': -1,
            'Cl': -1, 'Br': -1, 'I': -1, 'F': -1,
        }

        for ion, charge in ION_CHARGES.items():
            if ion in species:
                return charge

        # 默认返回 None（非离子）
        return None

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

    def _parse_lammps_log(self, log_file: Path) -> Dict[str, Any]:
        """从 LAMMPS 日志文件解析能量、温度、密度、盒子尺寸等"""
        import re

        result = {}

        try:
            content = log_file.read_text()

            # 查找所有 thermo 输出块
            # 通常格式: Step Temp Press E_total KinEng PotEng Density Lx Ly Lz ...
            thermo_pattern = r'Step\s+Temp\s+.*?\n([\s\S]*?)Loop time'
            matches = re.findall(thermo_pattern, content)

            if matches:
                # 尝试解析第一个 thermo 块的第一行（初始状态）
                first_thermo = matches[0].strip().split('\n')
                if first_thermo:
                    first_line = first_thermo[0].split()
                    if len(first_line) >= 7:
                        try:
                            density_val = float(first_line[6])
                            if 0.5 < density_val < 3.0:  # 合理的密度范围
                                result['initial_density'] = density_val
                            # 尝试提取盒子尺寸 (Lx, Ly, Lz 通常在第 7, 8, 9 列)
                            if len(first_line) >= 10:
                                result['initial_box_x'] = float(first_line[7])
                                result['initial_box_y'] = float(first_line[8])
                                result['initial_box_z'] = float(first_line[9])
                        except (ValueError, IndexError):
                            pass

                # 解析最后一个 thermo 块的最后一行（最终状态）
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
                                density_val = float(last_line[6])
                                if 0.5 < density_val < 3.0:
                                    result['final_density'] = density_val
                            # 提取最终盒子尺寸
                            if len(last_line) >= 10:
                                result['box_x'] = float(last_line[7])
                                result['box_y'] = float(last_line[8])
                                result['box_z'] = float(last_line[9])
                        except (ValueError, IndexError):
                            pass

            # 如果没有从 thermo 找到盒子尺寸，尝试从其他地方解析
            if 'box_x' not in result:
                # 尝试匹配 "orthogonal box = ..." 或类似模式
                box_pattern = r'orthogonal box = \(([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\) to \(([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\)'
                box_matches = re.findall(box_pattern, content)
                if box_matches:
                    # 使用最后一个匹配（最终状态）
                    last_box = box_matches[-1]
                    try:
                        result['box_x'] = float(last_box[3]) - float(last_box[0])
                        result['box_y'] = float(last_box[4]) - float(last_box[1])
                        result['box_z'] = float(last_box[5]) - float(last_box[2])
                    except (ValueError, IndexError):
                        pass

                    # 初始盒子
                    if len(box_matches) > 1:
                        first_box = box_matches[0]
                        try:
                            result['initial_box_x'] = float(first_box[3]) - float(first_box[0])
                            result['initial_box_y'] = float(first_box[4]) - float(first_box[1])
                            result['initial_box_z'] = float(first_box[5]) - float(first_box[2])
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

    def _upload_additional_files(self, job_id: int, file_paths: List[Path]) -> List[str]:
        """上传额外的文件（如生成的图片）到对象存储"""
        uploaded_files = []

        try:
            # 获取结果文件前缀
            if self.storage_type == 'cos':
                result_prefix = self.config['cos']['result_prefix']
            else:
                result_prefix = self.config['oss']['result_prefix']

            for file_path in file_paths:
                if not file_path.exists():
                    continue

                file_size = file_path.stat().st_size

                # 构建对象存储 Key
                object_key = f"{result_prefix}{job_id}/{file_path.name}"

                self.logger.info(f"上传额外文件: {file_path.name} ({file_size/1024/1024:.1f}MB)")

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
                self.logger.info(f"✅ 额外上传成功: {file_path.name}")

            return uploaded_files

        except Exception as e:
            self.logger.error(f"上传额外文件失败: {e}", exc_info=True)
            return []

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

