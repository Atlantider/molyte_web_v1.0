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
from typing import Dict, List, Optional, Any, Tuple
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

        # 恢复运行中的任务
        self._recover_running_jobs()

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
                max_jobs = self.config['worker']['max_concurrent_jobs']
                current_jobs = len(self.running_jobs)
                if current_jobs < max_jobs:
                    self._fetch_and_process_new_jobs()
                else:
                    self.logger.info(f"当前运行任务数 ({current_jobs}) 已达到最大并发数 ({max_jobs})，跳过拉取新任务")
                    self.logger.debug(f"运行中的任务: {list(self.running_jobs.keys())}")

                # 等待下一次轮询
                time.sleep(self.config['api']['poll_interval'])

            except KeyboardInterrupt:
                self.logger.info("收到中断信号，正在退出...")
                break
            except Exception as e:
                self.logger.error(f"主循环错误: {e}", exc_info=True)
                time.sleep(10)

    def _recover_running_jobs(self):
        """恢复运行中的任务（Worker 重启后恢复追踪）"""
        self.logger.info("正在恢复运行中的任务...")

        try:
            # 获取当前 Slurm 中运行的任务
            result = subprocess.run(
                ['squeue', '-u', os.environ.get('USER', 'xiaoji'),
                 '-o', '%i %j', '-h'],
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode != 0:
                self.logger.error(f"获取 Slurm 队列失败: {result.stderr}")
                return

            running_slurm_jobs = {}
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split()
                    if len(parts) >= 2:
                        slurm_job_id = parts[0]
                        job_name = parts[1]
                        running_slurm_jobs[job_name] = slurm_job_id

            self.logger.info(f"Slurm 中运行的任务: {running_slurm_jobs}")

            # 从 API 获取状态为 RUNNING 的 MD 任务
            response = requests.get(
                f"{self.api_base_url}/workers/jobs/running",
                headers=self.api_headers,
                timeout=self.config['api']['timeout']
            )

            if response.status_code != 200:
                self.logger.warning(f"获取运行中任务失败: {response.status_code}")
                # 尝试从本地工作目录恢复
                self._recover_from_local_dirs(running_slurm_jobs)
                return

            running_jobs = response.json()
            self.logger.info(f"API 返回运行中的任务: {len(running_jobs)} 个")

            # 用于收集需要检查状态的任务
            jobs_to_check = []

            for job in running_jobs:
                job_id = job['id']
                job_type = job.get('type', 'MD').lower()
                slurm_job_id = job.get('slurm_job_id')
                work_dir = job.get('work_dir')

                if not slurm_job_id or not work_dir:
                    continue

                # 检查 Slurm 任务是否仍在运行
                if slurm_job_id in running_slurm_jobs.values():
                    # 任务仍在运行，恢复追踪
                    self.running_jobs[job_id] = {
                        'type': job_type,
                        'slurm_job_id': slurm_job_id,
                        'work_dir': work_dir,
                        'start_time': time.time(),
                        'resp_cpu_hours': 0.0
                    }
                    self.logger.info(f"恢复追踪任务 {job_id} (Slurm: {slurm_job_id})")
                else:
                    # 任务不在 Slurm 队列中，可能已完成或失败，需要检查状态
                    self.logger.info(f"任务 {job_id} (Slurm: {slurm_job_id}) 不在队列中，检查状态...")
                    jobs_to_check.append({
                        'job_id': job_id,
                        'job_type': job_type,
                        'slurm_job_id': slurm_job_id,
                        'work_dir': work_dir
                    })

            self.logger.info(f"恢复了 {len(self.running_jobs)} 个运行中的任务")

            # 检查并更新已完成/失败的任务
            if jobs_to_check:
                self.logger.info(f"检查 {len(jobs_to_check)} 个可能已完成的任务...")
                self._check_and_update_completed_jobs(jobs_to_check)

        except Exception as e:
            self.logger.error(f"恢复运行中任务失败: {e}", exc_info=True)

    def _check_and_update_completed_jobs(self, jobs_to_check: List[Dict]):
        """检查并更新已完成/失败的任务状态"""
        for job_info in jobs_to_check:
            try:
                job_id = job_info['job_id']
                job_type = job_info['job_type']
                slurm_job_id = job_info['slurm_job_id']
                work_dir = job_info['work_dir']

                # 检查 Slurm 任务状态
                status = self._check_slurm_status(slurm_job_id)

                if status == 'COMPLETED':
                    self.logger.info(f"任务 {job_id} (Slurm: {slurm_job_id}) 已完成，开始处理结果...")
                    # 临时添加到 running_jobs 以便 _handle_job_completion 可以处理
                    temp_job_info = {
                        'type': job_type,
                        'slurm_job_id': slurm_job_id,
                        'work_dir': work_dir,
                        'start_time': time.time(),
                        'resp_cpu_hours': 0.0
                    }
                    self._handle_job_completion(job_id, temp_job_info)

                elif status == 'FAILED':
                    self.logger.error(f"任务 {job_id} (Slurm: {slurm_job_id}) 失败")
                    error_msg = self._get_job_failure_reason(slurm_job_id, work_dir)
                    self._update_job_status(job_id, 'FAILED', job_type, error_message=error_msg)

                elif status == 'CANCELLED':
                    self.logger.info(f"任务 {job_id} (Slurm: {slurm_job_id}) 已取消")
                    self._update_job_status(job_id, 'CANCELLED', job_type, error_message="任务被取消")

                else:
                    self.logger.warning(f"任务 {job_id} (Slurm: {slurm_job_id}) 状态未知: {status}")

            except Exception as e:
                self.logger.error(f"检查任务 {job_info.get('job_id')} 状态失败: {e}", exc_info=True)

    def _recover_from_local_dirs(self, running_slurm_jobs: Dict[str, str]):
        """从本地工作目录恢复任务追踪"""
        md_work_dir = Path(self.config['local']['work_base_path'])

        for job_name, slurm_job_id in running_slurm_jobs.items():
            # 尝试从任务名解析任务 ID
            # 任务名格式: MD-{job_id}-{system_name}
            if job_name.startswith('MD-'):
                parts = job_name.split('-')
                if len(parts) >= 2:
                    try:
                        # 从工作目录名中提取任务 ID
                        for work_dir in md_work_dir.iterdir():
                            if work_dir.is_dir() and job_name in work_dir.name:
                                # 尝试从 slurm 脚本中获取任务 ID
                                slurm_script = work_dir / "run_md.slurm"
                                if slurm_script.exists():
                                    # 从目录名解析任务 ID
                                    # 格式: MD-{date}-{job_id}-EL-...
                                    dir_parts = work_dir.name.split('-')
                                    if len(dir_parts) >= 3:
                                        job_id_str = dir_parts[2]  # 第三部分是任务 ID
                                        job_id = int(job_id_str)

                                        self.running_jobs[job_id] = {
                                            'type': 'md',
                                            'slurm_job_id': slurm_job_id,
                                            'work_dir': str(work_dir),
                                            'start_time': time.time(),
                                            'resp_cpu_hours': 0.0
                                        }
                                        self.logger.info(f"从本地目录恢复任务 {job_id} (Slurm: {slurm_job_id})")
                                        break
                    except (ValueError, IndexError):
                        continue
    
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

            # 获取待处理的后处理任务（包括去溶剂化能计算）
            postprocess_jobs = self._fetch_pending_jobs('postprocess')
            for job in postprocess_jobs:
                if len(self.running_jobs) >= self.config['worker']['max_concurrent_jobs']:
                    break
                self._process_postprocess_job(job)

            # 获取待处理的 Binding 分析任务
            binding_jobs = self._fetch_pending_jobs('binding')
            for job in binding_jobs:
                if len(self.running_jobs) >= self.config['worker']['max_concurrent_jobs']:
                    break
                self._process_binding_job(job)

            # 获取待处理的 Redox 热力学循环任务
            redox_jobs = self._fetch_pending_jobs('redox')
            for job in redox_jobs:
                if len(self.running_jobs) >= self.config['worker']['max_concurrent_jobs']:
                    break
                self._process_redox_job(job)

            # 获取待处理的重组能计算任务
            reorg_jobs = self._fetch_pending_jobs('reorg')
            for job in reorg_jobs:
                if len(self.running_jobs) >= self.config['worker']['max_concurrent_jobs']:
                    break
                self._process_reorg_energy_job(job)

            # 获取待处理的 Cluster 高级计算任务
            cluster_analysis_jobs = self._fetch_pending_jobs('cluster_analysis')
            for job in cluster_analysis_jobs:
                if len(self.running_jobs) >= self.config['worker']['max_concurrent_jobs']:
                    break
                self._process_cluster_analysis_job(job)

            # 检查等待 QC 任务完成的 Cluster 高级计算任务
            self._check_waiting_cluster_analysis_jobs()

            # 检查等待 QC 任务完成的去溶剂化任务（通过 API 调用后端处理）
            self._check_waiting_desolvation_jobs_via_api()

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

            # 2. 导入模块
            from app.workers.molyte_wrapper import MolyteWrapper
            from app.workers.resp_calculator import RESPCalculator
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

            job_data = job['config']
            charge_method = job_data.get("charge_method", "ligpargen")

            # 4. 检查是否需要 RESP 计算
            if charge_method == "resp":
                solvents = job_data.get("solvents", [])
                solvents_needing_resp = wrapper.get_solvents_needing_resp(solvents)

                if solvents_needing_resp:
                    # 需要先运行 RESP 计算
                    self.logger.info(f"MD 任务 {job_id} 需要 RESP 计算: {[s['name'] for s in solvents_needing_resp]}")
                    self._start_resp_calculations(job_id, job, solvents_needing_resp, wrapper)
                    return  # RESP 完成后会继续 MD

            # 5. 直接生成 LAMMPS 输入文件（无需 RESP 或 RESP 已完成）
            self._continue_md_job(job_id, job, wrapper)

        except Exception as e:
            self.logger.error(f"处理 MD 任务 {job_id} 失败: {e}", exc_info=True)
            self._update_job_status(job_id, 'FAILED', 'md', error_message=str(e))

    def _start_resp_calculations(self, job_id: int, job: Dict, solvents: List[Dict], wrapper):
        """
        启动 RESP 电荷计算

        Args:
            job_id: MD 任务 ID
            job: 任务配置
            solvents: 需要计算的溶剂列表
            wrapper: MolyteWrapper 实例
        """
        from app.workers.resp_calculator import RESPCalculator
        import shutil

        try:
            # 创建 RESP 工作目录
            resp_base_dir = Path(self.config['local']['work_base_path']) / f"RESP_{job['config']['name']}"
            resp_base_dir.mkdir(parents=True, exist_ok=True)

            # 初始化 RESP 计算器
            slurm_partition = self.config.get('slurm', {}).get('partition', 'cpu')
            resp_calculator = RESPCalculator(
                charge_save_path=Path(self.config['local']['charge_save_path']),
                slurm_partition=slurm_partition
            )

            resp_jobs = []

            for solvent in solvents:
                name = solvent['name']
                smiles = solvent.get('smiles', '')

                if not smiles:
                    self.logger.warning(f"No SMILES for solvent {name}, skipping RESP")
                    continue

                # 为每个溶剂创建工作目录
                solvent_dir = resp_base_dir / name
                solvent_dir.mkdir(parents=True, exist_ok=True)

                # 运行 LigParGen 生成 PDB 文件
                import subprocess
                ligpargen_cmd = f"{self.config['local']['ligpargen_path']}/ligpargen -s '{smiles}' -n {name} -r MOL -c 0 -o 0 -cgen CM1A"
                result = subprocess.run(
                    ligpargen_cmd, shell=True, cwd=str(solvent_dir),
                    capture_output=True, text=True
                )

                if result.returncode != 0:
                    self.logger.error(f"LigParGen failed for {name}: {result.stderr}")
                    continue

                # 生成 RESP Slurm 脚本
                pdb_file = f"{name}.charmm.pdb"
                script_path = resp_calculator.generate_resp_slurm_script(
                    work_dir=solvent_dir,
                    pdb_file=pdb_file,
                    molecule_name=name,
                    charge=0,
                    spin_multiplicity=1,
                    solvent="water",
                    cpus=16,
                    time_limit_hours=24
                )

                # 提交到 Slurm
                success, slurm_job_id, error = resp_calculator.submit_resp_job(solvent_dir, script_path)

                if success:
                    resp_jobs.append({
                        'molecule_name': name,
                        'slurm_job_id': slurm_job_id,
                        'work_dir': str(solvent_dir),
                        'status': 'RUNNING',
                        'cpu_hours': 0.0
                    })
                    self.logger.info(f"RESP 任务已提交: {name} (Slurm Job: {slurm_job_id})")
                else:
                    self.logger.error(f"Failed to submit RESP job for {name}: {error}")

            if resp_jobs:
                # 将 MD 任务添加到等待 RESP 的队列
                self.running_jobs[job_id] = {
                    'type': 'md_waiting_resp',
                    'job': job,
                    'wrapper': wrapper,
                    'resp_jobs': resp_jobs,
                    'resp_base_dir': str(resp_base_dir),
                    'start_time': time.time(),
                    'total_resp_cpu_hours': 0.0
                }
                self.logger.info(f"MD 任务 {job_id} 正在等待 {len(resp_jobs)} 个 RESP 计算完成")
            else:
                # 没有成功提交任何 RESP 任务，直接继续 MD（使用 LigParGen 电荷）
                self.logger.warning(f"No RESP jobs submitted for MD job {job_id}, falling back to LigParGen")
                job['config']['charge_method'] = 'ligpargen'
                self._continue_md_job(job_id, job, wrapper)

        except Exception as e:
            self.logger.error(f"启动 RESP 计算失败: {e}", exc_info=True)
            self._update_job_status(job_id, 'FAILED', 'md', error_message=f"RESP calculation failed: {e}")

    def _check_resp_jobs(self, job_id: int, job_info: Dict) -> Tuple[bool, bool]:
        """
        检查 RESP 任务状态

        Args:
            job_id: MD 任务 ID
            job_info: 任务信息

        Returns:
            (是否所有 RESP 任务都已完成, 是否有任何 RESP 任务失败)
        """
        from app.workers.resp_calculator import RESPCalculator

        resp_calculator = RESPCalculator(
            charge_save_path=Path(self.config['local']['charge_save_path'])
        )

        all_completed = True
        any_failed = False
        total_cpu_hours = 0.0

        for resp_job in job_info['resp_jobs']:
            if resp_job['status'] in ['COMPLETED', 'FAILED']:
                total_cpu_hours += resp_job.get('cpu_hours', 0.0)
                if resp_job['status'] == 'FAILED':
                    any_failed = True
                continue

            status = resp_calculator.check_job_status(resp_job['slurm_job_id'])

            if status in ['PENDING', 'RUNNING']:
                all_completed = False
            elif status == 'COMPLETED':
                resp_job['status'] = 'COMPLETED'
                resp_job['cpu_hours'] = resp_calculator.get_job_cpu_hours(resp_job['slurm_job_id'])
                total_cpu_hours += resp_job['cpu_hours']
                self.logger.info(f"RESP 计算完成: {resp_job['molecule_name']} (CPU hours: {resp_job['cpu_hours']:.2f})")
            else:
                resp_job['status'] = 'FAILED'
                any_failed = True
                self.logger.error(f"RESP 计算失败: {resp_job['molecule_name']} (status: {status})")

        job_info['total_resp_cpu_hours'] = total_cpu_hours

        if any_failed:
            # 有 RESP 任务失败，需要回退到 LigParGen 电荷
            self.logger.warning(f"Some RESP jobs failed for MD job {job_id}, will use LigParGen for failed molecules")

        return all_completed, any_failed

    def _continue_md_job(self, job_id: int, job: Dict, wrapper):
        """
        继续 MD 任务（RESP 完成后或不需要 RESP）

        Args:
            job_id: MD 任务 ID
            job: 任务配置
            wrapper: MolyteWrapper 实例
        """
        try:
            job_data = job['config']

            # 生成 LAMMPS 输入文件
            result = wrapper.generate_lammps_input(job_data)

            if not result['success']:
                raise Exception(result.get('error', 'Unknown error'))

            work_dir = result['work_dir']

            # 提交到 Slurm
            slurm_result = wrapper.submit_to_slurm(work_dir)

            if not slurm_result['success']:
                raise Exception(slurm_result.get('error', 'Failed to submit to Slurm'))

            slurm_job_id = slurm_result['slurm_job_id']

            # 获取 RESP CPU 核时数
            resp_cpu_hours = 0.0
            if job_id in self.running_jobs and self.running_jobs[job_id].get('type') == 'md_waiting_resp':
                resp_cpu_hours = self.running_jobs[job_id].get('total_resp_cpu_hours', 0.0)

            # 更新任务状态为 RUNNING
            self._update_job_status(
                job_id, 'RUNNING', 'md',
                slurm_job_id=slurm_job_id,
                work_dir=str(work_dir)
            )

            # 更新运行中的任务列表
            self.running_jobs[job_id] = {
                'type': 'md',
                'slurm_job_id': slurm_job_id,
                'work_dir': work_dir,
                'start_time': time.time(),
                'resp_cpu_hours': resp_cpu_hours  # 保存 RESP 核时数
            }

            self.logger.info(f"MD 任务 {job_id} 已提交到 Slurm (Job ID: {slurm_job_id}), RESP CPU hours: {resp_cpu_hours:.2f}")

        except Exception as e:
            self.logger.error(f"继续 MD 任务 {job_id} 失败: {e}", exc_info=True)
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

            # 获取 xyz_content（来自 desolvation 任务创建的 QC 任务）
            xyz_content = config.get('xyz_content', '')

            # 检测是否是 desolvation 任务创建的 QC 任务
            desolvation_job_type = config.get('desolvation_job_type', '')
            is_desolvation_qc = bool(desolvation_job_type)

            # 获取溶剂配置（支持自定义溶剂参数）
            solvent_config = config.get('solvent_config', {})
            if solvent_config:
                solvent_model = solvent_config.get('model', solvent_model)
                solvent_name = solvent_config.get('solvent_name', solvent_name)

            # 3. 创建工作目录
            # 如果是 desolvation 任务的 cluster/ligand 计算，使用 cluster_work 目录
            if is_desolvation_qc:
                cluster_work_base = Path(self.config['local'].get('cluster_work_base_path',
                                         self.config['local']['qc_work_base_path']))
                work_dir = cluster_work_base / f"Cluster-{job_id}-{molecule_name}"
            else:
                qc_work_base = Path(self.config['local']['qc_work_base_path'])
                work_dir = qc_work_base / f"QC-{job_id}-{molecule_name}"
            work_dir.mkdir(parents=True, exist_ok=True)

            # 4. 生成安全的文件名
            safe_name = self._sanitize_filename(molecule_name)

            # 5. 生成 Gaussian 输入文件
            gjf_path = work_dir / f"{safe_name}.gjf"
            self._generate_gaussian_input(
                gjf_path, molecule_name, smiles,
                charge, spin_multiplicity,
                functional, basis_set,
                solvent_model, solvent_name,
                solvent_config,
                xyz_content=xyz_content,  # 传递 xyz 坐标
                nprocs=slurm_cpus  # 使用配置的 CPU 核心数
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

    def _process_postprocess_job(self, job: Dict):
        """处理后处理任务（包括去溶剂化能计算）"""
        job_id = job['id']
        config = job.get('config', {})
        job_type = config.get('job_type', 'UNKNOWN')

        self.logger.info(f"开始处理后处理任务 {job_id} (类型: {job_type})")

        try:
            # 1. 更新任务状态为 QUEUED
            self._update_job_status(job_id, 'QUEUED', 'postprocess')

            # 2. 根据任务类型分发处理
            if job_type == 'DESOLVATION_ENERGY':
                self._process_desolvation_energy_job(job_id, config)
            else:
                raise ValueError(f"Unknown postprocess job type: {job_type}")

        except Exception as e:
            self.logger.error(f"处理后处理任务 {job_id} 失败: {e}", exc_info=True)
            self._update_job_status(job_id, 'FAILED', 'postprocess', error_message=str(e))

    def _process_desolvation_energy_job(self, job_id: int, config: Dict):
        """处理去溶剂化能计算任务（通过 API 调用后端处理）"""
        self.logger.info(f"开始去溶剂化能计算任务 {job_id}")

        try:
            # 调用后端 API 处理任务
            endpoint = f"{self.api_base_url}/workers/jobs/{job_id}/process_desolvation"

            response = requests.post(
                endpoint,
                headers=self.api_headers,
                timeout=300  # 5分钟超时
            )

            if response.status_code == 200:
                result = response.json()
                # API 返回 {"status": "ok", ...} 表示成功
                if result.get('status') == 'ok':
                    self.logger.info(f"去溶剂化能任务 {job_id} 处理成功: {result}")
                    # 注意：任务状态由后端 API 更新，worker 不需要再次更新
                    # self._update_job_status(job_id, 'COMPLETED', 'postprocess')
                else:
                    error_msg = result.get('error', 'Unknown error')
                    self.logger.error(f"去溶剂化能任务 {job_id} 处理失败: {error_msg}")
                    self._update_job_status(job_id, 'FAILED', 'postprocess', error_message=error_msg)
            else:
                error_msg = f"API returned status {response.status_code}: {response.text}"
                self.logger.error(f"去溶剂化能任务 {job_id} API调用失败: {error_msg}")
                self._update_job_status(job_id, 'FAILED', 'postprocess', error_message=error_msg)

        except Exception as e:
            self.logger.error(f"去溶剂化能任务 {job_id} 失败: {e}", exc_info=True)
            self._update_job_status(job_id, 'FAILED', 'postprocess', error_message=str(e))
            raise

    def _process_binding_job(self, job: Dict):
        """
        处理 Binding 分析任务

        简化版 Li-配体 binding energy 计算：
        E_bind_shell = E_cluster - (E_center_ion + Σ n_j × E_ligand_j)
        """
        job_id = job['id']
        config = job.get('config', {})
        md_job_id = config.get('md_job_id')

        self.logger.info(f"开始处理 Binding 分析任务 {job_id} (MD job: {md_job_id})")

        try:
            # 1. 更新状态为 RUNNING
            self._update_job_status(job_id, 'RUNNING', 'binding', progress=0)

            # 2. 获取可用的 cluster 信息
            endpoint = f"{self.api_base_url}/binding/available-clusters/{md_job_id}"
            response = requests.get(
                endpoint,
                headers=self.api_headers,
                timeout=60
            )

            if response.status_code != 200:
                raise Exception(f"获取 cluster 信息失败: {response.status_code}")

            clusters_info = response.json()
            composition_keys = config.get('composition_keys', [])

            if not composition_keys:
                composition_keys = clusters_info.get('composition_keys', [])

            if not composition_keys:
                raise Exception("没有找到可分析的 cluster")

            self.logger.info(f"Binding 任务 {job_id}: 分析 {len(composition_keys)} 个 cluster 类型")

            # 3. 获取已有的 QC 结果
            existing_qc = clusters_info.get('existing_qc_by_type', {})

            # 4. 调用后端 API 执行 binding 计算
            # 这里简化处理：直接调用后端 API 处理
            process_endpoint = f"{self.api_base_url}/binding/jobs/{job_id}/process"

            # 尝试调用后端处理（如果后端有这个接口）
            try:
                process_response = requests.post(
                    process_endpoint,
                    headers=self.api_headers,
                    json={'composition_keys': composition_keys},
                    timeout=300
                )

                if process_response.status_code == 200:
                    result = process_response.json()
                    self.logger.info(f"Binding 任务 {job_id} 后端处理成功")
                    return
                elif process_response.status_code == 404:
                    # 后端没有实现 process 接口，使用本地处理
                    self.logger.info(f"Binding 任务 {job_id}: 使用本地处理逻辑")
                else:
                    self.logger.warning(f"Binding 后端处理返回 {process_response.status_code}")
            except Exception as e:
                self.logger.warning(f"Binding 后端处理调用失败，使用本地逻辑: {e}")

            # 5. 本地处理：从现有 QC 结果计算 binding energy
            self._calculate_binding_from_existing_qc(job_id, md_job_id, composition_keys, config)

        except Exception as e:
            self.logger.error(f"处理 Binding 任务 {job_id} 失败: {e}", exc_info=True)
            self._update_job_status(job_id, 'FAILED', 'binding', error_message=str(e))

    def _calculate_binding_from_existing_qc(self, job_id: int, md_job_id: int,
                                            composition_keys: List[str], config: Dict):
        """从现有 QC 结果计算 binding energy"""
        self.logger.info(f"Binding 任务 {job_id}: 从现有 QC 结果计算")

        try:
            # 获取该 MD 任务下的 QC 结果
            endpoint = f"{self.api_base_url}/qc/jobs"
            params = {
                'md_job_id': md_job_id,
                'status': 'COMPLETED',
                'limit': 1000
            }

            response = requests.get(
                endpoint,
                headers=self.api_headers,
                params=params,
                timeout=60
            )

            if response.status_code != 200:
                raise Exception(f"获取 QC 结果失败: {response.status_code}")

            qc_data = response.json()
            qc_jobs = qc_data.get('items', [])

            if not qc_jobs:
                raise Exception("没有找到已完成的 QC 计算结果")

            # 按分子类型分组
            energy_by_type = {}  # {molecule_name: energy_au}
            cluster_energies = {}  # {composition_key: energy_au}

            for qc in qc_jobs:
                mol_name = qc.get('molecule_name', '')
                mol_type = qc.get('molecule_type', 'unknown')

                # 获取能量（需要从 results 中获取）
                # 这里简化处理，假设 config 中有 energy_au
                qc_config = qc.get('config', {}) or {}
                energy_au = qc_config.get('energy_au')

                if energy_au is None:
                    # 需要获取详细结果
                    detail_endpoint = f"{self.api_base_url}/qc/jobs/{qc['id']}"
                    detail_resp = requests.get(
                        detail_endpoint,
                        headers=self.api_headers,
                        timeout=30
                    )
                    if detail_resp.status_code == 200:
                        detail = detail_resp.json()
                        results = detail.get('results', [])
                        if results:
                            energy_au = results[0].get('energy_au')

                if energy_au is not None:
                    energy_by_type[mol_name] = energy_au

                    # 检查是否是 cluster
                    if 'cluster' in mol_name.lower() or mol_type == 'cluster':
                        # 提取 composition_key
                        for key in composition_keys:
                            if key in mol_name:
                                cluster_energies[key] = energy_au
                                break

            # 计算 binding energy（简化版）
            per_cluster_results = []
            qc_job_ids = []

            # 获取 Li+ 能量
            li_energy = None
            for name, energy in energy_by_type.items():
                if 'li' in name.lower() and 'cluster' not in name.lower():
                    li_energy = energy
                    break

            total_bindings = []

            for comp_key, cluster_energy in cluster_energies.items():
                if li_energy is None:
                    continue

                # 简化计算：只有 cluster 和 Li+
                # 实际应该解析 ligand 类型和数量
                binding_au = cluster_energy - li_energy
                binding_kcal = binding_au * 627.509  # Hartree to kcal/mol

                per_cluster_results.append({
                    'composition_key': comp_key,
                    'cluster_energy_au': cluster_energy,
                    'center_ion_energy_au': li_energy,
                    'binding_energy_au': binding_au,
                    'binding_energy_kcal': binding_kcal,
                    'converged': True,
                    'warnings': ['简化计算：未减去配体能量']
                })
                total_bindings.append(binding_kcal)

            # 生成汇总
            import statistics
            summary = {
                'total_clusters': len(composition_keys),
                'completed_clusters': len(per_cluster_results),
                'failed_clusters': len(composition_keys) - len(per_cluster_results),
                'warnings': []
            }

            if total_bindings:
                summary['mean_total_binding_kcal'] = statistics.mean(total_bindings)
                if len(total_bindings) > 1:
                    summary['std_total_binding_kcal'] = statistics.stdev(total_bindings)

            result = {
                'per_cluster_results': per_cluster_results,
                'summary': summary
            }

            # 更新任务状态为完成
            self._update_binding_result(job_id, result, qc_job_ids)

        except Exception as e:
            self.logger.error(f"计算 binding energy 失败: {e}", exc_info=True)
            raise

    def _update_binding_result(self, job_id: int, result: Dict, qc_job_ids: List[int]):
        """更新 Binding 任务结果"""
        try:
            endpoint = f"{self.api_base_url}/workers/jobs/{job_id}/status"

            data = {
                'status': 'COMPLETED',
                'job_type': 'BINDING',
                'worker_name': self.worker_name,
                'progress': 100.0,
                'result': result,
                'qc_job_ids': qc_job_ids
            }

            response = requests.put(
                endpoint,
                headers=self.api_headers,
                json=data,
                timeout=30
            )

            if response.status_code == 200:
                self.logger.info(f"Binding 任务 {job_id} 完成，已更新结果")
            else:
                self.logger.error(f"更新 Binding 任务 {job_id} 状态失败: {response.status_code}")

        except Exception as e:
            self.logger.error(f"更新 Binding 任务结果失败: {e}")

    # =========================================================================
    # Redox 热力学循环计算
    # =========================================================================

    def _get_qc_job_structure(self, qc_job_id: int) -> Optional[Dict]:
        """
        从 API 获取 QC 任务的结构信息

        返回:
            {
                'xyz_content': str,  # XYZ 格式的结构
                'smiles': str,
                'charge': int,
                'multiplicity': int,
                'molecule_name': str
            }
        """
        try:
            endpoint = f"{self.api_base_url}/qc/jobs/{qc_job_id}"
            response = requests.get(
                endpoint,
                headers=self.api_headers,
                timeout=self.config['api']['timeout']
            )

            if response.status_code == 200:
                job_data = response.json()

                # 尝试从 COS 获取优化后的结构
                xyz_content = None
                result = job_data.get('result', {})
                if result and result.get('optimized_xyz_url'):
                    try:
                        xyz_response = requests.get(result['optimized_xyz_url'], timeout=30)
                        if xyz_response.status_code == 200:
                            xyz_content = xyz_response.text
                    except Exception as e:
                        self.logger.warning(f"获取优化结构失败: {e}")

                # 如果没有优化结构，尝试获取输入结构
                if not xyz_content:
                    input_data = job_data.get('input_data', {})
                    xyz_content = input_data.get('xyz_content', '')

                return {
                    'xyz_content': xyz_content,
                    'smiles': job_data.get('smiles', ''),
                    'charge': job_data.get('charge', 0),
                    'multiplicity': job_data.get('spin_multiplicity', 1),
                    'molecule_name': job_data.get('molecule_name', 'unknown')
                }
            else:
                self.logger.error(f"获取 QC 任务 {qc_job_id} 失败: {response.status_code}")
                return None

        except Exception as e:
            self.logger.error(f"获取 QC 任务结构失败: {e}")
            return None

    def _process_redox_job(self, job: Dict):
        """
        处理 Redox 热力学循环任务

        ⚠️ 高风险警告：
        - 结果对方法/基组/溶剂模型高度敏感
        - 计算量大，经常不收敛
        - 仅供研究参考

        热力学循环：
        ΔG°(sol) = ΔG°(gas) + ΔG_solv(Red) - ΔG_solv(Ox)
        E° = -ΔG°(sol) / nF
        """
        job_id = job['id']
        config = job.get('config', {})

        self.logger.info(f"开始处理 Redox 热力学循环任务 {job_id}")
        self.logger.warning(f"⚠️ Redox 任务 {job_id}: 高风险功能，结果可能不准确")

        try:
            # 1. 更新状态为 RUNNING
            self._update_job_status(job_id, 'RUNNING', 'redox', progress=0)

            species_list = config.get('species_list', [])
            mode = config.get('mode', 'cheap')
            functional = config.get('functional', 'B3LYP')
            basis_set = config.get('basis_set', '6-31G*')
            solvent_model = config.get('solvent_model', 'SMD')
            solvent = config.get('solvent', 'water')
            use_dispersion = config.get('use_dispersion', True)
            li_ref_potential = config.get('li_reference_potential', -3.04)

            if not species_list:
                raise Exception("物种列表为空")

            if len(species_list) > 20:
                raise Exception(f"物种数量 ({len(species_list)}) 超过限制 (20)")

            self.logger.info(f"Redox 任务 {job_id}: {len(species_list)} 个物种, 模式={mode}")

            # 2. 为每个物种创建 QC 任务
            species_results = []
            all_qc_job_ids = []
            global_warnings = []

            for i, species in enumerate(species_list):
                progress = (i / len(species_list)) * 80
                self._update_job_status(job_id, 'RUNNING', 'redox', progress=progress)

                try:
                    species_result = self._calculate_species_redox(
                        job_id, species, mode, functional, basis_set,
                        solvent_model, solvent, use_dispersion, li_ref_potential
                    )
                    species_results.append(species_result)
                    all_qc_job_ids.extend(species_result.get('qc_job_ids', []))
                except Exception as e:
                    self.logger.error(f"Redox 物种 {species.get('name')} 计算失败: {e}")
                    species_results.append({
                        'name': species.get('name', 'unknown'),
                        'redox_type': species.get('redox_type', 'oxidation'),
                        'converged': False,
                        'warnings': [str(e)],
                        'qc_job_ids': []
                    })
                    global_warnings.append(f"{species.get('name')}: {str(e)}")

            # 3. 计算电化学窗口
            ox_potentials = [s['e_vs_li_v'] for s in species_results
                           if s.get('converged') and s.get('redox_type') == 'oxidation' and s.get('e_vs_li_v')]
            red_potentials = [s['e_vs_li_v'] for s in species_results
                            if s.get('converged') and s.get('redox_type') == 'reduction' and s.get('e_vs_li_v')]

            result = {
                'species_results': species_results,
                'oxidation_potentials_v': ox_potentials,
                'reduction_potentials_v': red_potentials,
                'global_warnings': global_warnings,
                'calculation_mode': mode,
                'reference_note': f'所有电位相对于 Li+/Li ({li_ref_potential} V vs SHE)'
            }

            # 计算窗口（如果有足够数据）
            if ox_potentials:
                import statistics
                result['oxidation_limit_v'] = sorted(ox_potentials)[max(0, len(ox_potentials)//20)]  # 5% 分位数
            if red_potentials:
                result['reduction_limit_v'] = sorted(red_potentials, reverse=True)[max(0, len(red_potentials)//20)]  # 95% 分位数
            if result.get('oxidation_limit_v') and result.get('reduction_limit_v'):
                result['electrochemical_window_v'] = result['oxidation_limit_v'] - result['reduction_limit_v']

            # 4. 更新结果
            self._update_redox_result(job_id, result, all_qc_job_ids)

        except Exception as e:
            self.logger.error(f"处理 Redox 任务 {job_id} 失败: {e}", exc_info=True)
            self._update_job_status(job_id, 'FAILED', 'redox', error_message=str(e))

    def _calculate_species_redox(self, job_id: int, species: Dict, mode: str,
                                  functional: str, basis_set: str,
                                  solvent_model: str, solvent: str,
                                  use_dispersion: bool, li_ref_potential: float) -> Dict:
        """计算单个物种的氧化还原电位"""
        name = species.get('name', 'unknown')
        qc_job_id = species.get('qc_job_id')  # 优先使用已有 QC 任务
        smiles = species.get('smiles', '')
        xyz_content = species.get('xyz_content', '')
        charge = species.get('charge', 0)
        multiplicity = species.get('multiplicity', 1)
        redox_type = species.get('redox_type', 'oxidation')

        # 如果提供了 qc_job_id，从已有 QC 任务获取结构
        if qc_job_id:
            self.logger.info(f"从 QC 任务 {qc_job_id} 获取物种 {name} 的结构")
            qc_structure = self._get_qc_job_structure(qc_job_id)
            if qc_structure:
                xyz_content = qc_structure.get('xyz_content', xyz_content)
                smiles = qc_structure.get('smiles', smiles)
                charge = qc_structure.get('charge', charge)
                multiplicity = qc_structure.get('multiplicity', multiplicity)
                self.logger.info(f"成功获取结构: charge={charge}, mult={multiplicity}")
            else:
                self.logger.warning(f"无法从 QC 任务 {qc_job_id} 获取结构，使用默认值")

        self.logger.info(f"计算物种 {name} 的 {redox_type} 电位")

        # 根据模式决定计算步骤
        # cheap: 垂直近似，只做单点
        # standard: gas 优化 + freq，溶剂单点
        # heavy: gas 和溶剂都优化 + freq

        qc_job_ids = []
        warnings = []

        # 物理常数
        HARTREE_TO_KCAL = 627.509
        FARADAY_KCAL = 23.061  # kcal/(mol·V)

        # 确定氧化态/还原态的电荷
        if redox_type == 'oxidation':
            charged_charge = charge + 1
            # 氧化：失去电子，自由基可能增加
            charged_multiplicity = multiplicity + 1 if multiplicity == 1 else multiplicity - 1
        else:  # reduction
            charged_charge = charge - 1
            charged_multiplicity = multiplicity + 1 if multiplicity == 1 else multiplicity - 1

        try:
            if mode == 'cheap':
                # 垂直近似：直接在初始构型上做单点
                # gas 单点
                e_neutral_gas = self._run_redox_qc_calculation(
                    name, smiles, xyz_content, charge, multiplicity,
                    functional, basis_set, 'SP', None, use_dispersion, job_id
                )
                qc_job_ids.append(e_neutral_gas.get('qc_job_id'))

                e_charged_gas = self._run_redox_qc_calculation(
                    name + '_charged', smiles, xyz_content, charged_charge, charged_multiplicity,
                    functional, basis_set, 'SP', None, use_dispersion, job_id
                )
                qc_job_ids.append(e_charged_gas.get('qc_job_id'))

                # 溶剂单点
                e_neutral_sol = self._run_redox_qc_calculation(
                    name + '_sol', smiles, xyz_content, charge, multiplicity,
                    functional, basis_set, 'SP', solvent_model, use_dispersion, job_id,
                    solvent=solvent
                )
                qc_job_ids.append(e_neutral_sol.get('qc_job_id'))

                e_charged_sol = self._run_redox_qc_calculation(
                    name + '_charged_sol', smiles, xyz_content, charged_charge, charged_multiplicity,
                    functional, basis_set, 'SP', solvent_model, use_dispersion, job_id,
                    solvent=solvent
                )
                qc_job_ids.append(e_charged_sol.get('qc_job_id'))

                warnings.append('使用垂直近似（cheap模式），结果仅供参考')

            else:
                # standard/heavy: 需要优化
                # 这里简化处理，实际应该等待 QC 任务完成
                warnings.append(f'{mode}模式需要几何优化，当前使用简化流程')

                # 先做 gas 优化
                opt_result = self._run_redox_qc_calculation(
                    name + '_opt', smiles, xyz_content, charge, multiplicity,
                    functional, basis_set, 'OPT', None, use_dispersion, job_id
                )
                qc_job_ids.append(opt_result.get('qc_job_id'))
                e_neutral_gas = opt_result

                # 其他计算类似 cheap 模式
                e_charged_gas = self._run_redox_qc_calculation(
                    name + '_charged', smiles, xyz_content, charged_charge, charged_multiplicity,
                    functional, basis_set, 'SP', None, use_dispersion, job_id
                )
                qc_job_ids.append(e_charged_gas.get('qc_job_id'))

                e_neutral_sol = self._run_redox_qc_calculation(
                    name + '_sol', smiles, xyz_content, charge, multiplicity,
                    functional, basis_set, 'SP', solvent_model, use_dispersion, job_id,
                    solvent=solvent
                )
                qc_job_ids.append(e_neutral_sol.get('qc_job_id'))

                e_charged_sol = self._run_redox_qc_calculation(
                    name + '_charged_sol', smiles, xyz_content, charged_charge, charged_multiplicity,
                    functional, basis_set, 'SP', solvent_model, use_dispersion, job_id,
                    solvent=solvent
                )
                qc_job_ids.append(e_charged_sol.get('qc_job_id'))

            # 计算自由能变化
            # 注意：这里是简化计算，实际应该等待 QC 完成并获取结果
            # 暂时使用占位值
            dg_gas_kcal = None
            dg_solv_neutral_kcal = None
            dg_solv_charged_kcal = None
            dg_sol_kcal = None
            e_abs_v = None
            e_vs_li_v = None

            # 如果能获取到能量值，计算电位
            if all([e_neutral_gas.get('energy'), e_charged_gas.get('energy'),
                    e_neutral_sol.get('energy'), e_charged_sol.get('energy')]):

                E_n_g = e_neutral_gas['energy']
                E_c_g = e_charged_gas['energy']
                E_n_s = e_neutral_sol['energy']
                E_c_s = e_charged_sol['energy']

                # ΔG(gas) = E(charged) - E(neutral) (in Hartree)
                dg_gas_au = E_c_g - E_n_g
                dg_gas_kcal = dg_gas_au * HARTREE_TO_KCAL

                # ΔG_solv = E(sol) - E(gas)
                dg_solv_neutral_kcal = (E_n_s - E_n_g) * HARTREE_TO_KCAL
                dg_solv_charged_kcal = (E_c_s - E_c_g) * HARTREE_TO_KCAL

                # ΔG(sol) = ΔG(gas) + ΔG_solv(charged) - ΔG_solv(neutral)
                dg_sol_kcal = dg_gas_kcal + dg_solv_charged_kcal - dg_solv_neutral_kcal

                # E° = -ΔG / nF
                e_abs_v = -dg_sol_kcal / FARADAY_KCAL

                # vs Li+/Li
                e_vs_li_v = e_abs_v - li_ref_potential

            return {
                'name': name,
                'redox_type': redox_type,
                'e_neutral_gas': e_neutral_gas.get('energy'),
                'e_charged_gas': e_charged_gas.get('energy'),
                'e_neutral_sol': e_neutral_sol.get('energy'),
                'e_charged_sol': e_charged_sol.get('energy'),
                'dg_gas_kcal': dg_gas_kcal,
                'dg_solv_neutral_kcal': dg_solv_neutral_kcal,
                'dg_solv_charged_kcal': dg_solv_charged_kcal,
                'dg_sol_kcal': dg_sol_kcal,
                'e_abs_v': e_abs_v,
                'e_vs_li_v': e_vs_li_v,
                'converged': e_vs_li_v is not None,
                'warnings': warnings,
                'qc_job_ids': [jid for jid in qc_job_ids if jid]
            }

        except Exception as e:
            self.logger.error(f"物种 {name} 计算失败: {e}")
            return {
                'name': name,
                'redox_type': redox_type,
                'converged': False,
                'warnings': [str(e)],
                'qc_job_ids': [jid for jid in qc_job_ids if jid]
            }

    def _run_redox_qc_calculation(self, name: str, smiles: str, xyz_content: str,
                                   charge: int, multiplicity: int,
                                   functional: str, basis_set: str,
                                   calc_type: str, solvent_model: Optional[str],
                                   use_dispersion: bool, parent_job_id: int,
                                   solvent: str = None,
                                   max_wait_seconds: int = 3600) -> Dict:
        """
        运行单个 QC 计算（用于 Redox）

        返回：{'qc_job_id': int, 'energy': float or None}
        """
        self.logger.info(f"创建 Redox QC 任务: {name}, {calc_type}, solvent={solvent_model}")

        try:
            # 提取基础分子名称（去掉 _charged, _sol, _opt 等后缀）
            # 这样可以实现全局复用：Li-PF6_sol 和 Li-PF6 都用 Li-PF6 作为基础名称
            import re
            base_molecule_name = re.sub(r'(_charged|_sol|_opt|_sp_ox_at_neutral|_sp_neutral_at_ox|_neutral_opt|_oxidized_opt)+$', '', name)

            # 1. 构建 QC 任务创建请求
            config = {
                'job_type': calc_type,  # SP, OPT, or FREQ
                'use_dispersion': use_dispersion,
                'source': 'redox_thermodynamic_cycle',
                'parent_redox_job_id': parent_job_id,
                'original_task_name': name  # 保留原始名称用于日志
            }

            qc_job_data = {
                'smiles': smiles,
                'molecule_name': base_molecule_name,  # 使用基础名称，便于复用
                'molecule_type': 'redox_species',
                'charge': charge,
                'spin_multiplicity': multiplicity,
                'functional': functional,
                'basis_set': basis_set,
                'config': config
            }

            # 溶剂配置：作为顶层字段传递，用于后端查重复用
            if solvent_model and solvent_model != 'gas':
                qc_job_data['solvent_config'] = {
                    'model': solvent_model,
                    'solvent_name': solvent or 'water'
                }

            # 如果有 XYZ 坐标，添加到配置中
            if xyz_content:
                qc_job_data['config']['initial_xyz'] = xyz_content

            # 2. 调用 API 创建 QC 任务
            create_endpoint = f"{self.api_base_url}/qc/jobs"
            response = requests.post(
                create_endpoint,
                headers=self.api_headers,
                json=qc_job_data,
                timeout=60
            )

            if response.status_code not in [200, 201]:
                self.logger.error(f"创建 QC 任务失败: {response.status_code} - {response.text}")
                return {'qc_job_id': None, 'energy': None}

            qc_job = response.json()
            qc_job_id = qc_job['id']
            self.logger.info(f"创建 QC 任务成功: ID={qc_job_id}")

            # 3. 提交任务
            submit_endpoint = f"{self.api_base_url}/qc/jobs/{qc_job_id}/submit"
            submit_response = requests.post(
                submit_endpoint,
                headers=self.api_headers,
                timeout=30
            )

            if submit_response.status_code != 200:
                self.logger.error(f"提交 QC 任务失败: {submit_response.status_code}")
                return {'qc_job_id': qc_job_id, 'energy': None}

            # 4. 等待任务完成（轮询）
            # 注意：这里简化处理，实际的 QC 任务会被其他 Worker 处理
            # 我们只需要等待并检查状态
            start_time = time.time()
            poll_interval = 30  # 30秒轮询一次

            while time.time() - start_time < max_wait_seconds:
                status_endpoint = f"{self.api_base_url}/qc/jobs/{qc_job_id}"
                status_response = requests.get(
                    status_endpoint,
                    headers=self.api_headers,
                    timeout=30
                )

                if status_response.status_code != 200:
                    self.logger.warning(f"获取 QC 任务 {qc_job_id} 状态失败")
                    time.sleep(poll_interval)
                    continue

                job_status = status_response.json()
                current_status = job_status.get('status')

                if current_status == 'COMPLETED':
                    # 获取能量结果
                    results = job_status.get('results', [])
                    energy = None
                    if results:
                        energy = results[0].get('energy_au')

                    self.logger.info(f"QC 任务 {qc_job_id} 完成, 能量={energy}")
                    return {'qc_job_id': qc_job_id, 'energy': energy}

                elif current_status == 'FAILED':
                    error_msg = job_status.get('error_message', 'Unknown error')
                    self.logger.error(f"QC 任务 {qc_job_id} 失败: {error_msg}")
                    return {'qc_job_id': qc_job_id, 'energy': None}

                elif current_status == 'CANCELLED':
                    self.logger.warning(f"QC 任务 {qc_job_id} 被取消")
                    return {'qc_job_id': qc_job_id, 'energy': None}

                self.logger.debug(f"QC 任务 {qc_job_id} 状态: {current_status}, 等待中...")
                time.sleep(poll_interval)

            # 超时
            self.logger.error(f"QC 任务 {qc_job_id} 等待超时 ({max_wait_seconds}秒)")
            return {'qc_job_id': qc_job_id, 'energy': None}

        except Exception as e:
            self.logger.error(f"运行 Redox QC 计算失败: {e}", exc_info=True)
            return {'qc_job_id': None, 'energy': None}

    def _update_redox_result(self, job_id: int, result: Dict, qc_job_ids: List[int]):
        """更新 Redox 任务结果"""
        try:
            endpoint = f"{self.api_base_url}/workers/jobs/{job_id}/status"

            data = {
                'status': 'COMPLETED',
                'job_type': 'REDOX',
                'worker_name': self.worker_name,
                'progress': 100.0,
                'result': result,
                'qc_job_ids': qc_job_ids
            }

            response = requests.put(
                endpoint,
                headers=self.api_headers,
                json=data,
                timeout=30
            )

            if response.status_code == 200:
                self.logger.info(f"Redox 任务 {job_id} 完成，已更新结果")
            else:
                self.logger.error(f"更新 Redox 任务 {job_id} 状态失败: {response.status_code}")

        except Exception as e:
            self.logger.error(f"更新 Redox 任务结果失败: {e}")

    # =========================================================================
    # 重组能计算 (Marcus 理论)
    # =========================================================================

    def _process_reorg_energy_job(self, job: Dict):
        """
        处理重组能计算任务 (Marcus 理论)

        ⚠️ 极高风险警告：
        - 每个物种需要 4 次几何优化 + 4 次单点
        - 计算量极大，经常不收敛
        - 结果对初始构型极其敏感

        4点方案：
        λ_ox = E(q+1, R_q) - E(q+1, R_{q+1})
        λ_red = E(q, R_{q+1}) - E(q, R_q)
        λ_total = (λ_ox + λ_red) / 2
        """
        job_id = job['id']
        config = job.get('config', {})

        self.logger.info(f"开始处理重组能计算任务 {job_id}")
        self.logger.warning(f"⚠️⚠️ 重组能任务 {job_id}: 极高风险功能，计算量大，容易失败")

        try:
            # 1. 更新状态为 RUNNING
            self._update_job_status(job_id, 'RUNNING', 'reorg', progress=0)

            species_list = config.get('species_list', [])
            functional = config.get('functional', 'B3LYP')
            basis_set = config.get('basis_set', '6-31G*')
            use_dispersion = config.get('use_dispersion', True)

            if not species_list:
                raise Exception("物种列表为空")

            if len(species_list) > 5:
                raise Exception(f"物种数量 ({len(species_list)}) 超过限制 (5)")

            self.logger.info(f"重组能任务 {job_id}: {len(species_list)} 个物种")

            # 2. 为每个物种计算重组能
            species_results = []
            all_qc_job_ids = []
            global_warnings = []

            # 物理常数
            HARTREE_TO_EV = 27.2114

            for i, species in enumerate(species_list):
                progress = (i / len(species_list)) * 80
                self._update_job_status(job_id, 'RUNNING', 'reorg', progress=progress)

                try:
                    species_result = self._calculate_species_reorg(
                        job_id, species, functional, basis_set, use_dispersion
                    )
                    species_results.append(species_result)
                    all_qc_job_ids.extend(species_result.get('qc_job_ids', []))
                except Exception as e:
                    self.logger.error(f"重组能物种 {species.get('name')} 计算失败: {e}")
                    species_results.append({
                        'name': species.get('name', 'unknown'),
                        'converged': False,
                        'warnings': [str(e)],
                        'qc_job_ids': []
                    })
                    global_warnings.append(f"{species.get('name')}: {str(e)}")

            # 3. 计算平均值
            lambda_ox_values = [s['lambda_ox_ev'] for s in species_results if s.get('converged') and s.get('lambda_ox_ev')]
            lambda_red_values = [s['lambda_red_ev'] for s in species_results if s.get('converged') and s.get('lambda_red_ev')]
            lambda_total_values = [s['lambda_total_ev'] for s in species_results if s.get('converged') and s.get('lambda_total_ev')]

            import statistics
            result = {
                'species_results': species_results,
                'global_warnings': global_warnings,
                'reference_note': '重组能基于 Marcus 理论 4 点方案计算，结果敏感度高'
            }

            if lambda_ox_values:
                result['lambda_ox_mean_ev'] = statistics.mean(lambda_ox_values)
            if lambda_red_values:
                result['lambda_red_mean_ev'] = statistics.mean(lambda_red_values)
            if lambda_total_values:
                result['lambda_total_mean_ev'] = statistics.mean(lambda_total_values)

            # 4. 更新结果
            self._update_reorg_result(job_id, result, all_qc_job_ids)

        except Exception as e:
            self.logger.error(f"处理重组能任务 {job_id} 失败: {e}", exc_info=True)
            self._update_job_status(job_id, 'FAILED', 'reorg', error_message=str(e))

    def _calculate_species_reorg(self, job_id: int, species: Dict,
                                  functional: str, basis_set: str,
                                  use_dispersion: bool) -> Dict:
        """
        计算单个物种的重组能 (Marcus 理论 4 点方案)

        4点方案：
        1. 优化中性态几何 R_q → 得到 E(q, R_q)
        2. 优化氧化态几何 R_{q+1} → 得到 E(q+1, R_{q+1})
        3. 单点 E(q+1, R_q) - 在中性态几何上计算氧化态能量
        4. 单点 E(q, R_{q+1}) - 在氧化态几何上计算中性态能量

        λ_ox = E(q+1, R_q) - E(q+1, R_{q+1})
        λ_red = E(q, R_{q+1}) - E(q, R_q)
        λ_total = (λ_ox + λ_red) / 2
        """
        name = species.get('name', 'unknown')
        qc_job_id = species.get('qc_job_id')  # 优先使用已有 QC 任务
        smiles = species.get('smiles', '')
        xyz_content = species.get('xyz_content', '')
        charge_neutral = species.get('charge_neutral', 0)
        charge_oxidized = species.get('charge_oxidized', 1)
        mult_neutral = species.get('multiplicity_neutral', 1)
        mult_oxidized = species.get('multiplicity_oxidized', 2)

        # 如果提供了 qc_job_id，从已有 QC 任务获取结构
        if qc_job_id:
            self.logger.info(f"从 QC 任务 {qc_job_id} 获取物种 {name} 的结构")
            qc_structure = self._get_qc_job_structure(qc_job_id)
            if qc_structure:
                xyz_content = qc_structure.get('xyz_content', xyz_content)
                smiles = qc_structure.get('smiles', smiles)
                # 注意：重组能使用用户指定的电荷/多重度，不从 QC 任务获取
                self.logger.info(f"成功获取结构，使用用户指定的电荷: neutral={charge_neutral}, oxidized={charge_oxidized}")
            else:
                self.logger.warning(f"无法从 QC 任务 {qc_job_id} 获取结构，使用默认值")

        self.logger.info(f"计算物种 {name} 的重组能 (4点方案)")
        self.logger.warning(f"⚠️ 重组能计算极其耗时，物种: {name}")

        qc_job_ids = []
        warnings = []

        HARTREE_TO_EV = 27.2114

        try:
            # Step 1: 优化中性态几何 R_q
            self.logger.info(f"{name}: Step 1/4 - 优化中性态几何")
            opt_neutral = self._run_redox_qc_calculation(
                name=f"{name}_neutral_opt",
                smiles=smiles,
                xyz_content=xyz_content,
                charge=charge_neutral,
                multiplicity=mult_neutral,
                functional=functional,
                basis_set=basis_set,
                calc_type='OPT',
                solvent_model=None,
                use_dispersion=use_dispersion,
                parent_job_id=job_id,
                max_wait_seconds=7200  # 2小时超时
            )
            if opt_neutral.get('qc_job_id'):
                qc_job_ids.append(opt_neutral['qc_job_id'])

            E_q_Rq = opt_neutral.get('energy')  # E(q, R_q)

            if E_q_Rq is None:
                warnings.append('中性态优化失败')
                return {
                    'name': name,
                    'converged': False,
                    'warnings': warnings,
                    'qc_job_ids': qc_job_ids
                }

            # Step 2: 优化氧化态几何 R_{q+1}
            self.logger.info(f"{name}: Step 2/4 - 优化氧化态几何")
            opt_oxidized = self._run_redox_qc_calculation(
                name=f"{name}_oxidized_opt",
                smiles=smiles,
                xyz_content=xyz_content,
                charge=charge_oxidized,
                multiplicity=mult_oxidized,
                functional=functional,
                basis_set=basis_set,
                calc_type='OPT',
                solvent_model=None,
                use_dispersion=use_dispersion,
                parent_job_id=job_id,
                max_wait_seconds=7200
            )
            if opt_oxidized.get('qc_job_id'):
                qc_job_ids.append(opt_oxidized['qc_job_id'])

            E_qp1_Rqp1 = opt_oxidized.get('energy')  # E(q+1, R_{q+1})

            if E_qp1_Rqp1 is None:
                warnings.append('氧化态优化失败')
                return {
                    'name': name,
                    'E_q_Rq': E_q_Rq,
                    'converged': False,
                    'warnings': warnings,
                    'qc_job_ids': qc_job_ids
                }

            # Step 3: 单点 E(q+1, R_q) - 在中性态几何上计算氧化态能量
            # 需要获取中性态优化后的几何
            self.logger.info(f"{name}: Step 3/4 - 单点 E(q+1, R_q)")
            sp_oxidized_at_neutral = self._run_redox_qc_calculation(
                name=f"{name}_sp_ox_at_neutral",
                smiles=smiles,
                xyz_content=xyz_content,  # 理想情况下应该用优化后的几何
                charge=charge_oxidized,
                multiplicity=mult_oxidized,
                functional=functional,
                basis_set=basis_set,
                calc_type='SP',
                solvent_model=None,
                use_dispersion=use_dispersion,
                parent_job_id=job_id,
                max_wait_seconds=3600
            )
            if sp_oxidized_at_neutral.get('qc_job_id'):
                qc_job_ids.append(sp_oxidized_at_neutral['qc_job_id'])

            E_qp1_Rq = sp_oxidized_at_neutral.get('energy')  # E(q+1, R_q)

            # Step 4: 单点 E(q, R_{q+1}) - 在氧化态几何上计算中性态能量
            self.logger.info(f"{name}: Step 4/4 - 单点 E(q, R_{q+1})")
            sp_neutral_at_oxidized = self._run_redox_qc_calculation(
                name=f"{name}_sp_neutral_at_ox",
                smiles=smiles,
                xyz_content=xyz_content,  # 理想情况下应该用优化后的几何
                charge=charge_neutral,
                multiplicity=mult_neutral,
                functional=functional,
                basis_set=basis_set,
                calc_type='SP',
                solvent_model=None,
                use_dispersion=use_dispersion,
                parent_job_id=job_id,
                max_wait_seconds=3600
            )
            if sp_neutral_at_oxidized.get('qc_job_id'):
                qc_job_ids.append(sp_neutral_at_oxidized['qc_job_id'])

            E_q_Rqp1 = sp_neutral_at_oxidized.get('energy')  # E(q, R_{q+1})

            # 计算重组能
            lambda_ox_ev = None
            lambda_red_ev = None
            lambda_total_ev = None

            if E_qp1_Rq is not None and E_qp1_Rqp1 is not None:
                lambda_ox_au = E_qp1_Rq - E_qp1_Rqp1
                lambda_ox_ev = lambda_ox_au * HARTREE_TO_EV
                self.logger.info(f"{name}: λ_ox = {lambda_ox_ev:.4f} eV")

            if E_q_Rqp1 is not None and E_q_Rq is not None:
                lambda_red_au = E_q_Rqp1 - E_q_Rq
                lambda_red_ev = lambda_red_au * HARTREE_TO_EV
                self.logger.info(f"{name}: λ_red = {lambda_red_ev:.4f} eV")

            if lambda_ox_ev is not None and lambda_red_ev is not None:
                lambda_total_ev = (lambda_ox_ev + lambda_red_ev) / 2
                self.logger.info(f"{name}: λ_total = {lambda_total_ev:.4f} eV")

            converged = lambda_total_ev is not None

            if not converged:
                warnings.append('部分单点计算失败，无法计算完整重组能')

            return {
                'name': name,
                'E_q_Rq': E_q_Rq,
                'E_qp1_Rqp1': E_qp1_Rqp1,
                'E_qp1_Rq': E_qp1_Rq,
                'E_q_Rqp1': E_q_Rqp1,
                'lambda_ox_ev': lambda_ox_ev,
                'lambda_red_ev': lambda_red_ev,
                'lambda_total_ev': lambda_total_ev,
                'converged': converged,
                'warnings': warnings,
                'qc_job_ids': qc_job_ids
            }

        except Exception as e:
            self.logger.error(f"物种 {name} 重组能计算失败: {e}", exc_info=True)
            return {
                'name': name,
                'converged': False,
                'warnings': [str(e)],
                'qc_job_ids': qc_job_ids
            }

    def _update_reorg_result(self, job_id: int, result: Dict, qc_job_ids: List[int]):
        """更新重组能任务结果"""
        try:
            endpoint = f"{self.api_base_url}/workers/jobs/{job_id}/status"

            data = {
                'status': 'COMPLETED',
                'job_type': 'REORG',
                'worker_name': self.worker_name,
                'progress': 100.0,
                'result': result,
                'qc_job_ids': qc_job_ids
            }

            response = requests.put(
                endpoint,
                headers=self.api_headers,
                json=data,
                timeout=30
            )

            if response.status_code == 200:
                self.logger.info(f"重组能任务 {job_id} 完成，已更新结果")
            else:
                self.logger.error(f"更新重组能任务 {job_id} 状态失败: {response.status_code}")

        except Exception as e:
            self.logger.error(f"更新重组能任务结果失败: {e}")

    def _check_waiting_desolvation_jobs(self):
        """
        检查等待 QC 任务完成的去溶剂化任务

        去溶剂化能计算分两个阶段：
        - Phase 1: 创建 QC 任务（cluster, ligands, cluster_minus）
        - Phase 2: 等待 QC 完成后计算去溶剂化能

        此函数检查处于 POSTPROCESSING 状态且 phase=2 的任务，
        如果其 QC 任务已完成，则继续执行 Phase 2 计算。
        """
        try:
            import sys
            sys.path.insert(0, str(Path(__file__).parent.parent / 'backend'))

            from app.database import SessionLocal
            from app.models.job import PostprocessJob, JobStatus
            from app.tasks.desolvation import run_desolvation_job

            db = SessionLocal()
            try:
                # 查找处于 POSTPROCESSING 状态的去溶剂化任务
                waiting_jobs = db.query(PostprocessJob).filter(
                    PostprocessJob.job_type == 'DESOLVATION_ENERGY',
                    PostprocessJob.status == JobStatus.POSTPROCESSING
                ).all()

                if waiting_jobs:
                    self.logger.info(f"发现 {len(waiting_jobs)} 个等待中的去溶剂化任务")

                for job in waiting_jobs:
                    config = job.config or {}
                    phase = config.get('phase', 1)

                    if phase == 2:
                        # 获取需要等待的 QC 任务（不包括复用的）
                        qc_job_ids = config.get('qc_job_ids', [])
                        reused_qc_job_ids = config.get('reused_qc_job_ids', [])

                        # 需要等待完成的任务 = 所有任务 - 复用的任务
                        # 复用的任务已经是 COMPLETED 状态
                        pending_qc_job_ids = [jid for jid in qc_job_ids if jid not in reused_qc_job_ids]

                        if not qc_job_ids and not reused_qc_job_ids:
                            continue

                        # 检查 QC 任务状态
                        from app.models.qc import QCJob, QCJobStatus
                        all_completed = True
                        any_failed = False

                        for qc_job_id in pending_qc_job_ids:
                            qc_job = db.query(QCJob).filter(QCJob.id == qc_job_id).first()
                            if not qc_job:
                                all_completed = False
                                break
                            if qc_job.status == QCJobStatus.FAILED:
                                any_failed = True
                                break
                            if qc_job.status != QCJobStatus.COMPLETED:
                                all_completed = False
                                break

                        if any_failed:
                            self.logger.warning(f"去溶剂化任务 {job.id} 的 QC 任务失败")
                            job.status = JobStatus.FAILED
                            job.error_message = "One or more QC jobs failed"
                            db.commit()
                            continue

                        if all_completed:
                            reused_count = len(reused_qc_job_ids)
                            self.logger.info(f"去溶剂化任务 {job.id} 的所有 QC 任务已完成 (复用了 {reused_count} 个)，开始 Phase 2 计算")
                            try:
                                result = run_desolvation_job(job, db)
                                if result['success']:
                                    self.logger.info(f"去溶剂化任务 {job.id} Phase 2 完成")
                                else:
                                    self.logger.error(f"去溶剂化任务 {job.id} Phase 2 失败: {result.get('error')}")
                            except Exception as e:
                                self.logger.error(f"去溶剂化任务 {job.id} Phase 2 异常: {e}", exc_info=True)
                                job.status = JobStatus.FAILED
                                job.error_message = str(e)
                                db.commit()
                        else:
                            pending_count = len(pending_qc_job_ids)
                            self.logger.debug(f"去溶剂化任务 {job.id} 仍在等待 {pending_count} 个 QC 任务完成")

            finally:
                db.close()

        except Exception as e:
            self.logger.error(f"检查等待中的去溶剂化任务失败: {e}", exc_info=True)

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
                                   solvent_model: str, solvent_name: str,
                                   solvent_config: Dict = None,
                                   xyz_content: str = None,
                                   retry_level: int = 0,
                                   nprocs: int = 16):
        """
        生成 Gaussian 输入文件

        Args:
            gjf_path: 输出文件路径
            molecule_name: 分子名称
            smiles: SMILES 字符串（用于普通 QC 任务）
            charge: 电荷
            spin_multiplicity: 自旋多重度
            functional: 泛函
            basis_set: 基组
            solvent_model: 溶剂模型 (gas/pcm/smd/custom)
            solvent_name: 溶剂名称
            solvent_config: 自定义溶剂参数
            xyz_content: XYZ 坐标内容（用于 desolvation 创建的 QC 任务）
            retry_level: 重试级别 (0=首次, 1+=重试)
            nprocs: CPU 核心数
        """
        # 导入 QC 安全机制模块
        try:
            from qc_safety import check_structure, get_robust_keywords, SCF_CONVERGENCE_STRATEGIES
            has_safety_module = True
        except ImportError:
            has_safety_module = False
            self.logger.warning("未找到 qc_safety 模块，使用默认设置")

        # 构建计算关键字
        keywords = f"opt freq {functional}/{basis_set}"

        # 添加色散校正（对于DFT方法）
        if functional.upper() not in ["HF"]:
            keywords += " em=gd3bj"

        # 添加 SCF 收敛辅助关键词
        atom_count = 0
        if xyz_content:
            atom_count = len([l for l in xyz_content.strip().split('\n') if l.strip() and len(l.split()) >= 4])

        # 检查是否含有金属（用于选择更稳健的 SCF 策略）
        has_metal = False
        metals = {'Li', 'Na', 'K', 'Mg', 'Ca', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Mn'}
        if xyz_content:
            for line in xyz_content.strip().split('\n'):
                parts = line.split()
                if parts:
                    elem = ''.join(c for c in parts[0] if not c.isdigit())
                    if elem in metals:
                        has_metal = True
                        break

        # 根据体系特点和重试级别添加 SCF 设置
        if has_safety_module:
            scf_keywords = get_robust_keywords(atom_count, has_metal, retry_level)
        else:
            # 默认的稳健设置
            if retry_level == 0:
                if atom_count > 50 or has_metal:
                    scf_keywords = "scf=(maxcycle=200,xqc)"
                else:
                    scf_keywords = "scf=(maxcycle=150)"
            elif retry_level == 1:
                scf_keywords = "scf=(maxcycle=300,xqc,novaracc)"
            elif retry_level == 2:
                scf_keywords = "scf=(maxcycle=300,xqc,vshift=500)"
            else:
                scf_keywords = "scf=(maxcycle=400,qc,nodamping,novaracc)"

        if scf_keywords:
            keywords += f" {scf_keywords}"

        # 自定义溶剂参数块（用于custom模型）
        custom_solvent_params = ""

        # 添加溶剂效应
        if solvent_model and solvent_model.lower() != 'gas':
            if solvent_model.lower() == 'pcm':
                keywords += f" scrf=(pcm,solvent={solvent_name or 'water'})"
            elif solvent_model.lower() == 'smd':
                keywords += f" scrf=(smd,solvent={solvent_name or 'water'})"
            elif solvent_model.lower() == 'custom':
                # 自定义溶剂参数 - 使用 SMD 模型的 Generic 溶剂
                keywords += " scrf=(smd,solvent=generic,read)"
                if solvent_config:
                    eps = solvent_config.get('eps', 78.3553)
                    eps_inf = solvent_config.get('eps_inf', 1.778)
                    hbond_acidity = solvent_config.get('hbond_acidity', 0.82)
                    hbond_basicity = solvent_config.get('hbond_basicity', 0.35)
                    surface_tension = solvent_config.get('surface_tension', 71.99)
                    carbon_aromaticity = solvent_config.get('carbon_aromaticity', 0.0)
                    halogenicity = solvent_config.get('halogenicity', 0.0)

                    custom_solvent_params = (
                        f"eps={eps}\n"
                        f"epsinf={eps_inf}\n"
                        f"HBondAcidity={hbond_acidity}\n"
                        f"HBondBasicity={hbond_basicity}\n"
                        f"SurfaceTension={surface_tension}\n"
                        f"AromaticCarbon={carbon_aromaticity}\n"
                        f"Halogen={halogenicity}\n"
                    )

        safe_name = self._sanitize_filename(molecule_name)

        # 获取坐标：优先使用 xyz_content（来自 desolvation），否则从 SMILES 生成
        coords = None
        if xyz_content:
            # 解析 xyz_content 获取坐标
            coords = self._parse_xyz_content(xyz_content)
            if coords:
                self.logger.info(f"从 xyz_content 解析坐标: {len(coords)} 个原子")

        if not coords:
            # 降级：尝试从 SMILES 生成 3D 坐标
            coords = self._get_3d_coordinates(smiles, molecule_name)

        # 验证和纠正自旋多重度（只对有有效坐标的情况）
        if coords:
            corrected_charge, corrected_spin = self._validate_and_correct_spin(
                smiles, charge, spin_multiplicity, coords
            )

            if corrected_charge != charge or corrected_spin != spin_multiplicity:
                self.logger.warning(
                    f"自旋多重度已纠正: charge {charge}->{corrected_charge}, "
                    f"spin {spin_multiplicity}->{corrected_spin} for {molecule_name}"
                )
                charge = corrected_charge
                spin_multiplicity = corrected_spin

        # 根据 CPU 核心数动态设置内存（大约每核 0.5-1GB）
        mem_gb = max(4, min(nprocs, 32))  # 4-32GB 范围

        # 生成 gjf 内容
        gjf_content = f"""%nprocshared={nprocs}
%mem={mem_gb}GB
%chk={safe_name}.chk
# {keywords}

{molecule_name}

{charge} {spin_multiplicity}
"""

        if coords:
            for atom, x, y, z in coords:
                gjf_content += f" {atom:<2}  {x:>12.8f}  {y:>12.8f}  {z:>12.8f}\n"
        else:
            # 如果无法生成坐标，记录错误并使用注释
            self.logger.error(
                f"无法为 {molecule_name} (SMILES: {smiles}) 生成3D坐标。"
                f"请在 {self.config['local']['initial_salts_path']} 目录下添加对应的 PDB 文件。"
            )
            gjf_content += f"! SMILES: {smiles}\n"
            gjf_content += "! 错误: 无法生成3D坐标\n"
            gjf_content += "! 请手动添加分子坐标或提供 PDB 文件\n"

        gjf_content += "\n"  # 空行结尾

        # 添加自定义溶剂参数块（如果有的话）
        if custom_solvent_params:
            gjf_content += custom_solvent_params
            gjf_content += "\n"  # 参数块后的空行

        with open(gjf_path, 'w') as f:
            f.write(gjf_content)

    def _get_3d_coordinates(self, smiles: str, molecule_name: str):
        """从 SMILES 或 PDB 文件获取 3D 坐标"""
        # 首先尝试从 initial_salts 目录加载 PDB 文件
        initial_salts_path = Path(self.config['local']['initial_salts_path'])

        # 提取基础分子名称（去除计算参数）
        # 例如: "FSI-B3LYP-6-31ppgd,p-pcm-TetraHydroFuran" -> "FSI"
        import re
        base_name = molecule_name.split('-')[0] if '-' in molecule_name else molecule_name
        base_name = base_name.split('_')[0] if '_' in molecule_name else base_name

        # 清理名称（去除 +/- 符号）
        clean_name = molecule_name.replace("+", "").replace("-", "").strip()
        clean_base_name = base_name.replace("+", "").replace("-", "").strip()

        # 尝试多种可能的文件名
        possible_paths = [
            initial_salts_path / f"{base_name}.pdb",  # 基础名称（保留符号）
            initial_salts_path / f"{clean_base_name}.pdb",  # 基础名称（去除符号）
            initial_salts_path / f"{molecule_name}.pdb",  # 完整名称
            initial_salts_path / f"{clean_name}.pdb",  # 完整名称（去除符号）
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

                # 尝试多种方法生成3D坐标
                result = AllChem.EmbedMolecule(mol, randomSeed=42)

                if result == -1:
                    # 第一次失败，尝试使用随机坐标
                    self.logger.warning(f"第一次尝试失败，使用随机坐标方法: {smiles}")
                    result = AllChem.EmbedMolecule(mol, useRandomCoords=True, maxAttempts=100, randomSeed=42)

                if result == -1:
                    # 第二次失败，尝试使用ETKDGv3方法
                    self.logger.warning(f"第二次尝试失败，使用ETKDGv3方法: {smiles}")
                    params = AllChem.ETKDGv3()
                    params.randomSeed = 42
                    result = AllChem.EmbedMolecule(mol, params)

                if result == -1:
                    self.logger.error(f"无法为 {smiles} 生成3D坐标")
                    return None

                # 尝试优化几何结构
                try:
                    AllChem.MMFFOptimizeMolecule(mol)
                except Exception as opt_error:
                    self.logger.warning(f"MMFF优化失败: {opt_error}，使用未优化的坐标")

                coords = []
                conf = mol.GetConformer()
                for i, atom in enumerate(mol.GetAtoms()):
                    pos = conf.GetAtomPosition(i)
                    coords.append((atom.GetSymbol(), pos.x, pos.y, pos.z))

                self.logger.info(f"从 SMILES 生成 3D 坐标 (共 {len(coords)} 个原子)")
                return coords
        except ImportError:
            self.logger.warning("RDKit 未安装，无法从 SMILES 生成坐标")
        except Exception as e:
            self.logger.warning(f"从 SMILES 生成坐标失败: {e}")

        return None

    def _parse_xyz_content(self, xyz_content: str) -> list:
        """
        解析 XYZ 格式的坐标内容

        XYZ 格式:
        第一行: 原子数
        第二行: 注释（可选）
        后续行: 元素符号 x y z

        Args:
            xyz_content: XYZ 格式的字符串

        Returns:
            坐标列表 [(element, x, y, z), ...]
        """
        if not xyz_content or not xyz_content.strip():
            return None

        try:
            lines = xyz_content.strip().split('\n')
            if len(lines) < 3:
                self.logger.warning(f"XYZ 内容行数不足: {len(lines)}")
                return None

            # 第一行是原子数
            try:
                num_atoms = int(lines[0].strip())
            except ValueError:
                self.logger.warning(f"无法解析原子数: {lines[0]}")
                return None

            # 跳过第一行（原子数）和第二行（注释），从第三行开始解析坐标
            coords = []
            for i, line in enumerate(lines[2:], start=1):
                line = line.strip()
                if not line:
                    continue

                parts = line.split()
                if len(parts) < 4:
                    self.logger.warning(f"坐标行格式错误 (行 {i+2}): {line}")
                    continue

                element = parts[0]
                try:
                    x = float(parts[1])
                    y = float(parts[2])
                    z = float(parts[3])
                    coords.append((element, x, y, z))
                except ValueError as e:
                    self.logger.warning(f"坐标解析错误 (行 {i+2}): {e}")
                    continue

            if len(coords) != num_atoms:
                self.logger.warning(f"解析的原子数 ({len(coords)}) 与声明的原子数 ({num_atoms}) 不匹配")

            if coords:
                self.logger.info(f"成功解析 XYZ 坐标: {len(coords)} 个原子")
                return coords
            else:
                return None

        except Exception as e:
            self.logger.error(f"解析 XYZ 内容失败: {e}")
            return None

    def _validate_and_correct_spin(self, smiles: str, charge: int, spin_multiplicity: int,
                                     coords: list = None) -> tuple:
        """
        验证并纠正自旋多重度

        Args:
            smiles: SMILES 字符串
            charge: 电荷
            spin_multiplicity: 自旋多重度
            coords: 坐标列表 [(atom, x, y, z), ...]

        Returns:
            (corrected_charge, corrected_spin)
        """
        try:
            from rdkit import Chem

            # 如果有坐标，从坐标计算电子数
            if coords:
                total_electrons = 0
                for atom_symbol, x, y, z in coords:
                    # 原子序数映射
                    atomic_numbers = {
                        'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'P': 15, 'S': 16, 'Cl': 17,
                        'Li': 3, 'Na': 11, 'K': 19, 'Mg': 12, 'Ca': 20, 'Al': 13, 'Si': 14,
                        'B': 5, 'Br': 35, 'I': 53
                    }
                    total_electrons += atomic_numbers.get(atom_symbol, 0)
                total_electrons -= charge

                # 根据电子数判断自旋多重度
                if total_electrons % 2 == 0:
                    correct_spin = 1
                else:
                    correct_spin = 2

                if correct_spin != spin_multiplicity:
                    self.logger.warning(
                        f"从坐标验证: {len(coords)} 个原子, {total_electrons} 个电子, "
                        f"自旋多重度应为 {correct_spin} (当前为 {spin_multiplicity})"
                    )
                    return charge, correct_spin

            # 如果没有坐标，从 SMILES 计算
            else:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    return charge, spin_multiplicity

                # 添加氢原子以获得完整的电子数
                mol_with_h = Chem.AddHs(mol)

                # 计算总电荷
                total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())

                # 计算总电子数
                total_electrons = sum(atom.GetAtomicNum() for atom in mol_with_h.GetAtoms())
                total_electrons -= total_charge

                # 检查是否有显式自由基
                num_radical_electrons = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())

                if num_radical_electrons > 0:
                    correct_spin = num_radical_electrons + 1
                else:
                    if total_electrons % 2 == 0:
                        correct_spin = 1
                    else:
                        correct_spin = 2

                if correct_spin != spin_multiplicity or total_charge != charge:
                    self.logger.warning(
                        f"从 SMILES 验证: {total_electrons} 个电子, {num_radical_electrons} 个自由基电子, "
                        f"电荷应为 {total_charge} (当前为 {charge}), "
                        f"自旋多重度应为 {correct_spin} (当前为 {spin_multiplicity})"
                    )
                    return total_charge, correct_spin

        except Exception as e:
            self.logger.warning(f"验证自旋多重度失败: {e}")

        return charge, spin_multiplicity

    def _parse_pdb_coordinates(self, pdb_path: Path):
        """解析 PDB 文件坐标"""
        coords = []

        # 尝试不同的编码
        encodings = ['utf-8', 'latin1', 'gbk', 'gb2312']

        for encoding in encodings:
            try:
                with open(pdb_path, 'r', encoding=encoding) as f:
                    for line in f:
                        if line.startswith('ATOM') or line.startswith('HETATM'):
                            try:
                                # PDB 格式: ATOM/HETATM, atom_num, atom_name, ..., x, y, z
                                # 标准 PDB 格式的列位置
                                atom_name = line[12:16].strip()

                                # 提取元素符号
                                # 优先使用 76-78 列的元素符号（如果存在）
                                if len(line) > 76:
                                    element = line[76:78].strip()
                                    if element:
                                        # 去掉电荷符号（如 Mg2+）
                                        element = ''.join(c for c in element if c.isalpha())

                                # 如果没有元素符号列，从原子名称提取
                                if not element or not element.strip():
                                    element = ''.join(c for c in atom_name if c.isalpha())

                                if not element:
                                    element = atom_name[0] if atom_name else 'C'

                                # 提取坐标（标准 PDB 格式）
                                x = float(line[30:38].strip())
                                y = float(line[38:46].strip())
                                z = float(line[46:54].strip())

                                coords.append((element, x, y, z))
                            except (ValueError, IndexError) as e:
                                # 跳过无法解析的行
                                self.logger.debug(f"跳过无法解析的 PDB 行: {line.strip()}, 错误: {e}")
                                continue

                # 如果成功读取到坐标，返回结果
                if coords:
                    self.logger.debug(f"使用 {encoding} 编码成功解析 PDB 文件: {pdb_path}")
                    return coords

            except UnicodeDecodeError:
                # 尝试下一个编码
                continue
            except Exception as e:
                self.logger.warning(f"使用 {encoding} 编码解析 PDB 文件失败: {e}")
                continue

        # 所有编码都失败
        if not coords:
            self.logger.warning(f"无法解析 PDB 文件 {pdb_path}（尝试了所有编码）")

        return coords if coords else None

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

        for job_id, job_info in list(self.running_jobs.items()):
            try:
                job_type = job_info['type']

                # 处理等待 RESP 的 MD 任务
                if job_type == 'md_waiting_resp':
                    self._check_md_waiting_resp(job_id, job_info, completed_jobs)
                    continue

                slurm_job_id = job_info['slurm_job_id']

                # 检查任务是否被用户取消
                actual_type = 'md' if job_type in ['md', 'md_waiting_resp'] else job_type
                if self._check_if_cancelled(job_id, actual_type):
                    self.logger.info(f"任务 {job_id} 已被用户取消，执行 scancel")
                    self._cancel_slurm_job(slurm_job_id)
                    completed_jobs.append(job_id)
                    continue

                # 获取当前任务在 API 中的状态，防止已完成的任务被重新标记为 CANCELLED
                current_api_status = self._get_job_api_status(job_id, actual_type)
                if current_api_status in ['COMPLETED', 'FAILED']:
                    self.logger.info(f"任务 {job_id} 已是 {current_api_status} 状态，跳过状态检查")
                    completed_jobs.append(job_id)
                    continue

                status = self._check_slurm_status(slurm_job_id)

                if status == 'QUEUED':
                    # Slurm 任务还在排队，更新数据库状态为 QUEUED
                    self._update_job_status(job_id, 'QUEUED', actual_type, slurm_job_id=slurm_job_id)

                elif status == 'RUNNING':
                    # Slurm 任务正在运行，更新数据库状态为 RUNNING
                    self._update_job_status(job_id, 'RUNNING', actual_type, slurm_job_id=slurm_job_id)

                elif status == 'COMPLETED':
                    self.logger.info(f"任务 {job_id} (Slurm: {slurm_job_id}) 已完成")
                    self._handle_job_completion(job_id, job_info)
                    completed_jobs.append(job_id)

                elif status == 'FAILED':
                    self.logger.error(f"任务 {job_id} (Slurm: {slurm_job_id}) 失败")
                    # 获取失败原因
                    error_msg = self._get_job_failure_reason(slurm_job_id, job_info.get('work_dir'))

                    # 对于 QC 任务，尝试自动重试
                    if actual_type == 'qc':
                        retry_success = self._try_qc_retry(job_id, job_info, error_msg)
                        if retry_success:
                            # 重试已提交，不移除任务
                            continue

                    self._update_job_status(job_id, 'FAILED', actual_type, error_message=error_msg)
                    completed_jobs.append(job_id)

                elif status == 'CANCELLED':
                    self.logger.info(f"任务 {job_id} (Slurm: {slurm_job_id}) 已取消")
                    self._update_job_status(job_id, 'CANCELLED', actual_type, error_message="任务被取消")
                    completed_jobs.append(job_id)

            except Exception as e:
                self.logger.error(f"检查任务 {job_id} 状态失败: {e}")

        # 移除已完成的任务
        for job_id in completed_jobs:
            if job_id in self.running_jobs:
                del self.running_jobs[job_id]

    def _try_qc_retry(self, job_id: int, job_info: Dict, error_msg: str) -> bool:
        """
        尝试自动重试 QC 任务

        Args:
            job_id: 任务 ID
            job_info: 任务信息
            error_msg: 错误信息

        Returns:
            是否成功提交重试
        """
        try:
            # 导入 QC 安全机制模块
            try:
                from qc_safety import analyze_gaussian_error, generate_retry_gjf, QCRetryManager
                has_safety_module = True
            except ImportError:
                self.logger.warning("未找到 qc_safety 模块，无法自动重试")
                return False

            # 检查重试次数
            retry_count = job_info.get('retry_count', 0)
            max_retries = 3

            if retry_count >= max_retries:
                self.logger.warning(f"QC 任务 {job_id} 已达到最大重试次数 ({max_retries})，不再重试")
                return False

            work_dir = Path(job_info['work_dir'])

            # 读取 Gaussian 日志文件
            log_files = list(work_dir.glob("*_out.log")) + list(work_dir.glob("*.log"))
            if not log_files:
                self.logger.warning(f"QC 任务 {job_id} 未找到日志文件，无法分析错误")
                return False

            log_content = log_files[0].read_text(errors='ignore')

            # 分析错误
            error_analysis = analyze_gaussian_error(log_content)

            if not error_analysis.can_retry:
                self.logger.info(f"QC 任务 {job_id} 错误类型 {error_analysis.error_type.value} 不可自动重试")
                # 更新错误信息，添加建议
                if error_analysis.suggestions:
                    error_msg += f"\n建议: {'; '.join(error_analysis.suggestions)}"
                return False

            self.logger.info(f"QC 任务 {job_id} 准备第 {retry_count + 1} 次重试，错误类型: {error_analysis.error_type.value}")

            # 找到原始 GJF 文件
            gjf_files = list(work_dir.glob("*.gjf"))
            if not gjf_files:
                self.logger.warning(f"QC 任务 {job_id} 未找到 GJF 文件")
                return False

            original_gjf = gjf_files[0]

            # 生成重试 GJF
            new_gjf_content = generate_retry_gjf(original_gjf, retry_count + 1, error_analysis)

            # 备份原始 GJF
            backup_gjf = work_dir / f"{original_gjf.stem}_retry{retry_count}.gjf"
            import shutil
            shutil.copy(original_gjf, backup_gjf)

            # 写入新的 GJF
            with open(original_gjf, 'w') as f:
                f.write(new_gjf_content)

            self.logger.info(f"已更新 GJF 文件，使用策略: {error_analysis.retry_keywords or 'default'}")

            # 重新提交到 Slurm
            slurm_result = self._submit_to_slurm(work_dir)

            if not slurm_result['success']:
                self.logger.error(f"重试提交 Slurm 失败: {slurm_result.get('error')}")
                return False

            new_slurm_job_id = slurm_result['slurm_job_id']

            # 更新任务信息
            job_info['slurm_job_id'] = new_slurm_job_id
            job_info['retry_count'] = retry_count + 1
            job_info['last_error'] = error_analysis.error_type.value

            # 更新 API 状态
            self._update_job_status(
                job_id, 'RUNNING', 'qc',
                slurm_job_id=new_slurm_job_id,
                error_message=f"自动重试第 {retry_count + 1} 次 (原因: {error_analysis.error_message})"
            )

            self.logger.info(f"✅ QC 任务 {job_id} 重试已提交 (Slurm: {new_slurm_job_id})")
            return True

        except Exception as e:
            self.logger.error(f"QC 任务 {job_id} 重试失败: {e}", exc_info=True)
            return False

    def _check_md_waiting_resp(self, job_id: int, job_info: Dict, completed_jobs: List):
        """
        检查等待 RESP 完成的 MD 任务

        Args:
            job_id: MD 任务 ID
            job_info: 任务信息
            completed_jobs: 已完成任务列表
        """
        try:
            # 检查任务是否被取消
            if self._check_if_cancelled(job_id, 'md'):
                self.logger.info(f"MD 任务 {job_id} 已被用户取消，取消所有 RESP 任务")
                for resp_job in job_info.get('resp_jobs', []):
                    if resp_job['status'] == 'RUNNING':
                        self._cancel_slurm_job(resp_job['slurm_job_id'])
                completed_jobs.append(job_id)
                return

            # 检查所有 RESP 任务状态
            all_completed, any_failed = self._check_resp_jobs(job_id, job_info)

            if all_completed:
                self.logger.info(f"MD 任务 {job_id} 的所有 RESP 计算已完成，继续 MD 模拟")

                # 从 job_info 恢复 wrapper 和 job
                job = job_info['job']

                # 如果有 RESP 任务失败，回退到 LigParGen 电荷
                if any_failed:
                    self.logger.warning(f"MD 任务 {job_id} 有 RESP 计算失败，回退到 LigParGen 电荷")
                    job['config']['charge_method'] = 'ligpargen'

                # 重新初始化 wrapper
                from app.workers.molyte_wrapper import MolyteWrapper
                wrapper = MolyteWrapper(
                    work_base_path=Path(self.config['local']['work_base_path']),
                    initial_salts_path=Path(self.config['local']['initial_salts_path']),
                    ligpargen_path=Path(self.config['local']['ligpargen_path']),
                    packmol_path=Path(self.config['local']['packmol_path']),
                    ltemplify_path=Path(self.config['local']['ltemplify_path']),
                    moltemplate_path=Path(self.config['local']['moltemplate_path']),
                    charge_save_path=Path(self.config['local']['charge_save_path']),
                )

                # 继续 MD 任务
                self._continue_md_job(job_id, job, wrapper)
                # 注意：_continue_md_job 会更新 running_jobs[job_id]，所以不需要删除

        except Exception as e:
            self.logger.error(f"检查 RESP 状态失败: {e}", exc_info=True)
            self._update_job_status(job_id, 'FAILED', 'md', error_message=f"RESP check failed: {e}")
            completed_jobs.append(job_id)

    def _get_job_api_status(self, job_id: int, job_type: str) -> str:
        """获取任务在 API 中的当前状态"""
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
                return data.get('status', 'UNKNOWN')

            return 'UNKNOWN'

        except Exception as e:
            self.logger.warning(f"获取任务 {job_id} API 状态失败: {e}")
            return 'UNKNOWN'

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

                                    # 将COS图片路径添加到QC结果中
                                    cos_image_paths = self._map_local_to_cos_paths(job_id, png_files, additional_uploaded)
                                    qc_result.update(cos_image_paths)
                                    self.logger.info(f"图片COS路径: {cos_image_paths}")
                        except Exception as e:
                            self.logger.warning(f"QC 可视化生成失败: {e}")

                    # 上传 QC 结果到数据库
                    try:
                        self._upload_qc_result(job_id, qc_result)
                    except Exception as e:
                        # 上传失败时，不继续标记为COMPLETED
                        self.logger.error(f"❌ QC结果上传失败，任务 {job_id} 已标记为FAILED")
                        return  # 提前返回，不执行后续的COMPLETED状态更新
                else:
                    error_msg = f"未能解析 Gaussian 输出文件，可能计算未正常完成"
                    self.logger.warning(f"⚠️ QC 任务 {job_id} {error_msg}")
                    # 标记为失败
                    self._update_job_status(job_id, 'FAILED', 'qc', error_message=error_msg)
                    return  # 提前返回

            elif job_type == 'md':
                # MD 任务：解析 RDF、MSD 等结果并上传
                md_results = self._parse_md_results(work_dir)
                if md_results:
                    self._upload_md_results(job_id, md_results)
                else:
                    self.logger.warning(f"MD 任务 {job_id} 未能解析结果数据")

            # 3. 获取核时数
            cpu_hours = None
            resp_cpu_hours = None

            if job_type in ['md', 'qc']:
                slurm_job_id = job_info.get('slurm_job_id')
                if slurm_job_id:
                    cpu_hours = self._get_job_cpu_hours(slurm_job_id)
                    self.logger.info(f"任务 {job_id} MD/QC CPU hours: {cpu_hours:.2f}")

                # 获取 RESP 核时数（如果有）
                resp_cpu_hours = job_info.get('resp_cpu_hours', 0.0)
                if resp_cpu_hours > 0:
                    self.logger.info(f"任务 {job_id} RESP CPU hours: {resp_cpu_hours:.2f}")

            # 4. 更新任务状态为 COMPLETED
            self._update_job_status(
                job_id, 'COMPLETED', job_type,
                result_files=uploaded_files,
                progress=100.0,
                cpu_hours=cpu_hours,
                resp_cpu_hours=resp_cpu_hours
            )

            self.logger.info(f"任务 {job_id} 处理完成，上传了 {len(uploaded_files)} 个文件")

        except Exception as e:
            self.logger.error(f"处理任务 {job_id} 完成失败: {e}", exc_info=True)
            self._update_job_status(job_id, 'FAILED', job_info['type'], error_message=str(e))

    def _get_job_cpu_hours(self, slurm_job_id: str) -> float:
        """
        获取 Slurm 任务的 CPU 核时数

        Args:
            slurm_job_id: Slurm 任务 ID

        Returns:
            CPU 核时数
        """
        try:
            # sacct -j JOB_ID -o CPUTimeRAW -n -X
            result = subprocess.run(
                ['sacct', '-j', slurm_job_id, '-o', 'CPUTimeRAW', '-n', '-X'],
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode == 0 and result.stdout.strip():
                cpu_time_seconds = int(result.stdout.strip().split()[0])
                cpu_hours = cpu_time_seconds / 3600.0
                return cpu_hours

            return 0.0

        except Exception as e:
            self.logger.error(f"Failed to get CPU hours for job {slurm_job_id}: {e}")
            return 0.0

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
                self.logger.info(f"✅ QC 结果上传成功: job_id={job_id}, result_id={data.get('result_id')}")
            else:
                error_msg = f"QC 结果上传失败: HTTP {response.status_code} - {response.text}"
                self.logger.error(f"❌ {error_msg}")
                # 更新任务状态为失败，并记录详细错误信息
                self._update_job_status(
                    job_id, 'FAILED', 'qc',
                    error_message=f"结果上传失败: {error_msg}"
                )
                raise Exception(error_msg)

        except requests.exceptions.RequestException as e:
            error_msg = f"网络错误，无法上传QC结果: {str(e)}"
            self.logger.error(f"❌ {error_msg}", exc_info=True)
            # 更新任务状态为失败
            self._update_job_status(
                job_id, 'FAILED', 'qc',
                error_message=error_msg
            )
            raise
        except Exception as e:
            error_msg = f"上传 QC 结果时发生异常: {str(e)}"
            self.logger.error(f"❌ {error_msg}", exc_info=True)
            # 更新任务状态为失败
            self._update_job_status(
                job_id, 'FAILED', 'qc',
                error_message=error_msg
            )
            raise

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

                # 尝试提取盒子体积和离子数量（用于电导率计算）
                box_volume = None
                ion_counts = None
                try:
                    from app.tasks.msd_processor import extract_box_volume_and_ion_counts
                    box_volume, ion_counts = extract_box_volume_and_ion_counts(work_dir)
                    self.logger.info(f"提取到 box_volume={box_volume}, ion_counts={ion_counts}")
                except Exception as e:
                    self.logger.warning(f"无法提取 box_volume 和 ion_counts: {e}")

                for msd_file in msd_files:
                    msd_data = self._parse_lammps_msd_simple(msd_file, box_volume, ion_counts)
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

            # 5. 提取系统结构（混合方案：优先从dump文件，降级到后端服务）
            try:
                system_xyz = self._extract_system_structure_from_dump(work_dir)
                if system_xyz:
                    results['system_xyz_content'] = system_xyz
                    self.logger.info(f"从dump文件提取系统结构成功")
                else:
                    # 降级：尝试使用后端的 get_system_structure
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
                            self.logger.info(f"从后端服务提取系统结构成功: {system_result.get('atom_count', 0)} 原子")
                        else:
                            self.logger.warning(f"后端服务返回错误: {system_result.get('error', 'Unknown error')}")
                    except ImportError as e:
                        self.logger.warning(f"无法导入 get_system_structure: {e}")
                    except Exception as e:
                        self.logger.warning(f"后端服务提取系统结构失败: {e}")
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

    def _parse_lammps_msd_simple(self, msd_file: Path, box_volume: Optional[float] = None,
                                  ion_counts: Optional[Dict[str, int]] = None) -> Optional[Dict[str, Any]]:
        """
        简化版 MSD 解析（当无法导入后端模块时使用）
        文件名格式: out_Li_msd.dat, out_FSI_msd.dat 等

        Args:
            msd_file: MSD 文件路径
            box_volume: 模拟盒子体积 (Å³)，用于计算电导率
            ion_counts: 各离子数量，用于计算电导率
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
                        # 单位转换：Å²/fs -> cm²/s
                        # 1 Å²/fs = 1e-1 cm²/s
                        diffusion_coeff = slope / 6.0 * 1e-1  # Å²/fs -> cm²/s

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
                # D 单位是 cm²/s，需要转换为 m²/s
                D_m2s = diffusion_coeff * 1e-4
                mu_m2Vs = abs(ion_charge) * ELEMENTARY_CHARGE * D_m2s / (BOLTZMANN * temperature)
                # 转换为 cm²/(V·s)
                mobility = mu_m2Vs * 1e4

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
                    # 计算数密度 n (ions/cm³)
                    n = ion_count / (box_volume * 1e-24)

                    # 扩散系数单位转换: cm²/s -> m²/s
                    D_m2s = diffusion_coeff * 1e-4

                    # 计算电导率
                    q = abs(ion_charge) * 1.602176634e-19
                    sigma_Sm = n * 1e6 * q * q * D_m2s / (1.380649e-23 * 298.15)

                    # 转换为 S/cm
                    ionic_conductivity = sigma_Sm / 100

            return {
                'species': species,
                't_values': t_values,
                'msd_x_values': msd_x,
                'msd_y_values': msd_y,
                'msd_z_values': msd_z,
                'msd_total_values': msd_total,
                'diffusion_coefficient': diffusion_coeff,
                'mobility': mobility,
                'ionic_conductivity': ionic_conductivity,
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

    def _extract_system_structure_from_dump(self, work_dir: Path) -> Optional[str]:
        """
        从LAMMPS dump文件提取最后一帧的系统结构，转换为XYZ格式

        这是一个轻量级的替代方案，不依赖轨迹文件的完整解析
        """
        try:
            # 查找dump文件
            dump_files = list(work_dir.glob("*.dump"))
            if not dump_files:
                self.logger.debug(f"未找到dump文件: {work_dir}")
                return None

            dump_file = dump_files[0]
            self.logger.debug(f"从dump文件提取系统结构: {dump_file.name}")

            # 加载atom_mapping以获取类型到元素的映射
            type_to_element = self._load_atom_type_mapping(work_dir)

            # 读取最后一帧
            atoms = {}
            current_frame = []
            in_frame = False
            atom_section = False

            with open(dump_file, 'r') as f:
                for line in f:
                    if line.startswith('ITEM: TIMESTEP'):
                        # 新的一帧开始，保存前一帧
                        if current_frame and atom_section:
                            atoms = self._parse_dump_frame(current_frame, type_to_element)
                        current_frame = [line]
                        in_frame = True
                        atom_section = False
                    elif in_frame:
                        current_frame.append(line)
                        if line.startswith('ITEM: ATOMS'):
                            atom_section = True

                # 处理最后一帧
                if current_frame and atom_section:
                    atoms = self._parse_dump_frame(current_frame, type_to_element)

            if not atoms:
                self.logger.warning("未找到有效的原子数据")
                return None

            # 转换为XYZ格式
            xyz_lines = [str(len(atoms)), "System structure (final frame)"]
            for atom_id in sorted(atoms.keys()):
                atom = atoms[atom_id]
                element = atom.get('element', 'C')
                x = atom.get('x', 0.0)
                y = atom.get('y', 0.0)
                z = atom.get('z', 0.0)
                xyz_lines.append(f"{element} {x:.4f} {y:.4f} {z:.4f}")

            xyz_content = "\n".join(xyz_lines)
            self.logger.info(f"系统结构XYZ生成成功: {len(atoms)} 原子")
            return xyz_content

        except Exception as e:
            self.logger.warning(f"从dump文件提取系统结构失败: {e}")
            return None

    def _load_atom_type_mapping(self, work_dir: Path) -> Dict[int, str]:
        """
        从atom_mapping.json加载原子类型到元素的映射

        Returns:
            type_id -> element_symbol 的字典
        """
        type_to_element = {}

        try:
            atom_mapping_file = work_dir / "atom_mapping.json"
            if atom_mapping_file.exists():
                with open(atom_mapping_file) as f:
                    atom_mapping = json.load(f)

                # 从 atom_types 提取映射
                if 'atom_types' in atom_mapping:
                    for at in atom_mapping['atom_types']:
                        type_id = at.get('type_id')
                        element = at.get('element', 'C')
                        if type_id is not None:
                            type_to_element[type_id] = element

                self.logger.debug(f"加载atom_mapping成功: {len(type_to_element)} 个类型")
        except Exception as e:
            self.logger.warning(f"加载atom_mapping失败: {e}")

        return type_to_element

    def _parse_dump_frame(self, frame_lines: List[str], type_to_element: Dict[int, str] = None) -> Dict[int, Dict[str, Any]]:
        """
        解析LAMMPS dump文件的单个帧

        Args:
            frame_lines: 帧的所有行
            type_to_element: 原子类型到元素的映射字典

        Returns:
            原子字典，key为原子ID，value为原子信息（element, x, y, z）
        """
        if type_to_element is None:
            type_to_element = {}

        atoms = {}
        atom_section = False
        atom_format = None  # 'id type x y z' 或其他格式

        for line in frame_lines:
            if line.startswith('ITEM: ATOMS'):
                atom_section = True
                # 解析原子列格式
                parts = line.split()
                if len(parts) > 2:
                    atom_format = parts[2:]  # ['id', 'type', 'x', 'y', 'z', ...]
                continue

            if atom_section and not line.startswith('ITEM:'):
                parts = line.split()
                if len(parts) < 5:
                    continue

                try:
                    # 标准格式: id type x y z
                    atom_id = int(parts[0])
                    atom_type = int(parts[1])
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])

                    # 从映射中获取元素，如果没有则使用默认推断
                    if atom_type in type_to_element:
                        element = type_to_element[atom_type]
                    else:
                        element = self._get_element_from_type(atom_type)

                    atoms[atom_id] = {
                        'element': element,
                        'x': x,
                        'y': y,
                        'z': z,
                        'type': atom_type,
                    }
                except (ValueError, IndexError):
                    continue

        return atoms

    def _get_element_from_type(self, atom_type: int) -> str:
        """
        根据LAMMPS原子类型推断元素符号

        这是一个简化版本，假设：
        - 类型1-2: Li (锂)
        - 类型3-4: F (氟)
        - 类型5-6: S (硫)
        - 类型7-8: C (碳)
        - 类型9-10: O (氧)
        - 类型11-12: H (氢)
        - 其他: C (默认碳)
        """
        type_to_element = {
            1: 'Li', 2: 'Li',
            3: 'F', 4: 'F',
            5: 'S', 6: 'S',
            7: 'C', 8: 'C',
            9: 'O', 10: 'O',
            11: 'H', 12: 'H',
        }
        return type_to_element.get(atom_type, 'C')

    def _parse_lammps_log(self, log_file: Path) -> Dict[str, Any]:
        """从 LAMMPS 日志文件解析能量、温度、密度、盒子尺寸等

        使用与后端API相同的解析逻辑，确保准确提取初始和最终状态
        """
        result = {}

        try:
            with open(log_file, 'r') as f:
                lines = f.readlines()

            # 辅助函数：从一行数据中提取密度和盒子尺寸
            # LAMMPS thermo 输出格式: Step CPU CPULeft Temp Density Lx Ly Lz ...
            # 列索引:                  0    1   2       3    4       5  6  7
            def parse_thermo_line(parts):
                density_val = None
                box_vals = []
                try:
                    # 根据列位置提取密度（第5列，索引4）
                    if len(parts) > 4:
                        density_val = float(parts[4])

                    # 根据列位置提取盒子尺寸（第6-8列，索引5-7）
                    if len(parts) > 7:
                        box_vals = [float(parts[5]), float(parts[6]), float(parts[7])]
                except (ValueError, IndexError):
                    # 如果按列位置提取失败，返回 None
                    density_val = None
                    box_vals = []

                return density_val, box_vals

            # 判断是否为有效数据行
            skip_prefixes = ('Loop', 'Performance', 'MPI', 'Section',
                             'Total', 'Pair', 'Bond', 'Kspace', 'Neigh',
                             'Comm', 'Output', 'Modify', 'Other', 'Nlocal',
                             'Nghost', 'Neighs', 'Ave', 'Neighbor', 'Dangerous',
                             'System', 'PPPM', 'WARNING', 'G vector', 'grid',
                             'stencil', 'estimated', 'using', '3d grid', '-',
                             'Step', 'Per', 'run', 'thermo', 'fix', 'dump',
                             'Memory', 'units', 'atom', 'pair', 'kspace',
                             'Lattice', 'Created', 'Reading', 'Replicate')

            def is_data_line(line_str):
                line_str = line_str.strip()
                if not line_str:
                    return False
                if line_str.startswith(skip_prefixes):
                    return False
                parts = line_str.split()
                if len(parts) < 8:
                    return False
                try:
                    int(parts[0])  # 第一列应该是时间步（整数）
                    return True
                except ValueError:
                    return False

            # 找到第一个 "Step" 表头行的位置（能量最小化阶段的开始，即真正的初始状态）
            first_step_header_idx = -1
            for i, line in enumerate(lines):
                if line.strip().startswith('Step'):
                    first_step_header_idx = i
                    break  # 只找第一个

            # 找到第一个 "Step" 表头之后的第一行数据（真正的初始状态）
            first_data_line = None
            if first_step_header_idx >= 0:
                for line in lines[first_step_header_idx + 1:]:
                    if is_data_line(line):
                        first_data_line = line.strip().split()
                        break

            # 找到最后一行数据（最终状态）
            last_data_line = None
            for line in reversed(lines[-300:]):
                if is_data_line(line):
                    last_data_line = line.strip().split()
                    break

            # 解析初始状态
            if first_data_line:
                init_density, init_box = parse_thermo_line(first_data_line)
                if init_density:
                    result['initial_density'] = round(init_density, 4)
                if len(init_box) >= 3:
                    result['initial_box_x'] = round(init_box[0], 2)
                    result['initial_box_y'] = round(init_box[1], 2)
                    result['initial_box_z'] = round(init_box[2], 2)

            # 解析最终状态
            if last_data_line:
                final_density, final_box = parse_thermo_line(last_data_line)
                if final_density:
                    result['final_density'] = round(final_density, 4)
                if len(final_box) >= 3:
                    result['box_x'] = round(final_box[0], 2)
                    result['box_y'] = round(final_box[1], 2)
                    result['box_z'] = round(final_box[2], 2)

                # 提取能量和温度（从最后一行）
                # LAMMPS thermo 输出格式: Step CPU CPULeft Temp Density Lx Ly Lz TotEng KinEng PotEng ...
                # 列索引:                  0    1   2       3    4       5  6  7  8      9       10
                if len(last_data_line) >= 11:
                    try:
                        result['final_temperature'] = float(last_data_line[3])   # Temp
                        result['total_energy'] = float(last_data_line[8])        # TotEng
                        result['kinetic_energy'] = float(last_data_line[9])      # KinEng
                        result['potential_energy'] = float(last_data_line[10])   # PotEng
                    except (ValueError, IndexError):
                        pass

            self.logger.info(
                f"解析LAMMPS日志: "
                f"初始密度={result.get('initial_density')}, "
                f"最终密度={result.get('final_density')}, "
                f"初始盒子={result.get('initial_box_x')}×{result.get('initial_box_y')}×{result.get('initial_box_z')}, "
                f"最终盒子={result.get('box_x')}×{result.get('box_y')}×{result.get('box_z')}"
            )

        except Exception as e:
            self.logger.error(f"解析 LAMMPS 日志失败: {e}", exc_info=True)

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

    def _map_local_to_cos_paths(self, job_id: int, local_files: List[Path], cos_keys: List[str]) -> Dict[str, str]:
        """将本地图片文件映射到COS路径，用于数据库存储"""
        result = {}

        # 创建文件名到COS key的映射
        filename_to_cos = {}
        for cos_key in cos_keys:
            filename = cos_key.split('/')[-1]  # 从 results/105/ESP.png 提取 ESP.png
            filename_to_cos[filename] = cos_key

        # 映射特定的图片类型
        for local_file in local_files:
            filename = local_file.name
            cos_key = filename_to_cos.get(filename)

            if cos_key:
                if filename.upper().startswith('ESP'):
                    result['esp_image_path'] = cos_key
                elif filename.upper().startswith('HOMO'):
                    result['homo_image_path'] = cos_key
                elif filename.upper().startswith('LUMO'):
                    result['lumo_image_path'] = cos_key

        return result

    def _update_job_status(
        self,
        job_id: int,
        status: str,
        job_type: str,
        slurm_job_id: Optional[str] = None,
        work_dir: Optional[str] = None,
        error_message: Optional[str] = None,
        result_files: Optional[List[str]] = None,
        progress: Optional[float] = None,
        cpu_hours: Optional[float] = None,
        resp_cpu_hours: Optional[float] = None,
        qc_job_ids: Optional[List[int]] = None,
        result: Optional[Dict] = None
    ):
        """更新任务状态到阿里云"""
        try:
            endpoint = f"{self.api_base_url}/workers/jobs/{job_id}/status"

            # 确保状态值有效（包括 Cluster Analysis 特有的状态）
            valid_statuses = ["CREATED", "SUBMITTED", "QUEUED", "RUNNING", "POSTPROCESSING",
                            "COMPLETED", "FAILED", "CANCELLED", "WAITING_QC", "CALCULATING"]
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
            if cpu_hours is not None:
                data['cpu_hours'] = cpu_hours
            if resp_cpu_hours is not None:
                data['resp_cpu_hours'] = resp_cpu_hours
            if qc_job_ids is not None:
                data['qc_job_ids'] = qc_job_ids
            if result is not None:
                data['result'] = result

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

    def _process_cluster_analysis_job(self, job: Dict):
        """
        处理 Cluster 高级计算任务

        统一处理多种计算类型：
        - BINDING_TOTAL: 总 Binding Energy
        - BINDING_PAIRWISE: 分子-Li Binding
        - DESOLVATION_STEPWISE: 逐级去溶剂化
        - DESOLVATION_FULL: 完全去溶剂化
        - REDOX: 氧化还原电位
        - REORGANIZATION: Marcus 重组能
        """
        job_id = job['id']
        config = job.get('config', {})
        md_job_id = config.get('md_job_id')
        calc_types = config.get('calc_types', [])
        qc_task_plan = config.get('qc_task_plan', {})

        self.logger.info(f"开始处理 Cluster 高级计算任务 {job_id} (MD job: {md_job_id})")
        self.logger.info(f"计算类型: {calc_types}")

        try:
            # 1. 更新状态为 RUNNING
            self._update_job_status(job_id, 'RUNNING', 'cluster_analysis', progress=0)

            # 2. 获取需要创建的 QC 任务
            planned_tasks = qc_task_plan.get('planned_qc_tasks', [])
            new_tasks = [t for t in planned_tasks if t.get('status') == 'new']
            reused_tasks = [t for t in planned_tasks if t.get('status') == 'reused']

            self.logger.info(f"Cluster 任务 {job_id}: 需要新建 {len(new_tasks)} 个 QC 任务, "
                           f"复用 {len(reused_tasks)} 个已有任务")

            # 3. 创建新的 QC 任务
            created_qc_job_ids = []
            for i, task in enumerate(new_tasks):
                try:
                    qc_job_id = self._create_qc_job_for_cluster_analysis(
                        job_id, md_job_id, task, config.get('qc_config', {})
                    )
                    if qc_job_id:
                        created_qc_job_ids.append(qc_job_id)
                        progress = (i + 1) / len(new_tasks) * 30  # 创建 QC 任务占 30%
                        self._update_job_status(job_id, 'RUNNING', 'cluster_analysis',
                                              progress=progress)
                except Exception as e:
                    self.logger.error(f"创建 QC 任务失败: {e}")

            # 4. 更新状态为 WAITING_QC
            self._update_job_status(
                job_id, 'WAITING_QC', 'cluster_analysis',
                progress=30,
                qc_job_ids=created_qc_job_ids
            )

            self.logger.info(f"Cluster 任务 {job_id}: 已创建 {len(created_qc_job_ids)} 个 QC 任务，"
                           f"等待 QC 计算完成")

        except Exception as e:
            self.logger.error(f"处理 Cluster 高级计算任务 {job_id} 失败: {e}", exc_info=True)
            self._update_job_status(job_id, 'FAILED', 'cluster_analysis', error_message=str(e))

    def _check_waiting_cluster_analysis_jobs(self):
        """检查等待 QC 任务完成的 Cluster 高级计算任务"""
        try:
            # 获取 WAITING_QC 状态的任务
            endpoint = f"{self.api_base_url}/workers/jobs/pending"
            params = {
                'job_type': 'CLUSTER_ANALYSIS',
                'status': 'WAITING_QC',
                'limit': 10
            }

            response = requests.get(
                endpoint,
                headers=self.api_headers,
                params=params,
                timeout=self.config['api']['timeout']
            )

            if response.status_code != 200:
                return

            waiting_jobs = response.json()

            for job in waiting_jobs:
                job_id = job['id']

                # 首先检查是否有缺少的 QC 任务需要创建
                missing_created = self._create_missing_cluster_qc_jobs(job)
                if missing_created > 0:
                    self.logger.info(f"Cluster 任务 {job_id}: 补充创建了 {missing_created} 个缺少的 QC 任务")

                # 检查 QC 任务状态
                qc_status = self._check_cluster_analysis_qc_status(job_id)

                if qc_status.get('all_completed'):
                    self.logger.info(f"Cluster 任务 {job_id}: 所有 QC 任务已完成，开始计算结果")
                    self._calculate_cluster_analysis_results(job)
                elif qc_status.get('failed', 0) > 0:
                    self.logger.error(f"Cluster 任务 {job_id}: 有 QC 任务失败")
                    self._update_job_status(
                        job_id, 'FAILED', 'cluster_analysis',
                        error_message=f"{qc_status.get('failed')} 个 QC 任务失败"
                    )
                else:
                    # 更新进度
                    total = qc_status.get('total_qc_jobs', 1)
                    completed = qc_status.get('completed', 0)
                    progress = 30 + (completed / total) * 60  # QC 完成占 30%-90%
                    self._update_job_status(job_id, 'WAITING_QC', 'cluster_analysis', progress=progress)

        except Exception as e:
            self.logger.error(f"检查 WAITING_QC 任务失败: {e}", exc_info=True)

    def _create_missing_cluster_qc_jobs(self, job: Dict) -> int:
        """
        检查并创建缺少的 QC 任务

        当 Worker 在创建 QC 任务时中断（如网络错误、超时等），
        可能只创建了部分任务。此方法检查计划任务和已创建任务的差异，
        并补充创建缺少的任务。

        Returns:
            创建的任务数量
        """
        job_id = job['id']
        md_job_id = job.get('md_job_id')
        config = job.get('config', {})
        qc_task_plan = job.get('qc_task_plan', {})

        if not qc_task_plan:
            return 0

        planned_tasks = qc_task_plan.get('planned_qc_tasks', [])
        if not planned_tasks:
            return 0

        try:
            # 获取已创建的 QC 任务类型
            endpoint = f"{self.api_base_url}/cluster-analysis/jobs/{job_id}/qc-jobs"
            response = requests.get(
                endpoint,
                headers=self.api_headers,
                timeout=self.config['api']['timeout']
            )

            if response.status_code != 200:
                self.logger.warning(f"获取 Cluster 任务 {job_id} 的 QC 任务列表失败: {response.status_code}")
                return 0

            existing_qc_jobs = response.json()
            existing_types = set(qc.get('task_type') for qc in existing_qc_jobs if qc.get('task_type'))

            # 找出缺少的任务
            missing_tasks = [
                t for t in planned_tasks
                if t.get('status') == 'new' and t.get('task_type') not in existing_types
            ]

            if not missing_tasks:
                return 0

            self.logger.info(f"Cluster 任务 {job_id}: 发现 {len(missing_tasks)} 个缺少的 QC 任务，开始补充创建")

            # 创建缺少的任务
            created_count = 0
            for task in missing_tasks:
                try:
                    qc_job_id = self._create_qc_job_for_cluster_analysis(
                        job_id, md_job_id, task, config.get('qc_config', {})
                    )
                    if qc_job_id:
                        created_count += 1
                except Exception as e:
                    self.logger.error(f"补充创建 QC 任务失败 (type={task.get('task_type')}): {e}")

            return created_count

        except Exception as e:
            self.logger.error(f"检查缺少的 QC 任务失败: {e}", exc_info=True)
            return 0

    def _check_waiting_desolvation_jobs_via_api(self):
        """
        通过 API 检查等待 QC 任务完成的去溶剂化任务

        对于处于 POSTPROCESSING 状态（Phase 2 等待中）的任务，
        检查其 QC 任务是否已完成，如果完成则继续计算去溶剂化能。
        """
        try:
            # 1. 获取等待中的去溶剂化任务
            endpoint = f"{self.api_base_url}/workers/jobs/waiting_desolvation"
            response = requests.get(
                endpoint,
                headers=self.api_headers,
                timeout=self.config['api']['timeout']
            )

            if response.status_code != 200:
                self.logger.warning(f"获取等待中的去溶剂化任务失败: {response.status_code}")
                return

            result = response.json()
            waiting_jobs = result.get('jobs', [])

            if not waiting_jobs:
                return

            self.logger.info(f"发现 {len(waiting_jobs)} 个等待中的去溶剂化任务")

            # 2. 对每个任务调用检查接口
            for job in waiting_jobs:
                job_id = job['id']
                try:
                    check_endpoint = f"{self.api_base_url}/workers/jobs/{job_id}/check_waiting_desolvation"
                    check_response = requests.post(
                        check_endpoint,
                        headers=self.api_headers,
                        timeout=120  # 计算可能需要较长时间
                    )

                    if check_response.status_code == 200:
                        check_result = check_response.json()
                        status = check_result.get('status')
                        if status == 'completed':
                            self.logger.info(f"去溶剂化任务 {job_id}: 计算完成")
                        elif status == 'waiting':
                            self.logger.debug(f"去溶剂化任务 {job_id}: 仍在等待 QC 完成")
                        else:
                            self.logger.info(f"去溶剂化任务 {job_id}: 状态 {status}")
                    else:
                        self.logger.warning(f"检查去溶剂化任务 {job_id} 失败: {check_response.status_code} - {check_response.text[:200]}")

                except Exception as e:
                    self.logger.error(f"检查去溶剂化任务 {job_id} 异常: {e}")

        except Exception as e:
            self.logger.error(f"检查等待中的去溶剂化任务失败: {e}", exc_info=True)

    def _check_cluster_analysis_qc_status(self, job_id: int) -> Dict:
        """检查 Cluster 高级计算的 QC 任务状态"""
        try:
            endpoint = f"{self.api_base_url}/cluster-analysis/jobs/{job_id}/qc-status"
            response = requests.get(
                endpoint,
                headers=self.api_headers,
                timeout=self.config['api']['timeout']
            )

            if response.status_code == 200:
                return response.json()
            else:
                self.logger.warning(f"获取 QC 状态失败: {response.status_code}")
                return {}

        except Exception as e:
            self.logger.error(f"检查 QC 状态异常: {e}")
            return {}

    def _calculate_cluster_analysis_results(self, job: Dict):
        """计算 Cluster 高级计算结果"""
        job_id = job['id']

        try:
            # 更新状态为 CALCULATING
            self._update_job_status(job_id, 'CALCULATING', 'cluster_analysis', progress=90)

            # 调用 API 计算结果
            endpoint = f"{self.api_base_url}/cluster-analysis/jobs/{job_id}/calculate"
            response = requests.post(
                endpoint,
                headers=self.api_headers,
                timeout=120  # 计算可能需要较长时间
            )

            if response.status_code == 200:
                result = response.json()
                self.logger.info(f"Cluster 任务 {job_id}: 结果计算完成")
                self._update_job_status(
                    job_id, 'COMPLETED', 'cluster_analysis',
                    progress=100,
                    result=result.get('results', {})
                )
            else:
                self.logger.error(f"计算结果失败: {response.status_code} - {response.text}")
                self._update_job_status(
                    job_id, 'FAILED', 'cluster_analysis',
                    error_message=f"计算结果失败: {response.text[:200]}"
                )

        except Exception as e:
            self.logger.error(f"计算 Cluster 任务 {job_id} 结果失败: {e}", exc_info=True)
            self._update_job_status(job_id, 'FAILED', 'cluster_analysis', error_message=str(e))

    def _create_qc_job_for_cluster_analysis(self, cluster_job_id: int, md_job_id: int,
                                            task: Dict, qc_config: Dict) -> Optional[int]:
        """为 Cluster 高级计算创建 QC 任务"""
        task_type = task.get('task_type', '')
        smiles = task.get('smiles') or None  # 空字符串转为 None，让后端使用 XYZ 结构
        structure_id = task.get('structure_id')
        charge = task.get('charge', 0)
        multiplicity = task.get('multiplicity', 1)

        # 根据 task_type 决定溶剂配置
        # 规则：
        # 1. *_gas 后缀 → 气相（无溶剂），用于 redox 热力学循环中的气相计算
        # 2. *_sol 后缀 → 使用配置的溶剂模型，用于 redox 热力学循环中的溶剂相计算
        # 3. 其他任务（cluster, ligand, dimer, ion, cluster_minus, reorg_* 等）→ 使用默认配置
        if task_type.endswith('_gas'):
            # 气相计算，不使用溶剂
            solvent_model = 'gas'
            solvent_name = None
        elif task_type.endswith('_sol'):
            # 溶剂相计算，使用配置的溶剂模型
            solvent_model = qc_config.get('solvent_model') or 'pcm'
            solvent_name = qc_config.get('solvent') or qc_config.get('solvent_name')
        else:
            # 其他任务使用默认配置
            solvent_model = qc_config.get('solvent_model') or 'gas'
            solvent_name = qc_config.get('solvent') or qc_config.get('solvent_name')

        # 构建溶剂配置对象（符合 QCJobCreate schema）
        solvent_config = None
        if solvent_model and solvent_model.lower() != 'gas':
            solvent_config = {
                'model': solvent_model.lower(),
                'solvent_name': solvent_name
            }

        # 从 task_type 提取分子名称，用于全局复用匹配
        # task_type 格式：
        #   - ion -> Li+
        #   - ligand_{name} -> {name} (如 ligand_PF6 -> PF6)
        #   - dimer_{name} -> {name} (如 dimer_Li-PF6 -> Li-PF6)
        #   - redox_mol_{name}_{state} -> {name} (如 redox_mol_PF6_neutral_gas -> PF6)
        #   - redox_dimer_{name}_{state} -> {name} (如 redox_dimer_Li-PF6_neutral_sol -> Li-PF6)
        #   - cluster / cluster_minus_* -> 使用 task_type（这些依赖具体结构，不复用）
        import re
        if task_type == 'ion':
            molecule_name = 'Li+'
        elif task_type.startswith('ligand_'):
            molecule_name = task_type.replace('ligand_', '')
        elif task_type.startswith('dimer_'):
            molecule_name = task_type.replace('dimer_', '')
        elif task_type.startswith('redox_mol_'):
            # redox_mol_{name}_{state} -> 提取 name 部分
            # 状态后缀：neutral_gas, charged_gas, neutral_sol, charged_sol
            match = re.match(r'redox_mol_(.+)_(neutral|charged)_(gas|sol)$', task_type)
            if match:
                molecule_name = match.group(1)
            else:
                molecule_name = task_type.replace('redox_mol_', '').rsplit('_', 2)[0]
        elif task_type.startswith('redox_dimer_'):
            # redox_dimer_{name}_{state} -> 提取 name 部分
            match = re.match(r'redox_dimer_(.+)_(neutral|charged)_(gas|sol)$', task_type)
            if match:
                molecule_name = match.group(1)
            else:
                molecule_name = task_type.replace('redox_dimer_', '').rsplit('_', 2)[0]
        else:
            # cluster, cluster_minus_* 等依赖具体结构，使用完整 task_type
            molecule_name = task_type

        # 构建 QC 任务配置
        qc_job_config = {
            'molecule_name': molecule_name,
            'charge': charge,
            'spin_multiplicity': multiplicity,
            'functional': qc_config.get('functional', 'B3LYP'),
            'basis_set': qc_config.get('basis_set', '6-31G(d)'),
            'solvent_config': solvent_config,
            'md_job_id': md_job_id,
            # Cluster Analysis 关联字段
            'cluster_analysis_job_id': cluster_job_id,
            'task_type': task_type,
            'solvation_structure_id': structure_id,
            # Slurm 资源配置
            'slurm_partition': qc_config.get('slurm_partition', 'cpu'),
            'slurm_cpus': qc_config.get('slurm_cpus', 16),
            'slurm_time': qc_config.get('slurm_time', 7200),
        }

        # 只有当 smiles 有值时才添加（cluster 类型任务没有 SMILES，使用 XYZ 结构）
        if smiles:
            qc_job_config['smiles'] = smiles

        # 调用 API 创建 QC 任务
        try:
            endpoint = f"{self.api_base_url}/qc/jobs"
            response = requests.post(
                endpoint,
                headers=self.api_headers,
                json=qc_job_config,
                timeout=30
            )

            if response.status_code in [200, 201]:
                result = response.json()
                qc_job_id = result.get('id')
                self.logger.info(f"创建 QC 任务成功: {qc_job_id} (type={task_type})")

                # 创建成功后立即提交任务，将状态从 CREATED 改为 SUBMITTED
                submit_success = self._submit_qc_job(qc_job_id)
                if not submit_success:
                    self.logger.warning(f"QC 任务 {qc_job_id} 创建成功但提交失败，需要手动提交")

                return qc_job_id
            else:
                self.logger.error(f"创建 QC 任务失败: {response.status_code} - {response.text}")
                return None

        except Exception as e:
            self.logger.error(f"创建 QC 任务异常: {e}")
            return None

    def _submit_qc_job(self, qc_job_id: int) -> bool:
        """
        提交 QC 任务，将状态从 CREATED 改为 SUBMITTED

        Args:
            qc_job_id: QC 任务 ID

        Returns:
            bool: 提交是否成功
        """
        try:
            endpoint = f"{self.api_base_url}/qc/jobs/{qc_job_id}/submit"
            response = requests.post(
                endpoint,
                headers=self.api_headers,
                timeout=30
            )

            if response.status_code == 200:
                self.logger.info(f"QC 任务 {qc_job_id} 提交成功，状态已改为 SUBMITTED")
                return True
            else:
                self.logger.error(f"QC 任务 {qc_job_id} 提交失败: {response.status_code} - {response.text}")
                return False

        except Exception as e:
            self.logger.error(f"QC 任务 {qc_job_id} 提交异常: {e}")
            return False


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

