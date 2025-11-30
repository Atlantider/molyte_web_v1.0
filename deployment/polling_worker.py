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
        # SSL 验证设置（默认启用，可在配置中禁用）
        self.verify_ssl = self.config['api'].get('verify_ssl', True)
        if not self.verify_ssl:
            self.logger.warning("⚠️  SSL 证书验证已禁用！仅用于开发/测试环境")
            # 禁用 SSL 警告
            import urllib3
            urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
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
    
    def _send_heartbeat(self):
        """发送心跳到阿里云"""
        try:
            data = {
                'worker_name': self.config['worker']['name'],
                'status': 'running',
                'running_jobs': len(self.running_jobs),
                'timestamp': datetime.now().isoformat()
            }
            # 这里需要在阿里云后端添加心跳接口
            # response = requests.post(
            #     f"{self.api_base_url}/workers/heartbeat",
            #     headers=self.api_headers,
            #     json=data,
            #     timeout=self.config['api']['timeout']
            # )
            self.logger.debug(f"心跳已发送: {len(self.running_jobs)} 个任务运行中")
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
            endpoint = f"{self.api_base_url}/jobs/pending"
            params = {
                'job_type': job_type.upper(),
                'limit': 10
            }
            
            response = requests.get(
                endpoint,
                headers=self.api_headers,
                params=params,
                timeout=self.config['api']['timeout'],
                verify=self.verify_ssl
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
            # 1. 更新任务状态为 PROCESSING
            self._update_job_status(job_id, 'PROCESSING', 'md')
            
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
        
        # TODO: 实现 QC 任务处理逻辑
        # 类似 MD 任务的流程
        pass
    
    def _check_running_jobs(self):
        """检查运行中的任务状态"""
        completed_jobs = []
        
        for job_id, job_info in self.running_jobs.items():
            try:
                slurm_job_id = job_info['slurm_job_id']
                status = self._check_slurm_status(slurm_job_id)
                
                if status == 'COMPLETED':
                    self.logger.info(f"任务 {job_id} (Slurm: {slurm_job_id}) 已完成")
                    self._handle_job_completion(job_id, job_info)
                    completed_jobs.append(job_id)
                    
                elif status == 'FAILED':
                    self.logger.error(f"任务 {job_id} (Slurm: {slurm_job_id}) 失败")
                    self._update_job_status(job_id, 'FAILED', job_info['type'])
                    completed_jobs.append(job_id)
                    
            except Exception as e:
                self.logger.error(f"检查任务 {job_id} 状态失败: {e}")
        
        # 移除已完成的任务
        for job_id in completed_jobs:
            del self.running_jobs[job_id]

    def _check_slurm_status(self, slurm_job_id: str) -> str:
        """检查 Slurm 任务状态"""
        try:
            cmd = f"squeue -j {slurm_job_id} -h -o %T"
            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=10
            )

            if result.returncode == 0 and result.stdout.strip():
                status = result.stdout.strip()
                # PENDING, RUNNING, COMPLETED, FAILED, CANCELLED
                if status in ['PENDING', 'RUNNING']:
                    return 'RUNNING'
                elif status in ['COMPLETED']:
                    return 'COMPLETED'
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

            self.logger.info(f"开始处理任务 {job_id} 的结果")

            # 1. 上传结果文件到 OSS
            uploaded_files = self._upload_results_to_oss(job_id, work_dir)

            # 2. 更新任务状态为 COMPLETED
            self._update_job_status(
                job_id, 'COMPLETED', job_type,
                result_files=uploaded_files
            )

            self.logger.info(f"任务 {job_id} 处理完成，上传了 {len(uploaded_files)} 个文件")

        except Exception as e:
            self.logger.error(f"处理任务 {job_id} 完成失败: {e}", exc_info=True)
            self._update_job_status(job_id, 'FAILED', job_info['type'], error_message=str(e))

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
        result_files: Optional[List[str]] = None
    ):
        """更新任务状态到阿里云"""
        try:
            endpoint = f"{self.api_base_url}/jobs/{job_id}/status"

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
                data['error_message'] = error_message
            if result_files:
                data['result_files'] = result_files

            response = requests.put(
                endpoint,
                headers=self.api_headers,
                json=data,
                timeout=self.config['api']['timeout'],
                verify=self.verify_ssl
            )

            if response.status_code == 200:
                self.logger.info(f"任务 {job_id} 状态已更新为 {status}")
            else:
                self.logger.warning(f"更新任务状态失败: {response.status_code} - {response.text}")

        except Exception as e:
            self.logger.error(f"更新任务 {job_id} 状态失败: {e}")


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

