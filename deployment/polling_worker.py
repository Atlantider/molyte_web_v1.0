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
    import oss2
except ImportError:
    print("请安装 OSS SDK: pip install oss2")
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
        """初始化 OSS 客户端"""
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
        self.logger.info("OSS 客户端已初始化")
    
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

    def _upload_results_to_oss(self, job_id: int, work_dir: Path) -> List[str]:
        """上传结果文件到 OSS"""
        uploaded_files = []

        try:
            # 获取需要上传的文件
            result_patterns = self.config['upload']['result_files']
            max_size = self.config['upload']['max_file_size'] * 1024 * 1024  # MB to bytes

            for pattern in result_patterns:
                for file_path in work_dir.glob(pattern):
                    if not file_path.is_file():
                        continue

                    # 检查文件大小
                    file_size = file_path.stat().st_size
                    if file_size > max_size:
                        self.logger.warning(f"文件 {file_path.name} 超过大小限制，跳过")
                        continue

                    # 上传到 OSS
                    oss_key = f"{self.config['oss']['result_prefix']}{job_id}/{file_path.name}"

                    self.logger.info(f"上传文件: {file_path.name} -> {oss_key}")
                    self.oss_bucket.put_object_from_file(oss_key, str(file_path))

                    uploaded_files.append(oss_key)

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
                timeout=self.config['api']['timeout']
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

