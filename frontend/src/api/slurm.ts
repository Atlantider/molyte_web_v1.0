/**
 * Slurm 集群管理 API
 */
import client from './client';

/**
 * 分区信息
 */
export interface PartitionInfo {
  name: string;
  state: string;
  total_nodes: number;
  available_nodes: number;
  total_cpus: number;
  available_cpus: number;
  max_time?: string;
}

/**
 * 资源推荐
 */
export interface SlurmSuggestion {
  partition: string;
  ntasks: number;
  cpus_per_task: number;
  reason: string;
}

/**
 * 集群状态
 */
export interface ClusterStatus {
  cluster_status: 'online' | 'offline';
  partition_count: number;
  total_nodes: number;
  total_cpus: number;
  available_cpus: number;
  cpu_utilization: number;
  partitions: Array<{
    name: string;
    state: string;
    available_cpus: number;
    total_cpus: number;
  }>;
}

/**
 * 获取所有 Slurm 分区信息
 */
export const getPartitions = async (): Promise<PartitionInfo[]> => {
  const response = await client.get('/slurm/partitions');
  return response.data;
};

/**
 * 获取 Slurm 资源推荐配置
 */
export const getSlurmSuggestion = async (params?: {
  job_type?: string;
  expected_runtime_hours?: number;
  system_size?: number;
}): Promise<SlurmSuggestion> => {
  const response = await client.get('/slurm/suggestion', { params });
  return response.data;
};

/**
 * 获取 Slurm 集群整体状态
 */
export const getClusterStatus = async (): Promise<ClusterStatus> => {
  const response = await client.get('/slurm/status');
  return response.data;
};

