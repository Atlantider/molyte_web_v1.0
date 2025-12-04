/**
 * Desolvation energy calculation API client
 * 去溶剂化能计算 API 客户端
 */
import client from './client';
import type {
  DesolvationJobCreate,
  DesolvationJobResponse,
  BatchDesolvationJobCreate,
  BatchDesolvationJobResponse,
  DesolvationOverviewResponse
} from '../types/desolvation';

/**
 * 创建去溶剂化能任务
 */
export const createDesolvationJob = async (
  data: DesolvationJobCreate
): Promise<DesolvationJobResponse> => {
  const response = await client.post('/desolvation/jobs', data);
  return response.data;
};

/**
 * 批量创建去溶剂化能任务
 */
export const batchCreateDesolvationJobs = async (
  data: BatchDesolvationJobCreate
): Promise<BatchDesolvationJobResponse> => {
  const response = await client.post('/desolvation/batch', data);
  return response.data;
};

/**
 * 获取去溶剂化能任务详情
 */
export const getDesolvationJob = async (
  jobId: number
): Promise<DesolvationJobResponse> => {
  const response = await client.get(`/desolvation/jobs/${jobId}`);
  return response.data;
};

/**
 * 获取某个 cluster 的所有去溶剂化能任务
 */
export const listClusterDesolvationJobs = async (
  clusterId: number
): Promise<DesolvationJobResponse[]> => {
  const response = await client.get(`/desolvation/cluster/${clusterId}/jobs`);
  return response.data;
};

/**
 * 获取某个 MD 任务下所有去溶剂化计算的总览
 */
export const getDesolvationOverview = async (
  mdJobId: number
): Promise<DesolvationOverviewResponse> => {
  const response = await client.get(`/desolvation/md/${mdJobId}/overview`);
  return response.data;
};

