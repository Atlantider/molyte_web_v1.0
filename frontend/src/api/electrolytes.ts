/**
 * 电解质体系管理 API
 */
import client from './client';
import type { ElectrolyteSystem, ElectrolyteSystemCreate, ElectrolyteSystemUpdate, AvailableIons } from '../types';

/**
 * 获取可用的离子列表
 */
export const getAvailableIons = async (): Promise<AvailableIons> => {
  const response = await client.get('/electrolytes/available-ions');
  return response.data;
};

// 验证 SMILES
export const validateSmiles = async (smiles: string): Promise<any> => {
  const response = await client.post('/electrolytes/validate-smiles', { smiles });
  return response.data;
};

// 刷新离子缓存
export const refreshIons = async (): Promise<AvailableIons> => {
  const response = await client.post('/electrolytes/refresh-ions');
  return response.data;
};

/**
 * 获取所有电解质体系
 */
export const getElectrolytes = async (): Promise<ElectrolyteSystem[]> => {
  const response = await client.get('/electrolytes/');
  return response.data;
};

/**
 * 获取单个电解质体系
 */
export const getElectrolyte = async (id: number): Promise<ElectrolyteSystem> => {
  const response = await client.get(`/electrolytes/${id}`);
  return response.data;
};

/**
 * 创建电解质体系（旧格式）
 */
export const createElectrolyte = async (data: ElectrolyteSystemCreate): Promise<ElectrolyteSystem> => {
  const response = await client.post('/electrolytes/', data);
  return response.data;
};

/**
 * 创建电解质体系（新格式 - 使用浓度）
 */
export const createElectrolyteNew = async (data: any): Promise<ElectrolyteSystem> => {
  const response = await client.post('/electrolytes/new', data);
  return response.data;
};

/**
 * 获取可编辑格式的电解质体系（新格式 - 浓度）
 */
export const getElectrolyteEditable = async (id: number): Promise<any> => {
  const response = await client.get(`/electrolytes/${id}/editable`);
  return response.data;
};

/**
 * 更新电解质体系（旧格式）
 */
export const updateElectrolyte = async (id: number, data: ElectrolyteSystemUpdate): Promise<ElectrolyteSystem> => {
  const response = await client.put(`/electrolytes/${id}`, data);
  return response.data;
};

/**
 * 更新电解质体系（新格式 - 使用浓度）
 */
export const updateElectrolyteNew = async (id: number, data: any): Promise<ElectrolyteSystem> => {
  const response = await client.put(`/electrolytes/${id}`, data);
  return response.data;
};

/**
 * 删除电解质体系
 */
export const deleteElectrolyte = async (id: number): Promise<void> => {
  await client.delete(`/electrolytes/${id}`);
};

/**
 * 批量删除电解质配方
 */
export const batchDeleteElectrolytes = async (ids: number[]): Promise<{
  deleted_count: number;
  failed_ids: number[];
  failed_reasons?: Record<number, string>;
  message: string;
}> => {
  const response = await client.delete('/electrolytes/batch/delete', { data: { ids } });
  return response.data;
};

/**
 * 批量更改电解质配方的项目归属
 */
export const batchUpdateProject = async (ids: number[], projectId: number): Promise<{
  updated_count: number;
  failed_ids: number[];
  message: string;
}> => {
  const response = await client.put('/electrolytes/batch/project', { ids, project_id: projectId });
  return response.data;
};

