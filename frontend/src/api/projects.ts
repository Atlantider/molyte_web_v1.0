/**
 * 项目管理 API
 */
import client from './client';
import type { Project, ProjectCreate, ProjectUpdate } from '../types';

/**
 * 获取所有项目
 */
export const getProjects = async (): Promise<Project[]> => {
  const response = await client.get('/projects/');
  return response.data;
};

/**
 * 获取单个项目
 */
export const getProject = async (id: number): Promise<Project> => {
  const response = await client.get(`/projects/${id}`);
  return response.data;
};

/**
 * 创建项目
 */
export const createProject = async (data: ProjectCreate): Promise<Project> => {
  const response = await client.post('/projects/', data);
  return response.data;
};

/**
 * 更新项目
 */
export const updateProject = async (id: number, data: ProjectUpdate): Promise<Project> => {
  const response = await client.put(`/projects/${id}`, data);
  return response.data;
};

/**
 * 删除项目
 */
export const deleteProject = async (id: number): Promise<void> => {
  await client.delete(`/projects/${id}`);
};

