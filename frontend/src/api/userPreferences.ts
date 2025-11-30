/**
 * 用户偏好设置 API
 * 管理用户自定义的常用溶剂组合、离子组合等
 */
import client from './client';

// ==================== 类型定义 ====================

export interface SolventItem {
  name: string;
  smiles: string;
  molar_ratio: number;
}

export interface CustomSolventCombination {
  name: string;
  description?: string;
  solvents: SolventItem[];
}

export interface CustomSolventCombinationResponse extends CustomSolventCombination {
  id: number;
  user_id: number;
  created_at: string;
  updated_at: string;
}

export interface IonCombination {
  name: string;
  charge: number;
  concentration: number;
}

export interface CustomIonCombination {
  name: string;
  description?: string;
  cations: IonCombination[];
  anions: IonCombination[];
}

export interface CustomIonCombinationResponse extends CustomIonCombination {
  id: number;
  user_id: number;
  created_at: string;
  updated_at: string;
}

// ==================== 溶剂组合 API ====================

/**
 * 获取用户的自定义溶剂组合列表
 */
export async function getUserSolventCombinations(): Promise<CustomSolventCombinationResponse[]> {
  const response = await client.get('/user-preferences/solvent-combinations');
  return response.data;
}

/**
 * 创建新的自定义溶剂组合
 */
export async function createSolventCombination(
  combination: CustomSolventCombination
): Promise<CustomSolventCombinationResponse> {
  const response = await client.post('/user-preferences/solvent-combinations', combination);
  return response.data;
}

/**
 * 更新自定义溶剂组合
 */
export async function updateSolventCombination(
  id: number,
  combination: CustomSolventCombination
): Promise<CustomSolventCombinationResponse> {
  const response = await client.put(`/user-preferences/solvent-combinations/${id}`, combination);
  return response.data;
}

/**
 * 删除自定义溶剂组合
 */
export async function deleteSolventCombination(id: number): Promise<void> {
  await client.delete(`/user-preferences/solvent-combinations/${id}`);
}

// ==================== 离子组合 API ====================

/**
 * 获取用户的自定义离子组合列表
 */
export async function getUserIonCombinations(): Promise<CustomIonCombinationResponse[]> {
  const response = await client.get('/user-preferences/ion-combinations');
  return response.data;
}

/**
 * 创建新的自定义离子组合
 */
export async function createIonCombination(
  combination: CustomIonCombination
): Promise<CustomIonCombinationResponse> {
  const response = await client.post('/user-preferences/ion-combinations', combination);
  return response.data;
}

/**
 * 删除自定义离子组合
 */
export async function deleteIonCombination(id: number): Promise<void> {
  await client.delete(`/user-preferences/ion-combinations/${id}`);
}

