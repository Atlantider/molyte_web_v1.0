/**
 * 电解液研发 API
 */
import client from './client';

export interface ElectrolyteSearchParams {
  cations?: string[];
  anions?: string[];
  solvents?: string[];
  solvent_smiles?: string;  // 按 SMILES 搜索溶剂
  temp_min?: number;
  temp_max?: number;
  skip?: number;
  limit?: number;
}

export interface ElectrolyteSearchResult {
  job_id: number;
  job_name?: string;
  user_note?: string;  // 用户备注
  system_id: number;
  system_name: string;
  cations: Array<{ name: string; number: number; smiles?: string; concentration?: number }>;
  anions: Array<{ name: string; number: number; smiles?: string; concentration?: number }>;
  solvents: Array<{ name: string; number: number; smiles: string; concentration?: number }>;
  temperature: number;
  pressure: number;
  density?: number;
  concentration?: number;
  // 计算方式
  charge_method?: 'ligpargen' | 'resp';  // 电荷计算方式
  qc_enabled?: boolean;  // 是否有QC计算
  created_at: string;
  finished_at: string;
  has_rdf: boolean;
  has_msd: boolean;
  has_solvation: boolean;
  rdf_count: number;
  msd_count: number;
  solvation_count: number;
}

export interface ElectrolyteSearchResponse {
  total: number;
  skip: number;
  limit: number;
  data: ElectrolyteSearchResult[];
}

/**
 * 构建搜索查询参数
 */
const buildSearchQueryParams = (params: ElectrolyteSearchParams): URLSearchParams => {
  const queryParams = new URLSearchParams();

  if (params.cations) {
    params.cations.forEach(c => queryParams.append('cations', c));
  }
  if (params.anions) {
    params.anions.forEach(a => queryParams.append('anions', a));
  }
  if (params.solvents) {
    params.solvents.forEach(s => queryParams.append('solvents', s));
  }
  if (params.solvent_smiles) {
    queryParams.append('solvent_smiles', params.solvent_smiles);
  }
  if (params.temp_min !== undefined) {
    queryParams.append('temp_min', params.temp_min.toString());
  }
  if (params.temp_max !== undefined) {
    queryParams.append('temp_max', params.temp_max.toString());
  }
  if (params.skip !== undefined) {
    queryParams.append('skip', params.skip.toString());
  }
  if (params.limit !== undefined) {
    queryParams.append('limit', params.limit.toString());
  }

  return queryParams;
};

/**
 * 搜索所有用户的电解液计算结果（公开搜索）
 */
export const searchElectrolytes = async (
  params: ElectrolyteSearchParams
): Promise<ElectrolyteSearchResponse> => {
  const queryParams = buildSearchQueryParams(params);
  const response = await client.get(`/research/search?${queryParams.toString()}`);
  return response.data;
};

/**
 * 搜索当前用户的电解液计算结果（工作台内使用）
 */
export const searchMyElectrolytes = async (
  params: ElectrolyteSearchParams
): Promise<ElectrolyteSearchResponse> => {
  const queryParams = buildSearchQueryParams(params);
  const response = await client.get(`/research/search/my?${queryParams.toString()}`);
  return response.data;
};

/**
 * 获取可用的搜索选项（从数据库中实际数据提取）
 */
export interface AvailableOptions {
  cations: string[];
  anions: string[];
  solvents: string[];
}

export const getAvailableSearchOptions = async (): Promise<AvailableOptions> => {
  const response = await client.get('/research/available-options');
  return response.data;
};
