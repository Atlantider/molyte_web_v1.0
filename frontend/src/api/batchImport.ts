/**
 * 批量导入API
 */
import client from './client';

export interface BatchImportResult {
  success: boolean;
  message: string;
  project_id?: number;
  project_name?: string;
  total_electrolytes: number;
  success_electrolytes: number;
  failed_electrolytes: number;
  total_md_jobs: number;
  success_md_jobs: number;
  failed_md_jobs: number;
  total_qc_jobs: number;
  success_qc_jobs: number;
  failed_qc_jobs: number;
  electrolyte_results: Array<{
    row: number;
    name: string;
    id: number;
    status: string;
  }>;
  md_job_results: Array<{
    row: number;
    electrolyte_name: string;
    job_id: number;
    status: string;
  }>;
  qc_job_results: Array<{
    row: number;
    molecule_name: string;
    job_id: number;
    status: string;
  }>;
  errors: Array<{
    row: number;
    type: string;
    message: string;
  }>;
}

/**
 * 下载批量导入模板
 */
export async function downloadTemplate(
  templateType: 'electrolyte' | 'qc' = 'electrolyte',
  includeExamples: boolean = true
): Promise<Blob> {
  const response = await client.get('/batch-import/template/download', {
    params: {
      template_type: templateType,
      include_examples: includeExamples,
    },
    responseType: 'blob',
  });
  return response.data;
}

/**
 * 批量导入配方和QC计算
 */
export async function batchImportUpload(
  file: File,
  projectId?: number,
  projectName?: string,
  projectDescription?: string,
  sheetName: string = '配方'
): Promise<BatchImportResult> {
  const formData = new FormData();
  formData.append('file', file);
  
  if (projectId) {
    formData.append('project_id', projectId.toString());
  }
  if (projectName) {
    formData.append('project_name', projectName);
  }
  if (projectDescription) {
    formData.append('project_description', projectDescription);
  }
  formData.append('sheet_name', sheetName);

  const response = await client.post<BatchImportResult>(
    '/batch-import/upload',
    formData,
    {
      headers: {
        'Content-Type': 'multipart/form-data',
      },
    }
  );
  return response.data;
}

