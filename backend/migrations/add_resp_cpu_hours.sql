-- 添加 RESP CPU 核时数字段到 md_jobs 表
-- 用于追踪 RESP 电荷计算消耗的核时数

-- 添加 resp_cpu_hours 字段到 md_jobs 表
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS resp_cpu_hours FLOAT DEFAULT 0.0;

-- 添加注释
COMMENT ON COLUMN md_jobs.resp_cpu_hours IS 'RESP 电荷计算消耗的 CPU 核时数';

-- 创建 resp_jobs 表（如果不存在）
CREATE TABLE IF NOT EXISTS resp_jobs (
    id SERIAL PRIMARY KEY,
    md_job_id INTEGER NOT NULL REFERENCES md_jobs(id) ON DELETE CASCADE,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    
    -- 分子信息
    molecule_name VARCHAR(255) NOT NULL,
    smiles TEXT,
    
    -- 任务状态
    status VARCHAR(20) NOT NULL DEFAULT 'CREATED',
    slurm_job_id VARCHAR(50),
    
    -- 工作目录和文件
    work_dir TEXT,
    charge_file TEXT,
    log_file TEXT,
    error_message TEXT,
    
    -- 核时数统计
    cpu_hours FLOAT DEFAULT 0.0,
    estimated_cpu_hours FLOAT,
    
    -- 配置
    config JSONB,
    
    -- 时间戳
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW() NOT NULL,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW() NOT NULL,
    started_at TIMESTAMP WITH TIME ZONE,
    finished_at TIMESTAMP WITH TIME ZONE
);

-- 创建索引
CREATE INDEX IF NOT EXISTS idx_resp_jobs_md_job_id ON resp_jobs(md_job_id);
CREATE INDEX IF NOT EXISTS idx_resp_jobs_user_id ON resp_jobs(user_id);
CREATE INDEX IF NOT EXISTS idx_resp_jobs_status ON resp_jobs(status);
CREATE INDEX IF NOT EXISTS idx_resp_jobs_slurm_job_id ON resp_jobs(slurm_job_id);

