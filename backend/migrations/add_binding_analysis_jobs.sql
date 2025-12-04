-- ============================================================================
-- 创建 Binding 分析任务表
-- Migration: Create binding_analysis_jobs table
-- Date: 2024-12-04
-- ============================================================================

-- 1. 创建 Binding 分析状态枚举类型
DO $$ BEGIN
    CREATE TYPE binding_analysis_status AS ENUM (
        'CREATED',
        'SUBMITTED',
        'RUNNING',
        'COMPLETED',
        'FAILED'
    );
EXCEPTION
    WHEN duplicate_object THEN null;
END $$;

-- 2. 创建 binding_analysis_jobs 表
CREATE TABLE IF NOT EXISTS binding_analysis_jobs (
    id SERIAL PRIMARY KEY,
    md_job_id INTEGER NOT NULL REFERENCES md_jobs(id) ON DELETE CASCADE,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    
    -- 任务状态
    status binding_analysis_status DEFAULT 'CREATED' NOT NULL,
    progress FLOAT DEFAULT 0.0,
    error_message TEXT,
    
    -- 配置和结果
    config JSONB DEFAULT '{}',
    result JSONB,
    qc_job_ids JSONB,  -- 关联的 QC 任务 ID 列表
    
    -- 时间戳
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP WITH TIME ZONE,
    finished_at TIMESTAMP WITH TIME ZONE
);

-- 3. 创建索引
CREATE INDEX IF NOT EXISTS idx_binding_analysis_jobs_md_job_id ON binding_analysis_jobs(md_job_id);
CREATE INDEX IF NOT EXISTS idx_binding_analysis_jobs_user_id ON binding_analysis_jobs(user_id);
CREATE INDEX IF NOT EXISTS idx_binding_analysis_jobs_status ON binding_analysis_jobs(status);
CREATE INDEX IF NOT EXISTS idx_binding_analysis_jobs_created_at ON binding_analysis_jobs(created_at);

-- 4. 创建触发器：自动更新 updated_at
CREATE OR REPLACE FUNCTION update_binding_analysis_jobs_updated_at()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

DROP TRIGGER IF EXISTS trigger_binding_analysis_jobs_updated_at ON binding_analysis_jobs;
CREATE TRIGGER trigger_binding_analysis_jobs_updated_at
    BEFORE UPDATE ON binding_analysis_jobs
    FOR EACH ROW
    EXECUTE FUNCTION update_binding_analysis_jobs_updated_at();

-- 5. 添加注释
COMMENT ON TABLE binding_analysis_jobs IS 'Li-配体 Binding Energy 分析任务表';
COMMENT ON COLUMN binding_analysis_jobs.config IS '任务配置: {composition_keys, functional, basis_set, solvent_model, reuse_existing_qc}';
COMMENT ON COLUMN binding_analysis_jobs.result IS '分析结果: {per_cluster_results, summary}';
COMMENT ON COLUMN binding_analysis_jobs.qc_job_ids IS '关联的 QC 任务 ID 列表';

-- ============================================================================
-- 完成信息
-- ============================================================================
\echo '========================================='
\echo 'Binding Analysis Jobs Migration Completed!'
\echo '========================================='
\echo 'Created table: binding_analysis_jobs'
\echo 'Created indexes:'
\echo '  - idx_binding_analysis_jobs_md_job_id'
\echo '  - idx_binding_analysis_jobs_user_id'
\echo '  - idx_binding_analysis_jobs_status'
\echo '  - idx_binding_analysis_jobs_created_at'
\echo '========================================='

\dt binding_analysis_jobs

