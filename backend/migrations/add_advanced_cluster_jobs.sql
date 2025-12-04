-- ============================================================================
-- 创建 Cluster 高级计算统一规划相关表
-- Migration: Add advanced_cluster_jobs table
-- Date: 2024-12-04
-- ============================================================================

-- 1. 创建 Cluster 高级计算任务状态枚举
DO $$ BEGIN
    CREATE TYPE advanced_cluster_job_status AS ENUM (
        'CREATED',
        'SUBMITTED',
        'RUNNING',
        'WAITING_QC',
        'CALCULATING',
        'COMPLETED',
        'FAILED',
        'CANCELLED'
    );
EXCEPTION
    WHEN duplicate_object THEN null;
END $$;

-- 2. 创建 Cluster 计算类型枚举
DO $$ BEGIN
    CREATE TYPE cluster_calc_type AS ENUM (
        'BINDING_TOTAL',
        'BINDING_PAIRWISE',
        'DESOLVATION_STEPWISE',
        'DESOLVATION_FULL',
        'REDOX',
        'REORGANIZATION'
    );
EXCEPTION
    WHEN duplicate_object THEN null;
END $$;

-- 3. 创建 Cluster 高级计算任务表
CREATE TABLE IF NOT EXISTS advanced_cluster_jobs (
    id SERIAL PRIMARY KEY,
    md_job_id INTEGER REFERENCES md_jobs(id) ON DELETE SET NULL,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    
    -- 任务状态
    status advanced_cluster_job_status DEFAULT 'CREATED' NOT NULL,
    progress FLOAT DEFAULT 0.0,
    error_message TEXT,
    
    -- 计算配置
    calc_types JSONB DEFAULT '[]' NOT NULL,           -- 选择的计算类型列表
    selected_structures JSONB DEFAULT '{}' NOT NULL,  -- 选中的溶剂化结构
    qc_config JSONB DEFAULT '{}',                     -- QC 计算配置
    
    -- QC 任务规划和跟踪
    qc_task_plan JSONB DEFAULT '{}',                  -- QC 任务规划详情
    
    -- 计算结果
    results JSONB DEFAULT '{}',                       -- 各类型的计算结果
    
    -- 时间戳
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP WITH TIME ZONE,
    finished_at TIMESTAMP WITH TIME ZONE
);

-- 4. 创建索引
CREATE INDEX IF NOT EXISTS idx_advanced_cluster_md_job_id ON advanced_cluster_jobs(md_job_id);
CREATE INDEX IF NOT EXISTS idx_advanced_cluster_user_id ON advanced_cluster_jobs(user_id);
CREATE INDEX IF NOT EXISTS idx_advanced_cluster_status ON advanced_cluster_jobs(status);
CREATE INDEX IF NOT EXISTS idx_advanced_cluster_created_at ON advanced_cluster_jobs(created_at);

-- 5. 创建触发器：自动更新 updated_at
CREATE OR REPLACE FUNCTION update_advanced_cluster_jobs_updated_at()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

DROP TRIGGER IF EXISTS trigger_advanced_cluster_jobs_updated_at ON advanced_cluster_jobs;
CREATE TRIGGER trigger_advanced_cluster_jobs_updated_at
    BEFORE UPDATE ON advanced_cluster_jobs
    FOR EACH ROW
    EXECUTE FUNCTION update_advanced_cluster_jobs_updated_at();

-- 6. 添加注释
COMMENT ON TABLE advanced_cluster_jobs IS 'Cluster 高级计算统一规划任务表';

COMMENT ON COLUMN advanced_cluster_jobs.calc_types IS '选择的计算类型列表: ["BINDING_TOTAL", "DESOLVATION_STEPWISE", ...]';
COMMENT ON COLUMN advanced_cluster_jobs.selected_structures IS '选中的溶剂化结构: {solvation_structure_ids: [...], count: N}';
COMMENT ON COLUMN advanced_cluster_jobs.qc_config IS 'QC 计算配置: {functional, basis_set, solvent_model, ...}';
COMMENT ON COLUMN advanced_cluster_jobs.qc_task_plan IS 'QC 任务规划: {planned_qc_tasks: [...], reused_qc_jobs: [...], new_qc_jobs: [...]}';
COMMENT ON COLUMN advanced_cluster_jobs.results IS '计算结果: {BINDING_TOTAL: {...}, DESOLVATION_STEPWISE: {...}, ...}';

-- ============================================================================
-- 完成信息
-- ============================================================================
\echo '========================================='
\echo 'Advanced Cluster Jobs Table Migration Completed!'
\echo '========================================='
\echo 'Created tables:'
\echo '  - advanced_cluster_jobs (Cluster 高级计算统一规划)'
\echo ''
\echo 'Supported calc_types:'
\echo '  - BINDING_TOTAL: 总 Binding Energy'
\echo '  - BINDING_PAIRWISE: 分子-Li Binding'
\echo '  - DESOLVATION_STEPWISE: 逐级去溶剂化'
\echo '  - DESOLVATION_FULL: 完全去溶剂化'
\echo '  - REDOX: 氧化还原电位'
\echo '  - REORGANIZATION: Marcus 重组能'
\echo '========================================='

