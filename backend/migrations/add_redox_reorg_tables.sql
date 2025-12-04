-- ============================================================================
-- 创建热力学循环和重组能计算相关表
-- Migration: Add redox_potential_jobs and reorganization_energy_jobs tables
-- Date: 2024-12-04
-- ============================================================================

-- 1. 创建热力学循环任务状态枚举
DO $$ BEGIN
    CREATE TYPE redox_job_status AS ENUM (
        'CREATED',
        'SUBMITTED',
        'RUNNING',
        'COMPLETED',
        'FAILED'
    );
EXCEPTION
    WHEN duplicate_object THEN null;
END $$;

-- 2. 创建重组能任务状态枚举
DO $$ BEGIN
    CREATE TYPE reorg_energy_job_status AS ENUM (
        'CREATED',
        'SUBMITTED',
        'RUNNING',
        'COMPLETED',
        'FAILED'
    );
EXCEPTION
    WHEN duplicate_object THEN null;
END $$;

-- 3. 创建热力学循环计算任务表
CREATE TABLE IF NOT EXISTS redox_potential_jobs (
    id SERIAL PRIMARY KEY,
    md_job_id INTEGER REFERENCES md_jobs(id) ON DELETE SET NULL,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    
    -- 任务状态
    status redox_job_status DEFAULT 'CREATED' NOT NULL,
    progress FLOAT DEFAULT 0.0,
    error_message TEXT,
    
    -- 配置和结果
    config JSONB DEFAULT '{}',
    result JSONB,
    qc_job_ids JSONB,
    
    -- 时间戳
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP WITH TIME ZONE,
    finished_at TIMESTAMP WITH TIME ZONE
);

-- 4. 创建重组能计算任务表
CREATE TABLE IF NOT EXISTS reorganization_energy_jobs (
    id SERIAL PRIMARY KEY,
    md_job_id INTEGER REFERENCES md_jobs(id) ON DELETE SET NULL,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    
    -- 任务状态
    status reorg_energy_job_status DEFAULT 'CREATED' NOT NULL,
    progress FLOAT DEFAULT 0.0,
    error_message TEXT,
    
    -- 配置和结果
    config JSONB DEFAULT '{}',
    result JSONB,
    qc_job_ids JSONB,
    
    -- 时间戳
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP WITH TIME ZONE,
    finished_at TIMESTAMP WITH TIME ZONE
);

-- 5. 创建索引
CREATE INDEX IF NOT EXISTS idx_redox_potential_md_job_id ON redox_potential_jobs(md_job_id);
CREATE INDEX IF NOT EXISTS idx_redox_potential_user_id ON redox_potential_jobs(user_id);
CREATE INDEX IF NOT EXISTS idx_redox_potential_status ON redox_potential_jobs(status);
CREATE INDEX IF NOT EXISTS idx_redox_potential_created_at ON redox_potential_jobs(created_at);

CREATE INDEX IF NOT EXISTS idx_reorg_energy_md_job_id ON reorganization_energy_jobs(md_job_id);
CREATE INDEX IF NOT EXISTS idx_reorg_energy_user_id ON reorganization_energy_jobs(user_id);
CREATE INDEX IF NOT EXISTS idx_reorg_energy_status ON reorganization_energy_jobs(status);
CREATE INDEX IF NOT EXISTS idx_reorg_energy_created_at ON reorganization_energy_jobs(created_at);

-- 6. 创建触发器：自动更新 updated_at
CREATE OR REPLACE FUNCTION update_redox_potential_jobs_updated_at()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

DROP TRIGGER IF EXISTS trigger_redox_potential_jobs_updated_at ON redox_potential_jobs;
CREATE TRIGGER trigger_redox_potential_jobs_updated_at
    BEFORE UPDATE ON redox_potential_jobs
    FOR EACH ROW
    EXECUTE FUNCTION update_redox_potential_jobs_updated_at();

CREATE OR REPLACE FUNCTION update_reorganization_energy_jobs_updated_at()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

DROP TRIGGER IF EXISTS trigger_reorganization_energy_jobs_updated_at ON reorganization_energy_jobs;
CREATE TRIGGER trigger_reorganization_energy_jobs_updated_at
    BEFORE UPDATE ON reorganization_energy_jobs
    FOR EACH ROW
    EXECUTE FUNCTION update_reorganization_energy_jobs_updated_at();

-- 7. 添加注释
COMMENT ON TABLE redox_potential_jobs IS '热力学循环计算氧化还原电位任务表 (高风险功能)';
COMMENT ON TABLE reorganization_energy_jobs IS '重组能计算任务表 (极高风险功能)';

COMMENT ON COLUMN redox_potential_jobs.config IS '任务配置: {species_list, mode, functional, basis_set, solvent_model, ...}';
COMMENT ON COLUMN redox_potential_jobs.result IS '计算结果: {species_results, oxidation_potentials_v, reduction_potentials_v, ...}';

COMMENT ON COLUMN reorganization_energy_jobs.config IS '任务配置: {species_list, functional, basis_set, ...}';
COMMENT ON COLUMN reorganization_energy_jobs.result IS '计算结果: {species_results, lambda_ox_mean_ev, lambda_red_mean_ev, ...}';

-- ============================================================================
-- 完成信息
-- ============================================================================
\echo '========================================='
\echo 'Redox & Reorganization Energy Tables Migration Completed!'
\echo '========================================='
\echo 'Created tables:'
\echo '  - redox_potential_jobs (热力学循环)'
\echo '  - reorganization_energy_jobs (重组能)'
\echo '========================================='

