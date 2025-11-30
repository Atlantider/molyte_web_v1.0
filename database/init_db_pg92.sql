-- ============================================================================
-- Molyte Web 平台数据库初始化脚本
-- 版本: 1.0 (MVP) - PostgreSQL 9.2 兼容版本
-- 创建日期: 2025-11-19
-- 说明: 创建核心数据表、索引和约束（兼容 PostgreSQL 9.2）
-- 主要修改: JSONB -> JSON, 移除部分 9.2 不支持的特性
-- ============================================================================

-- 设置客户端编码
SET client_encoding = 'UTF8';

-- ============================================================================
-- 1. 删除已存在的表（开发阶段使用，按依赖顺序删除）
-- ============================================================================

DROP TABLE IF EXISTS rdf_results CASCADE;
DROP TABLE IF EXISTS msd_results CASCADE;
DROP TABLE IF EXISTS solvation_structures CASCADE;
DROP TABLE IF EXISTS result_summary CASCADE;
DROP TABLE IF EXISTS postprocess_jobs CASCADE;
DROP TABLE IF EXISTS md_jobs CASCADE;
DROP TABLE IF EXISTS electrolyte_systems CASCADE;
DROP TABLE IF EXISTS projects CASCADE;
DROP TABLE IF EXISTS users CASCADE;

-- ============================================================================
-- 2. 创建用户表
-- ============================================================================

CREATE TABLE users (
    id SERIAL PRIMARY KEY,
    email VARCHAR(255) UNIQUE NOT NULL,
    username VARCHAR(100),
    password_hash VARCHAR(255) NOT NULL,
    role VARCHAR(20) DEFAULT 'user' CHECK (role IN ('user', 'admin')),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

COMMENT ON TABLE users IS '用户表';
COMMENT ON COLUMN users.role IS '用户角色: user=普通用户, admin=管理员';

-- ============================================================================
-- 3. 创建项目表
-- ============================================================================

CREATE TABLE projects (
    id SERIAL PRIMARY KEY,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    name VARCHAR(255) NOT NULL,
    description TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

COMMENT ON TABLE projects IS '项目表';

-- ============================================================================
-- 4. 创建电解液体系表
-- ============================================================================

CREATE TABLE electrolyte_systems (
    id SERIAL PRIMARY KEY,
    project_id INTEGER NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
    hash_key VARCHAR(64) UNIQUE NOT NULL,
    name VARCHAR(255),
    
    -- 配方信息（JSON 格式）
    cations JSON NOT NULL,
    anions JSON NOT NULL,
    solvents JSON,
    
    -- 模拟条件
    temperature FLOAT DEFAULT 298.15,
    pressure FLOAT DEFAULT 1.0,
    density FLOAT,
    concentration FLOAT,
    box_size FLOAT,
    
    -- 模拟参数
    nsteps_npt INTEGER DEFAULT 5000000,
    nsteps_nvt INTEGER DEFAULT 10000000,
    timestep FLOAT DEFAULT 1.0,
    force_field VARCHAR(50) DEFAULT 'OPLS',
    
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

COMMENT ON TABLE electrolyte_systems IS '电解液体系表（规范化配方）';
COMMENT ON COLUMN electrolyte_systems.hash_key IS '配方哈希值，用于去重';
COMMENT ON COLUMN electrolyte_systems.cations IS '阳离子列表 JSON: [{"name": "Li", "smiles": "[Li+]", "number": 50}]';
COMMENT ON COLUMN electrolyte_systems.anions IS '阴离子列表 JSON';
COMMENT ON COLUMN electrolyte_systems.solvents IS '溶剂列表 JSON';

-- ============================================================================
-- 5. 创建 MD 任务表
-- ============================================================================

CREATE TABLE md_jobs (
    id SERIAL PRIMARY KEY,
    system_id INTEGER NOT NULL REFERENCES electrolyte_systems(id) ON DELETE CASCADE,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    
    -- 任务状态
    status VARCHAR(20) NOT NULL DEFAULT 'CREATED' 
        CHECK (status IN ('CREATED', 'QUEUED', 'RUNNING', 'POSTPROCESSING', 'COMPLETED', 'FAILED')),
    
    -- Slurm 信息
    slurm_job_id VARCHAR(50),
    
    -- 进度和路径
    progress INTEGER DEFAULT 0 CHECK (progress >= 0 AND progress <= 100),
    work_dir TEXT,
    log_file TEXT,
    error_message TEXT,
    
    -- 时间戳
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP,
    finished_at TIMESTAMP
);

COMMENT ON TABLE md_jobs IS 'MD 任务表';
COMMENT ON COLUMN md_jobs.status IS '任务状态: CREATED/QUEUED/RUNNING/POSTPROCESSING/COMPLETED/FAILED';
COMMENT ON COLUMN md_jobs.progress IS '任务进度 0-100';

-- ============================================================================
-- 6. 创建后处理任务表
-- ============================================================================

CREATE TABLE postprocess_jobs (
    id SERIAL PRIMARY KEY,
    md_job_id INTEGER NOT NULL REFERENCES md_jobs(id) ON DELETE CASCADE,
    status VARCHAR(20) NOT NULL DEFAULT 'PENDING'
        CHECK (status IN ('PENDING', 'RUNNING', 'COMPLETED', 'FAILED')),
    required_results JSON DEFAULT '["rdf", "msd", "solvation"]'::json,
    error_message TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    finished_at TIMESTAMP
);

COMMENT ON TABLE postprocess_jobs IS '后处理任务表';
COMMENT ON COLUMN postprocess_jobs.required_results IS '需要计算的结果类型 JSON 数组';

-- ============================================================================
-- 7. 创建结果概览表
-- ============================================================================

CREATE TABLE result_summary (
    id SERIAL PRIMARY KEY,
    md_job_id INTEGER UNIQUE NOT NULL REFERENCES md_jobs(id) ON DELETE CASCADE,
    
    -- 基本物理量
    final_density FLOAT,
    avg_temperature FLOAT,
    avg_pressure FLOAT,
    avg_energy FLOAT,
    box_volume FLOAT,
    
    -- 可用结果类型
    available_results JSON DEFAULT '[]'::json,
    
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

COMMENT ON TABLE result_summary IS '结果概览表';
COMMENT ON COLUMN result_summary.available_results IS '可用结果类型 JSON 数组: ["rdf", "msd", "solvation"]';

-- ============================================================================
-- 8. 创建 RDF 结果表
-- ============================================================================

CREATE TABLE rdf_results (
    id SERIAL PRIMARY KEY,
    md_job_id INTEGER NOT NULL REFERENCES md_jobs(id) ON DELETE CASCADE,

    -- RDF 配置
    center_species VARCHAR(100) NOT NULL,
    shell_species VARCHAR(100) NOT NULL,

    -- RDF 数据
    r_values JSONB,
    g_r_values JSONB,
    coordination_number JSONB,

    -- 计算参数
    cutoff FLOAT DEFAULT 10.0,
    bin_width FLOAT DEFAULT 0.1,

    -- 文件路径
    file_path TEXT,

    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

COMMENT ON TABLE rdf_results IS 'RDF（径向分布函数）结果表';
COMMENT ON COLUMN rdf_results.center_species IS '中心原子/离子，如 "Li", "K"';
COMMENT ON COLUMN rdf_results.shell_species IS '壳层原子，如 "EC_O_carbonyl", "DMC_O_ether"';
COMMENT ON COLUMN rdf_results.r_values IS '距离数组 JSON';
COMMENT ON COLUMN rdf_results.g_r_values IS 'g(r) 值数组 JSON';
COMMENT ON COLUMN rdf_results.coordination_number IS '配位数数组 JSON';

-- ============================================================================
-- 9. 创建 MSD 结果表
-- ============================================================================

CREATE TABLE msd_results (
    id SERIAL PRIMARY KEY,
    md_job_id INTEGER NOT NULL REFERENCES md_jobs(id) ON DELETE CASCADE,

    -- MSD 配置
    species VARCHAR(100) NOT NULL,

    -- MSD 数据
    t_values JSONB,
    msd_values JSONB,

    -- 扩散系数
    diffusion_coeff FLOAT,
    fit_range JSONB,

    -- 文件路径
    file_path TEXT,

    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

COMMENT ON TABLE msd_results IS 'MSD（均方位移）结果表';
COMMENT ON COLUMN msd_results.species IS '粒子类型，如 "Li", "anion", "solvent"';
COMMENT ON COLUMN msd_results.t_values IS '时间数组 JSON';
COMMENT ON COLUMN msd_results.msd_values IS 'MSD 值数组 JSON';
COMMENT ON COLUMN msd_results.diffusion_coeff IS '扩散系数 (cm²/s)';
COMMENT ON COLUMN msd_results.fit_range IS '拟合时间范围 JSON: [t_start, t_end]';

-- ============================================================================
-- 10. 创建溶剂化结构表
-- ============================================================================

CREATE TABLE solvation_structures (
    id SERIAL PRIMARY KEY,
    md_job_id INTEGER NOT NULL REFERENCES md_jobs(id) ON DELETE CASCADE,

    -- 溶剂化结构信息
    center_ion VARCHAR(50) NOT NULL,
    structure_type VARCHAR(50) DEFAULT 'first_shell',
    coordination_num INTEGER,
    composition JSONB,

    -- 文件信息
    file_path TEXT,
    snapshot_frame INTEGER,
    description TEXT,

    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

COMMENT ON TABLE solvation_structures IS '溶剂化结构表';
COMMENT ON COLUMN solvation_structures.center_ion IS '中心离子，如 "Li", "K"';
COMMENT ON COLUMN solvation_structures.structure_type IS '结构类型，如 "first_shell", "typical_cluster"';
COMMENT ON COLUMN solvation_structures.composition IS '组成 JSON: {"EC": 3, "DMC": 1, "PF6": 0}';

-- ============================================================================
-- 11. 创建索引（优化查询性能）
-- ============================================================================

-- 用户表索引
CREATE INDEX idx_users_email ON users(email);

-- 项目表索引
CREATE INDEX idx_projects_user_id ON projects(user_id);

-- 体系表索引
CREATE INDEX idx_systems_hash_key ON electrolyte_systems(hash_key);
CREATE INDEX idx_systems_project_id ON electrolyte_systems(project_id);

-- 任务表索引
CREATE INDEX idx_jobs_user_id ON md_jobs(user_id);
CREATE INDEX idx_jobs_system_id ON md_jobs(system_id);
CREATE INDEX idx_jobs_status ON md_jobs(status);
CREATE INDEX idx_jobs_slurm_job_id ON md_jobs(slurm_job_id);
CREATE INDEX idx_jobs_created_at ON md_jobs(created_at DESC);

-- 后处理任务索引
CREATE INDEX idx_postprocess_md_job_id ON postprocess_jobs(md_job_id);
CREATE INDEX idx_postprocess_status ON postprocess_jobs(status);

-- 结果表索引
CREATE INDEX idx_result_summary_md_job_id ON result_summary(md_job_id);
CREATE INDEX idx_rdf_md_job_id ON rdf_results(md_job_id);
CREATE INDEX idx_rdf_species ON rdf_results(center_species, shell_species);
CREATE INDEX idx_msd_md_job_id ON msd_results(md_job_id);
CREATE INDEX idx_msd_species ON msd_results(species);
CREATE INDEX idx_solvation_md_job_id ON solvation_structures(md_job_id);

-- ============================================================================
-- 12. 创建更新时间戳触发器
-- ============================================================================

-- 创建更新时间戳函数
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- 为需要的表添加触发器
CREATE TRIGGER update_users_updated_at
    BEFORE UPDATE ON users
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_projects_updated_at
    BEFORE UPDATE ON projects
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_md_jobs_updated_at
    BEFORE UPDATE ON md_jobs
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- ============================================================================
-- 13. 插入测试数据（可选）
-- ============================================================================

-- 插入测试用户
INSERT INTO users (email, username, password_hash, role) VALUES
    ('admin@molyte.com', 'admin', '$2b$12$LQv3c1yqBWVHxkd0LHAkCOYz6TtxMQJqhN8/LewY5GyYzS8qB/Unu', 'admin'),
    ('test@molyte.com', 'testuser', '$2b$12$LQv3c1yqBWVHxkd0LHAkCOYz6TtxMQJqhN8/LewY5GyYzS8qB/Unu', 'user');

-- 插入测试项目
INSERT INTO projects (user_id, name, description) VALUES
    (1, 'Test Project 1', 'This is a test project for development'),
    (2, 'User Test Project', 'User test project');

-- 插入测试体系
INSERT INTO electrolyte_systems (
    project_id,
    hash_key,
    name,
    cations,
    anions,
    solvents,
    temperature,
    concentration
) VALUES (
    1,
    'test_hash_key_12345',
    '1M LiPF6 in EC:DMC (1:1)',
    '[{"name": "Li", "smiles": "[Li+]", "number": 50}]'::jsonb,
    '[{"name": "PF6", "smiles": "F[P-](F)(F)(F)(F)F", "number": 50}]'::jsonb,
    '[{"name": "EC", "smiles": "C1COC(=O)O1", "number": 100}, {"name": "DMC", "smiles": "COC(=O)OC", "number": 100}]'::jsonb,
    298.15,
    1.0
);

-- ============================================================================
-- 完成信息
-- ============================================================================

\echo '========================================='
\echo 'Database initialization completed!'
\echo '========================================='
\echo 'Created tables:'
\echo '  - users'
\echo '  - projects'
\echo '  - electrolyte_systems'
\echo '  - md_jobs'
\echo '  - postprocess_jobs'
\echo '  - result_summary'
\echo '  - rdf_results'
\echo '  - msd_results'
\echo '  - solvation_structures'
\echo ''
\echo 'Created indexes and triggers'
\echo 'Inserted test data'
\echo '========================================='

-- 显示表列表
\dt


