-- ============================================================================
-- QC (Quantum Chemistry) Tables Migration
-- 量子化学计算相关表
-- Created: 2024-01-XX
-- ============================================================================

-- ============================================================================
-- 1. QC任务表 (qc_jobs)
-- ============================================================================
CREATE TABLE IF NOT EXISTS qc_jobs (
    id SERIAL PRIMARY KEY,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    
    -- 关联信息（可选：独立任务时为NULL，MD附带时有值）
    md_job_id INTEGER REFERENCES md_jobs(id) ON DELETE SET NULL,
    
    -- 分子信息
    molecule_name VARCHAR(255) NOT NULL,
    smiles TEXT NOT NULL,
    molecule_type VARCHAR(20) DEFAULT 'custom',  -- 'solvent', 'cation', 'anion', 'custom'
    
    -- 计算参数
    basis_set VARCHAR(50) NOT NULL DEFAULT '6-31++g(d,p)',
    functional VARCHAR(50) DEFAULT 'B3LYP',
    charge INTEGER DEFAULT 0,
    spin_multiplicity INTEGER DEFAULT 1,
    
    -- 额外配置
    config JSONB DEFAULT '{}',
    
    -- 任务状态
    status VARCHAR(20) NOT NULL DEFAULT 'CREATED'
        CHECK (status IN ('CREATED', 'QUEUED', 'RUNNING', 'POSTPROCESSING', 'COMPLETED', 'FAILED', 'CANCELLED')),
    slurm_job_id VARCHAR(50),
    progress FLOAT DEFAULT 0 CHECK (progress >= 0 AND progress <= 100),
    work_dir TEXT,
    log_file TEXT,
    error_message TEXT,
    
    -- 时间戳
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP,
    finished_at TIMESTAMP
);

-- 索引
CREATE INDEX IF NOT EXISTS idx_qc_jobs_user_id ON qc_jobs(user_id);
CREATE INDEX IF NOT EXISTS idx_qc_jobs_md_job_id ON qc_jobs(md_job_id);
CREATE INDEX IF NOT EXISTS idx_qc_jobs_status ON qc_jobs(status);
CREATE INDEX IF NOT EXISTS idx_qc_jobs_smiles ON qc_jobs(smiles);
CREATE INDEX IF NOT EXISTS idx_qc_jobs_created_at ON qc_jobs(created_at DESC);

-- ============================================================================
-- 2. QC计算结果表 (qc_results)
-- ============================================================================
CREATE TABLE IF NOT EXISTS qc_results (
    id SERIAL PRIMARY KEY,
    qc_job_id INTEGER NOT NULL REFERENCES qc_jobs(id) ON DELETE CASCADE,
    smiles TEXT NOT NULL,  -- 冗余存储方便搜索
    
    -- 核心结果
    energy_au FLOAT,                    -- 能量 (A.U.)
    homo FLOAT,                         -- HOMO 能量 (Hartree)
    lumo FLOAT,                         -- LUMO 能量 (Hartree)
    homo_lumo_gap FLOAT,               -- HOMO-LUMO 能隙 (eV, 自动计算)
    
    -- ESP 相关
    esp_min_kcal FLOAT,                -- ESP最小值 (kcal/mol)
    esp_max_kcal FLOAT,                -- ESP最大值 (kcal/mol)
    
    -- 文件路径
    esp_image_path TEXT,               -- ESP可视化图像路径
    fchk_file_path TEXT,               -- .fchk 文件路径
    log_file_path TEXT,                -- Gaussian日志文件路径
    cube_density_path TEXT,            -- density.cub 文件路径
    cube_esp_path TEXT,                -- ESP.cub 文件路径
    
    -- 其他可扩展数据
    dipole_moment FLOAT,               -- 偶极矩
    polarizability FLOAT,              -- 极化率
    additional_properties JSONB DEFAULT '{}',
    
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- 索引
CREATE INDEX IF NOT EXISTS idx_qc_results_qc_job_id ON qc_results(qc_job_id);
CREATE INDEX IF NOT EXISTS idx_qc_results_smiles ON qc_results(smiles);

-- ============================================================================
-- 3. 分子QC缓存表 (molecule_qc_cache)
-- ============================================================================
CREATE TABLE IF NOT EXISTS molecule_qc_cache (
    id SERIAL PRIMARY KEY,
    smiles TEXT UNIQUE NOT NULL,
    molecule_name VARCHAR(255),
    
    -- 默认展示的QC结果
    preferred_qc_result_id INTEGER REFERENCES qc_results(id) ON DELETE SET NULL,
    basis_set VARCHAR(50),
    functional VARCHAR(50),
    
    -- 缓存的核心数据 (eV单位，方便展示)
    energy_au FLOAT,
    homo_ev FLOAT,
    lumo_ev FLOAT,
    homo_lumo_gap_ev FLOAT,
    esp_min_kcal FLOAT,
    esp_max_kcal FLOAT,
    
    -- ESP图像路径
    esp_image_path TEXT,
    
    -- 计算次数统计
    calculation_count INTEGER DEFAULT 1,
    
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- 索引
CREATE INDEX IF NOT EXISTS idx_molecule_qc_cache_smiles ON molecule_qc_cache(smiles);

-- ============================================================================
-- 4. 触发器：自动更新 updated_at
-- ============================================================================
CREATE OR REPLACE FUNCTION update_qc_jobs_updated_at()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER trigger_qc_jobs_updated_at
    BEFORE UPDATE ON qc_jobs
    FOR EACH ROW
    EXECUTE FUNCTION update_qc_jobs_updated_at();

CREATE OR REPLACE FUNCTION update_molecule_qc_cache_updated_at()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER trigger_molecule_qc_cache_updated_at
    BEFORE UPDATE ON molecule_qc_cache
    FOR EACH ROW
    EXECUTE FUNCTION update_molecule_qc_cache_updated_at();

-- ============================================================================
-- 完成信息
-- ============================================================================
\echo '========================================='
\echo 'QC Tables Migration Completed!'
\echo '========================================='
\echo 'Created tables:'
\echo '  - qc_jobs'
\echo '  - qc_results'
\echo '  - molecule_qc_cache'
\echo '========================================='

\dt qc_*
\dt molecule_qc_cache

