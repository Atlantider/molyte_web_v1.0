-- ============================================================================
-- 添加QC任务可见性管理字段
-- Migration: Add visibility management fields to qc_jobs table
-- Date: 2025-11-30
-- ============================================================================

-- 1. 添加可见性相关字段
ALTER TABLE qc_jobs ADD COLUMN IF NOT EXISTS visibility VARCHAR(20) DEFAULT 'PRIVATE';
ALTER TABLE qc_jobs ADD COLUMN IF NOT EXISTS visibility_delay_until TIMESTAMP WITH TIME ZONE;

-- 2. 添加统计字段
ALTER TABLE qc_jobs ADD COLUMN IF NOT EXISTS view_count INTEGER DEFAULT 0;
ALTER TABLE qc_jobs ADD COLUMN IF NOT EXISTS download_count INTEGER DEFAULT 0;

-- 3. 创建索引以提高查询性能
CREATE INDEX IF NOT EXISTS idx_qc_jobs_visibility ON qc_jobs(visibility);
CREATE INDEX IF NOT EXISTS idx_qc_jobs_user_visibility ON qc_jobs(user_id, visibility);
CREATE INDEX IF NOT EXISTS idx_qc_jobs_visibility_delay ON qc_jobs(visibility_delay_until);

-- 4. 添加注释
COMMENT ON COLUMN qc_jobs.visibility IS '数据可见性：PUBLIC（公开）、DELAYED（延期公开）、PRIVATE（私有）';
COMMENT ON COLUMN qc_jobs.visibility_delay_until IS '延期公开日期（仅当visibility=DELAYED时有效）';
COMMENT ON COLUMN qc_jobs.view_count IS '查看次数统计';
COMMENT ON COLUMN qc_jobs.download_count IS '下载次数统计';

-- ============================================================================
-- 完成信息
-- ============================================================================
\echo '========================================='
\echo 'QC Visibility Migration Completed!'
\echo '========================================='
\echo 'Added columns:'
\echo '  - visibility (VARCHAR(20), default: PRIVATE)'
\echo '  - visibility_delay_until (TIMESTAMP WITH TIME ZONE)'
\echo '  - view_count (INTEGER, default: 0)'
\echo '  - download_count (INTEGER, default: 0)'
\echo '========================================='

