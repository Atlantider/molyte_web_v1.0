-- 软删除机制迁移脚本
-- 为MD任务和QC任务添加软删除支持
-- 删除时只标记为已删除，不真正删除数据
-- 管理员可以看到所有数据，普通用户看不到自己已删除的任务

-- ============================================================================
-- 1. MD Jobs 添加软删除字段
-- ============================================================================

-- 添加软删除相关字段
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS is_deleted BOOLEAN DEFAULT FALSE;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS deleted_at TIMESTAMP WITH TIME ZONE;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS deleted_by INTEGER REFERENCES users(id) ON DELETE SET NULL;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS delete_reason VARCHAR(500);

-- 创建索引以加速查询（排除已删除的数据）
CREATE INDEX IF NOT EXISTS idx_md_jobs_is_deleted ON md_jobs(is_deleted);
CREATE INDEX IF NOT EXISTS idx_md_jobs_deleted_at ON md_jobs(deleted_at);

-- 创建部分索引优化查询（只索引未删除的数据）
CREATE INDEX IF NOT EXISTS idx_md_jobs_active ON md_jobs(id) WHERE is_deleted = FALSE;

COMMENT ON COLUMN md_jobs.is_deleted IS '是否已删除（软删除标记）';
COMMENT ON COLUMN md_jobs.deleted_at IS '删除时间';
COMMENT ON COLUMN md_jobs.deleted_by IS '删除操作者ID';
COMMENT ON COLUMN md_jobs.delete_reason IS '删除原因';

-- ============================================================================
-- 2. QC Jobs 添加软删除字段
-- ============================================================================

ALTER TABLE qc_jobs ADD COLUMN IF NOT EXISTS is_deleted BOOLEAN DEFAULT FALSE;
ALTER TABLE qc_jobs ADD COLUMN IF NOT EXISTS deleted_at TIMESTAMP WITH TIME ZONE;
ALTER TABLE qc_jobs ADD COLUMN IF NOT EXISTS deleted_by INTEGER REFERENCES users(id) ON DELETE SET NULL;
ALTER TABLE qc_jobs ADD COLUMN IF NOT EXISTS delete_reason VARCHAR(500);

-- 创建索引
CREATE INDEX IF NOT EXISTS idx_qc_jobs_is_deleted ON qc_jobs(is_deleted);
CREATE INDEX IF NOT EXISTS idx_qc_jobs_deleted_at ON qc_jobs(deleted_at);
CREATE INDEX IF NOT EXISTS idx_qc_jobs_active ON qc_jobs(id) WHERE is_deleted = FALSE;

COMMENT ON COLUMN qc_jobs.is_deleted IS '是否已删除（软删除标记）';
COMMENT ON COLUMN qc_jobs.deleted_at IS '删除时间';
COMMENT ON COLUMN qc_jobs.deleted_by IS '删除操作者ID';
COMMENT ON COLUMN qc_jobs.delete_reason IS '删除原因';

-- ============================================================================
-- 3. 电解质配方添加软删除字段
-- ============================================================================

ALTER TABLE electrolytes ADD COLUMN IF NOT EXISTS is_deleted BOOLEAN DEFAULT FALSE;
ALTER TABLE electrolytes ADD COLUMN IF NOT EXISTS deleted_at TIMESTAMP WITH TIME ZONE;
ALTER TABLE electrolytes ADD COLUMN IF NOT EXISTS deleted_by INTEGER REFERENCES users(id) ON DELETE SET NULL;
ALTER TABLE electrolytes ADD COLUMN IF NOT EXISTS delete_reason VARCHAR(500);

CREATE INDEX IF NOT EXISTS idx_electrolytes_is_deleted ON electrolytes(is_deleted);
CREATE INDEX IF NOT EXISTS idx_electrolytes_active ON electrolytes(id) WHERE is_deleted = FALSE;

COMMENT ON COLUMN electrolytes.is_deleted IS '是否已删除（软删除标记）';
COMMENT ON COLUMN electrolytes.deleted_at IS '删除时间';
COMMENT ON COLUMN electrolytes.deleted_by IS '删除操作者ID';
COMMENT ON COLUMN electrolytes.delete_reason IS '删除原因';

-- ============================================================================
-- 4. 项目添加软删除字段
-- ============================================================================

ALTER TABLE projects ADD COLUMN IF NOT EXISTS is_deleted BOOLEAN DEFAULT FALSE;
ALTER TABLE projects ADD COLUMN IF NOT EXISTS deleted_at TIMESTAMP WITH TIME ZONE;
ALTER TABLE projects ADD COLUMN IF NOT EXISTS deleted_by INTEGER REFERENCES users(id) ON DELETE SET NULL;
ALTER TABLE projects ADD COLUMN IF NOT EXISTS delete_reason VARCHAR(500);

CREATE INDEX IF NOT EXISTS idx_projects_is_deleted ON projects(is_deleted);
CREATE INDEX IF NOT EXISTS idx_projects_active ON projects(id) WHERE is_deleted = FALSE;

-- ============================================================================
-- 5. 创建已删除数据恢复视图（管理员使用）
-- ============================================================================

CREATE OR REPLACE VIEW deleted_md_jobs_view AS
SELECT 
    mj.id,
    mj.system_id,
    mj.user_id,
    u.username,
    mj.status,
    mj.deleted_at,
    mj.delete_reason,
    du.username as deleted_by_username,
    mj.config,
    mj.finished_at
FROM md_jobs mj
LEFT JOIN users u ON mj.user_id = u.id
LEFT JOIN users du ON mj.deleted_by = du.id
WHERE mj.is_deleted = TRUE
ORDER BY mj.deleted_at DESC;

CREATE OR REPLACE VIEW deleted_qc_jobs_view AS
SELECT 
    qj.id,
    qj.user_id,
    u.username,
    qj.molecule_name,
    qj.smiles,
    qj.status,
    qj.deleted_at,
    qj.delete_reason,
    du.username as deleted_by_username,
    qj.config,
    qj.finished_at
FROM qc_jobs qj
LEFT JOIN users u ON qj.user_id = u.id
LEFT JOIN users du ON qj.deleted_by = du.id
WHERE qj.is_deleted = TRUE
ORDER BY qj.deleted_at DESC;

-- ============================================================================
-- 6. 创建统计函数
-- ============================================================================

-- 统计已删除任务数
CREATE OR REPLACE FUNCTION count_deleted_jobs()
RETURNS TABLE(
    deleted_md_jobs BIGINT,
    deleted_qc_jobs BIGINT,
    deleted_electrolytes BIGINT,
    deleted_projects BIGINT
) AS $$
BEGIN
    RETURN QUERY
    SELECT 
        (SELECT COUNT(*) FROM md_jobs WHERE is_deleted = TRUE),
        (SELECT COUNT(*) FROM qc_jobs WHERE is_deleted = TRUE),
        (SELECT COUNT(*) FROM electrolytes WHERE is_deleted = TRUE),
        (SELECT COUNT(*) FROM projects WHERE is_deleted = TRUE);
END;
$$ LANGUAGE plpgsql;

COMMENT ON FUNCTION count_deleted_jobs() IS '统计所有表中软删除的记录数';

-- ============================================================================
-- 完成
-- ============================================================================

SELECT 'Soft delete migration completed successfully!' as status;

