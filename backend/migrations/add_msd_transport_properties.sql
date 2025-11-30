-- 添加 MSD 传输性质字段
-- 用于存储离子电导率、迁移率等计算结果

-- 添加离子电导率字段 (S/cm)
ALTER TABLE msd_results ADD COLUMN IF NOT EXISTS ionic_conductivity FLOAT;

-- 添加离子迁移率字段 (cm²/(V·s))
ALTER TABLE msd_results ADD COLUMN IF NOT EXISTS mobility FLOAT;

-- 添加离子电荷数字段
ALTER TABLE msd_results ADD COLUMN IF NOT EXISTS charge INTEGER;

-- 添加注释
COMMENT ON COLUMN msd_results.ionic_conductivity IS '离子电导率 (S/cm)，基于 Nernst-Einstein 方程计算';
COMMENT ON COLUMN msd_results.mobility IS '离子迁移率 (cm²/(V·s))';
COMMENT ON COLUMN msd_results.charge IS '离子电荷数';

