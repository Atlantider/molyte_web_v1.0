-- ============================================================================
-- 添加 VIP/VEA 相关字段到 qc_results 表
-- Migration: Add VIP/VEA fields for electrochemical window estimation
-- Date: 2024-12-04
-- ============================================================================

-- 1. 添加 VIP/VEA 相关字段到 qc_results 表
ALTER TABLE qc_results ADD COLUMN IF NOT EXISTS vip_ev FLOAT;
ALTER TABLE qc_results ADD COLUMN IF NOT EXISTS vea_ev FLOAT;
ALTER TABLE qc_results ADD COLUMN IF NOT EXISTS oxidation_potential_v FLOAT;
ALTER TABLE qc_results ADD COLUMN IF NOT EXISTS reduction_potential_v FLOAT;

-- 2. 添加注释
COMMENT ON COLUMN qc_results.vip_ev IS '垂直电离势 VIP (eV)';
COMMENT ON COLUMN qc_results.vea_ev IS '垂直电子亲和能 VEA (eV)';
COMMENT ON COLUMN qc_results.oxidation_potential_v IS '氧化电位 vs Li/Li+ (V)';
COMMENT ON COLUMN qc_results.reduction_potential_v IS '还原电位 vs Li/Li+ (V)';

-- ============================================================================
-- 完成信息
-- ============================================================================
\echo '========================================='
\echo 'VIP/VEA Fields Migration Completed!'
\echo '========================================='
\echo 'Added columns to qc_results:'
\echo '  - vip_ev (FLOAT)'
\echo '  - vea_ev (FLOAT)'
\echo '  - oxidation_potential_v (FLOAT)'
\echo '  - reduction_potential_v (FLOAT)'
\echo '========================================='

