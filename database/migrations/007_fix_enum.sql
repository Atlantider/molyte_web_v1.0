-- 修复 postprocesstype 枚举类型，添加 DESOLVATION_ENERGY
-- 这个脚本需要在 007_add_desolvation_energy.sql 之前或同时执行

-- 检查当前枚举值
-- SELECT enum_range(NULL::postprocesstype);

-- 添加新的枚举值（如果不存在）
DO $$
BEGIN
    -- 检查枚举值是否已存在
    IF NOT EXISTS (
        SELECT 1 FROM pg_enum 
        WHERE enumlabel = 'DESOLVATION_ENERGY' 
        AND enumtypid = (SELECT oid FROM pg_type WHERE typname = 'postprocesstype')
    ) THEN
        -- 添加新的枚举值
        ALTER TYPE postprocesstype ADD VALUE 'DESOLVATION_ENERGY';
        RAISE NOTICE 'Added DESOLVATION_ENERGY to postprocesstype enum';
    ELSE
        RAISE NOTICE 'DESOLVATION_ENERGY already exists in postprocesstype enum';
    END IF;
END$$;

-- 验证枚举值
SELECT enum_range(NULL::postprocesstype);

COMMENT ON TYPE postprocesstype IS 'Postprocess job types: RDF, MSD, SOLVATION, DESOLVATION_ENERGY';

