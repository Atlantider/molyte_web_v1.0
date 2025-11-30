-- 用户偏好设置表迁移脚本
-- 创建用户自定义溶剂组合和离子组合表

-- 1. 创建用户自定义溶剂组合表
CREATE TABLE IF NOT EXISTS user_solvent_combinations (
    id SERIAL PRIMARY KEY,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    name VARCHAR(100) NOT NULL,
    description TEXT,
    solvents JSONB NOT NULL,  -- [{"name": "EC", "smiles": "C1COC(=O)O1", "molar_ratio": 1.0}, ...]
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    
    -- 约束：同一用户的组合名称不能重复
    CONSTRAINT unique_user_solvent_combination UNIQUE (user_id, name)
);

-- 创建索引
CREATE INDEX idx_user_solvent_combinations_user_id ON user_solvent_combinations(user_id);
CREATE INDEX idx_user_solvent_combinations_created_at ON user_solvent_combinations(created_at DESC);

-- 添加注释
COMMENT ON TABLE user_solvent_combinations IS '用户自定义溶剂组合';
COMMENT ON COLUMN user_solvent_combinations.name IS '组合名称，如"我的EC/DMC"';
COMMENT ON COLUMN user_solvent_combinations.solvents IS '溶剂列表，JSON格式';


-- 2. 创建用户自定义离子组合表
CREATE TABLE IF NOT EXISTS user_ion_combinations (
    id SERIAL PRIMARY KEY,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    name VARCHAR(100) NOT NULL,
    description TEXT,
    cations JSONB NOT NULL,  -- [{"name": "Li", "charge": 1, "concentration": 1.0}, ...]
    anions JSONB NOT NULL,   -- [{"name": "PF6", "charge": -1, "concentration": 1.0}, ...]
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    
    -- 约束：同一用户的组合名称不能重复
    CONSTRAINT unique_user_ion_combination UNIQUE (user_id, name)
);

-- 创建索引
CREATE INDEX idx_user_ion_combinations_user_id ON user_ion_combinations(user_id);
CREATE INDEX idx_user_ion_combinations_created_at ON user_ion_combinations(created_at DESC);

-- 添加注释
COMMENT ON TABLE user_ion_combinations IS '用户自定义离子组合';
COMMENT ON COLUMN user_ion_combinations.name IS '组合名称，如"标准LiPF6"';
COMMENT ON COLUMN user_ion_combinations.cations IS '阳离子列表，JSON格式';
COMMENT ON COLUMN user_ion_combinations.anions IS '阴离子列表，JSON格式';


-- 3. 更新 users 表，添加关系（如果需要）
-- 注意：SQLAlchemy 的 relationship 不需要在数据库层面创建额外字段
-- 这里只是为了文档说明

-- 示例数据（可选）
-- INSERT INTO user_solvent_combinations (user_id, name, description, solvents) VALUES
-- (1, '我的EC/DMC', '常用的碳酸酯组合', '[{"name": "EC", "smiles": "C1COC(=O)O1", "molar_ratio": 1.0}, {"name": "DMC", "smiles": "COC(=O)OC", "molar_ratio": 1.0}]'),
-- (1, '高电压FEC/DMC', '高电压电池用', '[{"name": "FEC", "smiles": "C1C(OC(=O)O1)F", "molar_ratio": 1.0}, {"name": "DMC", "smiles": "COC(=O)OC", "molar_ratio": 2.0}]');

-- INSERT INTO user_ion_combinations (user_id, name, description, cations, anions) VALUES
-- (1, '标准LiPF6', '1M LiPF6', '[{"name": "Li", "charge": 1, "concentration": 1.0}]', '[{"name": "PF6", "charge": -1, "concentration": 1.0}]'),
-- (1, '高浓度LiTFSI', '2M LiTFSI', '[{"name": "Li", "charge": 1, "concentration": 2.0}]', '[{"name": "TFSI", "charge": -1, "concentration": 2.0}]');

