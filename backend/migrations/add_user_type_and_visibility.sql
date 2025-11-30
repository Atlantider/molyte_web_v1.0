-- 用户类型和数据展示控制迁移脚本
-- 执行: psql -U molyte_user -d molyte_db -f add_user_type_and_visibility.sql

-- 1. 添加用户类型枚举
DO $$ BEGIN
    CREATE TYPE usertype AS ENUM ('STUDENT', 'RESEARCHER', 'COMPANY');
EXCEPTION
    WHEN duplicate_object THEN null;
END $$;

-- 2. 添加数据展示状态枚举
DO $$ BEGIN
    CREATE TYPE datavisibility AS ENUM ('PRIVATE', 'DELAYED', 'PUBLIC', 'ADMIN_ONLY');
EXCEPTION
    WHEN duplicate_object THEN null;
END $$;

-- 3. 扩展用户表 - 用户类型和单位信息
ALTER TABLE users ADD COLUMN IF NOT EXISTS user_type usertype DEFAULT 'STUDENT' NOT NULL;
ALTER TABLE users ADD COLUMN IF NOT EXISTS organization VARCHAR(200);
ALTER TABLE users ADD COLUMN IF NOT EXISTS department VARCHAR(100);

-- 4. 扩展用户表 - 邮箱验证
ALTER TABLE users ADD COLUMN IF NOT EXISTS email_verified BOOLEAN DEFAULT FALSE NOT NULL;
ALTER TABLE users ADD COLUMN IF NOT EXISTS verification_token VARCHAR(100);
ALTER TABLE users ADD COLUMN IF NOT EXISTS verification_expires TIMESTAMP WITH TIME ZONE;

-- 5. 扩展用户表 - 贡献统计
ALTER TABLE users ADD COLUMN IF NOT EXISTS public_data_count INTEGER DEFAULT 0 NOT NULL;
ALTER TABLE users ADD COLUMN IF NOT EXISTS contribution_points FLOAT DEFAULT 0.0 NOT NULL;
ALTER TABLE users ADD COLUMN IF NOT EXISTS private_quota_used INTEGER DEFAULT 0 NOT NULL;
ALTER TABLE users ADD COLUMN IF NOT EXISTS free_cpu_hours_granted FLOAT DEFAULT 100.0 NOT NULL;

-- 6. 扩展任务表 - 数据展示控制
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS visibility datavisibility DEFAULT 'DELAYED' NOT NULL;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS visibility_delay_until TIMESTAMP WITH TIME ZONE;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS anonymous_public BOOLEAN DEFAULT FALSE;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS allow_download BOOLEAN DEFAULT TRUE;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS visibility_changed_by INTEGER REFERENCES users(id) ON DELETE SET NULL;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS visibility_changed_at TIMESTAMP WITH TIME ZONE;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS visibility_reason VARCHAR(500);

-- 7. 扩展任务表 - 贡献追踪
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS view_count INTEGER DEFAULT 0 NOT NULL;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS download_count INTEGER DEFAULT 0 NOT NULL;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS reward_claimed BOOLEAN DEFAULT FALSE;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS is_free_quota BOOLEAN DEFAULT TRUE;

-- 8. 创建索引
CREATE INDEX IF NOT EXISTS idx_jobs_visibility ON md_jobs(visibility);
CREATE INDEX IF NOT EXISTS idx_jobs_visibility_delay ON md_jobs(visibility_delay_until);
CREATE INDEX IF NOT EXISTS idx_users_user_type ON users(user_type);
CREATE INDEX IF NOT EXISTS idx_users_organization ON users(organization);

-- 9. 为现有已完成的任务设置默认延期公开时间（1年后）
UPDATE md_jobs 
SET visibility_delay_until = finished_at + INTERVAL '1 year'
WHERE status = 'COMPLETED' 
  AND visibility = 'DELAYED' 
  AND visibility_delay_until IS NULL 
  AND finished_at IS NOT NULL;

-- 10. 插入用户类型相关的系统配置
INSERT INTO system_configs (key, value, description) 
VALUES ('student_price_per_hour', '0.05', '学生用户每核时单价（元）')
ON CONFLICT (key) DO NOTHING;

INSERT INTO system_configs (key, value, description) 
VALUES ('researcher_price_per_hour', '0.08', '研究者用户每核时单价（元）')
ON CONFLICT (key) DO NOTHING;

INSERT INTO system_configs (key, value, description) 
VALUES ('company_price_per_hour', '0.15', '企业用户每核时单价（元）')
ON CONFLICT (key) DO NOTHING;

INSERT INTO system_configs (key, value, description) 
VALUES ('public_reward_cpu_hours', '10', '数据公开奖励核时')
ON CONFLICT (key) DO NOTHING;

INSERT INTO system_configs (key, value, description) 
VALUES ('view_reward_cpu_hours', '0.1', '数据被查看奖励核时（每次）')
ON CONFLICT (key) DO NOTHING;

INSERT INTO system_configs (key, value, description) 
VALUES ('daily_view_reward_limit', '5', '每日查看奖励上限核时')
ON CONFLICT (key) DO NOTHING;

INSERT INTO system_configs (key, value, description) 
VALUES ('download_reward_cpu_hours', '1', '数据被下载奖励核时（每次）')
ON CONFLICT (key) DO NOTHING;

-- 11. 更新现有用户的余额为100核时（免费赠送）
UPDATE users 
SET balance_cpu_hours = 100.0, free_cpu_hours_granted = 100.0
WHERE balance_cpu_hours < 100.0 AND role != 'ADMIN';

-- 完成
SELECT 'User type and visibility migration completed successfully!' as message;

