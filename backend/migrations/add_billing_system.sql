-- 计费系统数据库迁移脚本
-- 执行: psql -U molyte_user -d molyte_db -f add_billing_system.sql

-- 1. 添加支付方式枚举
DO $$ BEGIN
    CREATE TYPE paymentmethod AS ENUM ('wechat', 'alipay', 'admin', 'simulated');
EXCEPTION
    WHEN duplicate_object THEN null;
END $$;

-- 2. 添加支付状态枚举
DO $$ BEGIN
    CREATE TYPE paymentstatus AS ENUM ('pending', 'paid', 'failed', 'cancelled', 'refunded');
EXCEPTION
    WHEN duplicate_object THEN null;
END $$;

-- 3. 添加交易类型枚举
DO $$ BEGIN
    CREATE TYPE transactiontype AS ENUM ('recharge', 'consume', 'refund', 'admin_adjust', 'debt_repay');
EXCEPTION
    WHEN duplicate_object THEN null;
END $$;

-- 4. 创建系统配置表
CREATE TABLE IF NOT EXISTS system_configs (
    id SERIAL PRIMARY KEY,
    key VARCHAR(100) UNIQUE NOT NULL,
    value TEXT NOT NULL,
    description VARCHAR(500),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_by INTEGER REFERENCES users(id) ON DELETE SET NULL
);
CREATE INDEX IF NOT EXISTS idx_system_configs_key ON system_configs(key);

-- 5. 创建充值订单表
CREATE TABLE IF NOT EXISTS recharge_orders (
    id SERIAL PRIMARY KEY,
    order_no VARCHAR(64) UNIQUE NOT NULL,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    amount FLOAT NOT NULL,
    cpu_hours FLOAT NOT NULL,
    price_per_hour FLOAT NOT NULL,
    payment_method paymentmethod NOT NULL,
    payment_status paymentstatus DEFAULT 'pending' NOT NULL,
    transaction_id VARCHAR(128),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW() NOT NULL,
    paid_at TIMESTAMP WITH TIME ZONE,
    expired_at TIMESTAMP WITH TIME ZONE,
    remark VARCHAR(500)
);
CREATE INDEX IF NOT EXISTS idx_recharge_orders_order_no ON recharge_orders(order_no);
CREATE INDEX IF NOT EXISTS idx_recharge_orders_user_id ON recharge_orders(user_id);

-- 6. 创建配额变更流水表
CREATE TABLE IF NOT EXISTS quota_transactions (
    id SERIAL PRIMARY KEY,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    type transactiontype NOT NULL,
    amount FLOAT NOT NULL,
    balance_before FLOAT NOT NULL,
    balance_after FLOAT NOT NULL,
    reference_id INTEGER,
    reference_type VARCHAR(50),
    description VARCHAR(500),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW() NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_quota_transactions_user_id ON quota_transactions(user_id);

-- 7. 扩展用户表，添加余额字段
ALTER TABLE users ADD COLUMN IF NOT EXISTS balance_cpu_hours FLOAT DEFAULT 10.0 NOT NULL;
ALTER TABLE users ADD COLUMN IF NOT EXISTS frozen_cpu_hours FLOAT DEFAULT 0.0 NOT NULL;
ALTER TABLE users ADD COLUMN IF NOT EXISTS debt_cpu_hours FLOAT DEFAULT 0.0 NOT NULL;

-- 8. 扩展任务表，添加计费相关字段
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS cpu_cores INTEGER DEFAULT 1;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS estimated_cpu_hours FLOAT DEFAULT 0.0;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS actual_cpu_hours FLOAT DEFAULT 0.0;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS result_locked BOOLEAN DEFAULT FALSE;
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS locked_reason VARCHAR(200);
ALTER TABLE md_jobs ADD COLUMN IF NOT EXISTS billed BOOLEAN DEFAULT FALSE;

-- 9. 插入默认系统配置
INSERT INTO system_configs (key, value, description) 
VALUES ('cpu_hour_price', '0.1', '每核时单价（元）')
ON CONFLICT (key) DO NOTHING;

INSERT INTO system_configs (key, value, description) 
VALUES ('min_recharge_amount', '10', '最低充值金额（元）')
ON CONFLICT (key) DO NOTHING;

INSERT INTO system_configs (key, value, description) 
VALUES ('max_debt_cpu_hours', '100', '最大允许欠费机时')
ON CONFLICT (key) DO NOTHING;

-- 完成
SELECT 'Billing system migration completed successfully!' as message;

