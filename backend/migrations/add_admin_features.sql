-- Migration: Add admin features and user quotas
-- Date: 2025-11-27

-- 1. Add GUEST role to UserRole enum
ALTER TYPE userrole ADD VALUE IF NOT EXISTS 'GUEST';

-- 2. Add new columns to users table
ALTER TABLE users ADD COLUMN IF NOT EXISTS is_active BOOLEAN DEFAULT true NOT NULL;
ALTER TABLE users ADD COLUMN IF NOT EXISTS last_login_at TIMESTAMP WITH TIME ZONE;
ALTER TABLE users ADD COLUMN IF NOT EXISTS total_cpu_hours DOUBLE PRECISION DEFAULT 100.0 NOT NULL;
ALTER TABLE users ADD COLUMN IF NOT EXISTS daily_job_limit INTEGER DEFAULT 10 NOT NULL;
ALTER TABLE users ADD COLUMN IF NOT EXISTS concurrent_job_limit INTEGER DEFAULT 3 NOT NULL;
ALTER TABLE users ADD COLUMN IF NOT EXISTS storage_quota_gb DOUBLE PRECISION DEFAULT 10.0 NOT NULL;

-- 3. Create user_usage_stats table
CREATE TABLE IF NOT EXISTS user_usage_stats (
    id SERIAL PRIMARY KEY,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    date DATE NOT NULL,
    jobs_submitted INTEGER DEFAULT 0 NOT NULL,
    jobs_completed INTEGER DEFAULT 0 NOT NULL,
    jobs_failed INTEGER DEFAULT 0 NOT NULL,
    jobs_cancelled INTEGER DEFAULT 0 NOT NULL,
    cpu_hours_used DOUBLE PRECISION DEFAULT 0.0 NOT NULL,
    storage_used_gb DOUBLE PRECISION DEFAULT 0.0 NOT NULL,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP NOT NULL,
    CONSTRAINT uq_user_date UNIQUE (user_id, date)
);

-- Create indexes for user_usage_stats
CREATE INDEX IF NOT EXISTS idx_user_usage_stats_user_id ON user_usage_stats(user_id);
CREATE INDEX IF NOT EXISTS idx_user_usage_stats_date ON user_usage_stats(date);

-- 4. Create audit_logs table
CREATE TABLE IF NOT EXISTS audit_logs (
    id SERIAL PRIMARY KEY,
    user_id INTEGER REFERENCES users(id) ON DELETE SET NULL,
    action VARCHAR(100) NOT NULL,
    resource_type VARCHAR(50),
    resource_id INTEGER,
    details JSONB,
    ip_address VARCHAR(50),
    user_agent VARCHAR(500),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP NOT NULL
);

-- Create indexes for audit_logs
CREATE INDEX IF NOT EXISTS idx_audit_logs_user_id ON audit_logs(user_id);
CREATE INDEX IF NOT EXISTS idx_audit_logs_action ON audit_logs(action);
CREATE INDEX IF NOT EXISTS idx_audit_logs_created_at ON audit_logs(created_at);

-- 5. Update existing users with default quotas based on role
UPDATE users SET 
    total_cpu_hours = CASE 
        WHEN role = 'ADMIN' THEN 999999.0
        WHEN role = 'PREMIUM' THEN 500.0
        ELSE 100.0
    END,
    daily_job_limit = CASE 
        WHEN role = 'ADMIN' THEN 999999
        WHEN role = 'PREMIUM' THEN 50
        ELSE 10
    END,
    concurrent_job_limit = CASE 
        WHEN role = 'ADMIN' THEN 999999
        WHEN role = 'PREMIUM' THEN 10
        ELSE 3
    END,
    storage_quota_gb = CASE 
        WHEN role = 'ADMIN' THEN 999999.0
        WHEN role = 'PREMIUM' THEN 50.0
        ELSE 10.0
    END
WHERE total_cpu_hours = 100.0;  -- Only update if not already customized

-- 6. Create function to update updated_at timestamp
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ language 'plpgsql';

-- 7. Create trigger for user_usage_stats
DROP TRIGGER IF EXISTS update_user_usage_stats_updated_at ON user_usage_stats;
CREATE TRIGGER update_user_usage_stats_updated_at
    BEFORE UPDATE ON user_usage_stats
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- 8. Insert audit log for migration
INSERT INTO audit_logs (user_id, action, resource_type, details)
VALUES (NULL, 'database_migration', 'system', '{"migration": "add_admin_features", "date": "2025-11-27"}');

COMMIT;

