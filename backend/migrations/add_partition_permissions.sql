-- Add partition permissions to users table
-- This allows admins to control which queues/partitions each user can access

-- Add allowed_partitions column (JSON array of partition names)
ALTER TABLE users ADD COLUMN IF NOT EXISTS allowed_partitions JSON;

-- Set default values:
-- Admin users: NULL (can access all partitions)
-- Premium users: ["cpu", "gpu"] (can access cpu and gpu)
-- Regular users: ["cpu"] (can only access cpu)
-- Guest users: ["debug"] (can only access debug queue)

UPDATE users SET allowed_partitions = NULL WHERE role = 'ADMIN';
UPDATE users SET allowed_partitions = '["cpu", "gpu"]'::json WHERE role = 'PREMIUM' AND allowed_partitions IS NULL;
UPDATE users SET allowed_partitions = '["cpu"]'::json WHERE role = 'USER' AND allowed_partitions IS NULL;
UPDATE users SET allowed_partitions = '["debug"]'::json WHERE role = 'GUEST' AND allowed_partitions IS NULL;

-- Add comment
COMMENT ON COLUMN users.allowed_partitions IS 'JSON array of allowed Slurm partition names. NULL means all partitions (admin only).';

