-- ============================================================================
-- 数据库测试脚本
-- 用途: 验证数据库表结构和基本功能
-- ============================================================================

\echo '========================================='
\echo 'Testing Molyte Web Database'
\echo '========================================='

-- ============================================================================
-- 1. 测试表是否存在
-- ============================================================================

\echo ''
\echo '1. Checking if all tables exist...'
\echo ''

SELECT 
    table_name,
    CASE 
        WHEN table_name IN (
            'users', 'projects', 'electrolyte_systems', 'md_jobs',
            'postprocess_jobs', 'result_summary', 'rdf_results',
            'msd_results', 'solvation_structures'
        ) THEN '✓ EXISTS'
        ELSE '✗ MISSING'
    END as status
FROM information_schema.tables
WHERE table_schema = 'public'
ORDER BY table_name;

-- ============================================================================
-- 2. 测试数据插入
-- ============================================================================

\echo ''
\echo '2. Testing data insertion...'
\echo ''

-- 测试用户数据
SELECT COUNT(*) as user_count FROM users;
SELECT email, username, role FROM users;

-- 测试项目数据
SELECT COUNT(*) as project_count FROM projects;
SELECT p.name, u.username as owner FROM projects p
JOIN users u ON p.user_id = u.id;

-- 测试体系数据
SELECT COUNT(*) as system_count FROM electrolyte_systems;
SELECT name, temperature, concentration FROM electrolyte_systems;

-- ============================================================================
-- 3. 测试外键关系
-- ============================================================================

\echo ''
\echo '3. Testing foreign key relationships...'
\echo ''

-- 测试用户-项目关系
SELECT 
    u.username,
    COUNT(p.id) as project_count
FROM users u
LEFT JOIN projects p ON u.id = p.user_id
GROUP BY u.username;

-- 测试项目-体系关系
SELECT 
    p.name as project_name,
    COUNT(es.id) as system_count
FROM projects p
LEFT JOIN electrolyte_systems es ON p.id = es.project_id
GROUP BY p.name;

-- ============================================================================
-- 4. 测试 JSONB 字段
-- ============================================================================

\echo ''
\echo '4. Testing JSONB fields...'
\echo ''

-- 查询体系的阳离子信息
SELECT 
    name,
    cations->0->>'name' as cation_name,
    cations->0->>'number' as cation_number
FROM electrolyte_systems;

-- 查询体系的溶剂信息
SELECT 
    name,
    jsonb_array_length(solvents) as solvent_count
FROM electrolyte_systems
WHERE solvents IS NOT NULL;

-- ============================================================================
-- 5. 测试索引
-- ============================================================================

\echo ''
\echo '5. Checking indexes...'
\echo ''

SELECT 
    tablename,
    indexname,
    indexdef
FROM pg_indexes
WHERE schemaname = 'public'
ORDER BY tablename, indexname;

-- ============================================================================
-- 6. 测试触发器
-- ============================================================================

\echo ''
\echo '6. Testing triggers (updated_at)...'
\echo ''

-- 更新用户并检查 updated_at
UPDATE users SET username = 'admin_updated' WHERE email = 'admin@molyte.com';
SELECT email, username, created_at, updated_at FROM users WHERE email = 'admin@molyte.com';

-- 恢复原值
UPDATE users SET username = 'admin' WHERE email = 'admin@molyte.com';

-- ============================================================================
-- 7. 测试约束
-- ============================================================================

\echo ''
\echo '7. Testing constraints...'
\echo ''

-- 测试唯一约束（应该失败）
\echo 'Testing unique constraint on email (should fail)...'
DO $$
BEGIN
    INSERT INTO users (email, username, password_hash) 
    VALUES ('admin@molyte.com', 'duplicate', 'hash');
    RAISE NOTICE 'ERROR: Unique constraint not working!';
EXCEPTION
    WHEN unique_violation THEN
        RAISE NOTICE '✓ Unique constraint working correctly';
END $$;

-- 测试检查约束（应该失败）
\echo 'Testing check constraint on role (should fail)...'
DO $$
BEGIN
    INSERT INTO users (email, username, password_hash, role) 
    VALUES ('invalid@test.com', 'test', 'hash', 'invalid_role');
    RAISE NOTICE 'ERROR: Check constraint not working!';
EXCEPTION
    WHEN check_violation THEN
        RAISE NOTICE '✓ Check constraint working correctly';
END $$;

-- ============================================================================
-- 8. 测试级联删除
-- ============================================================================

\echo ''
\echo '8. Testing cascade delete...'
\echo ''

-- 创建测试数据
INSERT INTO users (email, username, password_hash) 
VALUES ('cascade_test@test.com', 'cascade_test', 'hash')
RETURNING id as test_user_id \gset

INSERT INTO projects (user_id, name) 
VALUES (:test_user_id, 'Test Cascade Project')
RETURNING id as test_project_id \gset

-- 删除用户（应该级联删除项目）
DELETE FROM users WHERE id = :test_user_id;

-- 检查项目是否被删除
SELECT COUNT(*) as remaining_projects 
FROM projects 
WHERE id = :test_project_id;

\echo '✓ Cascade delete working (should show 0 remaining projects)'

-- ============================================================================
-- 9. 性能测试（简单查询）
-- ============================================================================

\echo ''
\echo '9. Performance test (simple queries)...'
\echo ''

\timing on

-- 测试索引查询
SELECT * FROM users WHERE email = 'admin@molyte.com';
SELECT * FROM md_jobs WHERE status = 'RUNNING';
SELECT * FROM electrolyte_systems WHERE hash_key = 'test_hash_key_12345';

\timing off

-- ============================================================================
-- 10. 总结
-- ============================================================================

\echo ''
\echo '========================================='
\echo 'Database Test Summary'
\echo '========================================='
\echo ''
\echo 'All tests completed!'
\echo 'Please review the output above for any errors.'
\echo ''
\echo 'Next steps:'
\echo '  1. Review test results'
\echo '  2. Fix any issues found'
\echo '  3. Proceed to Step 2: Backend setup'
\echo ''
\echo '========================================='


