# Molyte Web 数据库文档

## 📋 数据库概览

- **数据库名称**: `molyte_web`
- **用户**: `molyte_user`
- **版本**: 1.0 (MVP)
- **表数量**: 9 个核心表

---

## 🗂️ 表结构

### 1. users - 用户表
存储用户账户信息

| 字段 | 类型 | 说明 |
|------|------|------|
| id | SERIAL | 主键 |
| email | VARCHAR(255) | 邮箱（唯一） |
| username | VARCHAR(100) | 用户名 |
| password_hash | VARCHAR(255) | 密码哈希 |
| role | VARCHAR(20) | 角色（user/admin） |
| created_at | TIMESTAMP | 创建时间 |
| updated_at | TIMESTAMP | 更新时间 |

### 2. projects - 项目表
存储用户创建的项目

| 字段 | 类型 | 说明 |
|------|------|------|
| id | SERIAL | 主键 |
| user_id | INTEGER | 用户 ID（外键） |
| name | VARCHAR(255) | 项目名称 |
| description | TEXT | 项目描述 |
| created_at | TIMESTAMP | 创建时间 |
| updated_at | TIMESTAMP | 更新时间 |

### 3. electrolyte_systems - 电解液体系表
存储规范化的电解液配方

| 字段 | 类型 | 说明 |
|------|------|------|
| id | SERIAL | 主键 |
| project_id | INTEGER | 项目 ID（外键） |
| hash_key | VARCHAR(64) | 配方哈希（唯一，用于去重） |
| name | VARCHAR(255) | 体系名称 |
| cations | JSONB | 阳离子列表 JSON |
| anions | JSONB | 阴离子列表 JSON |
| solvents | JSONB | 溶剂列表 JSON |
| temperature | FLOAT | 温度（K） |
| pressure | FLOAT | 压力（atm） |
| density | FLOAT | 密度（g/cm³） |
| concentration | FLOAT | 浓度（M） |
| nsteps_npt | INTEGER | NPT 步数 |
| nsteps_nvt | INTEGER | NVT 步数 |
| timestep | FLOAT | 时间步长（fs） |
| force_field | VARCHAR(50) | 力场类型 |

**JSON 格式示例**：
```json
{
  "cations": [{"name": "Li", "smiles": "[Li+]", "number": 50}],
  "anions": [{"name": "PF6", "smiles": "F[P-](F)(F)(F)(F)F", "number": 50}],
  "solvents": [
    {"name": "EC", "smiles": "C1COC(=O)O1", "number": 100},
    {"name": "DMC", "smiles": "COC(=O)OC", "number": 100}
  ]
}
```

### 4. md_jobs - MD 任务表
存储分子动力学模拟任务

| 字段 | 类型 | 说明 |
|------|------|------|
| id | SERIAL | 主键 |
| system_id | INTEGER | 体系 ID（外键） |
| user_id | INTEGER | 用户 ID（外键） |
| status | VARCHAR(20) | 任务状态 |
| slurm_job_id | VARCHAR(50) | Slurm 任务 ID |
| progress | INTEGER | 进度（0-100） |
| work_dir | TEXT | 工作目录路径 |
| log_file | TEXT | 日志文件路径 |
| error_message | TEXT | 错误信息 |
| created_at | TIMESTAMP | 创建时间 |
| updated_at | TIMESTAMP | 更新时间 |
| started_at | TIMESTAMP | 开始时间 |
| finished_at | TIMESTAMP | 完成时间 |

**状态值**：
- `CREATED` - 已创建
- `QUEUED` - 已提交到队列
- `RUNNING` - 运行中
- `POSTPROCESSING` - 后处理中
- `COMPLETED` - 已完成
- `FAILED` - 失败

### 5. postprocess_jobs - 后处理任务表

| 字段 | 类型 | 说明 |
|------|------|------|
| id | SERIAL | 主键 |
| md_job_id | INTEGER | MD 任务 ID（外键） |
| status | VARCHAR(20) | 状态 |
| required_results | JSONB | 需要计算的结果类型 |
| error_message | TEXT | 错误信息 |
| created_at | TIMESTAMP | 创建时间 |
| finished_at | TIMESTAMP | 完成时间 |

### 6. result_summary - 结果概览表

| 字段 | 类型 | 说明 |
|------|------|------|
| id | SERIAL | 主键 |
| md_job_id | INTEGER | MD 任务 ID（外键，唯一） |
| final_density | FLOAT | 最终密度 |
| avg_temperature | FLOAT | 平均温度 |
| avg_pressure | FLOAT | 平均压力 |
| avg_energy | FLOAT | 平均能量 |
| box_volume | FLOAT | 盒子体积 |
| available_results | JSONB | 可用结果类型 |

### 7. rdf_results - RDF 结果表

| 字段 | 类型 | 说明 |
|------|------|------|
| id | SERIAL | 主键 |
| md_job_id | INTEGER | MD 任务 ID（外键） |
| center_species | VARCHAR(100) | 中心原子 |
| shell_species | VARCHAR(100) | 壳层原子 |
| r_values | JSONB | 距离数组 |
| g_r_values | JSONB | g(r) 值数组 |
| coordination_number | JSONB | 配位数数组 |
| cutoff | FLOAT | 截断距离 |
| bin_width | FLOAT | 分箱宽度 |
| file_path | TEXT | 数据文件路径 |

### 8. msd_results - MSD 结果表

| 字段 | 类型 | 说明 |
|------|------|------|
| id | SERIAL | 主键 |
| md_job_id | INTEGER | MD 任务 ID（外键） |
| species | VARCHAR(100) | 粒子类型 |
| t_values | JSONB | 时间数组 |
| msd_values | JSONB | MSD 值数组 |
| diffusion_coeff | FLOAT | 扩散系数 |
| fit_range | JSONB | 拟合范围 |
| file_path | TEXT | 数据文件路径 |

### 9. solvation_structures - 溶剂化结构表

| 字段 | 类型 | 说明 |
|------|------|------|
| id | SERIAL | 主键 |
| md_job_id | INTEGER | MD 任务 ID（外键） |
| center_ion | VARCHAR(50) | 中心离子 |
| structure_type | VARCHAR(50) | 结构类型 |
| coordination_num | INTEGER | 配位数 |
| composition | JSONB | 组成 |
| file_path | TEXT | 结构文件路径 |
| snapshot_frame | INTEGER | 快照帧号 |
| description | TEXT | 描述 |

---

## 🔗 表关系

```
users (1) ──→ (N) projects
projects (1) ──→ (N) electrolyte_systems
electrolyte_systems (1) ──→ (N) md_jobs
users (1) ──→ (N) md_jobs

md_jobs (1) ──→ (1) postprocess_jobs
md_jobs (1) ──→ (1) result_summary
md_jobs (1) ──→ (N) rdf_results
md_jobs (1) ──→ (N) msd_results
md_jobs (1) ──→ (N) solvation_structures
```

---

## 🚀 使用方法

### 初始化数据库

```bash
# 1. 创建数据库和用户
sudo -u postgres psql
CREATE DATABASE molyte_web;
CREATE USER molyte_user WITH PASSWORD 'your_password';
GRANT ALL PRIVILEGES ON DATABASE molyte_web TO molyte_user;
\q

# 2. 执行初始化脚本
psql -U molyte_user -d molyte_web -f init_db.sql
```

### 常用查询

```sql
-- 查看所有表
\dt

-- 查看表结构
\d md_jobs

-- 查询用户的所有项目
SELECT p.* FROM projects p
JOIN users u ON p.user_id = u.id
WHERE u.email = 'test@molyte.com';

-- 查询某个任务的所有 RDF 结果
SELECT * FROM rdf_results WHERE md_job_id = 1;

-- 查询正在运行的任务
SELECT * FROM md_jobs WHERE status IN ('QUEUED', 'RUNNING');
```

---

## 📝 注意事项

1. **密码哈希**：测试数据中的密码哈希对应明文密码 `password123`
2. **JSONB 类型**：使用 JSONB 存储灵活的配方和结果数据
3. **级联删除**：删除用户会级联删除其所有项目和任务
4. **自动更新时间戳**：`updated_at` 字段通过触发器自动更新


PostgreSQL 15.14 on x86_64-pc-linux-gnu
Compiled by gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-44), 64-bit

数据库: molyte_web
用户: molyte_user
密码: molyte2025
主机: localhost
端口: 5432