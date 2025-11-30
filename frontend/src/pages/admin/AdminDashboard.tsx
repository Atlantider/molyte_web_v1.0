/**
 * Admin Dashboard Page
 */
import React, { useState, useEffect } from 'react';
import { Card, Row, Col, Statistic, Table, Tag, Progress, Spin, message } from 'antd';
import {
  UserOutlined,
  RocketOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  ThunderboltOutlined,
  TeamOutlined,
} from '@ant-design/icons';
import { Column, Pie } from '@ant-design/plots';
import { useNavigate } from 'react-router-dom';
import AdminNav from '../../components/AdminNav';
import {
  getGlobalStats,
  getUserUsageStats,
  getCpuUsageRanking,
  getJobCountRanking,
  GlobalStats,
  UserUsageStatsItem,
  UserRanking,
} from '../../api/admin';

const AdminDashboard: React.FC = () => {
  const navigate = useNavigate();
  const [loading, setLoading] = useState(false);
  const [stats, setStats] = useState<GlobalStats | null>(null);
  const [userStats, setUserStats] = useState<UserUsageStatsItem[]>([]);
  const [cpuRanking, setCpuRanking] = useState<UserRanking[]>([]);
  const [jobRanking, setJobRanking] = useState<UserRanking[]>([]);

  useEffect(() => {
    loadData();
  }, []);

  const loadData = async () => {
    setLoading(true);
    try {
      const [globalStats, usageStats, cpuRank, jobRank] = await Promise.all([
        getGlobalStats(),
        getUserUsageStats({ limit: 10, sort_by: 'usage_percentage' }),
        getCpuUsageRanking(5),
        getJobCountRanking(5),
      ]);

      setStats(globalStats);
      setUserStats(usageStats);
      setCpuRanking(cpuRank);
      setJobRanking(jobRank);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载数据失败');
    } finally {
      setLoading(false);
    }
  };

  if (loading || !stats) {
    return (
      <div style={{ textAlign: 'center', padding: '100px 0' }}>
        <Spin size="large" />
      </div>
    );
  }

  // Calculate usage percentage
  const cpuUsagePercentage = stats.total_cpu_hours_allocated > 0
    ? (stats.total_cpu_hours_used / stats.total_cpu_hours_allocated) * 100
    : 0;

  const storageUsagePercentage = stats.total_storage_allocated_gb > 0
    ? (stats.total_storage_used_gb / stats.total_storage_allocated_gb) * 100
    : 0;

  // Prepare job status distribution data
  const jobStatusData = [
    { type: '已完成', value: stats.completed_jobs },
    { type: '运行中', value: stats.running_jobs },
    { type: '排队中', value: stats.queued_jobs },
    { type: '失败', value: stats.failed_jobs },
  ].filter(item => item.value > 0);

  // Prepare CPU ranking chart data
  const cpuRankingData = cpuRanking.map(item => ({
    username: item.username,
    value: item.metric_value,
  }));

  return (
    <div style={{ padding: '24px', background: '#f5f7fb', minHeight: '100vh' }}>
      {/* Header */}
      <div style={{ marginBottom: '24px' }}>
        <h1 style={{ fontSize: '28px', fontWeight: 700, margin: 0 }}>管理员仪表盘</h1>
        <p style={{ color: '#8c8c8c', marginTop: '8px' }}>系统资源监控与用户管理</p>
      </div>

      {/* Admin Navigation Menu */}
      <AdminNav />

      {/* KPI Cards */}
      <Row gutter={[16, 16]} style={{ marginBottom: '24px' }}>
        <Col xs={24} sm={12} lg={6}>
          <Card
            bordered={false}
            style={{
              borderRadius: '12px',
              boxShadow: '0 10px 30px rgba(15, 100, 255, 0.08)',
            }}
          >
            <Statistic
              title="总用户数"
              value={stats.total_users}
              prefix={<TeamOutlined style={{ color: '#1677ff' }} />}
              suffix={
                <span style={{ fontSize: '14px', color: '#8c8c8c' }}>
                  / 活跃 {stats.active_users}
                </span>
              }
            />
          </Card>
        </Col>

        <Col xs={24} sm={12} lg={6}>
          <Card
            bordered={false}
            style={{
              borderRadius: '12px',
              boxShadow: '0 10px 30px rgba(15, 100, 255, 0.08)',
            }}
          >
            <Statistic
              title="总任务数"
              value={stats.total_jobs}
              prefix={<RocketOutlined style={{ color: '#9254de' }} />}
              suffix={
                <span style={{ fontSize: '14px', color: '#8c8c8c' }}>
                  / 运行中 {stats.running_jobs}
                </span>
              }
            />
          </Card>
        </Col>

        <Col xs={24} sm={12} lg={6}>
          <Card
            bordered={false}
            style={{
              borderRadius: '12px',
              boxShadow: '0 10px 30px rgba(15, 100, 255, 0.08)',
            }}
          >
            <Statistic
              title="已完成任务"
              value={stats.completed_jobs}
              prefix={<CheckCircleOutlined style={{ color: '#52c41a' }} />}
              suffix={
                <span style={{ fontSize: '14px', color: '#8c8c8c' }}>
                  / 失败 {stats.failed_jobs}
                </span>
              }
            />
          </Card>
        </Col>

        <Col xs={24} sm={12} lg={6}>
          <Card
            bordered={false}
            style={{
              borderRadius: '12px',
              boxShadow: '0 10px 30px rgba(15, 100, 255, 0.08)',
            }}
          >
            <div>
              <div style={{ marginBottom: '8px', color: '#8c8c8c', fontSize: '14px' }}>
                CPU 核时使用率
              </div>
              <div style={{ fontSize: '24px', fontWeight: 700, marginBottom: '12px' }}>
                {cpuUsagePercentage.toFixed(1)}%
              </div>
              <Progress
                percent={cpuUsagePercentage}
                strokeColor={{
                  '0%': '#fa8c16',
                  '100%': '#f5222d',
                }}
                showInfo={false}
              />
              <div style={{ marginTop: '8px', fontSize: '12px', color: '#8c8c8c' }}>
                {stats.total_cpu_hours_used.toFixed(1)} / {stats.total_cpu_hours_allocated.toFixed(1)} 小时
              </div>
            </div>
          </Card>
        </Col>
      </Row>

      {/* Charts Row */}
      <Row gutter={[16, 16]} style={{ marginBottom: '24px' }}>
        {/* CPU Usage Ranking */}
        <Col xs={24} lg={12}>
          <Card
            title="CPU 核时使用排行 Top 5"
            bordered={false}
            style={{
              borderRadius: '12px',
              boxShadow: '0 10px 30px rgba(15, 100, 255, 0.08)',
            }}
          >
            {cpuRankingData.length > 0 ? (
              <Column
                data={cpuRankingData}
                xField="username"
                yField="value"
                label={{
                  position: 'top',
                  style: {
                    fill: '#595959',
                    fontSize: 12,
                  },
                  content: (originData: any) => {
                    const val = originData?.value;
                    return val !== undefined ? `${val.toFixed(1)}h` : '';
                  },
                }}
                color="l(270) 0:#1677ff 1:#13c2c2"
                columnStyle={{
                  radius: [8, 8, 0, 0],
                }}
                xAxis={{
                  label: {
                    autoRotate: false,
                  },
                }}
                yAxis={{
                  title: {
                    text: 'CPU 核时 (h)',
                  },
                }}
                height={300}
              />
            ) : (
              <div style={{ textAlign: 'center', padding: '50px 0', color: '#8c8c8c' }}>
                暂无数据
              </div>
            )}
          </Card>
        </Col>

        {/* Job Status Distribution */}
        <Col xs={24} lg={12}>
          <Card
            title="任务状态分布"
            bordered={false}
            style={{
              borderRadius: '12px',
              boxShadow: '0 10px 30px rgba(15, 100, 255, 0.08)',
            }}
          >
            {jobStatusData.length > 0 ? (
              <Pie
                data={jobStatusData}
                angleField="value"
                colorField="type"
                radius={0.8}
                innerRadius={0.6}
                label={{
                  type: 'inner',
                  offset: '-30%',
                  content: '{value}',
                  style: {
                    fontSize: 14,
                    textAlign: 'center',
                  },
                }}
                legend={{
                  position: 'bottom',
                }}
                statistic={{
                  title: {
                    content: '总任务',
                  },
                  content: {
                    content: stats.total_jobs.toString(),
                  },
                }}
                height={300}
              />
            ) : (
              <div style={{ textAlign: 'center', padding: '50px 0', color: '#8c8c8c' }}>
                暂无数据
              </div>
            )}
          </Card>
        </Col>
      </Row>

      {/* User Usage Stats Table */}
      <Card
        title="用户使用情况 Top 10"
        bordered={false}
        extra={
          <a onClick={() => navigate('/workspace/admin/users')}>查看全部用户 →</a>
        }
        style={{
          borderRadius: '12px',
          boxShadow: '0 10px 30px rgba(15, 100, 255, 0.08)',
        }}
      >
        <Table
          dataSource={userStats}
          rowKey="user_id"
          pagination={false}
          columns={[
            {
              title: '用户名',
              dataIndex: 'username',
              key: 'username',
              render: (text: string, record: UserUsageStatsItem) => (
                <a onClick={() => navigate(`/workspace/admin/users/${record.user_id}`)}>
                  {text}
                </a>
              ),
            },
            {
              title: '角色',
              dataIndex: 'role',
              key: 'role',
              render: (role: string) => {
                const colorMap: any = {
                  ADMIN: 'red',
                  PREMIUM: 'gold',
                  USER: 'blue',
                  GUEST: 'default',
                };
                return <Tag color={colorMap[role] || 'default'}>{role}</Tag>;
              },
            },
            {
              title: 'CPU 核时使用',
              key: 'cpu_usage',
              render: (_: any, record: UserUsageStatsItem) => (
                <div>
                  <div style={{ marginBottom: '4px' }}>
                    {record.used_cpu_hours.toFixed(1)} / {record.total_cpu_hours.toFixed(1)} h
                  </div>
                  <Progress
                    percent={record.usage_percentage}
                    size="small"
                    strokeColor={
                      record.usage_percentage > 80
                        ? '#f5222d'
                        : record.usage_percentage > 50
                        ? '#fa8c16'
                        : '#52c41a'
                    }
                  />
                </div>
              ),
            },
            {
              title: '任务统计',
              key: 'jobs',
              render: (_: any, record: UserUsageStatsItem) => (
                <div>
                  总计: {record.total_jobs} | 运行中: {record.running_jobs} |
                  完成: {record.completed_jobs} | 失败: {record.failed_jobs}
                </div>
              ),
            },
            {
              title: '最后任务时间',
              dataIndex: 'last_job_at',
              key: 'last_job_at',
              render: (time: string | null) =>
                time ? new Date(time).toLocaleString('zh-CN') : '-',
            },
          ]}
        />
      </Card>
    </div>
  );
};

export default AdminDashboard;
