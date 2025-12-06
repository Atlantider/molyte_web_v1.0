/**
 * Admin Dashboard Page
 */
import React, { useState, useEffect } from 'react';
import { Card, Row, Col, Statistic, Table, Tag, Progress, Spin, message, Typography, Space, theme, Button, Tooltip } from 'antd';
import {
  UserOutlined,
  RocketOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  ThunderboltOutlined,
  TeamOutlined,
  ControlOutlined,
  DatabaseOutlined,
  ReloadOutlined,
} from '@ant-design/icons';
import { Column, Pie } from '@ant-design/plots';
import { useNavigate, useLocation } from 'react-router-dom';
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
import { useThemeStore } from '../../stores/themeStore';

const { Title, Text } = Typography;

const AdminDashboard: React.FC = () => {
  const navigate = useNavigate();
  const location = useLocation();
  const { token } = theme.useToken();
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';
  const [loading, setLoading] = useState(false);
  const [stats, setStats] = useState<GlobalStats | null>(null);
  const [userStats, setUserStats] = useState<UserUsageStatsItem[]>([]);
  const [cpuRanking, setCpuRanking] = useState<UserRanking[]>([]);
  const [jobRanking, setJobRanking] = useState<UserRanking[]>([]);

  // 每次导航到此页面时刷新数据（使用 key 强制刷新）
  useEffect(() => {
    loadData();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [location.key, location.pathname]);

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

  // Prepare CPU ranking chart data (格式化为1位小数避免显示过长)
  const cpuRankingData = cpuRanking.map(item => ({
    username: item.username,
    value: Math.round(item.metric_value * 10) / 10,
  }));

  return (
    <div style={{ padding: '24px', background: token.colorBgLayout, minHeight: 'calc(100vh - 64px)' }}>
      {/* 页面标题 */}
      <div style={{ marginBottom: 24 }}>
        <Title level={2} style={{ margin: 0, marginBottom: 8 }}>
          <ControlOutlined style={{ marginRight: 12, color: token.colorPrimary }} />
          管理面板
        </Title>
        <Text type="secondary">
          系统资源监控与用户管理
        </Text>
      </div>

      {/* Admin Navigation Menu */}
      <AdminNav />

      {/* KPI Cards */}
      <Row gutter={[16, 16]} style={{ marginBottom: 24 }}>
        <Col xs={24} sm={12} lg={6}>
          <Card
            bordered={false}
            style={{
              borderRadius: 12,
              background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
              boxShadow: '0 4px 12px rgba(102, 126, 234, 0.3)',
              height: '100%',
            }}
          >
            <Statistic
              title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>总用户数</span>}
              value={stats.total_users}
              valueStyle={{ color: '#fff', fontSize: 28 }}
              prefix={<TeamOutlined />}
            />
          </Card>
        </Col>

        <Col xs={24} sm={12} lg={6}>
          <Card
            bordered={false}
            style={{
              borderRadius: 12,
              background: 'linear-gradient(135deg, #f093fb 0%, #f5576c 100%)',
              boxShadow: '0 4px 12px rgba(240, 147, 251, 0.3)',
              height: '100%',
            }}
          >
            <Statistic
              title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>总任务数</span>}
              value={stats.total_jobs}
              valueStyle={{ color: '#fff', fontSize: 28 }}
              prefix={<RocketOutlined />}
            />
          </Card>
        </Col>

        <Col xs={24} sm={12} lg={6}>
          <Card
            bordered={false}
            style={{
              borderRadius: 12,
              background: 'linear-gradient(135deg, #11998e 0%, #38ef7d 100%)',
              boxShadow: '0 4px 12px rgba(17, 153, 142, 0.3)',
              height: '100%',
            }}
          >
            <Statistic
              title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>已完成任务</span>}
              value={stats.completed_jobs}
              valueStyle={{ color: '#fff', fontSize: 28 }}
              prefix={<CheckCircleOutlined />}
            />
          </Card>
        </Col>

        <Col xs={24} sm={12} lg={6}>
          <Card
            bordered={false}
            style={{
              borderRadius: 12,
              background: 'linear-gradient(135deg, #fa709a 0%, #fee140 100%)',
              boxShadow: '0 4px 12px rgba(250, 112, 154, 0.3)',
              height: '100%',
            }}
          >
            <Statistic
              title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>CPU 核时使用率</span>}
              value={cpuUsagePercentage.toFixed(1)}
              valueStyle={{ color: '#fff', fontSize: 28 }}
              prefix={<ThunderboltOutlined />}
              suffix={<span style={{ fontSize: 16, color: 'rgba(255,255,255,0.7)' }}>%</span>}
            />
          </Card>
        </Col>
      </Row>

      {/* Charts Row */}
      <Row gutter={[16, 16]} style={{ marginBottom: 24 }}>
        {/* CPU Usage Ranking */}
        <Col xs={24} lg={12}>
          <Card
            title={
              <Space>
                <ThunderboltOutlined style={{ color: token.colorPrimary }} />
                <span>CPU 核时使用排行 Top 5</span>
              </Space>
            }
            bordered={false}
            style={{
              borderRadius: 12,
              border: 'none',
              boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 2px 8px rgba(0,0,0,0.06)',
              background: token.colorBgContainer,
            }}
          >
            {cpuRankingData.length > 0 ? (
              <Column
                data={cpuRankingData}
                xField="username"
                yField="value"
                theme={isDark ? 'dark' : undefined}
                label={{
                  position: 'top',
                  style: {
                    fill: token.colorText,
                    fontSize: 12,
                  },
                  content: (originData: any) => {
                    const val = originData?.value;
                    return val !== undefined ? `${val.toFixed(1)}h` : '';
                  },
                }}
                tooltip={{
                  formatter: (datum: any) => ({
                    name: 'CPU核时',
                    value: `${datum.value?.toFixed(1)} h`,
                  }),
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
              <div style={{ textAlign: 'center', padding: '50px 0', color: token.colorTextSecondary }}>
                暂无数据
              </div>
            )}
          </Card>
        </Col>

        {/* Job Status Distribution */}
        <Col xs={24} lg={12}>
          <Card
            title={
              <Space>
                <DatabaseOutlined style={{ color: token.colorPrimary }} />
                <span>任务状态分布</span>
              </Space>
            }
            bordered={false}
            style={{
              borderRadius: 12,
              border: 'none',
              boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 2px 8px rgba(0,0,0,0.06)',
              background: token.colorBgContainer,
            }}
          >
            {jobStatusData.length > 0 ? (
              <Pie
                data={jobStatusData}
                angleField="value"
                colorField="type"
                radius={0.8}
                innerRadius={0.6}
                theme={isDark ? 'dark' : undefined}
                label={{
                  type: 'inner',
                  offset: '-30%',
                  content: (data: any) => data?.value ?? '',
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
              <div style={{ textAlign: 'center', padding: '50px 0', color: token.colorTextSecondary }}>
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
          <Space>
            <Button
              icon={<ReloadOutlined />}
              onClick={loadData}
              loading={loading}
            >
              刷新
            </Button>
            <a onClick={() => navigate('/workspace/admin/users')}>查看全部用户 →</a>
          </Space>
        }
        style={{
          borderRadius: '12px',
          boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 10px 30px rgba(15, 100, 255, 0.08)',
          background: token.colorBgContainer,
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
              render: (_: any, record: UserUsageStatsItem) => {
                const percent = Math.round(record.usage_percentage * 100) / 100;
                return (
                  <div>
                    <div style={{ marginBottom: '4px' }}>
                      {record.used_cpu_hours.toFixed(1)} / {record.total_cpu_hours.toFixed(1)} h
                    </div>
                    <Progress
                      percent={percent}
                      size="small"
                      format={(p) => `${p?.toFixed(1)}%`}
                      strokeColor={
                        percent > 80
                          ? '#f5222d'
                          : percent > 50
                          ? '#fa8c16'
                          : '#52c41a'
                      }
                    />
                  </div>
                );
              },
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
