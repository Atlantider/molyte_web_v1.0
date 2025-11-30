/**
 * 个人信息页面
 */
import { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Card,
  Row,
  Col,
  Typography,
  Descriptions,
  Tag,
  Statistic,
  Space,
  Spin,
  Alert,
  Button,
  Divider,
  Table,
  Tabs,
} from 'antd';
import {
  UserOutlined,
  MailOutlined,
  ClockCircleOutlined,
  RocketOutlined,
  DatabaseOutlined,
  SettingOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  SyncOutlined,
  CrownOutlined,
  TransactionOutlined,
  WalletOutlined,
} from '@ant-design/icons';
import { getUserProfile, type UserProfile } from '../api/auth';
import { getTransactions, getOrders, type Transaction, type RechargeOrder } from '../api/billing';

const { Title, Text } = Typography;

// 角色标签
const RoleTag = ({ role }: { role: string }) => {
  const roleConfig: Record<string, { color: string; label: string; icon: React.ReactNode }> = {
    ADMIN: { color: 'gold', label: '管理员', icon: <CrownOutlined /> },
    PREMIUM: { color: 'purple', label: '高级用户', icon: <CrownOutlined /> },
    USER: { color: 'blue', label: '普通用户', icon: <UserOutlined /> },
    GUEST: { color: 'default', label: '访客', icon: <UserOutlined /> },
  };
  const config = roleConfig[role] || roleConfig.USER;
  return <Tag color={config.color} icon={config.icon}>{config.label}</Tag>;
};

export default function Profile() {
  const navigate = useNavigate();
  const [loading, setLoading] = useState(true);
  const [profile, setProfile] = useState<UserProfile | null>(null);
  const [transactions, setTransactions] = useState<Transaction[]>([]);
  const [orders, setOrders] = useState<RechargeOrder[]>([]);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    const loadData = async () => {
      try {
        setLoading(true);
        const [profileData, transData, orderData] = await Promise.all([
          getUserProfile(),
          getTransactions(0, 10).catch(() => []),
          getOrders(0, 10).catch(() => []),
        ]);
        setProfile(profileData);
        setTransactions(transData);
        setOrders(orderData);
      } catch (err: any) {
        setError(err.response?.data?.detail || '获取用户资料失败');
      } finally {
        setLoading(false);
      }
    };
    loadData();
  }, []);

  if (loading) {
    return (
      <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '60vh' }}>
        <Spin size="large" tip="加载中..." />
      </div>
    );
  }

  if (error || !profile) {
    return (
      <Alert
        message="加载失败"
        description={error || '无法获取用户资料'}
        type="error"
        showIcon
        action={<Button onClick={() => window.location.reload()}>重试</Button>}
      />
    );
  }

  const formatDate = (dateStr: string | null) => {
    if (!dateStr) return '-';
    return new Date(dateStr).toLocaleString('zh-CN');
  };

  return (
    <div style={{ padding: '24px', background: '#f5f7fb', minHeight: 'calc(100vh - 64px)' }}>
      {/* 页面标题 */}
      <div style={{ marginBottom: 24 }}>
        <Title level={2} style={{ margin: 0, marginBottom: 8 }}>
          <UserOutlined style={{ marginRight: 12, color: '#1677ff' }} />
          个人信息
        </Title>
        <Text type="secondary">查看您的账户信息和资源使用情况</Text>
      </div>

      <Row gutter={[24, 24]}>
        {/* 基本信息卡片 */}
        <Col xs={24} lg={12}>
          <Card
            title={<Space><UserOutlined />基本信息</Space>}
            style={{
              borderRadius: 12,
              boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
              border: 'none'
            }}
          >
            <Descriptions column={1} bordered size="small">
              <Descriptions.Item label="用户名">
                <Text strong>{profile.username}</Text>
              </Descriptions.Item>
              <Descriptions.Item label="邮箱">
                <Space><MailOutlined />{profile.email}</Space>
              </Descriptions.Item>
              <Descriptions.Item label="角色">
                <RoleTag role={profile.role} />
              </Descriptions.Item>
              <Descriptions.Item label="账户状态">
                {profile.is_active ? (
                  <Tag color="success" icon={<CheckCircleOutlined />}>正常</Tag>
                ) : (
                  <Tag color="error" icon={<CloseCircleOutlined />}>已禁用</Tag>
                )}
              </Descriptions.Item>
              <Descriptions.Item label="注册时间">
                <Space><ClockCircleOutlined />{formatDate(profile.created_at)}</Space>
              </Descriptions.Item>
              <Descriptions.Item label="上次登录">
                {formatDate(profile.last_login_at)}
              </Descriptions.Item>
            </Descriptions>
            <div style={{ marginTop: 16, textAlign: 'right' }}>
              <Button icon={<SettingOutlined />} onClick={() => navigate('/workspace/change-password')}>
                修改密码
              </Button>
            </div>
          </Card>
        </Col>

        {/* 任务统计卡片 */}
        <Col xs={24} lg={12}>
          <Card
            title={<Space><RocketOutlined />任务统计</Space>}
            style={{
              borderRadius: 12,
              boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
              border: 'none'
            }}
          >
            <Row gutter={16}>
              <Col span={8}>
                <Statistic title="总任务数" value={profile.total_jobs} valueStyle={{ color: '#1677ff' }} />
              </Col>
              <Col span={8}>
                <Statistic title="已完成" value={profile.completed_jobs}
                  valueStyle={{ color: '#52c41a' }} prefix={<CheckCircleOutlined />} />
              </Col>
              <Col span={8}>
                <Statistic title="运行中" value={profile.running_jobs}
                  valueStyle={{ color: '#faad14' }} prefix={<SyncOutlined spin />} />
              </Col>
            </Row>
            <Divider />
            <Row gutter={16}>
              <Col span={12}>
                <Statistic title="今日提交" value={profile.today_jobs}
                  suffix={`/ ${profile.daily_job_limit}`} />
              </Col>
              <Col span={12}>
                <Statistic title="失败任务" value={profile.failed_jobs}
                  valueStyle={{ color: '#ff4d4f' }} prefix={<CloseCircleOutlined />} />
              </Col>
            </Row>
          </Card>
        </Col>

        {/* 核时余额 */}
        <Col xs={24} lg={12}>
          <Card
            title={<Space><ClockCircleOutlined />CPU 核时余额</Space>}
            style={{
              borderRadius: 12,
              boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
              border: 'none'
            }}
          >
            <Row gutter={16} style={{ marginBottom: 16 }}>
              <Col span={8}>
                <Statistic
                  title="可用余额"
                  value={profile.available_cpu_hours || profile.remaining_cpu_hours}
                  precision={1}
                  suffix="核时"
                  valueStyle={{ color: '#52c41a', fontSize: 20 }}
                />
              </Col>
              <Col span={8}>
                <Statistic
                  title="冻结中"
                  value={profile.frozen_cpu_hours || 0}
                  precision={1}
                  suffix="核时"
                  valueStyle={{ color: '#faad14', fontSize: 20 }}
                />
              </Col>
              <Col span={8}>
                <Statistic
                  title="账户余额"
                  value={profile.balance_cpu_hours}
                  precision={1}
                  suffix="核时"
                  valueStyle={{ color: '#1677ff', fontSize: 20 }}
                />
              </Col>
            </Row>
            <Divider style={{ margin: '12px 0' }} />
            <Descriptions column={2} size="small">
              <Descriptions.Item label="初始赠送">{profile.free_cpu_hours_granted || 100} 核时</Descriptions.Item>
              <Descriptions.Item label="总充值">{profile.total_recharged || 0} 核时</Descriptions.Item>
              <Descriptions.Item label="总消耗">{profile.total_consumed || profile.used_cpu_hours} 核时</Descriptions.Item>
              <Descriptions.Item label="贡献积分">{profile.contribution_points || 0} 分</Descriptions.Item>
            </Descriptions>
            {(profile.available_cpu_hours || profile.remaining_cpu_hours) < 10 && (
              <Alert
                message="核时余额不足，请充值或使用积分兑换"
                type="warning"
                showIcon
                style={{ marginTop: 12, borderRadius: 8 }}
              />
            )}
          </Card>
        </Col>

        {/* 资源权限 */}
        <Col xs={24} lg={12}>
          <Card
            title={<Space><DatabaseOutlined />资源权限</Space>}
            style={{
              borderRadius: 12,
              boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
              border: 'none'
            }}
          >
            <Descriptions column={1} bordered size="small">
              <Descriptions.Item label="每日任务限制">
                <Space>
                  <Progress
                    type="circle"
                    percent={profile.daily_jobs_usage_percent}
                    size={50}
                    format={() => `${profile.today_jobs}/${profile.daily_job_limit}`}
                  />
                  <Text type="secondary">今日已提交 {profile.today_jobs} 个任务</Text>
                </Space>
              </Descriptions.Item>
              <Descriptions.Item label="并发任务限制">
                <Space>
                  <Text strong>{profile.concurrent_job_limit}</Text>
                  <Text type="secondary">个任务可同时运行</Text>
                </Space>
              </Descriptions.Item>
              <Descriptions.Item label="存储配额">
                <Space>
                  <Text strong>{profile.storage_quota_gb}</Text>
                  <Text type="secondary">GB</Text>
                </Space>
              </Descriptions.Item>
              <Descriptions.Item label="可用队列">
                {profile.allowed_partitions && profile.allowed_partitions.length > 0 ? (
                  <Space wrap>
                    {profile.allowed_partitions.map(p => (
                      <Tag key={p} color="blue">{p}</Tag>
                    ))}
                  </Space>
                ) : (
                  <Tag color="gold">全部队列（管理员）</Tag>
                )}
              </Descriptions.Item>
            </Descriptions>
          </Card>
        </Col>

        {/* 核时变更记录 */}
        <Col xs={24}>
          <Card
            title={<Space><TransactionOutlined />核时变更记录</Space>}
            style={{
              borderRadius: 12,
              boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
              border: 'none'
            }}
          >
            <Tabs
              items={[
                {
                  key: 'transactions',
                  label: <span><WalletOutlined /> 交易流水</span>,
                  children: (
                    <Table
                      dataSource={transactions}
                      rowKey="id"
                      size="small"
                      pagination={false}
                      locale={{ emptyText: '暂无交易记录' }}
                      columns={[
                        {
                          title: '时间',
                          dataIndex: 'created_at',
                          key: 'created_at',
                          width: 180,
                          render: (val: string) => new Date(val).toLocaleString('zh-CN'),
                        },
                        {
                          title: '类型',
                          dataIndex: 'type',
                          key: 'type',
                          width: 100,
                          render: (type: string) => {
                            const typeMap: Record<string, { color: string; text: string }> = {
                              recharge: { color: 'green', text: '充值' },
                              consume: { color: 'red', text: '消耗' },
                              points_exchange: { color: 'blue', text: '积分兑换' },
                              admin_adjust: { color: 'orange', text: '管理员调整' },
                              bonus: { color: 'purple', text: '奖励' },
                              freeze: { color: 'gold', text: '冻结' },
                              unfreeze: { color: 'cyan', text: '解冻' },
                            };
                            const config = typeMap[type] || { color: 'default', text: type };
                            return <Tag color={config.color}>{config.text}</Tag>;
                          },
                        },
                        {
                          title: '变更核时',
                          dataIndex: 'amount',
                          key: 'amount',
                          width: 120,
                          render: (val: number) => (
                            <Text style={{ color: val >= 0 ? '#52c41a' : '#ff4d4f' }}>
                              {val >= 0 ? '+' : ''}{val.toFixed(2)}
                            </Text>
                          ),
                        },
                        {
                          title: '变更后余额',
                          dataIndex: 'balance_after',
                          key: 'balance_after',
                          width: 120,
                          render: (val: number) => `${val.toFixed(2)} 核时`,
                        },
                        {
                          title: '说明',
                          dataIndex: 'description',
                          key: 'description',
                          ellipsis: true,
                        },
                      ]}
                    />
                  ),
                },
                {
                  key: 'orders',
                  label: <span><CrownOutlined /> 充值订单</span>,
                  children: (
                    <Table
                      dataSource={orders}
                      rowKey="id"
                      size="small"
                      pagination={false}
                      locale={{ emptyText: '暂无充值订单' }}
                      columns={[
                        {
                          title: '订单号',
                          dataIndex: 'order_no',
                          key: 'order_no',
                          width: 200,
                        },
                        {
                          title: '金额',
                          dataIndex: 'amount',
                          key: 'amount',
                          width: 100,
                          render: (val: number) => `¥${val.toFixed(2)}`,
                        },
                        {
                          title: '核时',
                          dataIndex: 'cpu_hours',
                          key: 'cpu_hours',
                          width: 100,
                          render: (val: number) => `${val.toFixed(2)}`,
                        },
                        {
                          title: '状态',
                          dataIndex: 'payment_status',
                          key: 'payment_status',
                          width: 100,
                          render: (status: string) => {
                            const statusMap: Record<string, { color: string; text: string }> = {
                              pending: { color: 'gold', text: '待支付' },
                              paid: { color: 'green', text: '已支付' },
                              cancelled: { color: 'default', text: '已取消' },
                              expired: { color: 'red', text: '已过期' },
                            };
                            const config = statusMap[status] || { color: 'default', text: status };
                            return <Tag color={config.color}>{config.text}</Tag>;
                          },
                        },
                        {
                          title: '创建时间',
                          dataIndex: 'created_at',
                          key: 'created_at',
                          render: (val: string) => new Date(val).toLocaleString('zh-CN'),
                        },
                      ]}
                    />
                  ),
                },
              ]}
            />
          </Card>
        </Col>
      </Row>
    </div>
  );
}

