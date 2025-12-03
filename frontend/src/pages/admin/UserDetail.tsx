/**
 * User Detail Page (Admin)
 */
import React, { useState, useEffect } from 'react';
import {
  Card,
  Row,
  Col,
  Descriptions,
  Tag,
  Progress,
  Button,
  Modal,
  Form,
  InputNumber,
  Select,
  message,
  Spin,
  Table,
  theme,
} from 'antd';
import { EditOutlined, ArrowLeftOutlined } from '@ant-design/icons';
import { useParams, useNavigate } from 'react-router-dom';
import AdminNav from '../../components/AdminNav';
import {
  getUserDetail,
  updateUser,
  checkUserQuota,
  getAllJobs,
  UserDetail as UserDetailType,
  UserUpdate,
  QuotaCheckResponse,
  getAllPartitions,
  PartitionInfo,
} from '../../api/admin';
import { useThemeStore } from '../../stores/themeStore';

const UserDetail: React.FC = () => {
  const { id } = useParams<{ id: string }>();
  const navigate = useNavigate();
  const { token } = theme.useToken();
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';
  const [loading, setLoading] = useState(false);
  const [user, setUser] = useState<UserDetailType | null>(null);
  const [quota, setQuota] = useState<QuotaCheckResponse | null>(null);
  const [jobs, setJobs] = useState<any[]>([]);
  const [editModalVisible, setEditModalVisible] = useState(false);
  const [partitions, setPartitions] = useState<PartitionInfo[]>([]);
  const [form] = Form.useForm();

  useEffect(() => {
    if (id) {
      loadUserDetail();
      loadUserJobs();
    }
    loadPartitions();
  }, [id]);

  const loadUserDetail = async () => {
    setLoading(true);
    try {
      const [userData, quotaData] = await Promise.all([
        getUserDetail(Number(id)),
        checkUserQuota(Number(id)),
      ]);
      setUser(userData);
      setQuota(quotaData);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载用户详情失败');
    } finally {
      setLoading(false);
    }
  };

  const loadUserJobs = async () => {
    try {
      const jobsData = await getAllJobs({ user_id: Number(id), limit: 10 });
      setJobs(jobsData);
    } catch (error: any) {
      console.error('加载用户任务失败:', error);
    }
  };

  const loadPartitions = async () => {
    try {
      const data = await getAllPartitions();
      setPartitions(data);
    } catch (error: any) {
      console.error('Failed to load partitions:', error);
      setPartitions([]);
    }
  };

  const handleEdit = () => {
    if (user) {
      form.setFieldsValue({
        role: user.role,
        total_cpu_hours: user.total_cpu_hours,
        daily_job_limit: user.daily_job_limit,
        concurrent_job_limit: user.concurrent_job_limit,
        storage_quota_gb: user.storage_quota_gb,
      });
      setEditModalVisible(true);
    }
  };

  const handleUpdate = async () => {
    try {
      const values = await form.validateFields();
      await updateUser(Number(id), values as UserUpdate);
      message.success('用户信息更新成功');
      setEditModalVisible(false);
      loadUserDetail();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '更新失败');
    }
  };

  if (loading || !user) {
    return (
      <div style={{ textAlign: 'center', padding: '100px 0' }}>
        <Spin size="large" />
      </div>
    );
  }

  const usagePercentage = user.total_cpu_hours > 0
    ? (user.used_cpu_hours / user.total_cpu_hours) * 100
    : 0;

  const roleColorMap: any = {
    ADMIN: 'red',
    PREMIUM: 'gold',
    USER: 'blue',
    GUEST: 'default',
  };

  return (
    <div style={{ padding: '24px', background: token.colorBgLayout, minHeight: '100vh' }}>
      <AdminNav />

      {/* Header */}
      <div style={{ marginBottom: '24px' }}>
        <Button
          icon={<ArrowLeftOutlined />}
          onClick={() => navigate('/workspace/admin/users')}
          style={{ marginBottom: '16px' }}
        >
          返回用户列表
        </Button>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
          <h1 style={{ fontSize: '28px', fontWeight: 700, margin: 0 }}>
            用户详情 - {user.username}
          </h1>
          <Button type="primary" icon={<EditOutlined />} onClick={handleEdit}>
            编辑配额
          </Button>
        </div>
      </div>

      {/* Basic Info */}
      <Row gutter={[16, 16]} style={{ marginBottom: '24px' }}>
        <Col xs={24} lg={12}>
          <Card
            title="基本信息"
            bordered={false}
            style={{
              borderRadius: '12px',
              boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 10px 30px rgba(15, 100, 255, 0.08)',
              background: token.colorBgContainer,
            }}
          >
            <Descriptions column={1}>
              <Descriptions.Item label="用户 ID">{user.id}</Descriptions.Item>
              <Descriptions.Item label="用户名">{user.username}</Descriptions.Item>
              <Descriptions.Item label="邮箱">{user.email}</Descriptions.Item>
              <Descriptions.Item label="角色">
                <Tag color={roleColorMap[user.role]}>{user.role}</Tag>
              </Descriptions.Item>
              <Descriptions.Item label="状态">
                <Tag color={user.is_active ? 'success' : 'error'}>
                  {user.is_active ? '激活' : '禁用'}
                </Tag>
              </Descriptions.Item>
              <Descriptions.Item label="最后登录">
                {user.last_login_at
                  ? new Date(user.last_login_at).toLocaleString('zh-CN')
                  : '从未登录'}
              </Descriptions.Item>
              <Descriptions.Item label="创建时间">
                {new Date(user.created_at).toLocaleString('zh-CN')}
              </Descriptions.Item>
            </Descriptions>
          </Card>
        </Col>

        <Col xs={24} lg={12}>
          <Card
            title="资源配额"
            bordered={false}
            style={{
              borderRadius: '12px',
              boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 10px 30px rgba(15, 100, 255, 0.08)',
              background: token.colorBgContainer,
            }}
          >
            <Descriptions column={1}>
              <Descriptions.Item label="CPU 核时配额">
                {user.total_cpu_hours.toFixed(1)} 小时
              </Descriptions.Item>
              <Descriptions.Item label="每日任务限制">
                {user.daily_job_limit} 个
              </Descriptions.Item>
              <Descriptions.Item label="并发任务限制">
                {user.concurrent_job_limit} 个
              </Descriptions.Item>
              <Descriptions.Item label="存储配额">
                {user.storage_quota_gb.toFixed(1)} GB
              </Descriptions.Item>
              <Descriptions.Item label="可用队列">
                {!user.allowed_partitions || user.allowed_partitions.length === 0 ? (
                  <Tag color="red">全部队列 (管理员)</Tag>
                ) : (
                  <>
                    {user.allowed_partitions.map((p) => (
                      <Tag key={p} color="blue">
                        {p}
                      </Tag>
                    ))}
                  </>
                )}
              </Descriptions.Item>
            </Descriptions>
          </Card>
        </Col>
      </Row>

      {/* Usage Statistics */}
      <Row gutter={[16, 16]} style={{ marginBottom: '24px' }}>
        <Col xs={24} sm={12} lg={6}>
          <Card
            bordered={false}
            style={{
              borderRadius: '12px',
              boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 10px 30px rgba(15, 100, 255, 0.08)',
              background: token.colorBgContainer,
            }}
          >
            <div style={{ textAlign: 'center' }}>
              <div style={{ fontSize: '14px', color: token.colorTextSecondary, marginBottom: '8px' }}>
                CPU 核时使用
              </div>
              <Progress
                type="circle"
                percent={usagePercentage}
                format={() => `${usagePercentage.toFixed(1)}%`}
                strokeColor={
                  usagePercentage > 80
                    ? '#f5222d'
                    : usagePercentage > 50
                    ? '#fa8c16'
                    : '#52c41a'
                }
              />
              <div style={{ marginTop: '12px', fontSize: '12px', color: token.colorTextSecondary }}>
                {user.used_cpu_hours.toFixed(1)} / {user.total_cpu_hours.toFixed(1)} 小时
              </div>
            </div>
          </Card>
        </Col>

        <Col xs={24} sm={12} lg={6}>
          <Card
            bordered={false}
            style={{
              borderRadius: '12px',
              boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 10px 30px rgba(15, 100, 255, 0.08)',
              background: token.colorBgContainer,
            }}
          >
            <div style={{ textAlign: 'center' }}>
              <div style={{ fontSize: '14px', color: token.colorTextSecondary, marginBottom: '8px' }}>
                今日任务
              </div>
              <div style={{ fontSize: '32px', fontWeight: 700, color: token.colorPrimary }}>
                {user.today_jobs}
              </div>
              <div style={{ marginTop: '8px', fontSize: '12px', color: token.colorTextSecondary }}>
                / {user.daily_job_limit} 限制
              </div>
            </div>
          </Card>
        </Col>

        <Col xs={24} sm={12} lg={6}>
          <Card
            bordered={false}
            style={{
              borderRadius: '12px',
              boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 10px 30px rgba(15, 100, 255, 0.08)',
              background: token.colorBgContainer,
            }}
          >
            <div style={{ textAlign: 'center' }}>
              <div style={{ fontSize: '14px', color: token.colorTextSecondary, marginBottom: '8px' }}>
                运行中任务
              </div>
              <div style={{ fontSize: '32px', fontWeight: 700, color: '#9254de' }}>
                {user.running_jobs}
              </div>
              <div style={{ marginTop: '8px', fontSize: '12px', color: token.colorTextSecondary }}>
                / {user.concurrent_job_limit} 限制
              </div>
            </div>
          </Card>
        </Col>

        <Col xs={24} sm={12} lg={6}>
          <Card
            bordered={false}
            style={{
              borderRadius: '12px',
              boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 10px 30px rgba(15, 100, 255, 0.08)',
              background: token.colorBgContainer,
            }}
          >
            <div style={{ textAlign: 'center' }}>
              <div style={{ fontSize: '14px', color: token.colorTextSecondary, marginBottom: '8px' }}>
                总任务数
              </div>
              <div style={{ fontSize: '32px', fontWeight: 700, color: '#13c2c2' }}>
                {user.total_jobs}
              </div>
              <div style={{ marginTop: '8px', fontSize: '12px', color: token.colorTextSecondary }}>
                完成 {user.completed_jobs} / 失败 {user.failed_jobs}
              </div>
            </div>
          </Card>
        </Col>
      </Row>

      {/* Quota Status */}
      {quota && (
        <Card
          title="配额状态"
          bordered={false}
          style={{
            borderRadius: '12px',
            boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 10px 30px rgba(15, 100, 255, 0.08)',
            marginBottom: '24px',
            background: token.colorBgContainer,
          }}
        >
          <Tag color={quota.allowed ? 'success' : 'error'} style={{ fontSize: '16px', padding: '8px 16px' }}>
            {quota.allowed ? '✓ 可以提交新任务' : '✗ 无法提交新任务'}
          </Tag>
          {!quota.allowed && quota.reason && (
            <div style={{ marginTop: '12px', color: '#f5222d' }}>
              原因：{quota.reason}
            </div>
          )}
          {quota.details && (
            <div style={{ marginTop: '16px' }}>
              <pre style={{ background: isDark ? 'rgba(255,255,255,0.04)' : '#f5f5f5', padding: '12px', borderRadius: '4px' }}>
                {JSON.stringify(quota.details, null, 2)}
              </pre>
            </div>
          )}
        </Card>
      )}

      {/* Recent Jobs */}
      <Card
        title="最近任务"
        bordered={false}
        style={{
          borderRadius: '12px',
          boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 10px 30px rgba(15, 100, 255, 0.08)',
          background: token.colorBgContainer,
        }}
      >
        <Table
          dataSource={jobs}
          rowKey="id"
          pagination={false}
          columns={[
            {
              title: '任务 ID',
              dataIndex: 'id',
              key: 'id',
            },
            {
              title: '状态',
              dataIndex: 'status',
              key: 'status',
              render: (status: string) => {
                const colorMap: any = {
                  COMPLETED: 'success',
                  RUNNING: 'processing',
                  QUEUED: 'default',
                  FAILED: 'error',
                  CANCELLED: 'warning',
                };
                return <Tag color={colorMap[status]}>{status}</Tag>;
              },
            },
            {
              title: 'Slurm Job ID',
              dataIndex: 'slurm_job_id',
              key: 'slurm_job_id',
            },
            {
              title: '创建时间',
              dataIndex: 'created_at',
              key: 'created_at',
              render: (time: string) => new Date(time).toLocaleString('zh-CN'),
            },
          ]}
        />
      </Card>

      {/* Edit Modal */}
      <Modal
        title="编辑用户配额"
        open={editModalVisible}
        onOk={handleUpdate}
        onCancel={() => setEditModalVisible(false)}
        okText="确定"
        cancelText="取消"
      >
        <Form form={form} layout="vertical">
          <Form.Item name="role" label="角色" rules={[{ required: true }]}>
            <Select>
              <Select.Option value="ADMIN">管理员</Select.Option>
              <Select.Option value="PREMIUM">高级用户</Select.Option>
              <Select.Option value="USER">普通用户</Select.Option>
              <Select.Option value="GUEST">访客</Select.Option>
            </Select>
          </Form.Item>

          <Form.Item
            name="total_cpu_hours"
            label="CPU 核时配额 (小时)"
            rules={[{ required: true }]}
          >
            <InputNumber min={0} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="daily_job_limit"
            label="每日任务限制"
            rules={[{ required: true }]}
          >
            <InputNumber min={0} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="concurrent_job_limit"
            label="并发任务限制"
            rules={[{ required: true }]}
          >
            <InputNumber min={0} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="storage_quota_gb"
            label="存储配额 (GB)"
            rules={[{ required: true }]}
          >
            <InputNumber min={0} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="allowed_partitions"
            label="可用队列"
            tooltip="选择用户可以使用的队列。管理员留空表示可以使用所有队列"
          >
            <Select
              mode="multiple"
              placeholder="选择可用队列（留空表示全部队列）"
              allowClear
            >
              {partitions.length > 0 ? (
                partitions.map((p) => (
                  <Select.Option key={p.name} value={p.name}>
                    {p.name} ({p.state === 'up' ? '可用' : '不可用'})
                  </Select.Option>
                ))
              ) : (
                <>
                  <Select.Option value="cpu">cpu</Select.Option>
                  <Select.Option value="gpu">gpu</Select.Option>
                  <Select.Option value="debug">debug</Select.Option>
                </>
              )}
            </Select>
          </Form.Item>
        </Form>
      </Modal>
    </div>
  );
};

export default UserDetail;

