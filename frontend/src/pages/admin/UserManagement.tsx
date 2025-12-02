/**
 * User Management Page
 */
import React, { useState, useEffect } from 'react';
import {
  Card,
  Table,
  Button,
  Tag,
  Space,
  Modal,
  Form,
  Input,
  Select,
  InputNumber,
  message,
  Popconfirm,
  Switch,
  Row,
  Col,
  Statistic,
} from 'antd';
import {
  PlusOutlined,
  EditOutlined,
  DeleteOutlined,
  UserOutlined,
  SearchOutlined,
  FilterOutlined,
  ReloadOutlined,
  TeamOutlined,
  CrownOutlined,
  UserSwitchOutlined,
} from '@ant-design/icons';
import { useNavigate } from 'react-router-dom';
import AdminNav from '../../components/AdminNav';
import {
  getAllUsers,
  createUser,
  updateUser,
  deleteUser,
  updateUserStatus,
  UserListItem,
  UserCreate,
  UserUpdate,
  getAllPartitions,
  PartitionInfo,
} from '../../api/admin';

const UserManagement: React.FC = () => {
  const navigate = useNavigate();
  const [loading, setLoading] = useState(false);
  const [users, setUsers] = useState<UserListItem[]>([]);
  const [filteredUsers, setFilteredUsers] = useState<UserListItem[]>([]);
  const [modalVisible, setModalVisible] = useState(false);
  const [editingUser, setEditingUser] = useState<UserListItem | null>(null);
  const [partitions, setPartitions] = useState<PartitionInfo[]>([]);
  const [form] = Form.useForm();

  // 筛选和搜索状态
  const [searchText, setSearchText] = useState('');
  const [roleFilter, setRoleFilter] = useState<string | undefined>(undefined);
  const [userTypeFilter, setUserTypeFilter] = useState<string | undefined>(undefined);
  const [statusFilter, setStatusFilter] = useState<boolean | undefined>(undefined);
  const [organizationFilter, setOrganizationFilter] = useState('');

  useEffect(() => {
    loadUsers();
    loadPartitions();
  }, []);

  // 应用筛选
  useEffect(() => {
    applyFilters();
  }, [users, searchText, roleFilter, userTypeFilter, statusFilter, organizationFilter]);

  const loadUsers = async () => {
    setLoading(true);
    try {
      const data = await getAllUsers();
      setUsers(data);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载用户列表失败');
    } finally {
      setLoading(false);
    }
  };

  const applyFilters = () => {
    let filtered = [...users];

    // 搜索过滤
    if (searchText) {
      const search = searchText.toLowerCase();
      filtered = filtered.filter(
        (user) =>
          user.username.toLowerCase().includes(search) ||
          user.email.toLowerCase().includes(search) ||
          (user.organization && user.organization.toLowerCase().includes(search))
      );
    }

    // 角色过滤
    if (roleFilter) {
      filtered = filtered.filter((user) => user.role === roleFilter);
    }

    // 用户类型过滤
    if (userTypeFilter) {
      filtered = filtered.filter((user) => user.user_type === userTypeFilter);
    }

    // 状态过滤
    if (statusFilter !== undefined) {
      filtered = filtered.filter((user) => user.is_active === statusFilter);
    }

    // 组织过滤
    if (organizationFilter) {
      const org = organizationFilter.toLowerCase();
      filtered = filtered.filter(
        (user) => user.organization && user.organization.toLowerCase().includes(org)
      );
    }

    setFilteredUsers(filtered);
  };

  const handleResetFilters = () => {
    setSearchText('');
    setRoleFilter(undefined);
    setUserTypeFilter(undefined);
    setStatusFilter(undefined);
    setOrganizationFilter('');
  };

  const loadPartitions = async () => {
    try {
      const data = await getAllPartitions();
      setPartitions(data);
    } catch (error: any) {
      console.error('Failed to load partitions:', error);
      // 如果获取失败，使用默认值
      setPartitions([]);
    }
  };

  const handleCreate = () => {
    setEditingUser(null);
    form.resetFields();
    setModalVisible(true);
  };

  const handleEdit = (user: UserListItem) => {
    setEditingUser(user);
    form.setFieldsValue(user);
    setModalVisible(true);
  };

  const handleSubmit = async () => {
    try {
      const values = await form.validateFields();
      
      if (editingUser) {
        // Update user
        await updateUser(editingUser.id, values as UserUpdate);
        message.success('用户更新成功');
      } else {
        // Create user
        await createUser(values as UserCreate);
        message.success('用户创建成功');
      }
      
      setModalVisible(false);
      loadUsers();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '操作失败');
    }
  };

  const handleDelete = async (userId: number) => {
    try {
      await deleteUser(userId);
      message.success('用户删除成功');
      loadUsers();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '删除失败');
    }
  };

  const handleStatusChange = async (userId: number, isActive: boolean) => {
    try {
      await updateUserStatus(userId, isActive);
      message.success(`用户已${isActive ? '启用' : '禁用'}`);
      loadUsers();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '操作失败');
    }
  };

  const columns = [
    {
      title: 'ID',
      dataIndex: 'id',
      key: 'id',
      width: 80,
      fixed: 'left' as const,
      sorter: (a: UserListItem, b: UserListItem) => a.id - b.id,
    },
    {
      title: '用户名',
      dataIndex: 'username',
      key: 'username',
      width: 150,
      fixed: 'left' as const,
      sorter: (a: UserListItem, b: UserListItem) => a.username.localeCompare(b.username),
      render: (text: string, record: UserListItem) => (
        <a onClick={() => navigate(`/workspace/admin/users/${record.id}`)}>
          <UserOutlined /> {text}
        </a>
      ),
    },
    {
      title: '邮箱',
      dataIndex: 'email',
      key: 'email',
      width: 200,
      ellipsis: true,
      sorter: (a: UserListItem, b: UserListItem) => a.email.localeCompare(b.email),
    },
    {
      title: '角色',
      dataIndex: 'role',
      key: 'role',
      width: 100,
      filters: [
        { text: '管理员', value: 'ADMIN' },
        { text: '高级用户', value: 'PREMIUM' },
        { text: '普通用户', value: 'USER' },
        { text: '访客', value: 'GUEST' },
      ],
      onFilter: (value: any, record: UserListItem) => record.role === value,
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
      title: '状态',
      dataIndex: 'is_active',
      key: 'is_active',
      width: 80,
      filters: [
        { text: '激活', value: true },
        { text: '禁用', value: false },
      ],
      onFilter: (value: any, record: UserListItem) => record.is_active === value,
      render: (isActive: boolean, record: UserListItem) => (
        <Switch
          checked={isActive}
          onChange={(checked) => handleStatusChange(record.id, checked)}
        />
      ),
    },
    {
      title: 'CPU配额',
      dataIndex: 'total_cpu_hours',
      key: 'total_cpu_hours',
      width: 120,
      sorter: (a: UserListItem, b: UserListItem) => a.total_cpu_hours - b.total_cpu_hours,
      render: (hours: number) => `${hours.toFixed(1)} h`,
    },
    {
      title: '每日限制',
      dataIndex: 'daily_job_limit',
      key: 'daily_job_limit',
      width: 100,
      sorter: (a: UserListItem, b: UserListItem) => a.daily_job_limit - b.daily_job_limit,
    },
    {
      title: '并发限制',
      dataIndex: 'concurrent_job_limit',
      key: 'concurrent_job_limit',
      width: 100,
      sorter: (a: UserListItem, b: UserListItem) => a.concurrent_job_limit - b.concurrent_job_limit,
    },
    {
      title: '用户类型',
      dataIndex: 'user_type',
      key: 'user_type',
      width: 100,
      filters: [
        { text: '学生', value: 'STUDENT' },
        { text: '研究者', value: 'RESEARCHER' },
        { text: '企业', value: 'COMPANY' },
      ],
      onFilter: (value: any, record: UserListItem) => record.user_type === value,
      render: (type: string) => {
        const typeMap: any = {
          STUDENT: { text: '学生', color: 'cyan' },
          RESEARCHER: { text: '研究者', color: 'purple' },
          COMPANY: { text: '企业', color: 'orange' },
        };
        const config = typeMap[type] || { text: type, color: 'default' };
        return <Tag color={config.color}>{config.text}</Tag>;
      },
    },
    {
      title: '组织',
      dataIndex: 'organization',
      key: 'organization',
      width: 150,
      ellipsis: true,
    },
    {
      title: '可用队列',
      dataIndex: 'allowed_partitions',
      key: 'allowed_partitions',
      width: 200,
      render: (partitions: string[] | null) => {
        if (!partitions || partitions.length === 0) {
          return <Tag color="red">全部</Tag>;
        }
        return (
          <Space size={4} wrap>
            {partitions.map((p) => (
              <Tag key={p} color="blue" style={{ margin: 0 }}>
                {p}
              </Tag>
            ))}
          </Space>
        );
      },
    },
    {
      title: '创建时间',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 180,
      sorter: (a: UserListItem, b: UserListItem) =>
        new Date(a.created_at).getTime() - new Date(b.created_at).getTime(),
      defaultSortOrder: 'descend' as const,
      render: (time: string) => new Date(time).toLocaleString('zh-CN'),
    },
    {
      title: '操作',
      key: 'actions',
      width: 180,
      fixed: 'right' as const,
      render: (_: any, record: UserListItem) => (
        <Space size="small">
          <Button
            type="link"
            size="small"
            icon={<EditOutlined />}
            onClick={() => handleEdit(record)}
          >
            编辑
          </Button>
          <Popconfirm
            title="确定要删除这个用户吗？"
            onConfirm={() => handleDelete(record.id)}
            okText="确定"
            cancelText="取消"
          >
            <Button type="link" size="small" danger icon={<DeleteOutlined />}>
              删除
            </Button>
          </Popconfirm>
        </Space>
      ),
    },
  ];

  // 统计数据
  const totalUsers = users.length;
  const activeUsers = users.filter((u) => u.is_active).length;
  const adminUsers = users.filter((u) => u.role === 'ADMIN').length;
  const premiumUsers = users.filter((u) => u.role === 'PREMIUM').length;

  return (
    <div style={{ padding: '24px', background: '#f5f7fb', minHeight: '100vh' }}>
      <AdminNav />

      {/* 统计卡片 */}
      <Row gutter={16} style={{ marginBottom: 24 }}>
        <Col xs={24} sm={12} md={6}>
          <Card bordered={false} style={{ borderRadius: 12 }}>
            <Statistic
              title="总用户数"
              value={totalUsers}
              prefix={<TeamOutlined />}
              valueStyle={{ color: '#1677ff' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card bordered={false} style={{ borderRadius: 12 }}>
            <Statistic
              title="活跃用户"
              value={activeUsers}
              prefix={<UserSwitchOutlined />}
              valueStyle={{ color: '#52c41a' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card bordered={false} style={{ borderRadius: 12 }}>
            <Statistic
              title="管理员"
              value={adminUsers}
              prefix={<CrownOutlined />}
              valueStyle={{ color: '#ff4d4f' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card bordered={false} style={{ borderRadius: 12 }}>
            <Statistic
              title="高级用户"
              value={premiumUsers}
              prefix={<CrownOutlined />}
              valueStyle={{ color: '#faad14' }}
            />
          </Card>
        </Col>
      </Row>

      <Card
        title="用户管理"
        bordered={false}
        extra={
          <Button type="primary" icon={<PlusOutlined />} onClick={handleCreate}>
            创建用户
          </Button>
        }
        style={{
          borderRadius: '12px',
          boxShadow: '0 10px 30px rgba(15, 100, 255, 0.08)',
        }}
      >
        {/* 筛选栏 */}
        <div style={{ marginBottom: 16, padding: '16px', background: '#fafafa', borderRadius: 8 }}>
          <Row gutter={[16, 16]}>
            <Col xs={24} sm={12} md={6}>
              <Input
                placeholder="搜索用户名、邮箱、组织"
                prefix={<SearchOutlined />}
                value={searchText}
                onChange={(e) => setSearchText(e.target.value)}
                allowClear
              />
            </Col>
            <Col xs={24} sm={12} md={4}>
              <Select
                placeholder="角色"
                value={roleFilter}
                onChange={setRoleFilter}
                allowClear
                style={{ width: '100%' }}
              >
                <Select.Option value="ADMIN">管理员</Select.Option>
                <Select.Option value="PREMIUM">高级用户</Select.Option>
                <Select.Option value="USER">普通用户</Select.Option>
                <Select.Option value="GUEST">访客</Select.Option>
              </Select>
            </Col>
            <Col xs={24} sm={12} md={4}>
              <Select
                placeholder="用户类型"
                value={userTypeFilter}
                onChange={setUserTypeFilter}
                allowClear
                style={{ width: '100%' }}
              >
                <Select.Option value="STUDENT">学生</Select.Option>
                <Select.Option value="RESEARCHER">研究者</Select.Option>
                <Select.Option value="COMPANY">企业</Select.Option>
              </Select>
            </Col>
            <Col xs={24} sm={12} md={4}>
              <Select
                placeholder="状态"
                value={statusFilter}
                onChange={setStatusFilter}
                allowClear
                style={{ width: '100%' }}
              >
                <Select.Option value={true}>激活</Select.Option>
                <Select.Option value={false}>禁用</Select.Option>
              </Select>
            </Col>
            <Col xs={24} sm={12} md={6}>
              <Space>
                <Button icon={<FilterOutlined />} onClick={handleResetFilters}>
                  重置筛选
                </Button>
                <Button icon={<ReloadOutlined />} onClick={loadUsers}>
                  刷新
                </Button>
              </Space>
            </Col>
          </Row>
        </div>

        <Table
          dataSource={filteredUsers}
          columns={columns}
          rowKey="id"
          loading={loading}
          scroll={{ x: 1600 }}
          pagination={{
            pageSize: 20,
            showSizeChanger: true,
            showTotal: (total) => `共 ${total} 个用户`,
            pageSizeOptions: ['10', '20', '50', '100'],
          }}
        />
      </Card>

      {/* Create/Edit User Modal */}
      <Modal
        title={
          <Space>
            <UserOutlined style={{ color: '#1677ff' }} />
            <span style={{ fontWeight: 600 }}>{editingUser ? '编辑用户' : '创建用户'}</span>
          </Space>
        }
        open={modalVisible}
        onOk={handleSubmit}
        onCancel={() => setModalVisible(false)}
        width={600}
        centered
        okText="确定"
        cancelText="取消"
        styles={{
          body: { maxHeight: '70vh', overflowY: 'auto', padding: '24px' },
          header: { borderBottom: '1px solid #f0f0f0', paddingBottom: 16 },
        }}
      >
        <Form form={form} layout="vertical">
          <Form.Item
            name="username"
            label="用户名"
            rules={[
              { required: true, message: '请输入用户名' },
              { min: 3, message: '用户名至少 3 个字符' },
            ]}
          >
            <Input placeholder="请输入用户名" disabled={!!editingUser} />
          </Form.Item>

          <Form.Item
            name="email"
            label="邮箱"
            rules={[
              { required: true, message: '请输入邮箱' },
              { type: 'email', message: '请输入有效的邮箱地址' },
            ]}
          >
            <Input placeholder="请输入邮箱" />
          </Form.Item>

          {!editingUser && (
            <Form.Item
              name="password"
              label="密码"
              rules={[
                { required: true, message: '请输入密码' },
                { min: 6, message: '密码至少 6 个字符' },
              ]}
            >
              <Input.Password placeholder="请输入密码" />
            </Form.Item>
          )}

          <Form.Item
            name="role"
            label="角色"
            rules={[{ required: true, message: '请选择角色' }]}
            initialValue="USER"
          >
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
            rules={[{ required: true, message: '请输入 CPU 核时配额' }]}
            initialValue={100}
          >
            <InputNumber min={0} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="daily_job_limit"
            label="每日任务限制"
            rules={[{ required: true, message: '请输入每日任务限制' }]}
            initialValue={10}
          >
            <InputNumber min={0} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="concurrent_job_limit"
            label="并发任务限制"
            rules={[{ required: true, message: '请输入并发任务限制' }]}
            initialValue={3}
          >
            <InputNumber min={0} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="storage_quota_gb"
            label="存储配额 (GB)"
            rules={[{ required: true, message: '请输入存储配额' }]}
            initialValue={10}
          >
            <InputNumber min={0} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="allowed_partitions"
            label="可用队列"
            tooltip="选择用户可以使用的队列。管理员留空表示可以使用所有队列"
            initialValue={['cpu']}
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

export default UserManagement;

