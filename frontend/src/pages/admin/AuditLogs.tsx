/**
 * Audit Logs Page
 */
import React, { useState, useEffect } from 'react';
import { Card, Table, Tag, Input, Select, DatePicker, Space, Button, message, Typography, Row, Col, theme } from 'antd';
import { SearchOutlined, ReloadOutlined, FileTextOutlined } from '@ant-design/icons';
import AdminNav from '../../components/AdminNav';
import { getAuditLogs, AuditLogItem } from '../../api/admin';
import type { RangePickerProps } from 'antd/es/date-picker';
import dayjs from 'dayjs';
import { useThemeStore } from '../../stores/themeStore';

const { RangePicker } = DatePicker;
const { Title, Text } = Typography;

const AuditLogs: React.FC = () => {
  const { token } = theme.useToken();
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';
  const [loading, setLoading] = useState(false);
  const [logs, setLogs] = useState<AuditLogItem[]>([]);
  const [filters, setFilters] = useState({
    action: undefined as string | undefined,
    resource_type: undefined as string | undefined,
    start_date: undefined as string | undefined,
    end_date: undefined as string | undefined,
  });

  useEffect(() => {
    loadLogs();
  }, []);

  const loadLogs = async () => {
    setLoading(true);
    try {
      const data = await getAuditLogs({
        limit: 100,
        ...filters,
      });
      setLogs(data);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载审计日志失败');
    } finally {
      setLoading(false);
    }
  };

  const handleSearch = () => {
    loadLogs();
  };

  const handleReset = () => {
    setFilters({
      action: undefined,
      resource_type: undefined,
      start_date: undefined,
      end_date: undefined,
    });
    setTimeout(() => loadLogs(), 100);
  };

  const handleDateChange: RangePickerProps['onChange'] = (dates) => {
    if (dates && dates[0] && dates[1]) {
      setFilters({
        ...filters,
        start_date: dates[0].toISOString(),
        end_date: dates[1].toISOString(),
      });
    } else {
      setFilters({
        ...filters,
        start_date: undefined,
        end_date: undefined,
      });
    }
  };

  const actionColorMap: any = {
    create_user: 'green',
    update_user: 'blue',
    delete_user: 'red',
    update_user_quota: 'orange',
    enable_user: 'cyan',
    disable_user: 'volcano',
    cancel_job: 'magenta',
  };

  const columns = [
    {
      title: 'ID',
      dataIndex: 'id',
      key: 'id',
      width: 80,
      fixed: 'left' as const,
    },
    {
      title: '操作',
      dataIndex: 'action',
      key: 'action',
      fixed: 'left' as const,
      render: (action: string) => (
        <Tag color={actionColorMap[action] || 'default'}>{action}</Tag>
      ),
    },
    {
      title: '操作用户',
      dataIndex: 'username',
      key: 'username',
      render: (username: string | null) => username || '-',
    },
    {
      title: '资源类型',
      dataIndex: 'resource_type',
      key: 'resource_type',
      render: (type: string | null) => type || '-',
    },
    {
      title: '资源 ID',
      dataIndex: 'resource_id',
      key: 'resource_id',
      render: (id: number | null) => id || '-',
    },
    {
      title: 'IP 地址',
      dataIndex: 'ip_address',
      key: 'ip_address',
      render: (ip: string | null) => ip || '-',
    },
    {
      title: '时间',
      dataIndex: 'created_at',
      key: 'created_at',
      render: (time: string) => new Date(time).toLocaleString('zh-CN'),
    },
    {
      title: '详情',
      dataIndex: 'details',
      key: 'details',
      fixed: 'right' as const,
      render: (details: any) => (
        <div style={{ maxWidth: '300px', overflow: 'auto' }}>
          <pre style={{ margin: 0, fontSize: '12px' }}>
            {JSON.stringify(details, null, 2)}
          </pre>
        </div>
      ),
    },
  ];

  return (
    <div style={{ padding: '24px', background: token.colorBgLayout, minHeight: 'calc(100vh - 64px)' }}>
      {/* 页面标题 */}
      <div style={{ marginBottom: 24 }}>
        <Title level={2} style={{ margin: 0, marginBottom: 8 }}>
          <FileTextOutlined style={{ marginRight: 12, color: token.colorPrimary }} />
          审计日志
        </Title>
        <Text type="secondary">
          查看系统操作记录和审计信息
        </Text>
      </div>

      <AdminNav />

      {/* 筛选栏 */}
      <Card
        style={{
          marginBottom: 24,
          borderRadius: 12,
          border: 'none',
          boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 2px 8px rgba(0,0,0,0.06)',
          background: isDark ? 'rgba(255,255,255,0.04)' : '#fafafa',
        }}
        styles={{ body: { padding: '16px' } }}
      >
        <Row gutter={[16, 16]}>
          <Col xs={24} sm={12} md={6}>
            <Select
              placeholder="操作类型"
              style={{ width: '100%' }}
              allowClear
              value={filters.action}
              onChange={(value) => setFilters({ ...filters, action: value })}
            >
              <Select.Option value="create_user">创建用户</Select.Option>
              <Select.Option value="update_user">更新用户</Select.Option>
              <Select.Option value="delete_user">删除用户</Select.Option>
              <Select.Option value="update_user_quota">更新配额</Select.Option>
              <Select.Option value="enable_user">启用用户</Select.Option>
              <Select.Option value="disable_user">禁用用户</Select.Option>
              <Select.Option value="cancel_job">取消任务</Select.Option>
            </Select>
          </Col>

          <Col xs={24} sm={12} md={6}>
            <Select
              placeholder="资源类型"
              style={{ width: '100%' }}
              allowClear
              value={filters.resource_type}
              onChange={(value) => setFilters({ ...filters, resource_type: value })}
            >
              <Select.Option value="user">用户</Select.Option>
              <Select.Option value="job">任务</Select.Option>
              <Select.Option value="project">项目</Select.Option>
            </Select>
          </Col>

          <Col xs={24} sm={12} md={8}>
            <RangePicker
              onChange={handleDateChange}
              showTime
              format="YYYY-MM-DD HH:mm:ss"
              style={{ width: '100%' }}
            />
          </Col>

          <Col xs={24} sm={12} md={4}>
            <Space style={{ width: '100%' }}>
              <Button type="primary" icon={<SearchOutlined />} onClick={handleSearch}>
                搜索
              </Button>
              <Button icon={<ReloadOutlined />} onClick={handleReset}>
                重置
              </Button>
            </Space>
          </Col>
        </Row>
      </Card>

      {/* 表格 */}
      <Card
        title={
          <Space>
            <FileTextOutlined style={{ color: token.colorPrimary }} />
            <span>操作记录</span>
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
        <Table
          dataSource={logs}
          columns={columns}
          rowKey="id"
          loading={loading}
          pagination={{
            pageSize: 20,
            showSizeChanger: true,
            pageSizeOptions: ['10', '20', '50', '100'],
            showTotal: (total) => `共 ${total} 条记录`,
          }}
          scroll={{ x: 1200 }}
        />
      </Card>
    </div>
  );
};

export default AuditLogs;

