/**
 * Audit Logs Page
 */
import React, { useState, useEffect } from 'react';
import { Card, Table, Tag, Input, Select, DatePicker, Space, Button, message } from 'antd';
import { SearchOutlined, ReloadOutlined } from '@ant-design/icons';
import AdminNav from '../../components/AdminNav';
import { getAuditLogs, AuditLogItem } from '../../api/admin';
import type { RangePickerProps } from 'antd/es/date-picker';
import dayjs from 'dayjs';

const { RangePicker } = DatePicker;

const AuditLogs: React.FC = () => {
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
    },
    {
      title: '操作',
      dataIndex: 'action',
      key: 'action',
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
    <div style={{ padding: '24px', background: '#f5f7fb', minHeight: '100vh' }}>
      <AdminNav />

      <Card
        title="审计日志"
        bordered={false}
        style={{
          borderRadius: '12px',
          boxShadow: '0 10px 30px rgba(15, 100, 255, 0.08)',
        }}
      >
        {/* Filters */}
        <Space style={{ marginBottom: '16px' }} wrap>
          <Select
            placeholder="操作类型"
            style={{ width: 200 }}
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

          <Select
            placeholder="资源类型"
            style={{ width: 150 }}
            allowClear
            value={filters.resource_type}
            onChange={(value) => setFilters({ ...filters, resource_type: value })}
          >
            <Select.Option value="user">用户</Select.Option>
            <Select.Option value="job">任务</Select.Option>
            <Select.Option value="project">项目</Select.Option>
          </Select>

          <RangePicker
            onChange={handleDateChange}
            showTime
            format="YYYY-MM-DD HH:mm:ss"
          />

          <Button type="primary" icon={<SearchOutlined />} onClick={handleSearch}>
            搜索
          </Button>

          <Button icon={<ReloadOutlined />} onClick={handleReset}>
            重置
          </Button>
        </Space>

        <Table
          dataSource={logs}
          columns={columns}
          rowKey="id"
          loading={loading}
          pagination={{
            pageSize: 20,
            showSizeChanger: true,
            showTotal: (total) => `共 ${total} 条记录`,
          }}
          scroll={{ x: 1200 }}
        />
      </Card>
    </div>
  );
};

export default AuditLogs;

