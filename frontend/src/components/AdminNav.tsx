/**
 * Admin Navigation Component
 */
import React from 'react';
import { Menu, Card } from 'antd';
import { useNavigate, useLocation } from 'react-router-dom';
import {
  DashboardOutlined,
  TeamOutlined,
  AuditOutlined,
} from '@ant-design/icons';

const AdminNav: React.FC = () => {
  const navigate = useNavigate();
  const location = useLocation();

  const menuItems = [
    {
      key: '/workspace/admin',
      icon: <DashboardOutlined />,
      label: '管理员仪表盘',
    },
    {
      key: '/workspace/admin/users',
      icon: <TeamOutlined />,
      label: '用户管理',
    },
    {
      key: '/workspace/admin/logs',
      icon: <AuditOutlined />,
      label: '审计日志',
    },
  ];

  const handleMenuClick = (e: any) => {
    navigate(e.key);
  };

  // Determine selected key based on current path
  const getSelectedKey = () => {
    const path = location.pathname;
    if (path.startsWith('/workspace/admin/users')) return '/workspace/admin/users';
    if (path.startsWith('/workspace/admin/logs')) return '/workspace/admin/logs';
    return '/workspace/admin';
  };

  return (
    <Card
      bordered={false}
      style={{
        borderRadius: '12px',
        boxShadow: '0 10px 30px rgba(15, 100, 255, 0.08)',
        marginBottom: '24px',
      }}
    >
      <Menu
        mode="horizontal"
        selectedKeys={[getSelectedKey()]}
        items={menuItems}
        onClick={handleMenuClick}
      />
    </Card>
  );
};

export default AdminNav;

