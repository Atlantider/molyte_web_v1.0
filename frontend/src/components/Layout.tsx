/**
 * 主布局组件
 */
import { useState } from 'react';
import { Outlet, useNavigate, useLocation } from 'react-router-dom';
import { Layout as AntLayout, Menu, Dropdown, Avatar, Space, Button, Badge, Typography, Tooltip } from 'antd';
import {
  DashboardOutlined,
  ProjectOutlined,
  ExperimentOutlined,
  RocketOutlined,
  UserOutlined,
  LogoutOutlined,
  SettingOutlined,
  ControlOutlined,
  DatabaseOutlined,
  WalletOutlined,
  BellOutlined,
  ThunderboltOutlined,
  EyeOutlined,
} from '@ant-design/icons';
import { useAuthStore } from '../stores/authStore';
import type { MenuProps } from 'antd';

const { Header, Sider, Content } = AntLayout;
const { Text } = Typography;

export default function Layout() {
  const navigate = useNavigate();
  const location = useLocation();
  const { user, logout } = useAuthStore();
  const [collapsed, setCollapsed] = useState(false);

  // 根据当前路径确定选中的菜单项和打开的子菜单
  const getSelectedKey = () => {
    const path = location.pathname;
    // 精确匹配，优先匹配更具体的路径
    if (path.startsWith('/workspace/admin')) return '/workspace/admin';
    if (path.startsWith('/workspace/data-visibility')) return '/workspace/data-visibility';
    if (path.startsWith('/workspace/research')) return '/workspace/research';
    if (path.startsWith('/workspace/qc-jobs')) return '/workspace/qc-jobs';
    if (path.startsWith('/workspace/jobs')) return '/workspace/jobs';
    if (path.startsWith('/workspace/electrolytes')) return '/workspace/electrolytes';
    if (path.startsWith('/workspace/projects')) return '/workspace/projects';
    if (path.startsWith('/workspace/dashboard')) return '/workspace/dashboard';
    return '/workspace/dashboard';
  };

  // 获取打开的子菜单
  const getOpenKeys = () => {
    const path = location.pathname;
    if (path.startsWith('/workspace/admin')) {
      return ['/workspace/admin'];
    }
    return [];
  };

  // 侧边栏菜单项
  const menuItems: MenuProps['items'] = [
    {
      key: '/workspace/dashboard',
      icon: <DashboardOutlined />,
      label: '个人中心',
    },
    {
      key: '/workspace/projects',
      icon: <ProjectOutlined />,
      label: '项目管理',
    },
    {
      key: '/workspace/electrolytes',
      icon: <ExperimentOutlined />,
      label: '配方管理',
    },
    {
      key: '/workspace/jobs',
      icon: <RocketOutlined />,
      label: 'MD计算',
    },
    {
      key: '/workspace/qc-jobs',
      icon: <ThunderboltOutlined />,
      label: 'QC计算',
    },
    {
      key: '/workspace/research',
      icon: <DatabaseOutlined />,
      label: '数据管理',
    },
    {
      key: '/workspace/data-visibility',
      icon: <EyeOutlined />,
      label: '公开设置',
    },
    // Admin menu - only show for admin users
    ...(user?.role === 'ADMIN' ? [{
      key: '/workspace/admin',
      icon: <ControlOutlined />,
      label: '系统管理',
      children: [
        {
          key: '/workspace/admin',
          label: '管理面板',
        },
        {
          key: '/workspace/admin/users',
          label: '用户管理',
        },
        {
          key: '/workspace/admin/visibility',
          label: '数据可见性',
        },
        {
          key: '/workspace/admin/logs',
          label: '审计日志',
        },
      ],
    }] : []),
  ];

  // 用户下拉菜单
  const userMenuItems: MenuProps['items'] = [
    {
      key: 'profile',
      icon: <UserOutlined />,
      label: '个人信息',
    },
    {
      key: 'recharge',
      icon: <WalletOutlined />,
      label: '充值中心',
    },
    {
      key: 'settings',
      icon: <SettingOutlined />,
      label: '修改密码',
    },
    {
      type: 'divider',
    },
    {
      key: 'logout',
      icon: <LogoutOutlined />,
      label: '退出登录',
      danger: true,
    },
  ];

  const handleMenuClick = ({ key }: { key: string }) => {
    navigate(key);
  };

  const handleUserMenuClick = ({ key }: { key: string }) => {
    if (key === 'logout') {
      logout();
      navigate('/login');
    } else if (key === 'settings') {
      navigate('/workspace/change-password');
    } else if (key === 'profile') {
      navigate('/workspace/profile');
    } else if (key === 'recharge') {
      navigate('/workspace/recharge');
    }
  };

  return (
    <AntLayout style={{ minHeight: '100vh' }}>
      <Sider
        collapsible
        collapsed={collapsed}
        onCollapse={setCollapsed}
        theme="dark"
        style={{
          overflow: 'auto',
          height: '100vh',
          position: 'fixed',
          left: 0,
          top: 0,
          bottom: 0,
          background: 'linear-gradient(180deg, #1a1f36 0%, #0d1025 100%)',
          boxShadow: '2px 0 12px rgba(0, 0, 0, 0.15)',
        }}
      >
        {/* Logo 区域 - 点击返回首页 */}
        <Tooltip title="返回首页" placement="right">
          <div
            onClick={() => navigate('/')}
            style={{
              height: 64,
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              background: 'linear-gradient(135deg, rgba(24, 144, 255, 0.15) 0%, rgba(114, 46, 209, 0.15) 100%)',
              borderBottom: '1px solid rgba(255, 255, 255, 0.08)',
              margin: collapsed ? '12px 8px' : '12px 16px',
              borderRadius: 12,
              transition: 'all 0.3s',
              cursor: 'pointer',
            }}
          >
            <div style={{
              display: 'flex',
              alignItems: 'center',
              gap: 8,
            }}>
              <div style={{
                width: collapsed ? 32 : 36,
                height: collapsed ? 32 : 36,
                borderRadius: 10,
                background: 'linear-gradient(135deg, #1890ff 0%, #722ed1 100%)',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                boxShadow: '0 4px 12px rgba(24, 144, 255, 0.4)',
                transition: 'all 0.3s',
              }}>
                <ThunderboltOutlined style={{
                  color: '#fff',
                  fontSize: collapsed ? 18 : 20,
                  transition: 'all 0.3s',
                }} />
              </div>
              {!collapsed && (
                <span style={{
                  color: '#fff',
                  fontSize: 20,
                  fontWeight: 700,
                  letterSpacing: 2,
                  background: 'linear-gradient(135deg, #1890ff 0%, #722ed1 100%)',
                  WebkitBackgroundClip: 'text',
                  WebkitTextFillColor: 'transparent',
                }}>
                  Molyte
                </span>
              )}
            </div>
          </div>
        </Tooltip>

        <Menu
          theme="dark"
          mode="inline"
          selectedKeys={[getSelectedKey()]}
          defaultOpenKeys={getOpenKeys()}
          items={menuItems}
          onClick={handleMenuClick}
          style={{
            background: 'transparent',
            border: 'none',
            marginTop: 8,
          }}
        />

        {/* 底部装饰 */}
        {!collapsed && (
          <div style={{
            position: 'absolute',
            bottom: 60,
            left: 16,
            right: 16,
            padding: '12px 16px',
            background: 'linear-gradient(135deg, rgba(24, 144, 255, 0.1) 0%, rgba(114, 46, 209, 0.1) 100%)',
            borderRadius: 12,
            border: '1px solid rgba(255, 255, 255, 0.06)',
          }}>
            <Text style={{ color: 'rgba(255, 255, 255, 0.45)', fontSize: 11 }}>
              © 2025 Molyte Platform
            </Text>
          </div>
        )}
      </Sider>

      <AntLayout style={{ marginLeft: collapsed ? 80 : 200, transition: 'all 0.2s' }}>
        <Header
          style={{
            padding: '0 24px',
            background: '#ffffff',
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
            boxShadow: '0 1px 4px rgba(0, 21, 41, 0.08)',
            position: 'sticky',
            top: 0,
            zIndex: 100,
            height: 64,
          }}
        >
          {/* 左侧：页面标识 */}
          <div style={{ display: 'flex', alignItems: 'center' }}>
            <Text style={{
              fontSize: 15,
              color: '#8c8c8c',
              fontWeight: 400,
            }}>
              工作台
            </Text>
          </div>

          {/* 右侧：操作区 */}
          <Space size={12}>
            {/* 通知按钮 */}
            <Tooltip title="通知">
              <Badge count={0} size="small">
                <Button
                  type="text"
                  icon={<BellOutlined style={{ fontSize: 18, color: '#595959' }} />}
                  style={{
                    width: 36,
                    height: 36,
                    borderRadius: 8,
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    background: 'transparent',
                  }}
                />
              </Badge>
            </Tooltip>

            {/* 用户信息下拉 */}
            <Dropdown menu={{ items: userMenuItems, onClick: handleUserMenuClick }} placement="bottomRight">
              <div
                style={{
                  cursor: 'pointer',
                  display: 'flex',
                  alignItems: 'center',
                  gap: 10,
                  padding: '6px 12px 6px 6px',
                  borderRadius: 24,
                  background: '#f5f7fa',
                  transition: 'all 0.2s',
                }}
                onMouseEnter={(e) => {
                  e.currentTarget.style.background = '#e8ecf3';
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.background = '#f5f7fa';
                }}
              >
                <Avatar
                  size={32}
                  icon={<UserOutlined />}
                  style={{
                    background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
                  }}
                />
                <Text style={{
                  fontWeight: 500,
                  fontSize: 14,
                  color: '#262626',
                  maxWidth: 100,
                  overflow: 'hidden',
                  textOverflow: 'ellipsis',
                  whiteSpace: 'nowrap',
                }}>
                  {user?.username}
                </Text>
              </div>
            </Dropdown>
          </Space>
        </Header>

        <Content style={{
          margin: 0,
          padding: 0,
          background: '#f5f7fb',
          minHeight: 'calc(100vh - 64px)',
          overflow: 'auto',
        }}>
          <Outlet />
        </Content>
      </AntLayout>
    </AntLayout>
  );
}

