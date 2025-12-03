/**
 * 修改密码页面
 */
import { useState, useMemo } from 'react';
import { Form, Input, Button, Card, message, Typography, Progress, Alert, theme } from 'antd';
import {
  LockOutlined,
  SafetyCertificateOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  EyeInvisibleOutlined,
  EyeTwoTone,
  InfoCircleOutlined,
} from '@ant-design/icons';
import { useNavigate } from 'react-router-dom';
import * as authApi from '../api/auth';
import { useThemeStore } from '../stores/themeStore';

const { Title, Text } = Typography;

// 密码强度检查
const checkPasswordStrength = (password: string) => {
  if (!password) return { score: 0, level: '', color: '' };

  let score = 0;
  const checks = {
    length: password.length >= 8,
    lowercase: /[a-z]/.test(password),
    uppercase: /[A-Z]/.test(password),
    number: /[0-9]/.test(password),
    special: /[!@#$%^&*(),.?":{}|<>]/.test(password),
  };

  if (checks.length) score += 20;
  if (checks.lowercase) score += 20;
  if (checks.uppercase) score += 20;
  if (checks.number) score += 20;
  if (checks.special) score += 20;

  let level = '';
  let color = '';

  if (score <= 20) {
    level = '弱';
    color = '#ff4d4f';
  } else if (score <= 40) {
    level = '较弱';
    color = '#faad14';
  } else if (score <= 60) {
    level = '中等';
    color = '#fadb14';
  } else if (score <= 80) {
    level = '较强';
    color = '#a0d911';
  } else {
    level = '强';
    color = '#52c41a';
  }

  return { score, level, color, checks };
};

export default function ChangePassword() {
  const [loading, setLoading] = useState(false);
  const [newPassword, setNewPassword] = useState('');
  const [form] = Form.useForm();
  const navigate = useNavigate();
  const { mode } = useThemeStore();
  const { token } = theme.useToken();

  // 密码强度
  const passwordStrength = useMemo(() => checkPasswordStrength(newPassword), [newPassword]);

  const handleChangePassword = async (values: any) => {
    setLoading(true);
    try {
      await authApi.changePassword(values.oldPassword, values.newPassword);

      message.success('密码修改成功！请重新登录');

      // 清除 token，跳转到登录页
      authApi.logout();
      setTimeout(() => {
        navigate('/login');
      }, 1500);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '密码修改失败');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div style={{
      padding: '24px',
      background: token.colorBgLayout,
      minHeight: 'calc(100vh - 64px)',
      transition: 'background 0.3s',
    }}>
      {/* 页面标题 */}
      <div style={{ marginBottom: 24 }}>
        <Title level={2} style={{ margin: 0, marginBottom: 8 }}>
          <SafetyCertificateOutlined style={{ marginRight: 12, color: '#1677ff' }} />
          修改密码
        </Title>
        <Text type="secondary">更新您的账户密码以保护账户安全</Text>
      </div>

      <div style={{ maxWidth: 500 }}>
        <Alert
          message="密码安全提示"
          description="建议使用包含大小写字母、数字和特殊字符的强密码，长度至少8位。"
          type="info"
          showIcon
          icon={<InfoCircleOutlined />}
          style={{ marginBottom: 16, borderRadius: 8 }}
        />

        <Card
          style={{
            borderRadius: 12,
            boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
            border: 'none',
          }}
        >
          <Form
            form={form}
            onFinish={handleChangePassword}
            layout="vertical"
            autoComplete="off"
          >
            <Form.Item
              name="oldPassword"
              label="当前密码"
              rules={[
                { required: true, message: '请输入当前密码！' },
              ]}
            >
              <Input.Password
                prefix={<LockOutlined style={{ color: '#bfbfbf' }} />}
                placeholder="请输入当前密码"
                size="large"
                style={{ borderRadius: 8 }}
                iconRender={(visible) => (visible ? <EyeTwoTone /> : <EyeInvisibleOutlined />)}
              />
            </Form.Item>

            <Form.Item
              name="newPassword"
              label="新密码"
              rules={[
                { required: true, message: '请输入新密码！' },
                { min: 8, message: '密码至少8个字符！' },
                { max: 32, message: '密码最多32个字符！' },
              ]}
            >
              <Input.Password
                prefix={<LockOutlined style={{ color: '#bfbfbf' }} />}
                placeholder="请输入新密码（8-32位）"
                size="large"
                style={{ borderRadius: 8 }}
                onChange={(e) => setNewPassword(e.target.value)}
                iconRender={(visible) => (visible ? <EyeTwoTone /> : <EyeInvisibleOutlined />)}
              />
            </Form.Item>

            {/* 密码强度指示器 */}
            {newPassword && (
              <div style={{
                marginTop: -16,
                marginBottom: 16,
                padding: 12,
                background: '#f8f9fa',
                borderRadius: 8,
              }}>
                <div style={{ marginBottom: 8 }}>
                  <Progress
                    percent={passwordStrength.score}
                    showInfo={false}
                    strokeColor={passwordStrength.color}
                    size="small"
                  />
                </div>
                <div style={{ marginBottom: 8, fontSize: 13 }}>
                  <span style={{ color: passwordStrength.color, fontWeight: 500 }}>
                    密码强度：{passwordStrength.level}
                  </span>
                </div>
                <div style={{
                  display: 'grid',
                  gridTemplateColumns: '1fr 1fr',
                  gap: '4px 16px',
                  fontSize: 12,
                  color: '#8c8c8c',
                }}>
                  <div style={{
                    display: 'flex',
                    alignItems: 'center',
                    gap: 6,
                    color: passwordStrength.checks?.length ? '#52c41a' : '#8c8c8c',
                  }}>
                    {passwordStrength.checks?.length ?
                      <CheckCircleOutlined style={{ color: '#52c41a' }} /> :
                      <CloseCircleOutlined style={{ color: '#ff4d4f' }} />}
                    <span>至少8个字符</span>
                  </div>
                  <div style={{
                    display: 'flex',
                    alignItems: 'center',
                    gap: 6,
                    color: passwordStrength.checks?.lowercase ? '#52c41a' : '#8c8c8c',
                  }}>
                    {passwordStrength.checks?.lowercase ?
                      <CheckCircleOutlined style={{ color: '#52c41a' }} /> :
                      <CloseCircleOutlined style={{ color: '#ff4d4f' }} />}
                    <span>包含小写字母</span>
                  </div>
                  <div style={{
                    display: 'flex',
                    alignItems: 'center',
                    gap: 6,
                    color: passwordStrength.checks?.uppercase ? '#52c41a' : '#8c8c8c',
                  }}>
                    {passwordStrength.checks?.uppercase ?
                      <CheckCircleOutlined style={{ color: '#52c41a' }} /> :
                      <CloseCircleOutlined style={{ color: '#ff4d4f' }} />}
                    <span>包含大写字母</span>
                  </div>
                  <div style={{
                    display: 'flex',
                    alignItems: 'center',
                    gap: 6,
                    color: passwordStrength.checks?.number ? '#52c41a' : '#8c8c8c',
                  }}>
                    {passwordStrength.checks?.number ?
                      <CheckCircleOutlined style={{ color: '#52c41a' }} /> :
                      <CloseCircleOutlined style={{ color: '#ff4d4f' }} />}
                    <span>包含数字</span>
                  </div>
                  <div style={{
                    display: 'flex',
                    alignItems: 'center',
                    gap: 6,
                    color: passwordStrength.checks?.special ? '#52c41a' : '#8c8c8c',
                  }}>
                    {passwordStrength.checks?.special ?
                      <CheckCircleOutlined style={{ color: '#52c41a' }} /> :
                      <CloseCircleOutlined style={{ color: '#ff4d4f' }} />}
                    <span>包含特殊字符</span>
                  </div>
                </div>
              </div>
            )}

            <Form.Item
              name="confirmPassword"
              label="确认新密码"
              dependencies={['newPassword']}
              rules={[
                { required: true, message: '请确认新密码！' },
                ({ getFieldValue }) => ({
                  validator(_, value) {
                    if (!value || getFieldValue('newPassword') === value) {
                      return Promise.resolve();
                    }
                    return Promise.reject(new Error('两次输入的密码不一致！'));
                  },
                }),
              ]}
            >
              <Input.Password
                prefix={<LockOutlined style={{ color: '#bfbfbf' }} />}
                placeholder="请再次输入新密码"
                size="large"
                style={{ borderRadius: 8 }}
                iconRender={(visible) => (visible ? <EyeTwoTone /> : <EyeInvisibleOutlined />)}
              />
            </Form.Item>

            <Form.Item style={{ marginBottom: 0, marginTop: 24 }}>
              <Button
                type="primary"
                htmlType="submit"
                loading={loading}
                block
                size="large"
                style={{
                  borderRadius: 8,
                  height: 44,
                  boxShadow: '0 2px 8px rgba(22, 119, 255, 0.3)',
                }}
              >
                确认修改
              </Button>
              <Button
                style={{
                  marginTop: 12,
                  borderRadius: 8,
                  height: 44,
                }}
                block
                size="large"
                onClick={() => navigate('/workspace/dashboard')}
              >
                取消
              </Button>
            </Form.Item>
          </Form>
        </Card>
      </div>
    </div>
  );
}

