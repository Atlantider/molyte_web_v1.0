/**
 * 登录页面
 */
import { useState, useEffect, useMemo } from 'react';
import { useNavigate, useLocation } from 'react-router-dom';
import { Form, Input, Button, Card, message, Tabs, Modal, Alert, Progress, Checkbox, Tooltip, Select, Divider } from 'antd';
import {
  UserOutlined,
  LockOutlined,
  MailOutlined,
  HomeOutlined,
  ThunderboltOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  ExclamationCircleOutlined,
  EyeInvisibleOutlined,
  EyeTwoTone,
  BankOutlined,
  TeamOutlined,
  GiftOutlined,
} from '@ant-design/icons';
import { useAuthStore } from '../stores/authStore';
import { useThemeStore } from '../stores/themeStore';
import './Login.css';

// 用户类型选项
const userTypeOptions = [
  { value: 'STUDENT', label: '🎓 学生', desc: '在校学生，享受学生优惠价格' },
  { value: 'RESEARCHER', label: '🔬 研究者', desc: '高校/研究机构研究人员' },
  { value: 'COMPANY', label: '🏢 企业用户', desc: '企业研发人员' },
];

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

export default function Login() {
  const navigate = useNavigate();
  const location = useLocation();
  const { login, register, isAuthenticated, error, clearError } = useAuthStore();
  const { isDark } = useThemeStore();
  const [loading, setLoading] = useState(false);
  const [registerPassword, setRegisterPassword] = useState('');
  const [agreeTerms, setAgreeTerms] = useState(false);

  // 从 URL 参数获取初始 tab（支持 ?tab=register）
  const searchParams = new URLSearchParams(location.search);
  const initialTab = searchParams.get('tab') === 'register' ? 'register' : 'login';
  const [activeTab, setActiveTab] = useState(initialTab);

  const [resetPasswordModalVisible, setResetPasswordModalVisible] = useState(false);
  const [resetPasswordForm] = Form.useForm();

  // 密码强度
  const passwordStrength = useMemo(() => checkPasswordStrength(registerPassword), [registerPassword]);

  // 获取登录后需要跳转的页面
  const from = (location.state as any)?.from || '/workspace/dashboard';

  useEffect(() => {
    if (isAuthenticated) {
      navigate(from);
    }
  }, [isAuthenticated, navigate, from]);

  useEffect(() => {
    if (error) {
      message.error(error);
      clearError();
    }
  }, [error, clearError]);

  const handleLogin = async (values: any) => {
    setLoading(true);
    try {
      await login(values.username, values.password);
      message.success('登录成功！');
      navigate(from);
    } catch (err) {
      // Error handled by store
    } finally {
      setLoading(false);
    }
  };

  const handleRegister = async (values: any) => {
    setLoading(true);
    try {
      await register(
        values.email,
        values.username,
        values.password,
        values.user_type,
        values.organization,
        values.department
      );
      message.success('🎉 注册成功！已赠送 100 核时免费计算配额');
      navigate(from);
    } catch (err) {
      // Error handled by store
    } finally {
      setLoading(false);
    }
  };

  const handleResetPassword = async (values: any) => {
    setLoading(true);
    try {
      // TODO: 实现重置密码 API 调用
      message.info('重置密码功能开发中，请联系管理员重置密码');
      setResetPasswordModalVisible(false);
      resetPasswordForm.resetFields();
    } catch (err) {
      message.error('重置密码失败');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className={`login-container ${isDark ? 'dark' : ''}`}>
      {/* 返回首页按钮 */}
      <Button
        type="text"
        icon={<HomeOutlined />}
        onClick={() => navigate('/')}
        className="back-home-btn"
      >
        返回首页
      </Button>

      <Card className="login-card" title={
        <div className="login-card-header">
          <div className="login-logo">
            <div className="login-logo-icon">
              <ThunderboltOutlined />
            </div>
            <span className="login-logo-text">Molyte</span>
          </div>
          <p className="login-subtitle">电解液研发模拟平台</p>
        </div>
      }>
        
        <Tabs
          activeKey={activeTab}
          onChange={setActiveTab}
          centered
          items={[
            {
              key: 'login',
              label: '登录',
              children: (
                <>
                  <Alert
                    message="提示：可以使用邮箱或用户名登录"
                    type="info"
                    showIcon
                    style={{ marginBottom: 16 }}
                  />
                  <Form
                    name="login"
                    onFinish={handleLogin}
                    autoComplete="off"
                    size="large"
                  >
                    <Form.Item
                      name="username"
                      rules={[{ required: true, message: '请输入用户名或邮箱！' }]}
                    >
                      <Input
                        prefix={<UserOutlined />}
                        placeholder="用户名或邮箱"
                      />
                    </Form.Item>

                    <Form.Item
                      name="password"
                      rules={[{ required: true, message: '请输入密码！' }]}
                    >
                      <Input.Password
                        prefix={<LockOutlined />}
                        placeholder="密码"
                      />
                    </Form.Item>

                    <Form.Item>
                      <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: 16 }}>
                        <a onClick={() => setResetPasswordModalVisible(true)}>忘记密码？</a>
                      </div>
                      <Button type="primary" htmlType="submit" loading={loading} block>
                        登录
                      </Button>
                    </Form.Item>
                  </Form>
                </>
              ),
            },
            {
              key: 'register',
              label: '注册',
              children: (
                <Form
                  name="register"
                  onFinish={handleRegister}
                  autoComplete="off"
                  size="large"
                  initialValues={{ user_type: 'STUDENT' }}
                >
                  {/* 注册福利提示 */}
                  <Alert
                    message={
                      <span>
                        <GiftOutlined style={{ marginRight: 8 }} />
                        注册即送 <strong>100 核时</strong> 免费计算配额
                      </span>
                    }
                    description="使用免费核时计算的数据将在1年后自动公开，促进科研共享"
                    type="success"
                    showIcon={false}
                    style={{ marginBottom: 16, borderRadius: 8 }}
                  />

                  <Divider plain style={{ margin: '8px 0 16px' }}>基本信息</Divider>

                  <Form.Item
                    name="email"
                    rules={[
                      { required: true, message: '请输入单位邮箱！' },
                      { type: 'email', message: '请输入有效的邮箱地址！' },
                    ]}
                  >
                    <Input
                      prefix={<MailOutlined />}
                      placeholder="请输入单位邮箱（用于登录和找回密码）"
                    />
                  </Form.Item>

                  <Form.Item
                    name="username"
                    rules={[
                      { required: true, message: '请输入用户名！' },
                      { min: 3, message: '用户名至少3个字符！' },
                      { max: 20, message: '用户名最多20个字符！' },
                      { pattern: /^[a-zA-Z0-9_]+$/, message: '用户名只能包含字母、数字和下划线！' },
                    ]}
                  >
                    <Input
                      prefix={<UserOutlined />}
                      placeholder="请输入用户名（3-20位，字母数字下划线）"
                    />
                  </Form.Item>

                  <Divider plain style={{ margin: '8px 0 16px' }}>单位信息</Divider>

                  <Form.Item
                    name="user_type"
                    rules={[{ required: true, message: '请选择您的身份！' }]}
                  >
                    <Select
                      placeholder="请选择您的身份"
                      options={userTypeOptions.map(opt => ({
                        value: opt.value,
                        label: (
                          <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
                            <span>{opt.label}</span>
                            <span style={{ fontSize: 12, color: '#999' }}>{opt.desc}</span>
                          </div>
                        ),
                      }))}
                      suffixIcon={<TeamOutlined />}
                    />
                  </Form.Item>

                  <Form.Item
                    name="organization"
                    rules={[
                      { required: true, message: '请输入单位名称！' },
                      { min: 2, message: '单位名称至少2个字符！' },
                    ]}
                  >
                    <Input
                      prefix={<BankOutlined />}
                      placeholder="请输入单位名称（如：清华大学）"
                    />
                  </Form.Item>

                  <Form.Item
                    name="department"
                  >
                    <Input
                      prefix={<TeamOutlined />}
                      placeholder="部门（选填，如：化学系）"
                    />
                  </Form.Item>

                  <Divider plain style={{ margin: '8px 0 16px' }}>账户安全</Divider>

                  <Form.Item
                    name="password"
                    rules={[
                      { required: true, message: '请输入密码！' },
                      { min: 8, message: '密码至少8个字符！' },
                      { max: 32, message: '密码最多32个字符！' },
                    ]}
                  >
                    <Input.Password
                      prefix={<LockOutlined />}
                      placeholder="请输入密码（8-32位）"
                      onChange={(e) => setRegisterPassword(e.target.value)}
                      iconRender={(visible) => (visible ? <EyeTwoTone /> : <EyeInvisibleOutlined />)}
                    />
                  </Form.Item>

                  {/* 密码强度指示器 */}
                  {registerPassword && (
                    <div className="password-strength">
                      <div className="password-strength-bar">
                        <Progress
                          percent={passwordStrength.score}
                          showInfo={false}
                          strokeColor={passwordStrength.color}
                          size="small"
                        />
                      </div>
                      <div className="password-strength-info">
                        <span style={{ color: passwordStrength.color, fontWeight: 500 }}>
                          密码强度：{passwordStrength.level}
                        </span>
                      </div>
                      <div className="password-requirements">
                        <div className={`requirement ${passwordStrength.checks?.length ? 'met' : ''}`}>
                          {passwordStrength.checks?.length ? <CheckCircleOutlined /> : <CloseCircleOutlined />}
                          <span>至少8个字符</span>
                        </div>
                        <div className={`requirement ${passwordStrength.checks?.lowercase ? 'met' : ''}`}>
                          {passwordStrength.checks?.lowercase ? <CheckCircleOutlined /> : <CloseCircleOutlined />}
                          <span>包含小写字母</span>
                        </div>
                        <div className={`requirement ${passwordStrength.checks?.uppercase ? 'met' : ''}`}>
                          {passwordStrength.checks?.uppercase ? <CheckCircleOutlined /> : <CloseCircleOutlined />}
                          <span>包含大写字母</span>
                        </div>
                        <div className={`requirement ${passwordStrength.checks?.number ? 'met' : ''}`}>
                          {passwordStrength.checks?.number ? <CheckCircleOutlined /> : <CloseCircleOutlined />}
                          <span>包含数字</span>
                        </div>
                        <div className={`requirement ${passwordStrength.checks?.special ? 'met' : ''}`}>
                          {passwordStrength.checks?.special ? <CheckCircleOutlined /> : <CloseCircleOutlined />}
                          <span>包含特殊字符</span>
                        </div>
                      </div>
                    </div>
                  )}

                  <Form.Item
                    name="confirmPassword"
                    dependencies={['password']}
                    rules={[
                      { required: true, message: '请再次输入密码！' },
                      ({ getFieldValue }) => ({
                        validator(_, value) {
                          if (!value || getFieldValue('password') === value) {
                            return Promise.resolve();
                          }
                          return Promise.reject(new Error('两次输入的密码不一致！'));
                        },
                      }),
                    ]}
                  >
                    <Input.Password
                      prefix={<LockOutlined />}
                      placeholder="请再次输入密码确认"
                      iconRender={(visible) => (visible ? <EyeTwoTone /> : <EyeInvisibleOutlined />)}
                    />
                  </Form.Item>

                  <Form.Item
                    name="agreement"
                    valuePropName="checked"
                    rules={[
                      {
                        validator: (_, value) =>
                          value ? Promise.resolve() : Promise.reject(new Error('请阅读并同意用户协议和隐私政策')),
                      },
                    ]}
                  >
                    <Checkbox onChange={(e) => setAgreeTerms(e.target.checked)}>
                      我已阅读并同意{' '}
                      <a onClick={(e) => { e.preventDefault(); message.info('用户协议页面开发中'); }}>
                        《用户协议》
                      </a>
                      {' '}和{' '}
                      <a onClick={(e) => { e.preventDefault(); message.info('隐私政策页面开发中'); }}>
                        《隐私政策》
                      </a>
                    </Checkbox>
                  </Form.Item>

                  <Form.Item>
                    <Button
                      type="primary"
                      htmlType="submit"
                      loading={loading}
                      block
                    >
                      注册
                    </Button>
                  </Form.Item>
                </Form>
              ),
            },
          ]}
        />
      </Card>

      {/* 忘记密码模态框 */}
      <Modal
        title="重置密码"
        open={resetPasswordModalVisible}
        onCancel={() => {
          setResetPasswordModalVisible(false);
          resetPasswordForm.resetFields();
        }}
        footer={null}
      >
        <Alert
          message="请输入您的注册邮箱，我们将发送重置密码链接到您的邮箱"
          type="info"
          showIcon
          style={{ marginBottom: 16 }}
        />
        <Form
          form={resetPasswordForm}
          onFinish={handleResetPassword}
          layout="vertical"
        >
          <Form.Item
            name="email"
            label="邮箱地址"
            rules={[
              { required: true, message: '请输入邮箱！' },
              { type: 'email', message: '请输入有效的邮箱地址！' },
            ]}
          >
            <Input
              prefix={<MailOutlined />}
              placeholder="请输入注册时使用的邮箱"
            />
          </Form.Item>

          <Form.Item>
            <Button type="primary" htmlType="submit" loading={loading} block>
              发送重置链接
            </Button>
          </Form.Item>
        </Form>
      </Modal>
    </div>
  );
}

