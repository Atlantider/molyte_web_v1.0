/**
 * 项目管理页面
 */
import { useState, useEffect } from 'react';
import { useNavigate, useLocation } from 'react-router-dom';
import { Button, Input, Space, message, Modal, Form, Row, Col, Spin, Empty, Typography, Card, theme } from 'antd';
import { PlusOutlined, SearchOutlined, ProjectOutlined, FolderOpenOutlined, FolderAddOutlined } from '@ant-design/icons';
import ProjectCard from '../components/ProjectCard';
import { getProjects, createProject, updateProject, deleteProject } from '../api/projects';
import type { Project, ProjectCreate } from '../types';
import { useThemeStore } from '../stores/themeStore';

const { TextArea } = Input;
const { Title, Text } = Typography;

export default function Projects() {
  const navigate = useNavigate();
  const location = useLocation();
  const { mode } = useThemeStore();
  const { token } = theme.useToken();
  const isDark = mode === 'dark';
  const [projects, setProjects] = useState<Project[]>([]);
  const [filteredProjects, setFilteredProjects] = useState<Project[]>([]);
  const [loading, setLoading] = useState(false);
  const [searchText, setSearchText] = useState('');
  const [modalVisible, setModalVisible] = useState(false);
  const [editingProject, setEditingProject] = useState<Project | null>(null);
  const [form] = Form.useForm();

  // 加载项目列表
  const loadProjects = async () => {
    setLoading(true);
    try {
      const data = await getProjects();
      setProjects(data);
      setFilteredProjects(data);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载项目列表失败');
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    loadProjects();
  }, []);

  // 搜索过滤
  useEffect(() => {
    if (searchText) {
      const filtered = projects.filter(
        (p) =>
          p.name.toLowerCase().includes(searchText.toLowerCase()) ||
          p.description?.toLowerCase().includes(searchText.toLowerCase())
      );
      setFilteredProjects(filtered);
    } else {
      setFilteredProjects(projects);
    }
  }, [searchText, projects]);

  // 检查是否需要自动打开创建/编辑对话框
  useEffect(() => {
    const searchParams = new URLSearchParams(location.search);
    if (location.state?.openCreateModal || searchParams.get('action') === 'create') {
      setModalVisible(true);
      // 清除 state 和 URL 参数，避免刷新时重复打开
      window.history.replaceState({}, document.title, location.pathname);
    } else if (location.state?.editProject) {
      // 处理编辑项目
      handleOpenModal(location.state.editProject);
      // 清除 state，避免刷新时重复打开
      window.history.replaceState({}, document.title, location.pathname);
    }
  }, [location]);

  // 打开创建/编辑对话框
  const handleOpenModal = (project?: Project) => {
    if (project) {
      setEditingProject(project);
      form.setFieldsValue({
        name: project.name,
        description: project.description,
      });
    } else {
      setEditingProject(null);
      form.resetFields();
    }
    setModalVisible(true);
  };

  // 关闭对话框
  const handleCloseModal = () => {
    setModalVisible(false);
    setEditingProject(null);
    form.resetFields();
  };

  // 提交表单
  const handleSubmit = async () => {
    try {
      const values = await form.validateFields();
      if (editingProject) {
        // 更新项目
        await updateProject(editingProject.id, values);
        message.success('项目更新成功');
      } else {
        // 创建项目
        await createProject(values as ProjectCreate);
        message.success('项目创建成功');
      }
      handleCloseModal();
      loadProjects();
    } catch (error: any) {
      if (error.response) {
        message.error(error.response?.data?.detail || '操作失败');
      }
    }
  };

  // 删除项目
  const handleDelete = async (id: number) => {
    try {
      await deleteProject(id);
      message.success('项目删除成功');
      loadProjects();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '删除失败');
    }
  };

  // 查看项目详情
  const handleView = (project: Project) => {
    navigate(`/workspace/projects/${project.id}`);
  };

  return (
    <div style={{
      padding: '24px',
      background: token.colorBgLayout,
      minHeight: 'calc(100vh - 64px)',
      transition: 'background 0.3s',
    }}>
      {/* 页面标题区域 */}
      <div style={{ marginBottom: 24 }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
          <div>
            <Title level={2} style={{ margin: 0, marginBottom: 8 }}>
              <ProjectOutlined style={{ marginRight: 12, color: token.colorPrimary }} />
              项目管理
            </Title>
            <Text type="secondary">管理您的研究项目，组织电解质配方和计算任务</Text>
          </div>
          <Button
            type="primary"
            icon={<PlusOutlined />}
            onClick={() => handleOpenModal()}
            size="large"
            style={{
              borderRadius: 8,
              boxShadow: isDark ? '0 2px 8px rgba(107, 154, 255, 0.3)' : '0 2px 8px rgba(91, 141, 239, 0.3)',
            }}
          >
            创建新项目
          </Button>
        </div>
      </div>

      {/* 搜索和统计卡片 */}
      <Card
        style={{
          marginBottom: 24,
          borderRadius: 12,
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none',
        }}
      >
        <Row gutter={24} align="middle">
          <Col flex="auto">
            <Input
              placeholder="搜索项目名称或描述..."
              prefix={<SearchOutlined style={{ color: '#bfbfbf' }} />}
              value={searchText}
              onChange={(e) => setSearchText(e.target.value)}
              style={{ maxWidth: 400, borderRadius: 8 }}
              allowClear
              size="large"
            />
          </Col>
          <Col>
            <Space size={24}>
              <div style={{ textAlign: 'center' }}>
                <div style={{
                  fontSize: 28,
                  fontWeight: 700,
                  color: '#1677ff',
                  lineHeight: 1.2
                }}>
                  {projects.length}
                </div>
                <Text type="secondary" style={{ fontSize: 12 }}>全部项目</Text>
              </div>
              <div style={{ textAlign: 'center' }}>
                <div style={{
                  fontSize: 28,
                  fontWeight: 700,
                  color: '#52c41a',
                  lineHeight: 1.2
                }}>
                  {filteredProjects.length}
                </div>
                <Text type="secondary" style={{ fontSize: 12 }}>当前显示</Text>
              </div>
            </Space>
          </Col>
        </Row>
      </Card>

      {/* 项目列表 */}
      <Spin spinning={loading}>
        {filteredProjects.length === 0 ? (
          <Card
            style={{
              borderRadius: 12,
              boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
              border: 'none',
            }}
          >
            <Empty
              image={<FolderOpenOutlined style={{ fontSize: 64, color: '#d9d9d9' }} />}
              description={
                <Space direction="vertical" size={8}>
                  <Text type="secondary" style={{ fontSize: 16 }}>
                    {searchText ? '没有找到匹配的项目' : '还没有项目'}
                  </Text>
                  {!searchText && (
                    <Text type="secondary">点击上方按钮创建第一个项目</Text>
                  )}
                </Space>
              }
              style={{ padding: '60px 0' }}
            >
              {!searchText && (
                <Button
                  type="primary"
                  icon={<PlusOutlined />}
                  onClick={() => handleOpenModal()}
                >
                  创建新项目
                </Button>
              )}
            </Empty>
          </Card>
        ) : (
          <Row gutter={[16, 16]}>
            {filteredProjects.map((project) => (
              <Col xs={24} sm={24} md={12} lg={8} key={project.id}>
                <ProjectCard
                  project={project}
                  onEdit={handleOpenModal}
                  onDelete={handleDelete}
                  onView={handleView}
                />
              </Col>
            ))}
          </Row>
        )}
      </Spin>

      {/* 创建/编辑对话框 */}
      <Modal
        title={
          <Space>
            {editingProject ? (
              <ProjectOutlined style={{ color: '#1677ff' }} />
            ) : (
              <FolderAddOutlined style={{ color: '#1677ff' }} />
            )}
            <span style={{ fontWeight: 600 }}>
              {editingProject ? '编辑项目' : '创建新项目'}
            </span>
          </Space>
        }
        open={modalVisible}
        onOk={handleSubmit}
        onCancel={handleCloseModal}
        okText="保存"
        cancelText="取消"
        width={600}
        centered
        styles={{
          body: { maxHeight: '70vh', overflowY: 'auto', padding: '24px' },
          header: { borderBottom: '1px solid #f0f0f0', paddingBottom: 16 },
        }}
      >
        <Form form={form} layout="vertical" style={{ marginTop: 24 }}>
          <Form.Item
            name="name"
            label="项目名称"
            rules={[
              { required: true, message: '请输入项目名称' },
              { max: 100, message: '项目名称不能超过100个字符' },
            ]}
          >
            <Input placeholder="例如: 锂电池电解液优化项目" size="large" />
          </Form.Item>
          <Form.Item
            name="description"
            label="项目描述"
            rules={[{ max: 500, message: '描述不能超过500个字符' }]}
          >
            <TextArea
              rows={4}
              placeholder="简要描述项目的研究目标和内容..."
              showCount
              maxLength={500}
            />
          </Form.Item>
        </Form>
      </Modal>
    </div>
  );
}

