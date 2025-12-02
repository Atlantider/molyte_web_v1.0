/**
 * 项目详情页面
 */
import { useState, useEffect } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import {
  Card,
  Descriptions,
  Button,
  Space,
  message,
  Spin,
  Empty,
  Tabs,
  Row,
  Col,
  Statistic,
  Modal,
  Form,
  Input,
} from 'antd';
import {
  ArrowLeftOutlined,
  EditOutlined,
  ExperimentOutlined,
  RocketOutlined,
} from '@ant-design/icons';
import { getProject, updateProject } from '../api/projects';
import { getElectrolytes } from '../api/electrolytes';
import { getMDJobs } from '../api/jobs';
import type { Project, ElectrolyteSystem, MDJob } from '../types';
import ElectrolyteCard from '../components/ElectrolyteCard';
import JobCard from '../components/JobCard';
import dayjs from 'dayjs';

const { TextArea } = Input;

export default function ProjectDetail() {
  const { id } = useParams<{ id: string }>();
  const navigate = useNavigate();
  const [loading, setLoading] = useState(true);
  const [project, setProject] = useState<Project | null>(null);
  const [electrolytes, setElectrolytes] = useState<ElectrolyteSystem[]>([]);
  const [jobs, setJobs] = useState<MDJob[]>([]);
  const [editModalVisible, setEditModalVisible] = useState(false);
  const [form] = Form.useForm();

  useEffect(() => {
    loadData();
  }, [id]);

  const loadData = async () => {
    if (!id) return;
    
    setLoading(true);
    try {
      // 加载项目信息
      const projectData = await getProject(parseInt(id));
      setProject(projectData);

      // 加载该项目的所有电解质配方
      const allElectrolytes = await getElectrolytes();
      const projectElectrolytes = allElectrolytes.filter(
        (e) => e.project_id === parseInt(id)
      );
      setElectrolytes(projectElectrolytes);

      // 加载该项目的所有任务
      const allJobs = await getMDJobs();
      const electrolyteIds = projectElectrolytes.map((e) => e.id);
      const projectJobs = allJobs.filter((j) =>
        electrolyteIds.includes(j.system_id)
      );
      setJobs(projectJobs);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载项目详情失败');
      navigate('/workspace/projects');
    } finally {
      setLoading(false);
    }
  };

  const handleEdit = () => {
    if (project) {
      form.setFieldsValue({
        name: project.name,
        description: project.description,
      });
      setEditModalVisible(true);
    }
  };

  const handleEditSubmit = async () => {
    if (!project) return;

    try {
      const values = await form.validateFields();
      await updateProject(project.id, values);
      message.success('项目更新成功');
      setEditModalVisible(false);
      loadData(); // 重新加载数据
    } catch (error: any) {
      if (error.response) {
        message.error(error.response?.data?.detail || '更新失败');
      }
    }
  };

  const handleDeleteElectrolyte = async (electrolyteId: number) => {
    // 这里应该调用删除 API，但为了简化，我们只是重新加载数据
    message.success('删除功能将在配方管理页面中实现');
    loadData();
  };

  const handleViewJob = (job: MDJob) => {
    message.info(`查看任务详情功能开发中: 任务 #${job.id}`);
  };

  const handleCancelJob = async (jobId: number) => {
    message.info('取消任务功能将在计算任务管理页面中实现');
  };

  if (loading) {
    return (
      <div style={{
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        height: 'calc(100vh - 64px)',
        background: '#f5f7fb',
      }}>
        <Spin size="large" />
      </div>
    );
  }

  if (!project) {
    return (
      <div style={{
        padding: 24,
        background: '#f5f7fb',
        minHeight: 'calc(100vh - 64px)'
      }}>
        <Empty description="项目不存在" />
      </div>
    );
  }

  return (
    <div style={{ padding: 24, background: '#f5f7fb', minHeight: 'calc(100vh - 64px)' }}>
      {/* 返回按钮和操作栏 */}
      <div style={{ marginBottom: 24 }}>
        <Space style={{ width: '100%', justifyContent: 'space-between' }}>
          <Button
            icon={<ArrowLeftOutlined />}
            onClick={() => navigate('/workspace/projects')}
            style={{ borderRadius: 8 }}
          >
            返回项目列表
          </Button>
          <Button
            type="primary"
            icon={<EditOutlined />}
            onClick={handleEdit}
            style={{ borderRadius: 8 }}
          >
            编辑项目
          </Button>
        </Space>
      </div>

      {/* 项目基本信息 */}
      <Card
        title="项目信息"
        style={{
          marginBottom: 24,
          borderRadius: 12,
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none'
        }}
      >
        <Descriptions column={2}>
          <Descriptions.Item label="项目名称">{project.name}</Descriptions.Item>
          <Descriptions.Item label="项目 ID">{project.id}</Descriptions.Item>
          <Descriptions.Item label="创建时间">
            {dayjs(project.created_at).format('YYYY-MM-DD HH:mm:ss')}
          </Descriptions.Item>
          <Descriptions.Item label="更新时间">
            {dayjs(project.updated_at).format('YYYY-MM-DD HH:mm:ss')}
          </Descriptions.Item>
          <Descriptions.Item label="描述" span={2}>
            {project.description || '无'}
          </Descriptions.Item>
        </Descriptions>
      </Card>

      {/* 统计信息 */}
      <Row gutter={16} style={{ marginBottom: 24 }}>
        <Col xs={24} sm={12}>
          <Card style={{
            borderRadius: 12,
            boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
            border: 'none'
          }}>
            <Statistic
              title="电解质配方"
              value={electrolytes.length}
              prefix={<ExperimentOutlined />}
              suffix="个"
              valueStyle={{ color: '#1677ff' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12}>
          <Card style={{
            borderRadius: 12,
            boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
            border: 'none'
          }}>
            <Statistic
              title="计算任务"
              value={jobs.length}
              prefix={<RocketOutlined />}
              suffix="个"
              valueStyle={{ color: '#52c41a' }}
            />
          </Card>
        </Col>
      </Row>

      {/* 标签页：电解质配方和计算任务 */}
      <Tabs
        defaultActiveKey="electrolytes"
        items={[
          {
            key: 'electrolytes',
            label: `电解质配方 (${electrolytes.length})`,
            children: (
              <div>
                {electrolytes.length === 0 ? (
                  <Empty
                    description="该项目还没有电解质配方"
                    style={{ marginTop: 40 }}
                  >
                    <Button
                      type="primary"
                      onClick={() => navigate('/workspace/electrolytes', { state: { openCreateModal: true } })}
                    >
                      创建电解质配方
                    </Button>
                  </Empty>
                ) : (
                  <div>
                    <div style={{ marginBottom: 16, display: 'flex', justifyContent: 'flex-end' }}>
                      <Button
                        type="primary"
                        icon={<ExperimentOutlined />}
                        onClick={() => navigate('/workspace/electrolytes', { state: { openCreateModal: true, projectId: parseInt(id!) } })}
                      >
                        新建配方
                      </Button>
                    </div>
                    <Row gutter={[16, 16]}>
                      {electrolytes.map((electrolyte) => (
                        <Col xs={24} sm={24} md={12} lg={8} key={electrolyte.id}>
                          <ElectrolyteCard
                            electrolyte={electrolyte}
                            jobs={jobs}
                            onEdit={() => navigate('/workspace/electrolytes')}
                            onCopy={() => navigate('/workspace/electrolytes')}
                            onDelete={handleDeleteElectrolyte}
                            onCreateJob={() => navigate('/workspace/jobs')}
                          />
                        </Col>
                      ))}
                    </Row>
                  </div>
                )}
              </div>
            ),
          },
          {
            key: 'jobs',
            label: `计算任务 (${jobs.length})`,
            children: (
              <div>
                {jobs.length === 0 ? (
                  <Empty
                    description="该项目还没有计算任务"
                    style={{ marginTop: 40 }}
                  >
                    <Button
                      type="primary"
                      onClick={() => navigate('/workspace/jobs', { state: { openCreateModal: true } })}
                      disabled={electrolytes.length === 0}
                    >
                      创建计算任务
                    </Button>
                    {electrolytes.length === 0 && (
                      <div style={{ marginTop: 12, color: '#999', fontSize: 12 }}>
                        请先创建电解质配方
                      </div>
                    )}
                  </Empty>
                ) : (
                  <div>
                    {jobs.map((job) => {
                      // 找到对应的电解质配方
                      const electrolyte = electrolytes.find(e => e.id === job.system_id);
                      return (
                        <JobCard
                          key={job.id}
                          job={job}
                          electrolyte={electrolyte}
                          onCancel={handleCancelJob}
                        />
                      );
                    })}
                  </div>
                )}
              </div>
            ),
          },
        ]}
      />

      {/* 编辑项目对话框 */}
      <Modal
        title="编辑项目"
        open={editModalVisible}
        onOk={handleEditSubmit}
        onCancel={() => setEditModalVisible(false)}
        okText="保存"
        cancelText="取消"
        width={600}
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


