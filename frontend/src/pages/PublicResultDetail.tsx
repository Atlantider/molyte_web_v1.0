/**
 * 公开结果展示页面 - 展示电解液模拟结果
 */
import { useState, useEffect } from 'react';
import { useParams, useNavigate, useLocation } from 'react-router-dom';
import {
  Card,
  Tabs,
  Button,
  Space,
  message,
  Spin,
  Typography,
  Row,
  Col,
  Tag,
  Descriptions,
  theme,
} from 'antd';
import {
  ArrowLeftOutlined,
  ExperimentOutlined,
  LineChartOutlined,
  FileTextOutlined,
  LoginOutlined,
} from '@ant-design/icons';
import type { MDJob, ElectrolyteSystem } from '../types';
import { JobStatus } from '../types';
import MoleculeViewer from '../components/MoleculeViewer';
import RDFCalculatorNature from '../components/RDFCalculatorNature';
import MSDCalculatorNature from '../components/MSDCalculatorNature';
import SolvationStructureNature from '../components/SolvationStructureNature';
import { getMDJob } from '../api/jobs';
import { getElectrolyte } from '../api/electrolytes';
import { useAuthStore } from '../stores/authStore';
import { useThemeStore } from '../stores/themeStore';
import dayjs from 'dayjs';

const { Title, Text } = Typography;

export default function PublicResultDetail() {
  const { jobId } = useParams<{ jobId: string }>();
  const navigate = useNavigate();
  const location = useLocation();
  const { isAuthenticated } = useAuthStore();
  const { mode } = useThemeStore();
  const { token } = theme.useToken();

  const [loading, setLoading] = useState(true);
  const [job, setJob] = useState<MDJob | null>(null);
  const [electrolyte, setElectrolyte] = useState<ElectrolyteSystem | null>(null);

  useEffect(() => {
    if (!isAuthenticated) {
      message.warning('请先登录后查看详情');
      navigate('/login', { state: { from: location.pathname } });
      return;
    }
    loadData();
  }, [jobId, isAuthenticated]);

  const loadData = async () => {
    if (!jobId) return;
    setLoading(true);
    try {
      const jobData = await getMDJob(parseInt(jobId));
      setJob(jobData);
      
      if (jobData.system_id) {
        try {
          const electrolyteData = await getElectrolyte(jobData.system_id);
          setElectrolyte(electrolyteData);
        } catch (e) {
          console.warn('Failed to load electrolyte data');
        }
      }
    } catch (error: any) {
      if (error.response?.status === 401) {
        message.warning('请先登录');
        navigate('/login', { state: { from: location.pathname } });
      } else {
        message.error('加载数据失败');
        navigate('/research');
      }
    } finally {
      setLoading(false);
    }
  };

  if (loading) {
    return (
      <div style={{ minHeight: '100vh', display: 'flex', alignItems: 'center', justifyContent: 'center', background: token.colorBgLayout }}>
        <Spin size="large" />
      </div>
    );
  }

  if (!job) {
    return (
      <div style={{ minHeight: '100vh', display: 'flex', alignItems: 'center', justifyContent: 'center', background: token.colorBgLayout }}>
        <Text>任务不存在</Text>
      </div>
    );
  }

  // 获取组成信息
  const cations = electrolyte?.cations || [];
  const anions = electrolyte?.anions || [];
  const solvents = electrolyte?.solvents || [];

  return (
    <div style={{ minHeight: '100vh', background: token.colorBgLayout, transition: 'background 0.3s' }}>
      {/* 顶部导航 */}
      <div style={{
        background: 'linear-gradient(135deg, #1a1f36 0%, #0d1025 100%)',
        padding: '16px 24px',
        boxShadow: '0 2px 8px rgba(0,0,0,0.15)'
      }}>
        <Row justify="space-between" align="middle">
          <Col>
            <Space>
              <Button
                icon={<ArrowLeftOutlined />}
                onClick={() => navigate('/research')}
                style={{
                  borderRadius: 8,
                  background: 'rgba(255,255,255,0.1)',
                  borderColor: 'rgba(255,255,255,0.3)',
                  color: '#fff',
                }}
              >
                返回搜索
              </Button>
              <Title level={4} style={{ margin: 0, color: '#fff' }}>
                <ExperimentOutlined /> {electrolyte?.name || `任务 #${job.id}`}
              </Title>
            </Space>
          </Col>
          <Col>
            <Button
              type="primary"
              onClick={() => navigate('/workspace/dashboard')}
              style={{ borderRadius: 8 }}
            >
              进入工作台
            </Button>
          </Col>
        </Row>
      </div>

      {/* 主内容 */}
      <div style={{ padding: 24, maxWidth: 1400, margin: '0 auto' }}>
        {/* 基本信息卡片 */}
        <Card style={{
          marginBottom: 24,
          borderRadius: 12,
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none'
        }}>
          <Descriptions title="电解液配方信息" bordered column={{ xs: 1, sm: 2, md: 3 }}>
            <Descriptions.Item label="配方名称">{electrolyte?.name || '-'}</Descriptions.Item>
            <Descriptions.Item label="温度">{electrolyte?.temperature} K</Descriptions.Item>
            <Descriptions.Item label="压力">{electrolyte?.pressure} atm</Descriptions.Item>
            <Descriptions.Item label="完成时间">
              {job.finished_at ? dayjs(job.finished_at).format('YYYY-MM-DD HH:mm') : '-'}
            </Descriptions.Item>
            <Descriptions.Item label="阳离子" span={2}>
              <Space wrap>
                {cations.map((c: any, i: number) => (
                  <Tag key={i} color="red">{c.name} ({c.number}) - {c.smiles}</Tag>
                ))}
              </Space>
            </Descriptions.Item>
            <Descriptions.Item label="阴离子" span={2}>
              <Space wrap>
                {anions.map((a: any, i: number) => (
                  <Tag key={i} color="blue">{a.name} ({a.number}) - {a.smiles}</Tag>
                ))}
              </Space>
            </Descriptions.Item>
            <Descriptions.Item label="溶剂" span={3}>
              <Space wrap>
                {solvents.map((s: any, i: number) => (
                  <Tag key={i} color="green">{s.name} ({s.number}) - {s.smiles}</Tag>
                ))}
              </Space>
            </Descriptions.Item>
          </Descriptions>
        </Card>

        {/* 结果展示 Tabs */}
        <Card style={{
          borderRadius: 12,
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none'
        }}>
          <Tabs
            defaultActiveKey="molecule_structure"
            style={{ minHeight: '500px' }}
            items={[
              {
                key: 'molecule_structure',
                label: (
                  <span style={{ fontSize: 14, fontWeight: 500 }}>
                    <ExperimentOutlined style={{ marginRight: 6 }} />
                    分子结构
                  </span>
                ),
                children: (
                  <div style={{ padding: 0 }}>
                    <MoleculeViewer jobId={job.id} />
                  </div>
                ),
              },
              {
                key: 'rdf',
                label: (
                  <span style={{ fontSize: 14, fontWeight: 500 }}>
                    <LineChartOutlined style={{ marginRight: 6 }} />
                    径向分布函数 (RDF)
                  </span>
                ),
                children: (
                  <div style={{ padding: 0 }}>
                    <RDFCalculatorNature jobId={job.id} />
                  </div>
                ),
              },
              {
                key: 'msd',
                label: (
                  <span style={{ fontSize: 14, fontWeight: 500 }}>
                    <LineChartOutlined style={{ marginRight: 6 }} />
                    均方位移 (MSD)
                  </span>
                ),
                children: (
                  <div style={{ padding: 0 }}>
                    <MSDCalculatorNature jobId={job.id} />
                  </div>
                ),
              },
              {
                key: 'solvation',
                label: (
                  <span style={{ fontSize: 14, fontWeight: 500 }}>
                    <ExperimentOutlined style={{ marginRight: 6 }} />
                    溶剂化结构
                  </span>
                ),
                children: (
                  <div style={{ padding: 0 }}>
                    <SolvationStructureNature jobId={job.id} />
                  </div>
                ),
              },
            ]}
          />
        </Card>
      </div>
    </div>
  );
}

