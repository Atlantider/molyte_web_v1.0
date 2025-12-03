/**
 * QC计算结果展示组件
 */
import React, { useState, useEffect } from 'react';
import {
  Card,
  Descriptions,
  Statistic,
  Row,
  Col,
  Image,
  Typography,
  Divider,
  Tag,
  Tooltip,
  Space,
  Button,
  message,
  Spin,
  Tabs,
  theme,
} from 'antd';
import {
  ThunderboltOutlined,
  ExperimentOutlined,
  FireOutlined,
  DownloadOutlined,
  BgColorsOutlined,
  ApartmentOutlined,
  InfoCircleOutlined,
} from '@ant-design/icons';
import {
  getESPImage, downloadESPImage,
  getHOMOImage, downloadHOMOImage,
  getLUMOImage, downloadLUMOImage
} from '../api/qc';
import type { QCResult, QCJob } from '../types/qc';
import { useThemeStore } from '../stores/themeStore';

const { Text, Title } = Typography;

// Hartree to eV conversion
const HARTREE_TO_EV = 27.2114;

interface QCResultsPanelProps {
  results: QCResult[];
  compact?: boolean;  // 紧凑模式，用于分子查看器侧边栏
  job?: QCJob;  // QC任务信息，用于显示计算参数
}

export default function QCResultsPanel({ results, compact = false, job }: QCResultsPanelProps) {
  const { mode } = useThemeStore();
  const { token } = theme.useToken();
  const [espImageUrl, setEspImageUrl] = useState<string | null>(null);
  const [homoImageUrl, setHomoImageUrl] = useState<string | null>(null);
  const [lumoImageUrl, setLumoImageUrl] = useState<string | null>(null);
  const [espLoading, setEspLoading] = useState(false);
  const [homoLoading, setHomoLoading] = useState(false);
  const [lumoLoading, setLumoLoading] = useState(false);
  const [downloading, setDownloading] = useState<string | null>(null);

  // 取第一个结果（通常只有一个）
  const result = results && results.length > 0 ? results[0] : null;

  // 加载所有图片
  useEffect(() => {
    if (!result?.id) return;

    // 加载ESP图片
    if (result?.esp_image_path) {
      setEspLoading(true);
      getESPImage(result.id)
        .then((blob) => {
          const url = URL.createObjectURL(blob);
          setEspImageUrl(url);
        })
        .catch((err) => {
          console.error('加载ESP图片失败:', err);
        })
        .finally(() => {
          setEspLoading(false);
        });
    }

    // 加载HOMO图片
    if (result?.homo_image_path) {
      setHomoLoading(true);
      getHOMOImage(result.id)
        .then((blob) => {
          const url = URL.createObjectURL(blob);
          setHomoImageUrl(url);
        })
        .catch((err) => {
          console.error('加载HOMO图片失败:', err);
        })
        .finally(() => {
          setHomoLoading(false);
        });
    }

    // 加载LUMO图片
    if (result?.lumo_image_path) {
      setLumoLoading(true);
      getLUMOImage(result.id)
        .then((blob) => {
          const url = URL.createObjectURL(blob);
          setLumoImageUrl(url);
        })
        .catch((err) => {
          console.error('加载LUMO图片失败:', err);
        })
        .finally(() => {
          setLumoLoading(false);
        });
    }

    // 清理URL
    return () => {
      if (espImageUrl) URL.revokeObjectURL(espImageUrl);
      if (homoImageUrl) URL.revokeObjectURL(homoImageUrl);
      if (lumoImageUrl) URL.revokeObjectURL(lumoImageUrl);
    };
  }, [result?.id, result?.esp_image_path, result?.homo_image_path, result?.lumo_image_path]);

  // 下载图片
  const handleDownload = async (type: 'ESP' | 'HOMO' | 'LUMO') => {
    if (!result?.id) return;
    setDownloading(type);
    try {
      if (type === 'ESP') {
        await downloadESPImage(result.id, `ESP_${result.id}.png`);
      } else if (type === 'HOMO') {
        await downloadHOMOImage(result.id, `HOMO_${result.id}.png`);
      } else {
        await downloadLUMOImage(result.id, `LUMO_${result.id}.png`);
      }
      message.success(`${type}图片下载成功`);
    } catch (error) {
      message.error('下载失败');
    } finally {
      setDownloading(null);
    }
  };

  if (!results || results.length === 0 || !result) {
    return (
      <Card title="QC计算结果">
        <div style={{ textAlign: 'center', padding: 20, color: '#999' }}>
          暂无QC计算结果
        </div>
      </Card>
    );
  }

  // 计算eV值
  const homoEv = result.homo ? result.homo * HARTREE_TO_EV : null;
  const lumoEv = result.lumo ? result.lumo * HARTREE_TO_EV : null;

  if (compact) {
    // 紧凑模式
    return (
      <div>
        <Title level={5}>
          <ExperimentOutlined /> QC性质
        </Title>
        <Descriptions column={1} size="small" bordered>
          {result.energy_au && (
            <Descriptions.Item label="能量">
              {result.energy_au.toFixed(6)} A.U.
            </Descriptions.Item>
          )}
          {homoEv && (
            <Descriptions.Item label="HOMO">
              <Tag color="blue">{homoEv.toFixed(3)} eV</Tag>
            </Descriptions.Item>
          )}
          {lumoEv && (
            <Descriptions.Item label="LUMO">
              <Tag color="green">{lumoEv.toFixed(3)} eV</Tag>
            </Descriptions.Item>
          )}
          {result.homo_lumo_gap && (
            <Descriptions.Item label="HOMO-LUMO Gap">
              <Tag color="orange">{result.homo_lumo_gap.toFixed(3)} eV</Tag>
            </Descriptions.Item>
          )}
          {result.esp_min_kcal && (
            <Descriptions.Item label="ESP Min">
              {result.esp_min_kcal.toFixed(2)} kcal/mol
            </Descriptions.Item>
          )}
          {result.esp_max_kcal && (
            <Descriptions.Item label="ESP Max">
              {result.esp_max_kcal.toFixed(2)} kcal/mol
            </Descriptions.Item>
          )}
        </Descriptions>

        {result.esp_image_path && (
          <div style={{ marginTop: 12 }}>
            <Space style={{ marginBottom: 8 }}>
              <Text type="secondary">ESP表面</Text>
              <Button
                size="small"
                icon={<DownloadOutlined />}
                onClick={() => handleDownload('ESP')}
                loading={downloading === 'ESP'}
              >
                下载
              </Button>
            </Space>
            {espLoading ? (
              <div style={{ textAlign: 'center', padding: 20 }}>
                <Spin tip="加载中..." />
              </div>
            ) : espImageUrl ? (
              <Image
                src={espImageUrl}
                alt="ESP Surface"
                style={{ width: '100%' }}
                placeholder
              />
            ) : (
              <div style={{ textAlign: 'center', padding: 20, color: '#999' }}>
                ESP图片加载失败
              </div>
            )}
          </div>
        )}
      </div>
    );
  }

  // 能量数据卡片
  const EnergyCard = () => (
    <Card size="small" style={{ marginBottom: 16 }}>
      <Row gutter={[16, 16]}>
        <Col span={12}>
          <div style={{
            background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
            borderRadius: 8,
            padding: '16px 20px',
            color: '#fff'
          }}>
            <Text style={{ color: 'rgba(255,255,255,0.8)', fontSize: 12 }}>总能量</Text>
            <div style={{ fontSize: 20, fontWeight: 600, marginTop: 4 }}>
              {result.energy_au?.toFixed(6) || '-'}
              <span style={{ fontSize: 12, fontWeight: 400, marginLeft: 4 }}>A.U.</span>
            </div>
          </div>
        </Col>
        <Col span={12}>
          <div style={{
            background: 'linear-gradient(135deg, #f093fb 0%, #f5576c 100%)',
            borderRadius: 8,
            padding: '16px 20px',
            color: '#fff'
          }}>
            <Text style={{ color: 'rgba(255,255,255,0.8)', fontSize: 12 }}>HOMO-LUMO Gap</Text>
            <div style={{ fontSize: 20, fontWeight: 600, marginTop: 4 }}>
              {result.homo_lumo_gap?.toFixed(3) || '-'}
              <span style={{ fontSize: 12, fontWeight: 400, marginLeft: 4 }}>eV</span>
            </div>
          </div>
        </Col>
      </Row>
      <Row gutter={[16, 16]} style={{ marginTop: 16 }}>
        <Col span={6}>
          <Statistic
            title={<><Tag color="blue" style={{ marginRight: 4 }}>HOMO</Tag></>}
            value={homoEv?.toFixed(3) || '-'}
            suffix="eV"
            valueStyle={{ fontSize: 16, color: '#1890ff' }}
          />
        </Col>
        <Col span={6}>
          <Statistic
            title={<><Tag color="green" style={{ marginRight: 4 }}>LUMO</Tag></>}
            value={lumoEv?.toFixed(3) || '-'}
            suffix="eV"
            valueStyle={{ fontSize: 16, color: '#52c41a' }}
          />
        </Col>
        <Col span={6}>
          <Statistic
            title={<Space size={4}><FireOutlined style={{ color: '#1890ff' }} />ESP Min</Space>}
            value={result.esp_min_kcal?.toFixed(2) || '-'}
            suffix="kcal/mol"
            valueStyle={{ fontSize: 14 }}
          />
        </Col>
        <Col span={6}>
          <Statistic
            title={<Space size={4}><FireOutlined style={{ color: '#f5222d' }} />ESP Max</Space>}
            value={result.esp_max_kcal?.toFixed(2) || '-'}
            suffix="kcal/mol"
            valueStyle={{ fontSize: 14 }}
          />
        </Col>
      </Row>
    </Card>
  );

  // 分子轨道标签页
  const OrbitalTab = () => (
    <div>
      <Row gutter={16}>
        {/* HOMO */}
        <Col span={12}>
          <Card
            size="small"
            title={
              <Space>
                <div style={{
                  width: 10,
                  height: 10,
                  borderRadius: '50%',
                  background: '#1890ff'
                }} />
                <Text strong>HOMO</Text>
                <Tag color="blue">{homoEv?.toFixed(3)} eV</Tag>
              </Space>
            }
            extra={
              result.homo_image_path && (
                <Button
                  size="small"
                  type="text"
                  icon={<DownloadOutlined />}
                  onClick={() => handleDownload('HOMO')}
                  loading={downloading === 'HOMO'}
                />
              )
            }
            styles={{ body: { padding: 12 } }}
          >
            {homoLoading ? (
              <div style={{ textAlign: 'center', padding: 60 }}>
                <Spin tip="加载中..." />
              </div>
            ) : homoImageUrl ? (
              <div style={{
                textAlign: 'center',
                background: token.colorBgContainer,
                borderRadius: 8,
                padding: 8
              }}>
                <Image
                  src={homoImageUrl}
                  alt="HOMO Orbital"
                  style={{ maxWidth: '100%', maxHeight: 280 }}
                  placeholder
                />
              </div>
            ) : (
              <div style={{
                textAlign: 'center',
                padding: 60,
                color: token.colorTextSecondary,
                background: token.colorBgContainer,
                borderRadius: 8
              }}>
                HOMO图片暂未生成
              </div>
            )}
          </Card>
        </Col>

        {/* LUMO */}
        <Col span={12}>
          <Card
            size="small"
            title={
              <Space>
                <div style={{
                  width: 10,
                  height: 10,
                  borderRadius: '50%',
                  background: '#52c41a'
                }} />
                <Text strong>LUMO</Text>
                <Tag color="green">{lumoEv?.toFixed(3)} eV</Tag>
              </Space>
            }
            extra={
              result.lumo_image_path && (
                <Button
                  size="small"
                  type="text"
                  icon={<DownloadOutlined />}
                  onClick={() => handleDownload('LUMO')}
                  loading={downloading === 'LUMO'}
                />
              )
            }
            styles={{ body: { padding: 12 } }}
          >
            {lumoLoading ? (
              <div style={{ textAlign: 'center', padding: 60 }}>
                <Spin tip="加载中..." />
              </div>
            ) : lumoImageUrl ? (
              <div style={{
                textAlign: 'center',
                background: token.colorBgContainer,
                borderRadius: 8,
                padding: 8
              }}>
                <Image
                  src={lumoImageUrl}
                  alt="LUMO Orbital"
                  style={{ maxWidth: '100%', maxHeight: 280 }}
                  placeholder
                />
              </div>
            ) : (
              <div style={{
                textAlign: 'center',
                padding: 60,
                color: token.colorTextSecondary,
                background: token.colorBgContainer,
                borderRadius: 8
              }}>
                LUMO图片暂未生成
              </div>
            )}
          </Card>
        </Col>
      </Row>
      <div style={{
        textAlign: 'center',
        marginTop: 12,
        padding: '8px 16px',
        background: '#f5f5f5',
        borderRadius: 4
      }}>
        <Space split={<Divider type="vertical" />}>
          <Space size={4}>
            <div style={{ width: 12, height: 12, background: '#1890ff', borderRadius: 2 }} />
            <Text type="secondary">正相位</Text>
          </Space>
          <Space size={4}>
            <div style={{ width: 12, height: 12, background: '#ff4d4f', borderRadius: 2 }} />
            <Text type="secondary">负相位</Text>
          </Space>
        </Space>
      </div>
    </div>
  );

  // ESP标签页
  const ESPTab = () => (
    <div>
      <div style={{
        display: 'flex',
        justifyContent: 'space-between',
        alignItems: 'center',
        marginBottom: 12,
        padding: '8px 12px',
        background: token.colorBgContainer,
        borderRadius: 8
      }}>
        <Space>
          <BgColorsOutlined style={{ color: token.colorPrimary }} />
          <Text>静电势表面映射</Text>
        </Space>
        <Button
          type="primary"
          size="small"
          icon={<DownloadOutlined />}
          onClick={() => handleDownload('ESP')}
          loading={downloading === 'ESP'}
        >
          下载图片
        </Button>
      </div>
      {espLoading ? (
        <div style={{ textAlign: 'center', padding: 80 }}>
          <Spin tip="加载ESP图片中..." size="large" />
        </div>
      ) : espImageUrl ? (
        <div style={{
          textAlign: 'center',
          background: token.colorBgContainer,
          borderRadius: 8,
          padding: 16
        }}>
          <Image
            src={espImageUrl}
            alt="ESP Surface"
            style={{ maxWidth: '100%', maxHeight: 450 }}
            placeholder
          />
        </div>
      ) : (
        <div style={{
          textAlign: 'center',
          padding: 80,
          color: token.colorTextSecondary,
          background: token.colorBgContainer,
          borderRadius: 8
        }}>
          ESP图片加载失败，请刷新页面重试
        </div>
      )}
      <div style={{
        textAlign: 'center',
        marginTop: 12,
        padding: '8px 16px',
        background: mode === 'dark' ? 'rgba(255,255,255,0.04)' : '#f5f5f5',
        borderRadius: 4
      }}>
        <Space split={<Divider type="vertical" />}>
          <Space size={4}>
            <div style={{
              width: 60,
              height: 12,
              background: 'linear-gradient(to right, #1890ff, #fff, #ff4d4f)',
              borderRadius: 2
            }} />
          </Space>
          <Text type="secondary">蓝色(负电势) → 白色(中性) → 红色(正电势)</Text>
        </Space>
      </div>
      <Row gutter={16} style={{ marginTop: 16 }}>
        <Col span={12}>
          <Card size="small" style={{ background: mode === 'dark' ? 'rgba(24, 144, 255, 0.15)' : '#e6f7ff', border: `1px solid ${token.colorPrimary}` }}>
            <Statistic
              title="ESP 最小值"
              value={result.esp_min_kcal?.toFixed(2) || '-'}
              suffix="kcal/mol"
              valueStyle={{ color: token.colorPrimary }}
              prefix={<FireOutlined />}
            />
          </Card>
        </Col>
        <Col span={12}>
          <Card size="small" style={{ background: mode === 'dark' ? 'rgba(255, 77, 79, 0.15)' : '#fff1f0', border: `1px solid ${token.colorError}` }}>
            <Statistic
              title="ESP 最大值"
              value={result.esp_max_kcal?.toFixed(2) || '-'}
              suffix="kcal/mol"
              valueStyle={{ color: token.colorError }}
              prefix={<FireOutlined />}
            />
          </Card>
        </Col>
      </Row>
    </div>
  );

  // 构建标签页
  const tabItems = [
    {
      key: 'orbital',
      label: (
        <Space>
          <ApartmentOutlined />
          分子轨道
        </Space>
      ),
      children: <OrbitalTab />,
    },
    {
      key: 'esp',
      label: (
        <Space>
          <BgColorsOutlined />
          静电势表面
        </Space>
      ),
      children: <ESPTab />,
    },
    ...(job ? [{
      key: 'params',
      label: (
        <Space>
          <ExperimentOutlined />
          计算参数
        </Space>
      ),
      children: (
        <Descriptions column={2} size="small" bordered>
          <Descriptions.Item label="泛函">
            <Tag color="blue">{job.functional || '-'}</Tag>
          </Descriptions.Item>
          <Descriptions.Item label="基组">
            <Tag color="green">{job.basis_set || '-'}</Tag>
          </Descriptions.Item>
          <Descriptions.Item label="溶剂模型">
            {(() => {
              const solventConfig = job.config?.solvent_config;
              if (!solventConfig) {
                return <Tag color="default">气相</Tag>;
              }
              const model = solventConfig.model || 'gas';
              const solventName = solventConfig.solvent_name;

              if (model === 'gas') {
                return <Tag color="default">气相</Tag>;
              } else if (model === 'pcm') {
                return <Tag color="cyan">PCM - {solventName || 'Water'}</Tag>;
              } else if (model === 'smd') {
                return <Tag color="blue">SMD - {solventName || 'Water'}</Tag>;
              } else if (model === 'custom') {
                return <Tag color="purple">自定义 - {solventName || '自定义溶剂'}</Tag>;
              }
              return <Tag color="default">{model}</Tag>;
            })()}
          </Descriptions.Item>
          <Descriptions.Item label="电荷">
            {job.charge || 0}
          </Descriptions.Item>
          <Descriptions.Item label="自旋多重度">
            {job.spin_multiplicity || 1}
          </Descriptions.Item>
          <Descriptions.Item label="计算分区">
            {job.config?.slurm_partition || '-'}
          </Descriptions.Item>
          <Descriptions.Item label="CPU核数">
            {job.config?.slurm_cpus || 16} 核
          </Descriptions.Item>
        </Descriptions>
      ),
    }] : []),
    {
      key: 'data',
      label: (
        <Space>
          <InfoCircleOutlined />
          原始数据
        </Space>
      ),
      children: (
        <Descriptions column={2} size="small" bordered>
          <Descriptions.Item label="总能量 (A.U.)">
            <Text copyable>{result.energy_au?.toFixed(8) || '-'}</Text>
          </Descriptions.Item>
          <Descriptions.Item label="HOMO-LUMO Gap (eV)">
            {result.homo_lumo_gap?.toFixed(6) || '-'}
          </Descriptions.Item>
          <Descriptions.Item label="HOMO (Hartree)">
            <Text copyable>{result.homo?.toFixed(8) || '-'}</Text>
          </Descriptions.Item>
          <Descriptions.Item label="LUMO (Hartree)">
            <Text copyable>{result.lumo?.toFixed(8) || '-'}</Text>
          </Descriptions.Item>
          <Descriptions.Item label="HOMO (eV)">
            {homoEv?.toFixed(6) || '-'}
          </Descriptions.Item>
          <Descriptions.Item label="LUMO (eV)">
            {lumoEv?.toFixed(6) || '-'}
          </Descriptions.Item>
          <Descriptions.Item label="ESP Min (kcal/mol)">
            {result.esp_min_kcal?.toFixed(4) || '-'}
          </Descriptions.Item>
          <Descriptions.Item label="ESP Max (kcal/mol)">
            {result.esp_max_kcal?.toFixed(4) || '-'}
          </Descriptions.Item>
          <Descriptions.Item label="计算时间">
            {result.created_at ? new Date(result.created_at).toLocaleString('zh-CN') : '-'}
          </Descriptions.Item>
          <Descriptions.Item label="结果ID">
            {result.id}
          </Descriptions.Item>
        </Descriptions>
      ),
    },
  ];

  // 完整模式
  return (
    <div>
      <EnergyCard />

      <Card
        size="small"
        styles={{ body: { padding: '12px 16px' } }}
      >
        <Tabs
          items={tabItems}
          defaultActiveKey="orbital"
          size="small"
        />
      </Card>
    </div>
  );
}

