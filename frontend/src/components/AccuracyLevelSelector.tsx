import React, { useEffect, useState } from 'react';
import { Card, Radio, Space, Typography, Tag, Spin, Alert, Descriptions } from 'antd';
import { RocketOutlined, ExperimentOutlined, AimOutlined, SettingOutlined } from '@ant-design/icons';
import axios from 'axios';

const { Title, Text, Paragraph } = Typography;

interface AccuracyConfig {
  name: string;
  description: string;
  charge_method: string;
  nsteps_npt: number;
  nsteps_nvt: number;
  freq_trj_npt: number;
  freq_trj_nvt: number;
  thermo_freq: number;
  estimated_time_hours: number;
  recommended_for: string;
  color: string;
  icon: string;
}

interface AccuracyLevels {
  fast: AccuracyConfig;
  standard: AccuracyConfig;
  accurate: AccuracyConfig;
}

interface Props {
  value?: string;
  onChange?: (value: string) => void;
}

const AccuracyLevelSelector: React.FC<Props> = ({ value = 'standard', onChange }) => {
  const [levels, setLevels] = useState<AccuracyLevels | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    fetchAccuracyLevels();
  }, []);

  const fetchAccuracyLevels = async () => {
    try {
      const token = localStorage.getItem('access_token');
      const response = await axios.get('/api/v1/jobs/accuracy-levels', {
        headers: { Authorization: `Bearer ${token}` }
      });
      setLevels(response.data);
      setLoading(false);
    } catch (err: any) {
      setError(err.response?.data?.detail || '加载精度等级配置失败');
      setLoading(false);
    }
  };

  const getIcon = (level: string) => {
    switch (level) {
      case 'fast':
        return <RocketOutlined style={{ fontSize: 24 }} />;
      case 'standard':
        return <ExperimentOutlined style={{ fontSize: 24 }} />;
      case 'accurate':
        return <AimOutlined style={{ fontSize: 24 }} />;
      case 'custom':
        return <SettingOutlined style={{ fontSize: 24 }} />;
      default:
        return null;
    }
  };

  const formatNumber = (num: number | null) => {
    if (num === null || num === undefined) {
      return '自定义';
    }
    return num.toLocaleString();
  };

  const formatTime = (steps: number | null, timestep: number = 1.0) => {
    if (steps === null || steps === undefined) {
      return '自定义';
    }

    const timeFs = steps * timestep; // 总时间（飞秒）

    // 根据时间大小选择合适的单位
    if (timeFs >= 1_000_000) {
      // >= 1 ns，使用 ns
      const timeNs = timeFs / 1_000_000;
      return `${timeNs.toFixed(1)} ns`;
    } else if (timeFs >= 1_000) {
      // >= 1 ps，使用 ps
      const timePs = timeFs / 1_000;
      return `${timePs.toFixed(1)} ps`;
    } else {
      // < 1 ps，使用 fs
      return `${timeFs.toFixed(0)} fs`;
    }
  };

  if (loading) {
    return <Spin tip="加载精度等级配置..." />;
  }

  if (error) {
    return <Alert message="错误" description={error} type="error" showIcon />;
  }

  if (!levels) {
    return null;
  }

  return (
    <div>
      <Title level={5}>计算精度等级</Title>
      <Paragraph type="secondary">
        选择合适的精度等级，系统会自动配置电荷计算方法和模拟参数
      </Paragraph>

      <Radio.Group
        value={value}
        onChange={(e) => onChange?.(e.target.value)}
        style={{ width: '100%' }}
      >
        <Space direction="vertical" style={{ width: '100%' }} size="middle">
          {Object.entries(levels).map(([key, config]) => (
            <Radio key={key} value={key} style={{ width: '100%' }}>
              <Card
                hoverable
                style={{
                  borderColor: value === key ? config.color : undefined,
                  borderWidth: value === key ? 2 : 1,
                }}
              >
                <Space direction="vertical" style={{ width: '100%' }}>
                  <Space>
                    {getIcon(key)}
                    <Title level={4} style={{ margin: 0 }}>
                      {config.icon} {config.name}
                    </Title>
                    {config.estimated_time_hours && (
                      <Tag color={config.color}>
                        预计 {config.estimated_time_hours} 小时
                      </Tag>
                    )}
                  </Space>

                  <Paragraph style={{ margin: 0 }}>
                    {config.description}
                  </Paragraph>

                  <Descriptions size="small" column={2}>
                    <Descriptions.Item label="电荷方法">
                      <Tag color={config.charge_method === 'resp' ? 'red' : 'blue'}>
                        {config.charge_method === 'resp' ? 'RESP (高精度)' : 'LigParGen CM1A'}
                      </Tag>
                    </Descriptions.Item>
                    <Descriptions.Item label="NPT 步数">
                      {formatNumber(config.nsteps_npt)} ({formatTime(config.nsteps_npt)})
                    </Descriptions.Item>
                    <Descriptions.Item label="NVT 步数">
                      {formatNumber(config.nsteps_nvt)} ({formatTime(config.nsteps_nvt)})
                    </Descriptions.Item>
                    <Descriptions.Item label="轨迹输出">
                      每 {formatNumber(config.freq_trj_nvt)} 步
                    </Descriptions.Item>
                  </Descriptions>

                  <Text type="secondary" style={{ fontSize: 12 }}>
                    💡 适用场景：{config.recommended_for}
                  </Text>
                </Space>
              </Card>
            </Radio>
          ))}
        </Space>
      </Radio.Group>
    </div>
  );
};

export default AccuracyLevelSelector;

