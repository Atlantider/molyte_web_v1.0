/**
 * 任务基本信息组件 - 统一设计风格
 */
import { Card, Descriptions, Row, Col, Tag, Space, Typography, Spin, Alert } from 'antd';
import {
  ExperimentOutlined,
  DatabaseOutlined,
  SettingOutlined,
  FundOutlined,
  ThunderboltOutlined,
} from '@ant-design/icons';
import type { MDJob, ElectrolyteSystem } from '../types';
import { JobStatus } from '../types';
import { getStructureInfo, type StructureInfo } from '../api/jobs';
import dayjs from 'dayjs';
import duration from 'dayjs/plugin/duration';
import { useEffect, useState } from 'react';

dayjs.extend(duration);

const { Text } = Typography;

// 统一的设计风格常量
const DASHBOARD_STYLES = {
  pageBackground: '#F5F7FB',
  cardBackground: '#FFFFFF',
  cardBorderRadius: 12,
  cardShadow: '0 4px 12px rgba(15, 23, 42, 0.08)',
  cardShadowHover: '0 8px 24px rgba(15, 23, 42, 0.12)',
  cardPadding: 24,
  gutter: 24,
  titleColor: '#111827',
  titleFontSize: 16,
  titleFontWeight: 600,
};

interface JobBasicInfoProps {
  job: MDJob;
  electrolyte: ElectrolyteSystem;
  slurmStatus?: any;
}

export default function JobBasicInfo({ job, electrolyte, slurmStatus }: JobBasicInfoProps) {
  const [structureInfo, setStructureInfo] = useState<StructureInfo | null>(null);
  const [loadingStructure, setLoadingStructure] = useState(false);

  // 调试信息
  useEffect(() => {
    console.log('=== JobBasicInfo Debug ===');
    console.log('job:', job);
    console.log('job.started_at:', job.started_at);
    console.log('job.finished_at:', job.finished_at);
    console.log('slurmStatus:', slurmStatus);
    console.log('electrolyte:', electrolyte);
  }, [job, slurmStatus, electrolyte]);

  // 加载结构信息
  useEffect(() => {
    const loadStructureInfo = async () => {
      // 检查任务是否完成
      if (job.status === JobStatus.COMPLETED || job.status === JobStatus.POSTPROCESSING) {
        setLoadingStructure(true);
        try {
          const info = await getStructureInfo(job.id);
          console.log('Structure info loaded:', info);
          setStructureInfo(info);
        } catch (error) {
          console.error('Failed to load structure info:', error);
        } finally {
          setLoadingStructure(false);
        }
      }
    };

    loadStructureInfo();
  }, [job.id, job.status]);
  // 计算运行时间
  const getRunningTime = () => {
    if (job.started_at) {
      const end = job.finished_at ? dayjs(job.finished_at) : dayjs();
      const start = dayjs(job.started_at);
      const diff = end.diff(start);
      const dur = dayjs.duration(diff);
      
      const days = Math.floor(dur.asDays());
      const hours = dur.hours();
      const minutes = dur.minutes();
      const seconds = dur.seconds();
      
      if (days > 0) {
        return `${days}天 ${hours}小时 ${minutes}分钟`;
      } else if (hours > 0) {
        return `${hours}小时 ${minutes}分钟 ${seconds}秒`;
      } else if (minutes > 0) {
        return `${minutes}分钟 ${seconds}秒`;
      } else {
        return `${seconds}秒`;
      }
    }
    return '-';
  };

  // 计算总模拟时间（ps）
  const getTotalSimulationTime = () => {
    const nptSteps = job.config?.nsteps_npt || electrolyte.nsteps_npt || 0;
    const nvtSteps = job.config?.nsteps_nvt || electrolyte.nsteps_nvt || 0;
    const timestep = job.config?.timestep || electrolyte.timestep || 1.0;
    
    const totalSteps = nptSteps + nvtSteps;
    const totalTime = (totalSteps * timestep) / 1000; // fs -> ps
    
    return totalTime.toLocaleString(undefined, { maximumFractionDigits: 0 });
  };

  // 计算CPU核时（core-hours）
  const getCoreHours = () => {
    const ntasks = job.config?.slurm_ntasks || 1;
    const cpusPerTask = job.config?.slurm_cpus_per_task || 1;
    const totalCores = ntasks * cpusPerTask;

    let hours = 0;

    // 优先使用 slurmStatus.elapsed
    if (slurmStatus?.elapsed) {
      // 解析 elapsed 时间（格式：DD-HH:MM:SS 或 HH:MM:SS）
      const elapsedParts = slurmStatus.elapsed.split(/[-:]/);

      if (elapsedParts.length === 4) {
        // DD-HH:MM:SS
        hours = parseInt(elapsedParts[0]) * 24 + parseInt(elapsedParts[1]) +
                parseInt(elapsedParts[2]) / 60 + parseInt(elapsedParts[3]) / 3600;
      } else if (elapsedParts.length === 3) {
        // HH:MM:SS
        hours = parseInt(elapsedParts[0]) + parseInt(elapsedParts[1]) / 60 +
                parseInt(elapsedParts[2]) / 3600;
      }
    } else if (job.started_at) {
      // 如果没有 slurmStatus，使用 started_at 和 finished_at 计算
      const end = job.finished_at ? dayjs(job.finished_at) : dayjs();
      const start = dayjs(job.started_at);
      hours = end.diff(start, 'hour', true); // 精确到小数
    } else {
      return '-';
    }

    const coreHours = totalCores * hours;
    return coreHours.toFixed(1);
  };

  const dashboardCardStyle = {
    background: DASHBOARD_STYLES.cardBackground,
    borderRadius: DASHBOARD_STYLES.cardBorderRadius,
    boxShadow: DASHBOARD_STYLES.cardShadow,
    border: '1px solid #e8e8e8',
  };

  return (
    <div style={{ background: DASHBOARD_STYLES.pageBackground, padding: DASHBOARD_STYLES.gutter }}>
      <Row gutter={[DASHBOARD_STYLES.gutter, DASHBOARD_STYLES.gutter]}>
        {/* 1. 任务信息（100%宽度） */}
        <Col xs={24}>
          <Card
            className="dashboard-card"
            style={dashboardCardStyle}
            title={
              <Space size={8}>
                <DatabaseOutlined style={{ color: DASHBOARD_STYLES.titleColor }} />
                <span style={{ fontSize: 14, fontWeight: 600, color: DASHBOARD_STYLES.titleColor }}>
                  任务信息
                </span>
              </Space>
            }
          >
            <Row gutter={16}>
              {/* 左侧：基本信息 */}
              <Col xs={24} lg={12}>
                <Descriptions column={1} size="small" bordered>
                  <Descriptions.Item label="任务 ID">
                    <Text strong>#{job.id}</Text>
                  </Descriptions.Item>
                  <Descriptions.Item label="任务名称">
                    <Text strong>{job.config?.job_name || '-'}</Text>
                  </Descriptions.Item>
                  <Descriptions.Item label="Slurm Job ID">
                    <Text code>{job.slurm_job_id || job.config?.slurm_job_id || '-'}</Text>
                  </Descriptions.Item>
                  <Descriptions.Item label="计算分区">
                    <Tag color="blue">{job.config?.slurm_partition || '-'}</Tag>
                  </Descriptions.Item>
                  <Descriptions.Item label="计算资源">
                    {job.config?.slurm_ntasks && job.config?.slurm_cpus_per_task ? (
                      <Space>
                        <Text>{job.config.slurm_ntasks} 任务 × {job.config.slurm_cpus_per_task} 核/任务 =</Text>
                        <Text strong style={{ color: '#52c41a', fontSize: 16 }}>
                          {job.config.slurm_ntasks * job.config.slurm_cpus_per_task} 核
                        </Text>
                      </Space>
                    ) : '-'}
                  </Descriptions.Item>
                </Descriptions>
              </Col>

              {/* 右侧：时间信息 */}
              <Col xs={24} lg={12}>
                <Descriptions column={1} size="small" bordered>
                  <Descriptions.Item label="创建时间">
                    {dayjs(job.created_at).format('YYYY-MM-DD HH:mm:ss')}
                  </Descriptions.Item>
                  <Descriptions.Item label="开始时间">
                    {job.started_at ? dayjs(job.started_at).format('YYYY-MM-DD HH:mm:ss') : '-'}
                  </Descriptions.Item>
                  <Descriptions.Item label="结束时间">
                    {job.finished_at ? dayjs(job.finished_at).format('YYYY-MM-DD HH:mm:ss') : '-'}
                  </Descriptions.Item>
                  <Descriptions.Item label="运行时长">
                    <Text strong style={{ color: '#1890ff', fontSize: 15 }}>{getRunningTime()}</Text>
                  </Descriptions.Item>
                  <Descriptions.Item label="CPU 核时">
                    <Space>
                      <Text strong style={{ color: '#fa8c16', fontSize: 16 }}>{getCoreHours()}</Text>
                      <Text type="secondary">小时</Text>
                    </Space>
                  </Descriptions.Item>
                </Descriptions>
              </Col>

              {/* 工作目录（全宽） */}
              <Col xs={24} style={{ marginTop: 16 }}>
                <div style={{
                  padding: '12px 16px',
                  background: '#fafafa',
                  borderRadius: 8,
                  border: '1px solid #e8e8e8'
                }}>
                  <Space direction="vertical" size={4} style={{ width: '100%' }}>
                    <Text type="secondary" style={{ fontSize: 12 }}>工作目录</Text>
                    <Text code style={{ fontSize: 12, wordBreak: 'break-all' }}>
                      {job.work_dir || '-'}
                    </Text>
                  </Space>
                </div>
              </Col>
            </Row>
          </Card>
        </Col>

        {/* 2. 配方信息（100%宽度） */}
        <Col xs={24}>
          <Card
            className="dashboard-card"
            style={dashboardCardStyle}
            title={
              <Space size={8}>
                <ExperimentOutlined style={{ color: DASHBOARD_STYLES.titleColor }} />
                <span style={{ fontSize: 14, fontWeight: 600, color: DASHBOARD_STYLES.titleColor }}>
                  溶液配方
                </span>
              </Space>
            }
          >
            <Row gutter={16}>
              {/* 左侧：基本参数 */}
              <Col xs={24} lg={8}>
                <div style={{
                  padding: '16px',
                  background: 'linear-gradient(135deg, #f0f9ff 0%, #e6f7ff 100%)',
                  borderRadius: 8,
                  border: '1px solid #91d5ff',
                  height: '100%'
                }}>
                  <Space direction="vertical" size={12} style={{ width: '100%' }}>
                    <div>
                      <Text type="secondary" style={{ fontSize: 12 }}>配方名称</Text>
                      <div style={{ marginTop: 4 }}>
                        <Text strong style={{ fontSize: 15 }}>{electrolyte.name}</Text>
                      </div>
                    </div>
                    <div>
                      <Text type="secondary" style={{ fontSize: 12 }}>温度</Text>
                      <div style={{ marginTop: 4 }}>
                        <Text strong style={{ fontSize: 15, color: '#fa8c16' }}>
                          {job.config?.temperature || electrolyte.temperature} K
                        </Text>
                      </div>
                    </div>
                    <div>
                      <Text type="secondary" style={{ fontSize: 12 }}>压力</Text>
                      <div style={{ marginTop: 4 }}>
                        <Text strong style={{ fontSize: 15, color: '#1890ff' }}>
                          {job.config?.pressure || electrolyte.pressure} atm
                        </Text>
                      </div>
                    </div>
                  </Space>
                </div>
              </Col>

              {/* 右侧：组分详情 */}
              <Col xs={24} lg={16}>
                <div style={{
                  padding: '16px',
                  background: '#fafafa',
                  borderRadius: 8,
                  border: '1px solid #e8e8e8',
                  height: '100%'
                }}>
                  <Text strong style={{ fontSize: 13, marginBottom: 12, display: 'block', color: '#595959' }}>
                    组分详情
                  </Text>
                  <Row gutter={[8, 8]}>
                    {electrolyte.cations.map((cation, idx) => (
                      <Col xs={24} sm={12} md={8} key={`cation-${idx}`}>
                        <div style={{
                          padding: '10px 12px',
                          background: 'linear-gradient(135deg, #fff1f0 0%, #ffccc7 100%)',
                          borderRadius: 6,
                          border: '1px solid #ffa39e',
                          boxShadow: '0 2px 4px rgba(255, 77, 79, 0.1)'
                        }}>
                          <Space direction="vertical" size={2} style={{ width: '100%' }}>
                            <Tag color="red" style={{ margin: 0 }}>{cation.name}</Tag>
                            <Text style={{ fontSize: 12, color: '#595959' }}>
                              数量: <Text strong style={{ fontSize: 14, color: '#cf1322' }}>{cation.number}</Text>
                            </Text>
                          </Space>
                        </div>
                      </Col>
                    ))}
                    {electrolyte.anions.map((anion, idx) => (
                      <Col xs={24} sm={12} md={8} key={`anion-${idx}`}>
                        <div style={{
                          padding: '10px 12px',
                          background: 'linear-gradient(135deg, #e6f7ff 0%, #bae7ff 100%)',
                          borderRadius: 6,
                          border: '1px solid #91d5ff',
                          boxShadow: '0 2px 4px rgba(24, 144, 255, 0.1)'
                        }}>
                          <Space direction="vertical" size={2} style={{ width: '100%' }}>
                            <Tag color="blue" style={{ margin: 0 }}>{anion.name}</Tag>
                            <Text style={{ fontSize: 12, color: '#595959' }}>
                              数量: <Text strong style={{ fontSize: 14, color: '#096dd9' }}>{anion.number}</Text>
                            </Text>
                          </Space>
                        </div>
                      </Col>
                    ))}
                    {electrolyte.solvents && electrolyte.solvents.map((solvent, idx) => (
                      <Col xs={24} sm={12} md={8} key={`solvent-${idx}`}>
                        <div style={{
                          padding: '10px 12px',
                          background: 'linear-gradient(135deg, #f6ffed 0%, #d9f7be 100%)',
                          borderRadius: 6,
                          border: '1px solid #b7eb8f',
                          boxShadow: '0 2px 4px rgba(82, 196, 26, 0.1)'
                        }}>
                          <Space direction="vertical" size={2} style={{ width: '100%' }}>
                            <Tag color="green" style={{ margin: 0 }}>{solvent.name}</Tag>
                            <Text style={{ fontSize: 12, color: '#595959' }}>
                              数量: <Text strong style={{ fontSize: 14, color: '#389e0d' }}>{solvent.number}</Text>
                            </Text>
                          </Space>
                        </div>
                      </Col>
                    ))}
                  </Row>
                </div>
              </Col>
            </Row>
          </Card>
        </Col>

        {/* 3. 浓度对比（100%宽度） */}
        <Col xs={24}>
          <Card
            className="dashboard-card"
            style={dashboardCardStyle}
            title={
              <Space size={8}>
                <FundOutlined style={{ color: DASHBOARD_STYLES.titleColor }} />
                <span style={{ fontSize: 14, fontWeight: 600, color: DASHBOARD_STYLES.titleColor }}>
                  浓度与盒子尺寸对比
                </span>
              </Space>
            }
          >
            {loadingStructure ? (
              <div style={{ textAlign: 'center', padding: '60px 0' }}>
                <Spin size="large" />
                <div style={{ marginTop: 16, color: '#8c8c8c' }}>加载计算结果...</div>
              </div>
            ) : structureInfo?.available ? (
              <Row gutter={16}>
                {/* 左侧：初始值 */}
                <Col xs={24} lg={12}>
                  <div style={{
                    padding: '20px',
                    background: 'linear-gradient(135deg, #fafafa 0%, #f0f0f0 100%)',
                    borderRadius: 8,
                    border: '1px solid #d9d9d9',
                    height: '100%'
                  }}>
                    <div style={{
                      marginBottom: 16,
                      paddingBottom: 12,
                      borderBottom: '2px solid #d9d9d9'
                    }}>
                      <Text strong style={{ fontSize: 15, color: '#595959' }}>
                        📋 初始设置
                      </Text>
                    </div>
                    <Space direction="vertical" size={16} style={{ width: '100%' }}>
                      <div>
                        <Text type="secondary" style={{ fontSize: 12 }}>初始浓度</Text>
                        <div style={{ marginTop: 6 }}>
                          <Text strong style={{ fontSize: 20, color: '#8c8c8c' }}>
                            {structureInfo.initial_concentration?.toFixed(4) || '-'}
                          </Text>
                          <Text type="secondary" style={{ marginLeft: 8, fontSize: 13 }}>mol/L</Text>
                        </div>
                      </div>
                      <div>
                        <Text type="secondary" style={{ fontSize: 12 }}>初始盒子尺寸</Text>
                        <div style={{ marginTop: 6 }}>
                          <Text code style={{ fontSize: 14, color: '#595959' }}>
                            {structureInfo.initial_box_dimensions || '-'}
                          </Text>
                          <Text type="secondary" style={{ marginLeft: 8, fontSize: 12 }}>Å</Text>
                        </div>
                      </div>
                    </Space>
                  </div>
                </Col>

                {/* 右侧：计算结果 */}
                <Col xs={24} lg={12}>
                  <div style={{
                    padding: '20px',
                    background: 'linear-gradient(135deg, #e6fffb 0%, #b5f5ec 100%)',
                    borderRadius: 8,
                    border: '2px solid #5cdbd3',
                    boxShadow: '0 4px 12px rgba(19, 194, 194, 0.15)',
                    height: '100%'
                  }}>
                    <div style={{
                      marginBottom: 16,
                      paddingBottom: 12,
                      borderBottom: '2px solid #5cdbd3'
                    }}>
                      <Text strong style={{ fontSize: 15, color: '#006d75' }}>
                        ✨ 计算结果
                      </Text>
                    </div>
                    <Space direction="vertical" size={16} style={{ width: '100%' }}>
                      <div>
                        <Text type="secondary" style={{ fontSize: 12 }}>计算浓度</Text>
                        <div style={{ marginTop: 6 }}>
                          <Text strong style={{ fontSize: 24, color: '#13c2c2' }}>
                            {structureInfo.concentration?.toFixed(4) || '-'}
                          </Text>
                          <Text type="secondary" style={{ marginLeft: 8, fontSize: 13 }}>mol/L</Text>
                        </div>
                        {structureInfo.initial_concentration && structureInfo.concentration && (
                          <div style={{ marginTop: 6 }}>
                            <Tag color={
                              Math.abs((structureInfo.concentration - structureInfo.initial_concentration) / structureInfo.initial_concentration * 100) < 5
                                ? 'success'
                                : 'warning'
                            }>
                              偏差: {((structureInfo.concentration - structureInfo.initial_concentration) / structureInfo.initial_concentration * 100).toFixed(2)}%
                            </Tag>
                          </div>
                        )}
                      </div>
                      <div>
                        <Text type="secondary" style={{ fontSize: 12 }}>最终盒子尺寸</Text>
                        <div style={{ marginTop: 6 }}>
                          <Text code strong style={{ fontSize: 15, color: '#006d75' }}>
                            {structureInfo.box_dimensions || '-'}
                          </Text>
                          <Text type="secondary" style={{ marginLeft: 8, fontSize: 12 }}>Å</Text>
                        </div>
                      </div>
                    </Space>
                  </div>
                </Col>
              </Row>
            ) : (
              <Alert
                message="计算结果未就绪"
                description={
                  <div>
                    <p>任务完成后将显示浓度和盒子尺寸的计算结果</p>
                    <p style={{ marginTop: 8, fontSize: 12 }}>
                      当前状态: <Tag>{job.status}</Tag>
                    </p>
                  </div>
                }
                type="info"
                showIcon
                style={{ margin: '20px 0' }}
              />
            )}
          </Card>
        </Col>

        {/* 4. 密度对比（100%宽度） */}
        <Col xs={24}>
          <Card
            className="dashboard-card"
            style={dashboardCardStyle}
            title={
              <Space size={8}>
                <FundOutlined style={{ color: DASHBOARD_STYLES.titleColor }} />
                <span style={{ fontSize: 14, fontWeight: 600, color: DASHBOARD_STYLES.titleColor }}>
                  密度对比
                </span>
              </Space>
            }
          >
            {loadingStructure ? (
              <div style={{ textAlign: 'center', padding: '60px 0' }}>
                <Spin size="large" />
                <div style={{ marginTop: 16, color: '#8c8c8c' }}>加载计算结果...</div>
              </div>
            ) : structureInfo?.available ? (
              <Row gutter={16}>
                {/* 左侧：初始值 */}
                <Col xs={24} lg={12}>
                  <div style={{
                    padding: '20px',
                    background: 'linear-gradient(135deg, #fafafa 0%, #f0f0f0 100%)',
                    borderRadius: 8,
                    border: '1px solid #d9d9d9',
                    height: '100%',
                    display: 'flex',
                    flexDirection: 'column',
                    justifyContent: 'center'
                  }}>
                    <div style={{
                      marginBottom: 16,
                      paddingBottom: 12,
                      borderBottom: '2px solid #d9d9d9'
                    }}>
                      <Text strong style={{ fontSize: 15, color: '#595959' }}>
                        📋 初始设置
                      </Text>
                    </div>
                    <div>
                      <Text type="secondary" style={{ fontSize: 12 }}>初始密度</Text>
                      <div style={{ marginTop: 6 }}>
                        <Text strong style={{ fontSize: 24, color: '#8c8c8c' }}>
                          {structureInfo.initial_density?.toFixed(4) || '-'}
                        </Text>
                        <Text type="secondary" style={{ marginLeft: 8, fontSize: 13 }}>g/cm³</Text>
                      </div>
                    </div>
                  </div>
                </Col>

                {/* 右侧：计算结果 */}
                <Col xs={24} lg={12}>
                  <div style={{
                    padding: '20px',
                    background: 'linear-gradient(135deg, #e6f7ff 0%, #bae7ff 100%)',
                    borderRadius: 8,
                    border: '2px solid #69c0ff',
                    boxShadow: '0 4px 12px rgba(24, 144, 255, 0.15)',
                    height: '100%',
                    display: 'flex',
                    flexDirection: 'column',
                    justifyContent: 'center'
                  }}>
                    <div style={{
                      marginBottom: 16,
                      paddingBottom: 12,
                      borderBottom: '2px solid #69c0ff'
                    }}>
                      <Text strong style={{ fontSize: 15, color: '#0050b3' }}>
                        ✨ 计算结果
                      </Text>
                    </div>
                    <div>
                      <Text type="secondary" style={{ fontSize: 12 }}>计算密度</Text>
                      <div style={{ marginTop: 6 }}>
                        <Text strong style={{ fontSize: 28, color: '#1890ff' }}>
                          {structureInfo.density?.toFixed(4) || '-'}
                        </Text>
                        <Text type="secondary" style={{ marginLeft: 8, fontSize: 13 }}>g/cm³</Text>
                      </div>
                      {structureInfo.initial_density && structureInfo.density && (
                        <div style={{ marginTop: 8 }}>
                          <Tag color={
                            Math.abs((structureInfo.density - structureInfo.initial_density) / structureInfo.initial_density * 100) < 5
                              ? 'success'
                              : 'warning'
                          } style={{ fontSize: 13 }}>
                            偏差: {((structureInfo.density - structureInfo.initial_density) / structureInfo.initial_density * 100).toFixed(2)}%
                          </Tag>
                        </div>
                      )}
                    </div>
                  </div>
                </Col>
              </Row>
            ) : (
              <Alert
                message="计算结果未就绪"
                description={
                  <div>
                    <p>任务完成后将显示密度计算结果</p>
                    <p style={{ marginTop: 8, fontSize: 12 }}>
                      当前状态: <Tag>{job.status}</Tag>
                    </p>
                  </div>
                }
                type="info"
                showIcon
                style={{ margin: '20px 0' }}
              />
            )}
          </Card>
        </Col>

        {/* 5. 计算参数（100%） */}
        <Col xs={24}>
          <Card
            className="dashboard-card"
            style={dashboardCardStyle}
            title={
              <Space size={8}>
                <SettingOutlined style={{ color: DASHBOARD_STYLES.titleColor }} />
                <span style={{ fontSize: 14, fontWeight: 600, color: DASHBOARD_STYLES.titleColor }}>
                  计算参数
                </span>
              </Space>
            }
          >
            <Descriptions column={4} size="small" bordered>
              <Descriptions.Item label="盒子大小 (Å)">
                {electrolyte.box_size != null ? Number(electrolyte.box_size).toFixed(1) : '-'}
              </Descriptions.Item>
              <Descriptions.Item label="NPT 步数">
                {(job.config?.nsteps_npt || electrolyte.nsteps_npt)?.toLocaleString()}
              </Descriptions.Item>
              <Descriptions.Item label="NVT 步数">
                {(job.config?.nsteps_nvt || electrolyte.nsteps_nvt)?.toLocaleString()}
              </Descriptions.Item>
              <Descriptions.Item label="时间步长 (fs)">
                {job.config?.timestep || electrolyte.timestep}
              </Descriptions.Item>
              <Descriptions.Item label="力场" span={4}>
                {electrolyte.force_field || 'OPLS-AA'}
              </Descriptions.Item>
            </Descriptions>
          </Card>
        </Col>

        {/* 6. QC计算配置（仅当启用QC时显示） */}
        {job.config?.qc_enabled && (
          <Col xs={24}>
            <Card
              className="dashboard-card"
              style={{
                ...dashboardCardStyle,
                borderLeft: '4px solid #722ed1',
              }}
              title={
                <Space size={8}>
                  <ThunderboltOutlined style={{ color: '#722ed1' }} />
                  <span style={{ fontSize: 14, fontWeight: 600, color: DASHBOARD_STYLES.titleColor }}>
                    量子化学计算配置
                  </span>
                  <Tag color="purple">QC</Tag>
                </Space>
              }
            >
              {/* 全局配置 */}
              <Descriptions column={4} size="small" bordered>
                <Descriptions.Item label="精度等级">
                  <Tag color={
                    job.config.qc_accuracy_level === 'fast' ? 'green' :
                    job.config.qc_accuracy_level === 'standard' ? 'blue' :
                    job.config.qc_accuracy_level === 'accurate' ? 'orange' : 'purple'
                  }>
                    {job.config.qc_accuracy_level === 'fast' ? '快速' :
                     job.config.qc_accuracy_level === 'standard' ? '标准' :
                     job.config.qc_accuracy_level === 'accurate' ? '精确' : '自定义'}
                  </Tag>
                </Descriptions.Item>
                <Descriptions.Item label="泛函">
                  <Text code>{job.config.qc_functional || 'B3LYP'}</Text>
                </Descriptions.Item>
                <Descriptions.Item label="基组">
                  <Text code>{job.config.qc_basis_set || '6-31++G(d,p)'}</Text>
                </Descriptions.Item>
                <Descriptions.Item label="溶剂模型">
                  <Tag color={
                    job.config.qc_solvent_model === 'gas' ? 'default' :
                    job.config.qc_solvent_model === 'pcm' ? 'blue' : 'cyan'
                  }>
                    {job.config.qc_solvent_model === 'gas' ? '气相' :
                     job.config.qc_solvent_model === 'pcm' ? 'PCM' :
                     job.config.qc_solvent_model === 'smd' ? 'SMD' : job.config.qc_solvent_model}
                  </Tag>
                </Descriptions.Item>
                {job.config.qc_solvent_model !== 'gas' && job.config.qc_solvent_name && (
                  <Descriptions.Item label="隐式溶剂" span={2}>
                    <Text code>{job.config.qc_solvent_name}</Text>
                  </Descriptions.Item>
                )}
                <Descriptions.Item label="智能推荐" span={job.config.qc_solvent_model === 'gas' || !job.config.qc_solvent_name ? 3 : 1}>
                  <Tag color={job.config.qc_use_recommended_params !== false ? 'green' : 'default'}>
                    {job.config.qc_use_recommended_params !== false ? '已启用' : '未启用'}
                  </Tag>
                </Descriptions.Item>
              </Descriptions>

              {/* 分子详情列表 */}
              <div style={{ marginTop: 16 }}>
                <Text strong style={{ fontSize: 13, color: '#374151' }}>
                  📋 待计算分子列表 ({
                    (electrolyte.cations?.length || 0) +
                    (electrolyte.anions?.length || 0) +
                    (electrolyte.solvents?.length || 0)
                  } 个分子)
                </Text>
                <div style={{ marginTop: 8, display: 'flex', flexDirection: 'column', gap: 8 }}>
                  {/* 阳离子 */}
                  {electrolyte.cations?.map((mol, idx) => {
                    const charge = mol.smiles?.includes('+') ? 1 : 0;
                    const recommended = job.config?.qc_use_recommended_params !== false;
                    const recFunctional = recommended ? 'B3LYP' : (job.config?.qc_functional || 'B3LYP');
                    const recBasisSet = recommended ? '6-31+G(d,p)' : (job.config?.qc_basis_set || '6-31++G(d,p)');
                    return (
                      <div key={`cation-${idx}`} style={{
                        padding: '8px 12px',
                        background: '#fef2f2',
                        borderRadius: 6,
                        border: '1px solid #fecaca'
                      }}>
                        <Space size={12} wrap>
                          <Tag color="red">阳离子</Tag>
                          <Text strong>{mol.name}</Text>
                          <Text type="secondary" style={{ fontSize: 11 }}>SMILES: {mol.smiles}</Text>
                          <span style={{ fontSize: 11, color: '#666' }}>
                            泛函: <Text code style={{ fontSize: 10 }}>{recFunctional}</Text>
                          </span>
                          <span style={{ fontSize: 11, color: '#666' }}>
                            基组: <Text code style={{ fontSize: 10 }}>{recBasisSet}</Text>
                          </span>
                          <span style={{ fontSize: 11, color: '#666' }}>
                            电荷: <Text code style={{ fontSize: 10 }}>{charge}</Text>
                          </span>
                          {recommended && (
                            <Tag color="cyan" style={{ fontSize: 10 }}>智能推荐</Tag>
                          )}
                        </Space>
                      </div>
                    );
                  })}

                  {/* 阴离子 */}
                  {electrolyte.anions?.map((mol, idx) => {
                    const charge = mol.smiles?.includes('-') ? -1 : 0;
                    const recommended = job.config?.qc_use_recommended_params !== false;
                    const recFunctional = recommended ? 'B3LYP' : (job.config?.qc_functional || 'B3LYP');
                    const recBasisSet = recommended ? '6-31++G(d,p)' : (job.config?.qc_basis_set || '6-31++G(d,p)');
                    return (
                      <div key={`anion-${idx}`} style={{
                        padding: '8px 12px',
                        background: '#eff6ff',
                        borderRadius: 6,
                        border: '1px solid #bfdbfe'
                      }}>
                        <Space size={12} wrap>
                          <Tag color="blue">阴离子</Tag>
                          <Text strong>{mol.name}</Text>
                          <Text type="secondary" style={{ fontSize: 11 }}>SMILES: {mol.smiles}</Text>
                          <span style={{ fontSize: 11, color: '#666' }}>
                            泛函: <Text code style={{ fontSize: 10 }}>{recFunctional}</Text>
                          </span>
                          <span style={{ fontSize: 11, color: '#666' }}>
                            基组: <Text code style={{ fontSize: 10 }}>{recBasisSet}</Text>
                          </span>
                          <span style={{ fontSize: 11, color: '#666' }}>
                            电荷: <Text code style={{ fontSize: 10 }}>{charge}</Text>
                          </span>
                          {recommended && (
                            <Tag color="cyan" style={{ fontSize: 10 }}>智能推荐: 使用弥散函数(++)</Tag>
                          )}
                        </Space>
                      </div>
                    );
                  })}

                  {/* 溶剂 */}
                  {electrolyte.solvents?.map((mol, idx) => {
                    const recommended = job.config?.qc_use_recommended_params !== false;
                    const recFunctional = job.config?.qc_functional || 'B3LYP';
                    const recBasisSet = job.config?.qc_basis_set || '6-31++G(d,p)';
                    return (
                      <div key={`solvent-${idx}`} style={{
                        padding: '8px 12px',
                        background: '#f0fdf4',
                        borderRadius: 6,
                        border: '1px solid #bbf7d0'
                      }}>
                        <Space size={12} wrap>
                          <Tag color="green">溶剂</Tag>
                          <Text strong>{mol.name}</Text>
                          <Text type="secondary" style={{ fontSize: 11 }}>SMILES: {mol.smiles}</Text>
                          <span style={{ fontSize: 11, color: '#666' }}>
                            泛函: <Text code style={{ fontSize: 10 }}>{recFunctional}</Text>
                          </span>
                          <span style={{ fontSize: 11, color: '#666' }}>
                            基组: <Text code style={{ fontSize: 10 }}>{recBasisSet}</Text>
                          </span>
                          <span style={{ fontSize: 11, color: '#666' }}>
                            电荷: <Text code style={{ fontSize: 10 }}>0</Text>
                          </span>
                        </Space>
                      </div>
                    );
                  })}
                </div>
              </div>

              <div style={{ marginTop: 12 }}>
                <Text type="secondary" style={{ fontSize: 12 }}>
                  💡 智能推荐：阴离子使用弥散函数(++)描述扩展电子密度，阳离子使用极化函数(d,p)描述极化效应。
                </Text>
              </div>
            </Card>
          </Col>
        )}
      </Row>
    </div>
  );
}


