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

  // 优先使用任务创建时的配方快照，如果没有则使用当前配方数据
  // 这样可以避免修改配方后影响历史任务的显示
  const systemSnapshot = job.config?.system_snapshot;
  const displayElectrolyte = systemSnapshot ? {
    ...electrolyte,
    name: systemSnapshot.name || electrolyte.name,
    cations: systemSnapshot.cations || electrolyte.cations,
    anions: systemSnapshot.anions || electrolyte.anions,
    solvents: systemSnapshot.solvents || electrolyte.solvents,
    box_size: systemSnapshot.box_size !== undefined ? systemSnapshot.box_size : electrolyte.box_size,
    temperature: systemSnapshot.temperature || electrolyte.temperature,
    pressure: systemSnapshot.pressure || electrolyte.pressure,
    force_field: systemSnapshot.force_field || electrolyte.force_field,
  } : electrolyte;

  // 调试信息
  useEffect(() => {
    console.log('=== JobBasicInfo Debug ===');
    console.log('job:', job);
    console.log('job.started_at:', job.started_at);
    console.log('job.finished_at:', job.finished_at);
    console.log('slurmStatus:', slurmStatus);
    console.log('electrolyte (current):', electrolyte);
    console.log('systemSnapshot:', systemSnapshot);
    console.log('displayElectrolyte (used for display):', displayElectrolyte);
  }, [job, slurmStatus, electrolyte, systemSnapshot, displayElectrolyte]);

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
                  {job.config?.user_note && (
                    <Descriptions.Item label="备注">
                      <Text type="secondary">{job.config.user_note}</Text>
                    </Descriptions.Item>
                  )}
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

        {/* 2. 溶液配方（100%宽度） */}
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
            <Descriptions column={2} size="small" bordered>
              <Descriptions.Item label="配方名称">
                <Text strong>{electrolyte.name}</Text>
              </Descriptions.Item>
              <Descriptions.Item label="温度 (K)">
                <Text>{job.config?.temperature || electrolyte.temperature}</Text>
              </Descriptions.Item>
              <Descriptions.Item label="压力 (atm)">
                <Text>{job.config?.pressure || electrolyte.pressure}</Text>
              </Descriptions.Item>
              <Descriptions.Item label="盒子尺寸 (Å)">
                <Text>{electrolyte.box_size ? Number(electrolyte.box_size).toFixed(2) : '-'}</Text>
              </Descriptions.Item>
              <Descriptions.Item label="阳离子" span={2}>
                <Space size={[8, 8]} wrap>
                  {electrolyte.cations.map((cation, idx) => (
                    <Tag key={`cation-${idx}`} color="red">
                      {cation.name} × {cation.number}
                    </Tag>
                  ))}
                </Space>
              </Descriptions.Item>
              <Descriptions.Item label="阴离子" span={2}>
                <Space size={[8, 8]} wrap>
                  {electrolyte.anions.map((anion, idx) => (
                    <Tag key={`anion-${idx}`} color="blue">
                      {anion.name} × {anion.number}
                    </Tag>
                  ))}
                </Space>
              </Descriptions.Item>
              {electrolyte.solvents && electrolyte.solvents.length > 0 && (
                <Descriptions.Item label="溶剂" span={2}>
                  <Space size={[8, 8]} wrap>
                    {electrolyte.solvents.map((solvent, idx) => (
                      <Tag key={`solvent-${idx}`} color="green">
                        {solvent.name} × {solvent.number}
                      </Tag>
                    ))}
                  </Space>
                </Descriptions.Item>
              )}
            </Descriptions>
          </Card>
        </Col>

        {/* 3. 计算结果对比（100%宽度） */}
        <Col xs={24}>
          <Card
            className="dashboard-card"
            style={dashboardCardStyle}
            title={
              <Space size={8}>
                <FundOutlined style={{ color: DASHBOARD_STYLES.titleColor }} />
                <span style={{ fontSize: 14, fontWeight: 600, color: DASHBOARD_STYLES.titleColor }}>
                  计算结果对比
                </span>
              </Space>
            }
          >
            {loadingStructure ? (
              <div style={{ textAlign: 'center', padding: '40px 0' }}>
                <Spin />
              </div>
            ) : structureInfo?.available ? (
              <Row gutter={16}>
                {/* 左侧：初始值 */}
                <Col xs={24} lg={12}>
                  <Descriptions
                    column={1}
                    size="small"
                    bordered
                    title={<Text strong style={{ fontSize: 13 }}>初始设置</Text>}
                  >
                    <Descriptions.Item label="初始浓度 (mol/L)">
                      <Text strong style={{ fontSize: 14 }}>
                        {structureInfo.initial_concentration?.toFixed(4) || '-'}
                      </Text>
                    </Descriptions.Item>
                    <Descriptions.Item label="初始密度 (g/cm³)">
                      <Text strong style={{ fontSize: 14 }}>
                        {structureInfo.initial_density?.toFixed(4) || '-'}
                      </Text>
                    </Descriptions.Item>
                    <Descriptions.Item label="初始盒子尺寸 (Å)">
                      <Text code style={{ fontSize: 12 }}>
                        {structureInfo.initial_box_dimensions || '-'}
                      </Text>
                    </Descriptions.Item>
                  </Descriptions>
                </Col>

                {/* 右侧：计算结果 */}
                <Col xs={24} lg={12}>
                  <Descriptions
                    column={1}
                    size="small"
                    bordered
                    title={<Text strong style={{ fontSize: 13, color: '#1890ff' }}>计算结果</Text>}
                  >
                    <Descriptions.Item label="计算浓度 (mol/L)">
                      <Space direction="vertical" size={0}>
                        <Text strong style={{ fontSize: 16, color: '#52c41a' }}>
                          {structureInfo.concentration?.toFixed(4) || '-'}
                        </Text>
                        {structureInfo.initial_concentration && structureInfo.concentration && (
                          <Text type="secondary" style={{ fontSize: 12 }}>
                            偏差: {((structureInfo.concentration - structureInfo.initial_concentration) / structureInfo.initial_concentration * 100).toFixed(2)}%
                          </Text>
                        )}
                      </Space>
                    </Descriptions.Item>
                    <Descriptions.Item label="计算密度 (g/cm³)">
                      <Space direction="vertical" size={0}>
                        <Text strong style={{ fontSize: 16, color: '#1890ff' }}>
                          {structureInfo.density?.toFixed(4) || '-'}
                        </Text>
                        {structureInfo.initial_density && structureInfo.density && (
                          <Text type="secondary" style={{ fontSize: 12 }}>
                            偏差: {((structureInfo.density - structureInfo.initial_density) / structureInfo.initial_density * 100).toFixed(2)}%
                          </Text>
                        )}
                      </Space>
                    </Descriptions.Item>
                    <Descriptions.Item label="最终盒子尺寸 (Å)">
                      <Text code strong style={{ fontSize: 12 }}>
                        {structureInfo.box_dimensions || '-'}
                      </Text>
                    </Descriptions.Item>
                  </Descriptions>
                </Col>
              </Row>
            ) : (
              <Alert
                message="计算结果未就绪"
                description={
                  <div>
                    <p>任务完成后将显示计算结果</p>
                    <p style={{ marginTop: 8, fontSize: 12 }}>
                      当前状态: <Tag>{job.status}</Tag>
                    </p>
                  </div>
                }
                type="info"
                showIcon
              />
            )}
          </Card>
        </Col>

        {/* 4. 计算参数（100%宽度） */}
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
            <Row gutter={16}>
              {/* 左侧：模拟参数 */}
              <Col xs={24} lg={12}>
                <Descriptions column={1} size="small" bordered title={<Text strong style={{ fontSize: 13 }}>模拟参数</Text>}>
                  <Descriptions.Item label="NPT 步数">
                    <Text>{(job.config?.nsteps_npt || electrolyte.nsteps_npt)?.toLocaleString() || '-'}</Text>
                  </Descriptions.Item>
                  <Descriptions.Item label="NVT 步数">
                    <Text>{(job.config?.nsteps_nvt || electrolyte.nsteps_nvt)?.toLocaleString() || '-'}</Text>
                  </Descriptions.Item>
                  <Descriptions.Item label="时间步长 (fs)">
                    <Text>{job.config?.timestep || electrolyte.timestep || '-'}</Text>
                  </Descriptions.Item>
                  <Descriptions.Item label="总模拟时间 (ps)">
                    <Text strong style={{ color: '#1890ff' }}>{getTotalSimulationTime()}</Text>
                  </Descriptions.Item>
                </Descriptions>
              </Col>

              {/* 右侧：系统参数 */}
              <Col xs={24} lg={12}>
                <Descriptions column={1} size="small" bordered title={<Text strong style={{ fontSize: 13 }}>系统参数</Text>}>
                  <Descriptions.Item label="力场">
                    <Tag color="purple">{electrolyte.force_field || 'OPLS-AA'}</Tag>
                  </Descriptions.Item>
                  <Descriptions.Item label="温度 (K)">
                    <Text>{job.config?.temperature || electrolyte.temperature || '-'}</Text>
                  </Descriptions.Item>
                  <Descriptions.Item label="压力 (atm)">
                    <Text>{job.config?.pressure || electrolyte.pressure || '-'}</Text>
                  </Descriptions.Item>
                  <Descriptions.Item label="初始盒子尺寸 (Å)">
                    <Text>{electrolyte.box_size ? Number(electrolyte.box_size).toFixed(2) : '-'}</Text>
                  </Descriptions.Item>
                </Descriptions>
              </Col>
            </Row>
          </Card>
        </Col>

        {/* 5. QC计算配置（仅当启用QC时显示） */}
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
              <Descriptions column={2} size="small" bordered>
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
                <Descriptions.Item label="智能推荐">
                  <Tag color={job.config.qc_use_recommended_params !== false ? 'green' : 'default'}>
                    {job.config.qc_use_recommended_params !== false ? '已启用' : '未启用'}
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
                  <Descriptions.Item label="隐式溶剂">
                    <Text code>{job.config.qc_solvent_name}</Text>
                  </Descriptions.Item>
                )}
                <Descriptions.Item label="待计算分子" span={2}>
                  <Space size={[8, 8]} wrap>
                    {electrolyte.cations?.map((mol, idx) => (
                      <Tag key={`qc-cation-${idx}`} color="red">
                        {mol.name}
                      </Tag>
                    ))}
                    {electrolyte.anions?.map((mol, idx) => (
                      <Tag key={`qc-anion-${idx}`} color="blue">
                        {mol.name}
                      </Tag>
                    ))}
                    {electrolyte.solvents?.map((mol, idx) => (
                      <Tag key={`qc-solvent-${idx}`} color="green">
                        {mol.name}
                      </Tag>
                    ))}
                  </Space>
                </Descriptions.Item>
              </Descriptions>

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


