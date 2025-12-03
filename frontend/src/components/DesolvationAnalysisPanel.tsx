/**
 * Desolvation energy analysis panel component
 * 去溶剂化能分析面板组件
 */
import React, { useState, useEffect } from 'react';
import { Card, Button, Table, message, Tag, Space, Typography } from 'antd';
import { ThunderboltOutlined, ReloadOutlined } from '@ant-design/icons';
import {
  createDesolvationJob,
  getDesolvationJob,
  listClusterDesolvationJobs
} from '../api/desolvation';
import type { DesolvationJobResponse } from '../types/desolvation';
import DesolvationResultView from './DesolvationResultView';

const { Text } = Typography;

interface DesolvationAnalysisPanelProps {
  jobId: number;
  clusterId?: number;
}

// 任务状态标签
const JobStatusTag: React.FC<{ status: string }> = ({ status }) => {
  const statusConfig: Record<string, { color: string; text: string }> = {
    CREATED: { color: 'default', text: '已创建' },
    SUBMITTED: { color: 'blue', text: '已提交' },
    QUEUED: { color: 'cyan', text: '排队中' },
    RUNNING: { color: 'processing', text: '运行中' },
    POSTPROCESSING: { color: 'purple', text: '后处理' },
    COMPLETED: { color: 'success', text: '已完成' },
    FAILED: { color: 'error', text: '失败' },
    CANCELLED: { color: 'default', text: '已取消' },
  };

  const config = statusConfig[status] || { color: 'default', text: status };
  return <Tag color={config.color}>{config.text}</Tag>;
};

export default function DesolvationAnalysisPanel({
  jobId,
  clusterId
}: DesolvationAnalysisPanelProps) {
  const [jobs, setJobs] = useState<DesolvationJobResponse[]>([]);
  const [selectedJob, setSelectedJob] = useState<DesolvationJobResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [creating, setCreating] = useState(false);

  // 加载任务列表
  const loadJobs = async () => {
    if (!clusterId) return;

    setLoading(true);
    try {
      const data = await listClusterDesolvationJobs(clusterId);
      setJobs(data);
    } catch (error: any) {
      message.error(`加载任务列表失败: ${error.message || '未知错误'}`);
    } finally {
      setLoading(false);
    }
  };

  // 初始加载
  useEffect(() => {
    loadJobs();
  }, [clusterId]);

  // 创建新任务
  const handleCreateJob = async () => {
    if (!clusterId) {
      message.warning('请先选择一个溶剂化结构');
      return;
    }

    setCreating(true);
    try {
      const newJob = await createDesolvationJob({
        md_job_id: jobId,
        solvation_structure_id: clusterId,
        method_level: 'standard'  // 使用 QC 标准方法（B3LYP/6-31++G(d,p)）
      });
      message.success('去溶剂化能任务已创建');
      await loadJobs();
    } catch (error: any) {
      message.error(`创建任务失败: ${error.response?.data?.detail || error.message || '未知错误'}`);
    } finally {
      setCreating(false);
    }
  };

  // 查看结果
  const handleViewResult = async (job: DesolvationJobResponse) => {
    try {
      const detail = await getDesolvationJob(job.job_id);
      setSelectedJob(detail);
    } catch (error: any) {
      message.error(`获取任务详情失败: ${error.message || '未知错误'}`);
    }
  };

  return (
    <div>
      {/* 顶部：创建按钮 */}
      <Card size="small">
        <Space>
          <Button
            type="primary"
            icon={<ThunderboltOutlined />}
            onClick={handleCreateJob}
            loading={creating}
            disabled={!clusterId}
          >
            计算去溶剂化能
          </Button>
          <Button
            icon={<ReloadOutlined />}
            onClick={loadJobs}
            loading={loading}
            disabled={!clusterId}
          >
            刷新
          </Button>
          {!clusterId && (
            <Text type="secondary">请先在左侧选择一个溶剂化结构</Text>
          )}
        </Space>
      </Card>

      {/* 任务列表 */}
      <Card title="任务列表" style={{ marginTop: 16 }} size="small">
        <Table
          dataSource={jobs}
          rowKey="job_id"
          loading={loading}
          columns={[
            {
              title: 'ID',
              dataIndex: 'job_id',
              key: 'job_id',
              width: 80,
            },
            {
              title: '方法',
              dataIndex: 'method_level',
              key: 'method_level',
              width: 120,
            },
            {
              title: '状态',
              dataIndex: 'status',
              key: 'status',
              width: 100,
              render: (status: string) => <JobStatusTag status={status} />
            },
            {
              title: '创建时间',
              dataIndex: 'created_at',
              key: 'created_at',
              width: 180,
              render: (time: string) => new Date(time).toLocaleString('zh-CN')
            },
            {
              title: '机时 (秒)',
              dataIndex: 'elapsed_seconds',
              key: 'elapsed_seconds',
              width: 100,
              render: (seconds?: number) => seconds ? seconds.toFixed(1) : '-'
            },
            {
              title: '操作',
              key: 'action',
              width: 120,
              render: (_, record: DesolvationJobResponse) => (
                <Button
                  type="link"
                  size="small"
                  disabled={record.status !== 'COMPLETED'}
                  onClick={() => handleViewResult(record)}
                >
                  查看结果
                </Button>
              )
            },
          ]}
          pagination={false}
          size="small"
        />
      </Card>

      {/* 结果展示 */}
      {selectedJob && selectedJob.result && (
        <DesolvationResultView result={selectedJob.result} />
      )}
    </div>
  );
}

