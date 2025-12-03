/**
 * Desolvation energy result view component
 * 去溶剂化能结果展示组件
 */
import React from 'react';
import { Card, Descriptions, Table, Divider, Typography } from 'antd';
import type { DesolvationEnergyResult } from '../types/desolvation';

const { Title } = Typography;

interface DesolvationResultViewProps {
  result: DesolvationEnergyResult;
}

export default function DesolvationResultView({ result }: DesolvationResultViewProps) {
  return (
    <Card title="去溶剂化能结果" style={{ marginTop: 16 }}>
      {/* 总览信息 */}
      <Descriptions bordered column={3} size="small">
        <Descriptions.Item label="计算方法">
          {result.method_level}
        </Descriptions.Item>
        <Descriptions.Item label="簇能量">
          {result.e_cluster.toFixed(6)} A.U.
        </Descriptions.Item>
        <Descriptions.Item label="基组/泛函">
          {result.functional || 'N/A'} / {result.basis_set || 'N/A'}
        </Descriptions.Item>
      </Descriptions>

      {/* 按配体展示 */}
      <Divider orientation="left">按配体展示</Divider>
      <Table
        dataSource={result.per_ligand_results}
        rowKey="ligand_id"
        columns={[
          {
            title: '配体',
            dataIndex: 'ligand_label',
            key: 'ligand_label',
            width: 120,
          },
          {
            title: '类型',
            dataIndex: 'ligand_type',
            key: 'ligand_type',
            width: 100,
          },
          {
            title: 'ΔE (kcal/mol)',
            dataIndex: 'delta_e',
            key: 'delta_e',
            width: 150,
            render: (val: number) => (
              <span style={{ fontWeight: 500, color: val < 0 ? '#52c41a' : '#1890ff' }}>
                {val.toFixed(2)}
              </span>
            ),
            sorter: (a, b) => a.delta_e - b.delta_e,
          },
          {
            title: 'E_ligand (A.U.)',
            dataIndex: 'e_ligand',
            key: 'e_ligand',
            render: (val: number) => val.toFixed(6),
          },
          {
            title: 'E_cluster_minus (A.U.)',
            dataIndex: 'e_cluster_minus',
            key: 'e_cluster_minus',
            render: (val: number) => val.toFixed(6),
          },
        ]}
        pagination={false}
        size="small"
      />

      {/* 按类型汇总 */}
      <Divider orientation="left">按类型汇总</Divider>
      <Table
        dataSource={result.per_type_summary}
        rowKey="ligand_type"
        columns={[
          {
            title: '配体类型',
            dataIndex: 'ligand_type',
            key: 'ligand_type',
            width: 120,
          },
          {
            title: '数量',
            dataIndex: 'count',
            key: 'count',
            width: 80,
          },
          {
            title: '平均 ΔE (kcal/mol)',
            dataIndex: 'avg_delta_e',
            key: 'avg_delta_e',
            width: 150,
            render: (val: number) => (
              <span style={{ fontWeight: 500 }}>
                {val.toFixed(2)}
              </span>
            ),
            sorter: (a, b) => a.avg_delta_e - b.avg_delta_e,
          },
          {
            title: '标准差',
            dataIndex: 'std_delta_e',
            key: 'std_delta_e',
            width: 100,
            render: (val: number) => val.toFixed(2),
          },
          {
            title: '范围 (kcal/mol)',
            key: 'range',
            render: (_, record) => (
              <span>
                {record.min_delta_e.toFixed(2)} ~ {record.max_delta_e.toFixed(2)}
              </span>
            ),
          },
        ]}
        pagination={false}
        size="small"
      />
    </Card>
  );
}

