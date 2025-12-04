/**
 * Desolvation energy result view component
 * 去溶剂化能结果展示组件
 */
import React from 'react';
import { Card, Descriptions, Table, Divider, Typography, theme, Statistic, Row, Col } from 'antd';
import { ThunderboltOutlined } from '@ant-design/icons';
import type { DesolvationEnergyResult } from '../types/desolvation';
import { useThemeStore } from '../stores/themeStore';

const { Text } = Typography;

interface DesolvationResultViewProps {
  result: DesolvationEnergyResult;
  mode?: 'stepwise' | 'full';
}

export default function DesolvationResultView({ result, mode = 'stepwise' }: DesolvationResultViewProps) {
  const { token } = theme.useToken();
  const { isDark } = useThemeStore();

  // 检测是否为 full 模式（通过结果数据判断）
  const isFullMode = result.per_ligand_results.length === 1 &&
                     result.per_ligand_results[0]?.ligand_type === 'TOTAL';

  return (
    <Card
      title="去溶剂化能结果"
      style={{
        marginTop: 16,
        background: isDark ? token.colorBgContainer : undefined,
        borderColor: token.colorBorder,
      }}
    >
      {/* 总览信息 */}
      <Descriptions
        bordered
        column={3}
        size="small"
        style={{
          background: isDark ? 'transparent' : undefined,
        }}
        labelStyle={{
          background: isDark ? 'rgba(255,255,255,0.04)' : undefined,
          color: token.colorText,
        }}
        contentStyle={{
          background: isDark ? 'transparent' : undefined,
          color: token.colorText,
        }}
      >
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

      {isFullMode ? (
        // Full 模式：显示总去溶剂化能
        <>
          <Divider orientation="left">
            <Text style={{ color: token.colorText }}>总去溶剂化能</Text>
          </Divider>
          <Row gutter={16} style={{ marginBottom: 16 }}>
            <Col span={8}>
              <Card
                size="small"
                style={{
                  background: isDark ? 'rgba(255,255,255,0.04)' : '#f6ffed',
                  borderColor: isDark ? token.colorBorder : '#b7eb8f',
                }}
              >
                <Statistic
                  title={<Text style={{ color: token.colorTextSecondary }}>总去溶剂化能</Text>}
                  value={result.per_ligand_results[0]?.delta_e || 0}
                  precision={2}
                  suffix="kcal/mol"
                  prefix={<ThunderboltOutlined />}
                  valueStyle={{
                    color: (result.per_ligand_results[0]?.delta_e || 0) < 0 ? '#52c41a' : '#1890ff',
                    fontSize: 24,
                  }}
                />
              </Card>
            </Col>
            <Col span={8}>
              <Card
                size="small"
                style={{
                  background: isDark ? 'rgba(255,255,255,0.04)' : '#fff7e6',
                  borderColor: isDark ? token.colorBorder : '#ffd591',
                }}
              >
                <Statistic
                  title={<Text style={{ color: token.colorTextSecondary }}>中心离子能量</Text>}
                  value={result.per_ligand_results[0]?.e_ligand || 0}
                  precision={6}
                  suffix="A.U."
                  valueStyle={{ color: token.colorText, fontSize: 16 }}
                />
              </Card>
            </Col>
            <Col span={8}>
              <Card
                size="small"
                style={{
                  background: isDark ? 'rgba(255,255,255,0.04)' : '#e6f7ff',
                  borderColor: isDark ? token.colorBorder : '#91d5ff',
                }}
              >
                <Statistic
                  title={<Text style={{ color: token.colorTextSecondary }}>簇能量</Text>}
                  value={result.e_cluster}
                  precision={6}
                  suffix="A.U."
                  valueStyle={{ color: token.colorText, fontSize: 16 }}
                />
              </Card>
            </Col>
          </Row>
        </>
      ) : (
        // Stepwise 模式：显示每个配体的去溶剂化能
        <>
          <Divider orientation="left">
            <Text style={{ color: token.colorText }}>按配体展示</Text>
          </Divider>
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
          <Divider orientation="left">
            <Text style={{ color: token.colorText }}>按类型汇总</Text>
          </Divider>
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
        </>
      )}
    </Card>
  );
}

