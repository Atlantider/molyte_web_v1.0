/**
 * 原子映射查看组件
 */
import { useState, useEffect } from 'react';
import {
  Card,
  Table,
  message,
  Spin,
  Empty,
  Tag,
  Space,
  Typography,
  Statistic,
  Row,
  Col,
} from 'antd';
import { DatabaseOutlined, ExperimentOutlined } from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';

const { Title, Text } = Typography;

interface Atom {
  atom_id: number;
  atom_index: number;
  atom_name: string;
  label: string;
  element: string;
  type: string;
}

interface Molecule {
  molecule_id: number;
  molecule_name: string;
  molecule_type: string;
  atoms: Atom[];
}

interface AtomMapping {
  molecules: Molecule[];
}

interface AtomMappingViewerProps {
  jobId: number;
}

export default function AtomMappingViewer({ jobId }: AtomMappingViewerProps) {
  const [loading, setLoading] = useState(false);
  const [mapping, setMapping] = useState<AtomMapping | null>(null);

  // 加载原子映射
  const loadMapping = async () => {
    setLoading(true);
    try {
      const response = await fetch(`/api/v1/jobs/${jobId}/atom_mapping`, {
        headers: {
          'Authorization': `Bearer ${localStorage.getItem('access_token')}`,
        },
      });

      if (!response.ok) {
        throw new Error('Failed to load atom mapping');
      }

      const data = await response.json();
      setMapping(data);
    } catch (error: any) {
      message.error('加载原子映射失败: ' + error.message);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    loadMapping();
  }, [jobId]);

  if (loading) {
    return (
      <Card>
        <div style={{ textAlign: 'center', padding: '50px 0' }}>
          <Spin size="large" tip="加载原子映射..." />
        </div>
      </Card>
    );
  }

  if (!mapping) {
    return (
      <Card>
        <Empty description="未找到原子映射文件" />
      </Card>
    );
  }

  // 统计信息
  const totalMolecules = mapping.molecules.length;
  const totalAtoms = mapping.molecules.reduce((sum, mol) => sum + mol.atoms.length, 0);
  const moleculeTypes = {
    cation: mapping.molecules.filter((m) => m.molecule_type === 'cation').length,
    anion: mapping.molecules.filter((m) => m.molecule_type === 'anion').length,
    solvent: mapping.molecules.filter((m) => m.molecule_type === 'solvent').length,
  };

  // 按分子类型分组
  const groupedMolecules = {
    cation: mapping.molecules.filter((m) => m.molecule_type === 'cation'),
    anion: mapping.molecules.filter((m) => m.molecule_type === 'anion'),
    solvent: mapping.molecules.filter((m) => m.molecule_type === 'solvent'),
  };

  // 表格列定义
  const columns: ColumnsType<Atom> = [
    {
      title: 'Atom ID',
      dataIndex: 'atom_id',
      key: 'atom_id',
      width: 100,
    },
    {
      title: 'Index',
      dataIndex: 'atom_index',
      key: 'atom_index',
      width: 80,
    },
    {
      title: 'Name',
      dataIndex: 'atom_name',
      key: 'atom_name',
      width: 100,
    },
    {
      title: 'Label',
      dataIndex: 'label',
      key: 'label',
      render: (label: string) => <Tag color="blue">{label}</Tag>,
    },
    {
      title: 'Element',
      dataIndex: 'element',
      key: 'element',
      width: 100,
      render: (element: string) => <Tag color="green">{element}</Tag>,
    },
    {
      title: 'Type',
      dataIndex: 'type',
      key: 'type',
      render: (type: string) => <Tag>{type}</Tag>,
    },
  ];

  return (
    <Space direction="vertical" style={{ width: '100%' }} size="large">
      {/* 统计信息 */}
      <Card>
        <Row gutter={16}>
          <Col span={6}>
            <Statistic title="总分子数" value={totalMolecules} prefix={<ExperimentOutlined />} />
          </Col>
          <Col span={6}>
            <Statistic title="总原子数" value={totalAtoms} prefix={<DatabaseOutlined />} />
          </Col>
          <Col span={4}>
            <Statistic title="阳离子" value={moleculeTypes.cation} />
          </Col>
          <Col span={4}>
            <Statistic title="阴离子" value={moleculeTypes.anion} />
          </Col>
          <Col span={4}>
            <Statistic title="溶剂" value={moleculeTypes.solvent} />
          </Col>
        </Row>
      </Card>

      {/* 阳离子 */}
      {groupedMolecules.cation.length > 0 && (
        <Card title={<><Tag color="red">阳离子</Tag> ({groupedMolecules.cation.length} 个分子)</>}>
          {groupedMolecules.cation.map((molecule) => (
            <Card
              key={molecule.molecule_id}
              type="inner"
              title={`${molecule.molecule_name} #${molecule.molecule_id}`}
              size="small"
              style={{ marginBottom: 16 }}
            >
              <Table
                columns={columns}
                dataSource={molecule.atoms}
                rowKey="atom_id"
                pagination={false}
                size="small"
              />
            </Card>
          ))}
        </Card>
      )}

      {/* 阴离子 */}
      {groupedMolecules.anion.length > 0 && (
        <Card title={<><Tag color="blue">阴离子</Tag> ({groupedMolecules.anion.length} 个分子)</>}>
          {groupedMolecules.anion.map((molecule) => (
            <Card
              key={molecule.molecule_id}
              type="inner"
              title={`${molecule.molecule_name} #${molecule.molecule_id}`}
              size="small"
              style={{ marginBottom: 16 }}
            >
              <Table
                columns={columns}
                dataSource={molecule.atoms}
                rowKey="atom_id"
                pagination={false}
                size="small"
              />
            </Card>
          ))}
        </Card>
      )}

      {/* 溶剂 */}
      {groupedMolecules.solvent.length > 0 && (
        <Card title={<><Tag color="green">溶剂</Tag> ({groupedMolecules.solvent.length} 个分子)</>}>
          {groupedMolecules.solvent.slice(0, 5).map((molecule) => (
            <Card
              key={molecule.molecule_id}
              type="inner"
              title={`${molecule.molecule_name} #${molecule.molecule_id}`}
              size="small"
              style={{ marginBottom: 16 }}
            >
              <Table
                columns={columns}
                dataSource={molecule.atoms}
                rowKey="atom_id"
                pagination={false}
                size="small"
              />
            </Card>
          ))}
          {groupedMolecules.solvent.length > 5 && (
            <Text type="secondary">
              ... 还有 {groupedMolecules.solvent.length - 5} 个溶剂分子（仅显示前 5 个）
            </Text>
          )}
        </Card>
      )}
    </Space>
  );
}

