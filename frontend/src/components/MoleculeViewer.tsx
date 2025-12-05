/**
 * 分子3D结构查看器 - Nature 期刊风格
 * 使用 3Dmol.js 显示分子结构和原子电荷
 */
import { useState, useEffect, useRef } from 'react';
import {
  Card,
  message,
  Spin,
  Empty,
  Tag,
  Space,
  Typography,
  Row,
  Col,
  Tabs,
  Table,
  Statistic,
  Button,
  Descriptions,
  Divider,
  Image,
  theme,
} from 'antd';
import { ExperimentOutlined, ReloadOutlined, EyeOutlined, ThunderboltOutlined } from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';
import type { MoleculeQCCache } from '../types/qc';
import { useThemeStore } from '../stores/themeStore';

const { Text } = Typography;

// 3Dmol.js 类型声明
declare global {
  interface Window {
    $3Dmol: any;
  }
}

interface Atom {
  id: number;
  name: string;
  element: string;
  x: number;
  y: number;
  z: number;
  charge: number | null;
}

interface Molecule {
  name: string;
  type: string;
  pdb_content: string;
  atoms: Atom[];
  total_charge: number;
  charge_method: string;
  smiles?: string;  // 用于查询QC属性
}

interface MoleculeViewerProps {
  jobId: number;
}

// Dashboard 样式常量（与其他页面保持一致）
const DASHBOARD_STYLES = {
  cardBorderRadius: 12,
  cardPadding: 24,
  gutter: 24,
  titleFontSize: 16,
  titleFontWeight: 600,
};

// 响应式CSS样式 - 支持暗色模式
const getResponsiveStyles = (isDark: boolean) => `
  .mol-stats-grid {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 20px;
  }
  .mol-list-grid {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 16px;
  }
  @media (max-width: 1400px) {
    .mol-stats-grid {
      grid-template-columns: repeat(2, 1fr);
    }
    .mol-list-grid {
      grid-template-columns: repeat(3, 1fr);
    }
  }
  @media (max-width: 768px) {
    .mol-stats-grid {
      grid-template-columns: 1fr;
    }
    .mol-list-grid {
      grid-template-columns: repeat(2, 1fr);
    }
  }
  .dashboard-card {
    transition: all 0.3s ease;
  }
  .dashboard-card:hover {
    box-shadow: ${isDark ? '0 8px 24px rgba(0, 0, 0, 0.4)' : '0 8px 24px rgba(15, 23, 42, 0.12)'};
    transform: translateY(-2px);
  }
  .molecule-card {
    border: 1px solid ${isDark ? '#333' : '#e8e8e8'} !important;
    box-shadow: ${isDark ? '0 2px 8px rgba(0, 0, 0, 0.2)' : '0 2px 8px rgba(0, 0, 0, 0.06)'} !important;
    transition: all 0.3s ease !important;
  }
  .molecule-card:hover {
    box-shadow: ${isDark ? '0 4px 16px rgba(0, 0, 0, 0.3)' : '0 4px 16px rgba(0, 0, 0, 0.12)'} !important;
    border-color: #1890ff !important;
  }
  .molecule-card.selected {
    border-color: #1890ff !important;
    border-width: 2px !important;
    box-shadow: 0 4px 16px rgba(24, 144, 255, 0.2) !important;
  }
  .ant-table-small .ant-table-thead > tr > th {
    background: ${isDark ? 'rgba(255, 255, 255, 0.04)' : '#f8fafc'} !important;
    font-weight: 600 !important;
    font-size: 11px !important;
    padding: 8px 8px !important;
    color: ${isDark ? 'rgba(255, 255, 255, 0.85)' : '#475569'} !important;
  }
  .ant-table-small .ant-table-tbody > tr > td {
    padding: 6px 8px !important;
    font-size: 12px !important;
  }
  .ant-table-small .ant-table-tbody > tr:hover > td {
    background: ${isDark ? 'rgba(24, 144, 255, 0.1)' : '#f0f9ff'} !important;
  }
`;

export default function MoleculeViewer({ jobId }: MoleculeViewerProps) {
  const { mode } = useThemeStore();
  const { token } = theme.useToken();
  const [loading, setLoading] = useState(false);
  const [molecules, setMolecules] = useState<Molecule[]>([]);
  const [selectedMolecule, setSelectedMolecule] = useState<string | null>(null);
  const viewerRef = useRef<HTMLDivElement>(null);
  const viewerInstance = useRef<any>(null);
  const [qcCache, setQcCache] = useState<Record<string, MoleculeQCCache>>({});
  const [loadingQC, setLoadingQC] = useState(false);
  const isDark = mode === 'dark';

  // 卡片样式
  const dashboardCardStyle: React.CSSProperties = {
    background: token.colorBgContainer,
    borderRadius: DASHBOARD_STYLES.cardBorderRadius,
    boxShadow: isDark ? '0 4px 12px rgba(0, 0, 0, 0.3)' : '0 4px 12px rgba(15, 23, 42, 0.08)',
    border: `1px solid ${token.colorBorder}`,
    transition: 'all 0.3s ease',
  };

  // 加载分子模板
  const loadMolecules = async () => {
    setLoading(true);
    try {
      const response = await fetch(`/api/v1/jobs/${jobId}/molecule_templates`, {
        headers: {
          'Authorization': `Bearer ${localStorage.getItem('access_token')}`,
        },
      });

      if (!response.ok) {
        throw new Error('Failed to load molecule templates');
      }

      const data = await response.json();
      console.log('Loaded molecules:', data.molecules);
      console.log('Molecules count:', data.molecules.length);

      // 检查每个分子的数据
      data.molecules.forEach((mol: any) => {
        console.log(`Molecule ${mol.name}:`, {
          type: mol.type,
          atoms: mol.atoms?.length,
          hasPDB: !!mol.pdb_content,
          pdbLength: mol.pdb_content?.length,
          totalCharge: mol.total_charge
        });
      });

      setMolecules(data.molecules);

      // 默认选择第一个分子
      if (data.molecules.length > 0) {
        setSelectedMolecule(data.molecules[0].name);
      }
      // 加载QC缓存数据
      loadQCCache(data.molecules);
    } catch (error: any) {
      console.error('Error loading molecules:', error);
      message.error('加载分子结构失败: ' + error.message);
    } finally {
      setLoading(false);
    }
  };

  // 加载分子的QC缓存数据
  const loadQCCache = async (mols: Molecule[]) => {
    setLoadingQC(true);
    const cache: Record<string, MoleculeQCCache> = {};

    for (const mol of mols) {
      if (mol.smiles) {
        try {
          const response = await fetch(`/api/v1/qc/cache/${encodeURIComponent(mol.smiles)}`, {
            headers: {
              'Authorization': `Bearer ${localStorage.getItem('access_token')}`,
            },
          });

          if (response.ok) {
            const data = await response.json();
            cache[mol.smiles] = data;
          }
        } catch (error) {
          console.log(`No QC cache for ${mol.name}`);
        }
      }
    }

    setQcCache(cache);
    setLoadingQC(false);
  };

  // 初始化 3Dmol.js
  useEffect(() => {
    // 动态加载 3Dmol.js
    if (!window.$3Dmol) {
      const script = document.createElement('script');
      script.src = 'https://3Dmol.csb.pitt.edu/build/3Dmol-min.js';
      script.async = true;
      script.onload = () => {
        loadMolecules();
      };
      document.body.appendChild(script);
    } else {
      loadMolecules();
    }
  }, [jobId]);

  // 渲染分子结构
  useEffect(() => {
    if (!selectedMolecule || !viewerRef.current || !window.$3Dmol || molecules.length === 0) {
      return;
    }

    const molecule = molecules.find(m => m.name === selectedMolecule);
    if (!molecule) return;

    try {
      // 清空容器并重新创建 viewer
      viewerRef.current.innerHTML = '';

      // 创建新的 viewer - 根据主题设置背景色
      const config = { backgroundColor: isDark ? '#1a1a1a' : 'white' };
      const viewer = window.$3Dmol.createViewer(viewerRef.current, config);
      viewerInstance.current = viewer;

      console.log('Creating viewer for molecule:', molecule.name);
      console.log('Atoms count:', molecule.atoms.length);
      console.log('PDB content length:', molecule.pdb_content?.length);

      // 加载 PDB 数据
      if (molecule.pdb_content) {
        viewer.addModel(molecule.pdb_content, 'pdb');
      } else {
        console.error('No PDB content for molecule:', molecule.name);
        return;
      }

      // 设置样式：球棍模型
      // 对于单原子分子，使用更大的球体
      if (molecule.atoms.length === 1) {
        viewer.setStyle({}, {
          sphere: { scale: 1.0, colorscheme: 'Jmol' }
        });
      } else {
        viewer.setStyle({}, {
          stick: { radius: 0.15, colorscheme: 'Jmol' },
          sphere: { scale: 0.3, colorscheme: 'Jmol' }
        });
      }

      // 添加原子标签（显示原子名称和电荷）
      molecule.atoms.forEach((atom) => {
        const label = atom.charge !== null
          ? `${atom.element}\n${atom.charge.toFixed(3)}`
          : atom.element;

        viewer.addLabel(label, {
          position: { x: atom.x, y: atom.y, z: atom.z },
          backgroundColor: 'rgba(255, 255, 255, 0.8)',
          fontColor: 'black',
          fontSize: 12,
          showBackground: true,
        });
      });

      viewer.zoomTo();
      viewer.render();

      console.log('Viewer rendered successfully');
      console.log('Viewer instance saved:', viewerInstance.current);
    } catch (error) {
      console.error('Error rendering molecule:', error);
      message.error('渲染分子结构失败: ' + error);
    }

    // 不需要清理函数，因为我们在每次渲染时都会清空容器
  }, [selectedMolecule, molecules]);

  // 重置视角
  const resetView = (e?: React.MouseEvent) => {
    e?.preventDefault();
    e?.stopPropagation();

    console.log('=== Reset View Clicked ===');
    console.log('Viewer instance exists:', !!viewerInstance.current);
    console.log('Viewer instance:', viewerInstance.current);

    if (!viewerInstance.current) {
      console.error('No viewer instance!');
      message.warning('3D查看器未初始化');
      return;
    }

    try {
      const viewer = viewerInstance.current;
      console.log('Calling viewer methods...');

      // 方法1: 重置视图矩阵
      viewer.setView([0, 0, 0, 0, 0, 0, 0, 1]);

      // 方法2: 重新缩放到适合大小
      viewer.zoomTo();

      // 方法3: 重新渲染
      viewer.render();

      console.log('View reset completed');
      message.success('视角已重置');
    } catch (error) {
      console.error('Error resetting view:', error);
      message.error('重置视角失败: ' + String(error));
    }
  };

  if (loading) {
    return (
      <div style={{
        background: token.colorBgLayout,
        padding: DASHBOARD_STYLES.gutter,
        minHeight: '100%',
        borderRadius: 8,
      }}>
        <Card className="dashboard-card" style={dashboardCardStyle}>
          <div style={{ textAlign: 'center', padding: 60 }}>
            <Spin size="large" tip="加载分子结构中..." />
          </div>
        </Card>
      </div>
    );
  }

  if (molecules.length === 0) {
    return (
      <div style={{
        background: token.colorBgLayout,
        padding: DASHBOARD_STYLES.gutter,
        minHeight: '100%',
        borderRadius: 8,
      }}>
        <Card className="dashboard-card" style={dashboardCardStyle}>
          <Empty
            description="未找到分子结构文件"
            image={Empty.PRESENTED_IMAGE_SIMPLE}
          >
            <Button
              icon={<ReloadOutlined />}
              onClick={loadMolecules}
              style={{ borderRadius: 6 }}
            >
              重新加载
            </Button>
          </Empty>
        </Card>
      </div>
    );
  }

  const currentMolecule = molecules.find(m => m.name === selectedMolecule);

  // 按类型分组
  const cations = molecules.filter(m => m.type === 'cation');
  const anions = molecules.filter(m => m.type === 'anion');
  const solvents = molecules.filter(m => m.type === 'solvent');

  // 原子表格列定义
  const atomColumns: ColumnsType<Atom> = [
    {
      title: 'ID',
      dataIndex: 'id',
      key: 'id',
      width: 50,
      align: 'center',
    },
    {
      title: '元素',
      dataIndex: 'element',
      key: 'element',
      width: 70,
      align: 'center',
      render: (element: string) => (
        <Tag color="green" style={{ margin: 0 }}>{element}</Tag>
      ),
    },
    {
      title: '电荷 (e)',
      dataIndex: 'charge',
      key: 'charge',
      width: 90,
      align: 'center',
      render: (charge: number | null) => {
        if (charge === null) return <Text type="secondary">N/A</Text>;
        const color = charge > 0 ? 'red' : charge < 0 ? 'blue' : 'default';
        return <Tag color={color} style={{ margin: 0, fontWeight: 600 }}>{charge.toFixed(3)}</Tag>;
      },
    },
    {
      title: 'X (Å)',
      dataIndex: 'x',
      key: 'x',
      width: 80,
      align: 'right',
      render: (x: number) => (
        <span style={{ fontFamily: 'monospace', fontSize: 12, color: token.colorText }}>{x.toFixed(3)}</span>
      ),
    },
    {
      title: 'Y (Å)',
      dataIndex: 'y',
      key: 'y',
      width: 80,
      align: 'right',
      render: (y: number) => (
        <span style={{ fontFamily: 'monospace', fontSize: 12, color: token.colorText }}>{y.toFixed(3)}</span>
      ),
    },
    {
      title: 'Z (Å)',
      dataIndex: 'z',
      key: 'z',
      width: 80,
      align: 'right',
      render: (z: number) => (
        <span style={{ fontFamily: 'monospace', fontSize: 12, color: token.colorText }}>{z.toFixed(3)}</span>
      ),
    },
  ];

  // 渲染分子卡片
  const renderMoleculeCard = (molecule: Molecule) => {
    const isSelected = selectedMolecule === molecule.name;
    const typeColor = molecule.type === 'cation' ? '#ff4d4f' :
                      molecule.type === 'anion' ? '#1890ff' : '#52c41a';

    return (
      <Card
        key={molecule.name}
        size="small"
        hoverable
        className={`molecule-card ${isSelected ? 'selected' : ''}`}
        style={{
          cursor: 'pointer',
          borderRadius: 8,
        }}
        onClick={() => setSelectedMolecule(molecule.name)}
      >
        <Space direction="vertical" size={8} style={{ width: '100%' }}>
          <div style={{
            fontSize: 15,
            fontWeight: 600,
            color: token.colorText,
            display: 'flex',
            alignItems: 'center',
            gap: 8,
          }}>
            <ExperimentOutlined style={{ color: typeColor }} />
            {molecule.name}
          </div>
          <Space wrap size={[4, 4]}>
            <Tag
              color={molecule.type === 'cation' ? 'red' : molecule.type === 'anion' ? 'blue' : 'green'}
              style={{ borderRadius: 4, fontSize: 11 }}
            >
              {molecule.type === 'cation' ? '阳离子' :
               molecule.type === 'anion' ? '阴离子' :
               '溶剂'}
            </Tag>
            <Tag
              color={molecule.total_charge > 0 ? 'red' : molecule.total_charge < 0 ? 'blue' : 'default'}
              style={{ borderRadius: 4, fontSize: 11, fontWeight: 600 }}
            >
              {molecule.total_charge > 0 ? '+' : ''}{molecule.total_charge.toFixed(2)} e
            </Tag>
          </Space>
          <div style={{
            fontSize: 11,
            color: token.colorTextSecondary,
            display: 'flex',
            justifyContent: 'space-between',
          }}>
            <span>{molecule.atoms.length} 个原子</span>
            <span style={{ fontStyle: 'italic' }}>{molecule.charge_method || 'ligpargen'}</span>
          </div>
        </Space>
      </Card>
    );
  };

  return (
    <div style={{
      background: token.colorBgLayout,
      padding: DASHBOARD_STYLES.gutter,
      minHeight: '100%',
      borderRadius: 8,
    }}>
      <style>{getResponsiveStyles(isDark)}</style>

      <Space direction="vertical" size={DASHBOARD_STYLES.gutter} style={{ width: '100%' }}>
        {/* 统计卡片网格 */}
        <div className="mol-stats-grid">
          {/* 分子种类 */}
          <div className="dashboard-card" style={dashboardCardStyle}>
            <div style={{ padding: DASHBOARD_STYLES.cardPadding }}>
              <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
                <div style={{
                  width: 48, height: 48, borderRadius: 12,
                  background: 'linear-gradient(135deg, #1890ff 0%, #096dd9 100%)',
                  display: 'flex', alignItems: 'center', justifyContent: 'center',
                }}>
                  <ExperimentOutlined style={{ fontSize: 24, color: '#fff' }} />
                </div>
                <div>
                  <div style={{ fontSize: 12, color: '#6b7280' }}>分子种类 (Species)</div>
                  <div style={{ fontSize: 28, fontWeight: 700, color: '#1890ff' }}>
                    {molecules.length}
                  </div>
                  <div style={{ fontSize: 11, color: '#94a3b8' }}>电解液组分</div>
                </div>
              </div>
            </div>
          </div>

          {/* 阳离子 */}
          <div className="dashboard-card" style={dashboardCardStyle}>
            <div style={{ padding: DASHBOARD_STYLES.cardPadding }}>
              <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
                <div style={{
                  width: 48, height: 48, borderRadius: 12,
                  background: 'linear-gradient(135deg, #ff4d4f 0%, #cf1322 100%)',
                  display: 'flex', alignItems: 'center', justifyContent: 'center',
                }}>
                  <span style={{ fontSize: 18, color: '#fff', fontWeight: 700 }}>+</span>
                </div>
                <div>
                  <div style={{ fontSize: 12, color: '#6b7280' }}>阳离子 (Cations)</div>
                  <div style={{ fontSize: 28, fontWeight: 700, color: '#ff4d4f' }}>
                    {cations.length}
                  </div>
                  <div style={{ fontSize: 11, color: '#94a3b8' }}>
                    {cations.map(c => c.name).join(', ') || '-'}
                  </div>
                </div>
              </div>
            </div>
          </div>

          {/* 阴离子 */}
          <div className="dashboard-card" style={dashboardCardStyle}>
            <div style={{ padding: DASHBOARD_STYLES.cardPadding }}>
              <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
                <div style={{
                  width: 48, height: 48, borderRadius: 12,
                  background: 'linear-gradient(135deg, #1890ff 0%, #0050b3 100%)',
                  display: 'flex', alignItems: 'center', justifyContent: 'center',
                }}>
                  <span style={{ fontSize: 18, color: '#fff', fontWeight: 700 }}>−</span>
                </div>
                <div>
                  <div style={{ fontSize: 12, color: '#6b7280' }}>阴离子 (Anions)</div>
                  <div style={{ fontSize: 28, fontWeight: 700, color: '#1890ff' }}>
                    {anions.length}
                  </div>
                  <div style={{ fontSize: 11, color: '#94a3b8' }}>
                    {anions.map(a => a.name).join(', ') || '-'}
                  </div>
                </div>
              </div>
            </div>
          </div>

          {/* 溶剂 */}
          <div className="dashboard-card" style={dashboardCardStyle}>
            <div style={{ padding: DASHBOARD_STYLES.cardPadding }}>
              <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
                <div style={{
                  width: 48, height: 48, borderRadius: 12,
                  background: 'linear-gradient(135deg, #52c41a 0%, #389e0d 100%)',
                  display: 'flex', alignItems: 'center', justifyContent: 'center',
                }}>
                  <span style={{ fontSize: 16, color: '#fff', fontWeight: 700 }}>S</span>
                </div>
                <div>
                  <div style={{ fontSize: 12, color: '#6b7280' }}>溶剂 (Solvents)</div>
                  <div style={{ fontSize: 28, fontWeight: 700, color: '#52c41a' }}>
                    {solvents.length}
                  </div>
                  <div style={{ fontSize: 11, color: '#94a3b8' }}>
                    {solvents.map(s => s.name).join(', ') || '-'}
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* 阴阳离子 */}
        {(cations.length > 0 || anions.length > 0) && (
          <Card
            className="dashboard-card"
            style={dashboardCardStyle}
            title={
              <Space size={8}>
                <ExperimentOutlined style={{ color: token.colorText }} />
                <span style={{ fontSize: 14, fontWeight: 600, color: token.colorText }}>
                  阴阳离子 (Ion Species)
                </span>
                <Tag color="red" style={{ marginLeft: 8, borderRadius: 4, fontSize: 11 }}>{cations.length} 阳</Tag>
                <Tag color="blue" style={{ borderRadius: 4, fontSize: 11 }}>{anions.length} 阴</Tag>
              </Space>
            }
          >
            <div className="mol-list-grid">
              {cations.map((mol) => renderMoleculeCard(mol))}
              {anions.map((mol) => renderMoleculeCard(mol))}
            </div>
          </Card>
        )}

        {/* 溶剂 */}
        {solvents.length > 0 && (
          <Card
            className="dashboard-card"
            style={dashboardCardStyle}
            title={
              <Space size={8}>
                <ExperimentOutlined style={{ color: token.colorText }} />
                <span style={{ fontSize: 14, fontWeight: 600, color: token.colorText }}>
                  溶剂分子 (Solvent Molecules)
                </span>
                <Tag color="green" style={{ marginLeft: 8, borderRadius: 4, fontSize: 11 }}>{solvents.length} 种</Tag>
              </Space>
            }
          >
            <div className="mol-list-grid">
              {solvents.map((mol) => renderMoleculeCard(mol))}
            </div>
          </Card>
        )}

        {/* 分子详情 */}
        {currentMolecule && (
          <Card
            className="dashboard-card"
            style={dashboardCardStyle}
            title={
              <Space size={8}>
                <EyeOutlined style={{ color: token.colorText }} />
                <span style={{ fontSize: 14, fontWeight: 600, color: token.colorText }}>
                  {currentMolecule.name} 详细信息
                </span>
              </Space>
            }
          >
            <Row gutter={DASHBOARD_STYLES.gutter}>
              {/* 3D 结构显示 */}
              <Col xs={24} lg={12}>
                <Card
                  title={
                    <span style={{ fontSize: 13, fontWeight: 600, color: token.colorText }}>
                      3D 结构
                    </span>
                  }
                  size="small"
                  style={{
                    borderRadius: 8,
                    border: '1px solid #e8e8e8',
                  }}
                  extra={
                    <Space size={8}>
                      <Tag color="purple" style={{ borderRadius: 4, fontSize: 11 }}>
                        {currentMolecule.charge_method || 'ligpargen'}
                      </Tag>
                      <Button
                        size="small"
                        icon={<ReloadOutlined />}
                        onClick={(e) => {
                          console.log('Button clicked!');
                          resetView(e);
                        }}
                        style={{ borderRadius: 6 }}
                      >
                        重置视角
                      </Button>
                    </Space>
                  }
                >
                  <div
                    ref={viewerRef}
                    style={{
                      width: '100%',
                      height: '400px',
                      position: 'relative',
                      background: isDark ? '#1a1a1a' : 'white',
                      borderRadius: 6,
                    }}
                  />
                  <div style={{
                    marginTop: 12,
                    padding: 8,
                    textAlign: 'center',
                    background: isDark ? 'rgba(255,255,255,0.04)' : '#f8fafc',
                    borderRadius: 6,
                  }}>
                    <Text type="secondary" style={{ fontSize: 11 }}>
                      拖动旋转 | 滚轮缩放 | 右键平移
                    </Text>
                  </div>
                </Card>
              </Col>

              {/* 原子信息表格 */}
              <Col xs={24} lg={12}>
                <Card
                  title={
                    <span style={{ fontSize: 13, fontWeight: 600, color: token.colorText }}>
                      原子信息
                    </span>
                  }
                  size="small"
                  style={{
                    borderRadius: 8,
                    border: `1px solid ${token.colorBorder}`,
                  }}
                >
                  <div style={{
                    background: isDark ? 'rgba(255,255,255,0.02)' : 'linear-gradient(135deg, #fafbfc 0%, #f0f2f5 100%)',
                    border: `1px solid ${token.colorBorder}`,
                    borderRadius: 8,
                    padding: 8,
                  }}>
                    <Table
                      columns={atomColumns}
                      dataSource={currentMolecule.atoms}
                      rowKey="id"
                      pagination={false}
                      size="small"
                      scroll={{ y: 320 }}
                      style={{
                        borderRadius: 6,
                        background: token.colorBgContainer,
                      }}
                    />
                  </div>
                  <div style={{
                    marginTop: 12,
                    padding: 12,
                    background: isDark ? 'rgba(24, 144, 255, 0.1)' : 'linear-gradient(135deg, #f0f9ff 0%, #e0f2fe 100%)',
                    borderRadius: 6,
                    border: `1px solid ${isDark ? 'rgba(24, 144, 255, 0.3)' : '#bae6fd'}`,
                  }}>
                    <Row gutter={16}>
                      <Col span={12}>
                        <div style={{ fontSize: 11, color: token.colorTextSecondary }}>总电荷</div>
                        <div style={{
                          fontSize: 16,
                          fontWeight: 700,
                          color: currentMolecule.total_charge > 0 ? '#ff4d4f' :
                                 currentMolecule.total_charge < 0 ? '#1890ff' : '#52c41a',
                          marginTop: 4,
                        }}>
                          {currentMolecule.total_charge > 0 ? '+' : ''}{currentMolecule.total_charge.toFixed(3)} e
                        </div>
                      </Col>
                      <Col span={12}>
                        <div style={{ fontSize: 11, color: token.colorTextSecondary }}>原子数量</div>
                        <div style={{ fontSize: 16, fontWeight: 700, color: token.colorText, marginTop: 4 }}>
                          {currentMolecule.atoms.length}
                        </div>
                      </Col>
                    </Row>
                  </div>
                </Card>
              </Col>
            </Row>

            {/* QC量子化学属性 */}
            {currentMolecule.smiles && qcCache[currentMolecule.smiles] && (
              <>
                <Divider />
                <Card
                  title={
                    <Space size={8}>
                      <ThunderboltOutlined style={{ color: '#722ed1' }} />
                      <span style={{ fontSize: 13, fontWeight: 600, color: token.colorText }}>
                        量子化学性质
                      </span>
                    </Space>
                  }
                  size="small"
                  style={{
                    borderRadius: 8,
                    border: `1px solid ${token.colorBorder}`,
                  }}
                >
                  {(() => {
                    const qc = qcCache[currentMolecule.smiles!];
                    const hasImages = qc.esp_image_path || qc.homo_image_path || qc.lumo_image_path;
                    return (
                      <>
                        <Row gutter={24}>
                          <Col xs={24} md={hasImages ? 8 : 24}>
                            <Descriptions column={1} size="small" bordered>
                              <Descriptions.Item label="能量 (A.U.)">
                                {qc.energy_au?.toFixed(6) || '-'}
                              </Descriptions.Item>
                              <Descriptions.Item label="HOMO">
                                <Tag color="blue">{qc.homo_ev?.toFixed(3) || '-'} eV</Tag>
                              </Descriptions.Item>
                              <Descriptions.Item label="LUMO">
                                <Tag color="green">{qc.lumo_ev?.toFixed(3) || '-'} eV</Tag>
                              </Descriptions.Item>
                              <Descriptions.Item label="HOMO-LUMO Gap">
                                <Tag color="orange">{qc.homo_lumo_gap_ev?.toFixed(3) || '-'} eV</Tag>
                              </Descriptions.Item>
                              <Descriptions.Item label="ESP Min">
                                {qc.esp_min_kcal?.toFixed(2) || '-'} kcal/mol
                              </Descriptions.Item>
                              <Descriptions.Item label="ESP Max">
                                {qc.esp_max_kcal?.toFixed(2) || '-'} kcal/mol
                              </Descriptions.Item>
                              <Descriptions.Item label="基组/泛函">
                                {qc.basis_set || '-'} / {qc.functional || '-'}
                              </Descriptions.Item>
                            </Descriptions>
                          </Col>
                          {hasImages && (
                            <Col xs={24} md={16}>
                              <Row gutter={[8, 8]}>
                                {qc.esp_image_path && (
                                  <Col xs={24} sm={8}>
                                    <div style={{ textAlign: 'center' }}>
                                      <Text type="secondary" style={{ display: 'block', marginBottom: 4, fontSize: 12 }}>
                                        静电势 (ESP)
                                      </Text>
                                      <Image
                                        src={`/api/v1/qc/esp-image-by-smiles/${encodeURIComponent(currentMolecule.smiles!)}`}
                                        alt="ESP Surface"
                                        style={{ maxWidth: '100%', maxHeight: 180 }}
                                        placeholder
                                      />
                                    </div>
                                  </Col>
                                )}
                                {qc.homo_image_path && (
                                  <Col xs={24} sm={8}>
                                    <div style={{ textAlign: 'center' }}>
                                      <Text type="secondary" style={{ display: 'block', marginBottom: 4, fontSize: 12 }}>
                                        HOMO 轨道
                                      </Text>
                                      <Image
                                        src={`/api/v1/qc/homo-image-by-smiles/${encodeURIComponent(currentMolecule.smiles!)}`}
                                        alt="HOMO Orbital"
                                        style={{ maxWidth: '100%', maxHeight: 180 }}
                                        placeholder
                                      />
                                    </div>
                                  </Col>
                                )}
                                {qc.lumo_image_path && (
                                  <Col xs={24} sm={8}>
                                    <div style={{ textAlign: 'center' }}>
                                      <Text type="secondary" style={{ display: 'block', marginBottom: 4, fontSize: 12 }}>
                                        LUMO 轨道
                                      </Text>
                                      <Image
                                        src={`/api/v1/qc/lumo-image-by-smiles/${encodeURIComponent(currentMolecule.smiles!)}`}
                                        alt="LUMO Orbital"
                                        style={{ maxWidth: '100%', maxHeight: 180 }}
                                        placeholder
                                      />
                                    </div>
                                  </Col>
                                )}
                              </Row>
                            </Col>
                          )}
                        </Row>
                      </>
                    );
                  })()}
                </Card>
              </>
            )}

            {/* 无QC数据提示 */}
            {currentMolecule.smiles && !qcCache[currentMolecule.smiles] && !loadingQC && (
              <div style={{ marginTop: 16, textAlign: 'center' }}>
                <Text type="secondary">
                  暂无QC计算数据。可在创建MD任务时勾选"启用QC计算"，或前往
                  <a href="/workspace/liquid-electrolyte/qc" style={{ marginLeft: 4 }}>QC计算</a>
                  页面单独提交计算。
                </Text>
              </div>
            )}
          </Card>
        )}
      </Space>
    </div>
  );
}

