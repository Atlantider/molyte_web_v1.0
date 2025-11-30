/**
 * 用户指南页面
 */
import { useNavigate } from 'react-router-dom';
import {
  Typography,
  Card,
  Row,
  Col,
  Collapse,
  Button,
  Divider,
  Space,
} from 'antd';
import {
  HomeOutlined,
  ThunderboltOutlined,
  BookOutlined,
  ExperimentOutlined,
  RocketOutlined,
  LineChartOutlined,
  QuestionCircleOutlined,
  BulbOutlined,
  SafetyCertificateOutlined,
  AimOutlined,
  AreaChartOutlined,
  DotChartOutlined,
  ArrowRightOutlined,
  CheckCircleOutlined,
} from '@ant-design/icons';
import { useAuthStore } from '../stores/authStore';
import './Guide.css';

const { Title, Paragraph, Text } = Typography;

export default function Guide() {
  const navigate = useNavigate();
  const { isAuthenticated } = useAuthStore();

  return (
    <div className="guide-container">
      {/* 顶部导航栏 */}
      <header className="guide-header">
        <div className="guide-header-content">
          <div className="guide-logo" onClick={() => navigate('/')}>
            <div className="guide-logo-icon">
              <ThunderboltOutlined />
            </div>
            <span className="guide-logo-text">Molyte</span>
          </div>
          
          <nav className="guide-nav">
            <a href="#platform">平台介绍</a>
            <a href="#md-theory">MD原理</a>
            <a href="#workflow">使用流程</a>
            <a href="#results">结果解读</a>
            <a href="#faq">常见问题</a>
          </nav>
          
          <div className="guide-header-actions">
            <Button
              type="text"
              icon={<HomeOutlined />}
              onClick={() => navigate('/')}
              className="guide-back-btn"
            >
              返回首页
            </Button>
            {isAuthenticated ? (
              <Button
                type="primary"
                onClick={() => navigate('/workspace/dashboard')}
                className="guide-action-btn"
              >
                进入工作台
              </Button>
            ) : (
              <Button
                type="primary"
                onClick={() => navigate('/login')}
                className="guide-action-btn"
              >
                登录 / 注册
              </Button>
            )}
          </div>
        </div>
      </header>

      {/* 主内容区域 */}
      <main className="guide-main">
        {/* Hero Section */}
        <section className="guide-hero">
          <div className="guide-hero-content">
            <div className="guide-hero-icon">
              <BookOutlined />
            </div>
            <Title level={1} className="guide-hero-title">
              用户指南
            </Title>
            <Paragraph className="guide-hero-subtitle">
              全面了解 Molyte 电解液研发模拟平台的功能与使用方法
            </Paragraph>
          </div>
        </section>

        {/* 快速导航卡片 */}
        <section className="guide-quick-nav">
          <Row gutter={[24, 24]}>
            <Col xs={24} sm={12} lg={6}>
              <Card className="quick-nav-card" hoverable onClick={() => document.getElementById('platform')?.scrollIntoView({ behavior: 'smooth' })}>
                <div className="quick-nav-icon icon-blue"><RocketOutlined /></div>
                <Title level={4}>平台介绍</Title>
                <Text type="secondary">了解平台核心功能</Text>
              </Card>
            </Col>
            <Col xs={24} sm={12} lg={6}>
              <Card className="quick-nav-card" hoverable onClick={() => document.getElementById('md-theory')?.scrollIntoView({ behavior: 'smooth' })}>
                <div className="quick-nav-icon icon-purple"><ExperimentOutlined /></div>
                <Title level={4}>MD 原理</Title>
                <Text type="secondary">分子动力学基础知识</Text>
              </Card>
            </Col>
            <Col xs={24} sm={12} lg={6}>
              <Card className="quick-nav-card" hoverable onClick={() => document.getElementById('results')?.scrollIntoView({ behavior: 'smooth' })}>
                <div className="quick-nav-icon icon-green"><LineChartOutlined /></div>
                <Title level={4}>结果解读</Title>
                <Text type="secondary">RDF/MSD 结果分析</Text>
              </Card>
            </Col>
            <Col xs={24} sm={12} lg={6}>
              <Card className="quick-nav-card" hoverable onClick={() => document.getElementById('faq')?.scrollIntoView({ behavior: 'smooth' })}>
                <div className="quick-nav-icon icon-orange"><QuestionCircleOutlined /></div>
                <Title level={4}>常见问题</Title>
                <Text type="secondary">FAQ 与技术支持</Text>
              </Card>
            </Col>
          </Row>
        </section>

        {/* 平台介绍 */}
        <section id="platform" className="guide-section">
          <div className="section-header">
            <div className="section-icon"><RocketOutlined /></div>
            <Title level={2}>平台介绍</Title>
          </div>

          <Card className="guide-card">
            <Title level={3}>🎯 平台定位</Title>
            <Paragraph className="guide-paragraph">
              <Text strong>Molyte</Text> 是一个专注于<Text mark>电池电解液分子动力学模拟</Text>的一站式云计算平台。
              平台整合了高性能计算资源、自动化工作流和智能数据分析，帮助研究人员快速完成电解液配方的筛选与优化。
            </Paragraph>

            <Divider />

            <Title level={3}>✨ 核心功能</Title>
            <Row gutter={[24, 24]}>
              <Col xs={24} sm={12} lg={8}>
                <div className="feature-item">
                  <div className="feature-item-icon" style={{ background: 'linear-gradient(135deg, #667eea 0%, #5a67d8 100%)' }}>
                    <ExperimentOutlined />
                  </div>
                  <Title level={4}>配方管理</Title>
                  <Paragraph>
                    可视化配置电解液配方，支持多种阳离子（Li⁺、Na⁺、K⁺等）、阴离子（FSI⁻、TFSI⁻、PF₆⁻等）和溶剂（EC、DMC、EMC等）的组合。
                  </Paragraph>
                </div>
              </Col>
              <Col xs={24} sm={12} lg={8}>
                <div className="feature-item">
                  <div className="feature-item-icon" style={{ background: 'linear-gradient(135deg, #48bb78 0%, #38a169 100%)' }}>
                    <RocketOutlined />
                  </div>
                  <Title level={4}>自动化计算</Title>
                  <Paragraph>
                    一键提交分子动力学模拟任务，系统自动生成 LAMMPS 输入文件、调度 HPC 资源、执行计算并返回结果。
                  </Paragraph>
                </div>
              </Col>
              <Col xs={24} sm={12} lg={8}>
                <div className="feature-item">
                  <div className="feature-item-icon" style={{ background: 'linear-gradient(135deg, #ed8936 0%, #dd6b20 100%)' }}>
                    <LineChartOutlined />
                  </div>
                  <Title level={4}>结果分析</Title>
                  <Paragraph>
                    自动计算 RDF（径向分布函数）、MSD（均方位移）等关键参数，可视化展示离子溶剂化结构和传输性能。
                  </Paragraph>
                </div>
              </Col>
              <Col xs={24} sm={12} lg={8}>
                <div className="feature-item">
                  <div className="feature-item-icon" style={{ background: 'linear-gradient(135deg, #9f7aea 0%, #764ba2 100%)' }}>
                    <AimOutlined />
                  </div>
                  <Title level={4}>智能原子映射</Title>
                  <Paragraph>
                    自动识别分子中的原子类型和化学环境，精确区分不同分子中的同类原子（如 EC 的羰基氧与醚氧）。
                  </Paragraph>
                </div>
              </Col>
              <Col xs={24} sm={12} lg={8}>
                <div className="feature-item">
                  <div className="feature-item-icon" style={{ background: 'linear-gradient(135deg, #f56565 0%, #c53030 100%)' }}>
                    <SafetyCertificateOutlined />
                  </div>
                  <Title level={4}>数据复用</Title>
                  <Paragraph>
                    所有计算结果自动存入公共数据库，相同配方可直接查询历史结果，避免重复计算，节省资源。
                  </Paragraph>
                </div>
              </Col>
              <Col xs={24} sm={12} lg={8}>
                <div className="feature-item">
                  <div className="feature-item-icon" style={{ background: 'linear-gradient(135deg, #38b2ac 0%, #2c7a7b 100%)' }}>
                    <BulbOutlined />
                  </div>
                  <Title level={4}>多精度模式</Title>
                  <Paragraph>
                    支持快速（~1h）、标准（~12h）、精确（~36h）三种计算模式，满足不同场景的精度和时效需求。
                  </Paragraph>
                </div>
              </Col>
            </Row>
          </Card>
        </section>

        {/* 分子动力学原理 */}
        <section id="md-theory" className="guide-section">
          <div className="section-header">
            <div className="section-icon"><ExperimentOutlined /></div>
            <Title level={2}>分子动力学原理</Title>
          </div>

          <Card className="guide-card">
            <Title level={3}>📚 什么是分子动力学模拟？</Title>
            <Paragraph className="guide-paragraph">
              分子动力学（Molecular Dynamics, MD）模拟是一种通过数值求解牛顿运动方程来研究原子和分子运动的计算方法。
              通过追踪每个粒子的位置和速度随时间的演化，我们可以获得体系的微观结构和动力学性质。
            </Paragraph>

            <div className="theory-card">
              <Title level={4}>🔬 基本原理</Title>
              <Paragraph>
                在 MD 模拟中，每个原子被视为一个经典粒子，其运动遵循牛顿第二定律：
              </Paragraph>
              <div className="formula-box">
                F = ma = -∇U(r)
              </div>
              <Paragraph>
                其中 F 是原子受到的力，m 是原子质量，a 是加速度，U(r) 是势能函数。
                通过力场参数（如 OPLS-AA）描述原子间的相互作用，包括键合作用（键伸缩、键角弯曲、二面角扭转）和非键合作用（范德华力、静电作用）。
              </Paragraph>
            </div>

            <div className="theory-card">
              <Title level={4}>⚡ 电解液模拟流程</Title>
              <Paragraph>
                本平台采用经典的两阶段平衡策略：
              </Paragraph>
              <ol style={{ lineHeight: 2, color: '#4a5568' }}>
                <li><Text strong>NPT 平衡阶段</Text>：在恒温恒压条件下平衡体系密度（约 0.1-0.5 ns）</li>
                <li><Text strong>NVT 产出阶段</Text>：在恒温恒容条件下采集轨迹用于分析（约 0.5-5 ns）</li>
              </ol>
              <Paragraph>
                模拟使用 LAMMPS 软件包，力场参数由 LigParGen 生成（基于 OPLS-AA 力场），
                电荷通过 CM1A 或 CM1A-LBCC 方法计算。
              </Paragraph>
            </div>
          </Card>
        </section>

        {/* 使用流程 */}
        <section id="workflow" className="guide-section">
          <div className="section-header">
            <div className="section-icon"><ArrowRightOutlined /></div>
            <Title level={2}>使用流程</Title>
          </div>

          <Card className="guide-card">
            <div className="workflow-steps">
              <div className="workflow-step">
                <div className="workflow-step-number">1</div>
                <div className="workflow-step-content">
                  <Title level={4}>注册登录</Title>
                  <Paragraph>
                    创建账号并登录平台。新用户将获得一定的免费计算配额，可在个人中心查看剩余资源。
                  </Paragraph>
                </div>
              </div>
              <div className="workflow-step">
                <div className="workflow-step-number">2</div>
                <div className="workflow-step-content">
                  <Title level={4}>创建配方</Title>
                  <Paragraph>
                    在"配方管理"中创建电解液配方：选择阳离子、阴离子、溶剂组分，设置浓度和温度等参数。
                    系统支持多种常见的锂电池电解液组分。
                  </Paragraph>
                </div>
              </div>
              <div className="workflow-step">
                <div className="workflow-step-number">3</div>
                <div className="workflow-step-content">
                  <Title level={4}>提交计算</Title>
                  <Paragraph>
                    选择计算精度（快速/标准/精确），配置 Slurm 资源参数，一键提交任务。
                    系统自动生成输入文件并调度到 HPC 集群执行。
                  </Paragraph>
                </div>
              </div>
              <div className="workflow-step">
                <div className="workflow-step-number">4</div>
                <div className="workflow-step-content">
                  <Title level={4}>监控进度</Title>
                  <Paragraph>
                    在"计算任务"页面实时查看任务状态（排队中、运行中、已完成等），
                    支持查看 Slurm 日志和计算进度。
                  </Paragraph>
                </div>
              </div>
              <div className="workflow-step">
                <div className="workflow-step-number">5</div>
                <div className="workflow-step-content">
                  <Title level={4}>分析结果</Title>
                  <Paragraph>
                    任务完成后，进入详情页查看 RDF、MSD 等分析结果。
                    支持选择不同的原子对进行 RDF 计算，可视化展示离子配位结构。
                  </Paragraph>
                </div>
              </div>
            </div>
          </Card>
        </section>

        {/* 结果解读 */}
        <section id="results" className="guide-section">
          <div className="section-header">
            <div className="section-icon"><LineChartOutlined /></div>
            <Title level={2}>结果解读</Title>
          </div>

          <Row gutter={[24, 24]}>
            <Col xs={24} lg={12}>
              <Card className="result-card" styles={{ body: { padding: 0 } }}>
                <div className="result-card-header rdf">
                  <Title level={3}><AreaChartOutlined /> 径向分布函数 (RDF)</Title>
                  <Paragraph>Radial Distribution Function, g(r)</Paragraph>
                </div>
                <div className="result-card-body">
                  <Paragraph>
                    RDF 描述了以某一原子为中心，在距离 r 处找到另一原子的概率密度。
                    它是分析电解液微观结构的核心工具。
                  </Paragraph>

                  <Divider />

                  <div className="result-item">
                    <CheckCircleOutlined className="result-item-icon" />
                    <div className="result-item-content">
                      <Title level={5}>第一峰位置</Title>
                      <Paragraph>
                        表示最近邻原子的平均距离。例如 Li-O 的 RDF 第一峰约在 2.0 Å，
                        反映锂离子与氧原子的配位距离。
                      </Paragraph>
                    </div>
                  </div>

                  <div className="result-item">
                    <CheckCircleOutlined className="result-item-icon" />
                    <div className="result-item-content">
                      <Title level={5}>峰高与积分</Title>
                      <Paragraph>
                        峰高反映局部密度增强程度，对第一峰积分可得到配位数（Coordination Number），
                        即锂离子第一溶剂化壳层中的配位原子数目。
                      </Paragraph>
                    </div>
                  </div>

                  <div className="result-item">
                    <CheckCircleOutlined className="result-item-icon" />
                    <div className="result-item-content">
                      <Title level={5}>溶剂化结构</Title>
                      <Paragraph>
                        通过比较 Li 与不同氧原子（如 EC 羰基氧、DMC 醚氧、阴离子氧）的 RDF，
                        可以分析锂离子的优先溶剂化行为。
                      </Paragraph>
                    </div>
                  </div>
                </div>
              </Card>
            </Col>

            <Col xs={24} lg={12}>
              <Card className="result-card" styles={{ body: { padding: 0 } }}>
                <div className="result-card-header msd">
                  <Title level={3}><DotChartOutlined /> 均方位移 (MSD)</Title>
                  <Paragraph>Mean Square Displacement</Paragraph>
                </div>
                <div className="result-card-body">
                  <Paragraph>
                    MSD 描述粒子随时间的平均位移平方，是计算扩散系数的基础。
                    通过 MSD 可以评估电解液的离子传输性能。
                  </Paragraph>

                  <Divider />

                  <div className="result-item">
                    <CheckCircleOutlined className="result-item-icon" style={{ color: '#48bb78' }} />
                    <div className="result-item-content">
                      <Title level={5}>扩散系数</Title>
                      <Paragraph>
                        根据 Einstein 关系：D = lim(t→∞) MSD / 6t，
                        从 MSD 曲线的线性区斜率可计算扩散系数 D（单位：cm²/s）。
                      </Paragraph>
                    </div>
                  </div>

                  <div className="result-item">
                    <CheckCircleOutlined className="result-item-icon" style={{ color: '#48bb78' }} />
                    <div className="result-item-content">
                      <Title level={5}>离子电导率</Title>
                      <Paragraph>
                        通过 Nernst-Einstein 方程，可从扩散系数估算离子电导率：
                        σ = (F²/RT) × Σ(c_i × z_i² × D_i)，其中 c_i 是浓度，z_i 是电荷数。
                      </Paragraph>
                    </div>
                  </div>

                  <div className="result-item">
                    <CheckCircleOutlined className="result-item-icon" style={{ color: '#48bb78' }} />
                    <div className="result-item-content">
                      <Title level={5}>迁移数</Title>
                      <Paragraph>
                        比较阳离子和阴离子的扩散系数，可计算锂离子迁移数：
                        t⁺ = D_Li / (D_Li + D_anion)，迁移数越高表示锂离子传输效率越好。
                      </Paragraph>
                    </div>
                  </div>
                </div>
              </Card>
            </Col>
          </Row>
        </section>

        {/* 常见问题 */}
        <section id="faq" className="guide-section">
          <div className="section-header">
            <div className="section-icon"><QuestionCircleOutlined /></div>
            <Title level={2}>常见问题</Title>
          </div>

          <Collapse
            className="faq-collapse"
            expandIconPosition="end"
            items={[
              {
                key: '1',
                label: '如何选择合适的计算精度？',
                children: (
                  <div>
                    <Paragraph>
                      <Text strong>快速模式（~1小时）</Text>：适合快速测试配方、调试参数，结果仅供参考。
                    </Paragraph>
                    <Paragraph>
                      <Text strong>标准模式（~12小时）</Text>：适合一般研究使用，平衡精度和时间，满足大多数分析需求。
                    </Paragraph>
                    <Paragraph>
                      <Text strong>精确模式（~36小时）</Text>：适合论文发表级别的研究，提供高精度的 RDF 和扩散系数数据。
                    </Paragraph>
                  </div>
                ),
              },
              {
                key: '2',
                label: '任务一直排队怎么办？',
                children: (
                  <Paragraph>
                    任务排队时间取决于 HPC 集群的负载情况。如果长时间排队，可以尝试：
                    1) 选择负载较低的计算队列；
                    2) 减少请求的 CPU 核心数；
                    3) 在非高峰时段（如夜间）提交任务。
                  </Paragraph>
                ),
              },
              {
                key: '3',
                label: 'RDF 曲线为什么有振荡？',
                children: (
                  <Paragraph>
                    RDF 曲线的振荡反映了溶剂化壳层结构。第一个峰对应第一溶剂化壳层，
                    后续的峰表示更远处的溶剂化壳层。在液态体系中，长程结构逐渐消失，
                    g(r) 趋近于 1。如果振荡异常剧烈或不收敛到 1，可能需要增加模拟时间。
                  </Paragraph>
                ),
              },
              {
                key: '4',
                label: '如何解读负的 MSD 斜率？',
                children: (
                  <Paragraph>
                    物理上 MSD 应该是单调递增的。如果观察到负斜率，可能是：
                    1) 统计采样不足，需要更长的模拟时间；
                    2) 体系未充分平衡，NPT 阶段时间不够；
                    3) 分析的时间窗口选择不当。建议使用更长的 NVT 产出阶段。
                  </Paragraph>
                ),
              },
              {
                key: '5',
                label: '支持哪些溶剂和盐类？',
                children: (
                  <Paragraph>
                    目前支持常见的锂电池电解液组分：
                    <br /><Text strong>溶剂</Text>：EC、DMC、EMC、DEC、PC、FEC 等碳酸酯类溶剂
                    <br /><Text strong>阳离子</Text>：Li⁺、Na⁺、K⁺、Mg²⁺ 等
                    <br /><Text strong>阴离子</Text>：FSI⁻、TFSI⁻、PF₆⁻、BF₄⁻、ClO₄⁻ 等
                    <br />如需添加新的分子，请联系管理员。
                  </Paragraph>
                ),
              },
              {
                key: '6',
                label: '如何引用本平台？',
                children: (
                  <Paragraph>
                    如果您在研究中使用了本平台的计算结果，请在论文中引用：
                    <br /><br />
                    <Text code>
                      Molyte Platform: A Web-based Molecular Dynamics Simulation Platform for Battery Electrolytes.
                      https://molyte.example.com
                    </Text>
                  </Paragraph>
                ),
              },
            ]}
          />
        </section>

        {/* 底部 CTA */}
        <section className="guide-section">
          <Card className="guide-card" style={{ textAlign: 'center', background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)' }}>
            <Title level={3} style={{ color: '#fff', marginBottom: 16 }}>
              准备好开始您的研究了吗？
            </Title>
            <Paragraph style={{ color: 'rgba(255,255,255,0.9)', marginBottom: 24, fontSize: 16 }}>
              立即注册，获取免费计算配额，体验高效的电解液模拟平台
            </Paragraph>
            <Space size={16}>
              {isAuthenticated ? (
                <Button
                  size="large"
                  onClick={() => navigate('/workspace/dashboard')}
                  style={{ borderRadius: 8, fontWeight: 500 }}
                >
                  进入工作台
                </Button>
              ) : (
                <>
                  <Button
                    size="large"
                    onClick={() => navigate('/login?tab=register')}
                    style={{ borderRadius: 8, fontWeight: 500 }}
                  >
                    免费注册
                  </Button>
                  <Button
                    size="large"
                    type="default"
                    ghost
                    onClick={() => navigate('/login')}
                    style={{ borderRadius: 8, fontWeight: 500 }}
                  >
                    登录
                  </Button>
                </>
              )}
            </Space>
          </Card>
        </section>
      </main>

      {/* 页脚 */}
      <footer style={{
        background: '#1a1f36',
        padding: '32px 24px',
        textAlign: 'center',
        color: 'rgba(255,255,255,0.6)',
        fontSize: 14,
      }}>
        © 2025 Molyte Platform. All rights reserved.
      </footer>
    </div>
  );
}

