/**
 * 公开主页 - 项目介绍
 */
import { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { Button, Row, Col, Drawer } from 'antd';
import {
  ArrowRightOutlined,
  ThunderboltOutlined,
  ExperimentOutlined,
  LineChartOutlined,
  RocketOutlined,
  DatabaseOutlined,
  SafetyCertificateOutlined,
  MenuOutlined,
  CloseOutlined,
  MailOutlined,
  PhoneOutlined,
  EnvironmentOutlined,
} from '@ant-design/icons';
import { useAuthStore } from '../stores/authStore';
import './Home.css';

export default function Home() {
  const navigate = useNavigate();
  const { isAuthenticated } = useAuthStore();
  const [mobileMenuOpen, setMobileMenuOpen] = useState(false);

  const scrollToSection = (id: string) => {
    const element = document.getElementById(id);
    if (element) {
      element.scrollIntoView({ behavior: 'smooth' });
    }
    setMobileMenuOpen(false);
  };

  return (
    <div className="home-container">
      {/* 顶部导航栏 */}
      <header className="home-header">
        <div className="header-content">
          <div className="logo" onClick={() => scrollToSection('hero')}>
            <div className="logo-icon">
              <ThunderboltOutlined />
            </div>
            <span className="logo-text">Molyte</span>
          </div>

          {/* 桌面导航 */}
          <nav className="nav-menu desktop-nav">
            <a onClick={() => scrollToSection('hero')} className="nav-link">首页</a>
            <a onClick={() => scrollToSection('features')} className="nav-link">功能模块</a>
            <a onClick={() => scrollToSection('news')} className="nav-link">新闻动态</a>
            <a onClick={() => scrollToSection('about')} className="nav-link">关于我们</a>
            <a onClick={() => navigate('/research')} className="nav-link">数据库</a>
            <a onClick={() => navigate('/guide')} className="nav-link">用户指南</a>
          </nav>

          <div className="header-actions">
            {isAuthenticated ? (
              <Button
                type="primary"
                className="header-btn primary"
                icon={<RocketOutlined />}
                onClick={() => navigate('/workspace/dashboard')}
              >
                进入工作台
              </Button>
            ) : (
              <>
                <Button
                  type="text"
                  className="header-btn text"
                  onClick={() => navigate('/login')}
                >
                  登录
                </Button>
                <Button
                  type="primary"
                  className="header-btn primary"
                  onClick={() => navigate('/login?tab=register')}
                >
                  免费试用
                </Button>
              </>
            )}

            {/* 移动端菜单按钮 */}
            <Button
              className="mobile-menu-btn"
              type="text"
              icon={<MenuOutlined />}
              onClick={() => setMobileMenuOpen(true)}
            />
          </div>
        </div>
      </header>

      {/* 移动端抽屉菜单 */}
      <Drawer
        title={null}
        placement="right"
        onClose={() => setMobileMenuOpen(false)}
        open={mobileMenuOpen}
        width={280}
        className="mobile-drawer"
        closeIcon={<CloseOutlined />}
      >
        <div className="mobile-menu">
          <a onClick={() => scrollToSection('hero')} className="mobile-nav-link">首页</a>
          <a onClick={() => scrollToSection('features')} className="mobile-nav-link">功能模块</a>
          <a onClick={() => scrollToSection('news')} className="mobile-nav-link">新闻动态</a>
          <a onClick={() => scrollToSection('about')} className="mobile-nav-link">关于我们</a>
          <a onClick={() => { navigate('/research'); setMobileMenuOpen(false); }} className="mobile-nav-link">数据库</a>
          <a onClick={() => { navigate('/guide'); setMobileMenuOpen(false); }} className="mobile-nav-link">用户指南</a>
          <div className="mobile-menu-divider" />
          {isAuthenticated ? (
            <Button type="primary" block onClick={() => navigate('/workspace/dashboard')}>
              进入工作台
            </Button>
          ) : (
            <>
              <Button block onClick={() => navigate('/login')} style={{ marginBottom: 12 }}>登录</Button>
              <Button type="primary" block onClick={() => navigate('/login?tab=register')}>免费试用</Button>
            </>
          )}
        </div>
      </Drawer>

      {/* 主视觉区域 */}
      <section id="hero" className="hero-section">
        <div className="hero-background">
          <div className="bg-particles"></div>
          <div className="bg-gradient"></div>
          <div className="bg-grid"></div>
          <div className="floating-shapes">
            <div className="shape shape-1"></div>
            <div className="shape shape-2"></div>
            <div className="shape shape-3"></div>
          </div>
        </div>

        <div className="hero-content">
          <div className="hero-text">
            <h1 className="hero-title">
              <span className="title-highlight">新一代</span>电解液研发模拟平台
            </h1>
            <p className="hero-subtitle">
              基于分子动力学的高通量计算与AI深度预测，加速电池电解液研发创新
            </p>
            <div className="hero-actions">
              <Button
                type="primary"
                size="large"
                className="hero-btn primary"
                icon={<RocketOutlined />}
                onClick={() => isAuthenticated ? navigate('/workspace/dashboard') : navigate('/login')}
              >
                开始计算
              </Button>
              <Button
                size="large"
                className="hero-btn secondary"
                onClick={() => scrollToSection('features')}
              >
                了解更多
              </Button>
            </div>
            <div className="hero-stats">
              <div className="stat-item">
                <span className="stat-number">1000+</span>
                <span className="stat-label">配方数据</span>
              </div>
              <div className="stat-divider"></div>
              <div className="stat-item">
                <span className="stat-number">50+</span>
                <span className="stat-label">合作机构</span>
              </div>
              <div className="stat-divider"></div>
              <div className="stat-item">
                <span className="stat-number">99.9%</span>
                <span className="stat-label">计算准确率</span>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* 功能模块区域 */}
      <section id="features" className="features-section">
        <div className="section-content">
          <div className="section-header">
            <span className="section-tag">核心功能</span>
            <h2 className="section-title">强大的功能模块</h2>
            <p className="section-subtitle">集成高通量计算、数据分析、性能预测于一体的综合研发平台</p>
          </div>

          <Row gutter={[32, 32]}>
            <Col xs={24} md={8}>
              <div className="feature-card feature-card-1">
                <div className="feature-icon">
                  <ThunderboltOutlined />
                </div>
                <h3 className="feature-title">高通量计算</h3>
                <p className="feature-desc">
                  采用先进的分子动力学算法，支持批量电解液配方筛选，
                  快速获取RDF、MSD等关键物理化学参数
                </p>
                <Button
                  type="link"
                  className="feature-link"
                  onClick={() => isAuthenticated ? navigate('/workspace/electrolytes') : navigate('/login')}
                >
                  开始计算 <ArrowRightOutlined />
                </Button>
              </div>
            </Col>
            <Col xs={24} md={8}>
              <div className="feature-card feature-card-2">
                <div className="feature-icon">
                  <LineChartOutlined />
                </div>
                <h3 className="feature-title">数据分析</h3>
                <p className="feature-desc">
                  深度分析模拟结果，提供可视化图表展示，
                  支持RDF分布、扩散系数、溶剂化结构等多维度分析
                </p>
                <Button
                  type="link"
                  className="feature-link"
                  onClick={() => isAuthenticated ? navigate('/workspace/jobs') : navigate('/login')}
                >
                  查看分析 <ArrowRightOutlined />
                </Button>
              </div>
            </Col>
            <Col xs={24} md={8}>
              <div className="feature-card feature-card-3">
                <div className="feature-icon">
                  <ExperimentOutlined />
                </div>
                <h3 className="feature-title">性能预测</h3>
                <p className="feature-desc">
                  基于AI模型的电解液性能预测，结合历史数据和机器学习，
                  准确预测电导率、粘度等关键性能指标
                </p>
                <Button
                  type="link"
                  className="feature-link"
                  onClick={() => navigate('/research')}
                >
                  探索数据 <ArrowRightOutlined />
                </Button>
              </div>
            </Col>
          </Row>
        </div>
      </section>

      {/* 新闻动态 */}
      <section id="news" className="news-section">
        <div className="section-content">
          <div className="section-header">
            <span className="section-tag">最新动态</span>
            <h2 className="section-title">新闻与公告</h2>
            <p className="section-subtitle">了解平台最新进展和行业动态</p>
          </div>

          <Row gutter={[24, 24]}>
            <Col xs={24} md={8}>
              <div className="news-card">
                <div className="news-image">
                  <div className="news-image-placeholder">
                    <RocketOutlined />
                  </div>
                </div>
                <div className="news-content">
                  <div className="news-meta">
                    <span className="news-date">2025-02-23</span>
                    <span className="news-tag">平台动态</span>
                  </div>
                  <h3 className="news-title">Molyte平台正式上线</h3>
                  <p className="news-desc">
                    Molyte网站版正式上线，为电解液研发提供全新的技术支持平台...
                  </p>
                  <a href="#" className="news-link">
                    阅读更多 <ArrowRightOutlined />
                  </a>
                </div>
              </div>
            </Col>
            <Col xs={24} md={8}>
              <div className="news-card">
                <div className="news-image">
                  <div className="news-image-placeholder">
                    <SafetyCertificateOutlined />
                  </div>
                </div>
                <div className="news-content">
                  <div className="news-meta">
                    <span className="news-date">2025-02-20</span>
                    <span className="news-tag">合作动态</span>
                  </div>
                  <h3 className="news-title">与中科院研究机构达成战略合作</h3>
                  <p className="news-desc">
                    与中科院知名研究机构展开深度合作，共同推进电解液技术创新...
                  </p>
                  <a href="#" className="news-link">
                    阅读更多 <ArrowRightOutlined />
                  </a>
                </div>
              </div>
            </Col>
            <Col xs={24} md={8}>
              <div className="news-card">
                <div className="news-image">
                  <div className="news-image-placeholder">
                    <DatabaseOutlined />
                  </div>
                </div>
                <div className="news-content">
                  <div className="news-meta">
                    <span className="news-date">2025-02-15</span>
                    <span className="news-tag">功能更新</span>
                  </div>
                  <h3 className="news-title">平台功能全面升级</h3>
                  <p className="news-desc">
                    引入量子化学计算模块，数据分析能力大幅提升，用户体验优化...
                  </p>
                  <a href="#" className="news-link">
                    阅读更多 <ArrowRightOutlined />
                  </a>
                </div>
              </div>
            </Col>
          </Row>
        </div>
      </section>

      {/* 关于我们 */}
      <section id="about" className="about-section">
        <div className="section-content">
          <div className="section-header">
            <span className="section-tag">关于我们</span>
            <h2 className="section-title">专注电解液研发创新</h2>
            <p className="section-subtitle">以科技创新驱动新能源行业发展</p>
          </div>

          <Row gutter={[48, 48]} align="middle">
            <Col xs={24} md={12}>
              <div className="about-content">
                <p className="about-text">
                  Molyte是专注于电池电解液研发的高科技平台，致力于通过先进的分子动力学模拟
                  和人工智能技术，加速新型电解液的研发进程。
                </p>
                <p className="about-text">
                  我们的团队由来自材料科学、计算化学、人工智能等领域的专家组成，
                  拥有丰富的电解液研发经验和技术积累。
                </p>
                <div className="about-features">
                  <div className="about-feature">
                    <ThunderboltOutlined className="about-feature-icon" />
                    <span>高效计算引擎</span>
                  </div>
                  <div className="about-feature">
                    <SafetyCertificateOutlined className="about-feature-icon" />
                    <span>数据安全保障</span>
                  </div>
                  <div className="about-feature">
                    <ExperimentOutlined className="about-feature-icon" />
                    <span>专业技术支持</span>
                  </div>
                </div>
              </div>
            </Col>
            <Col xs={24} md={12}>
              <div className="about-stats-grid">
                <div className="about-stat-card">
                  <div className="about-stat-number">5+</div>
                  <div className="about-stat-label">年技术积累</div>
                </div>
                <div className="about-stat-card">
                  <div className="about-stat-number">100+</div>
                  <div className="about-stat-label">企业用户</div>
                </div>
                <div className="about-stat-card">
                  <div className="about-stat-number">10000+</div>
                  <div className="about-stat-label">模拟任务</div>
                </div>
                <div className="about-stat-card">
                  <div className="about-stat-number">24/7</div>
                  <div className="about-stat-label">技术支持</div>
                </div>
              </div>
            </Col>
          </Row>
        </div>
      </section>

      {/* 页脚 */}
      <footer className="home-footer">
        <div className="footer-content">
          <div className="footer-main">
            <Row gutter={[48, 48]}>
              <Col xs={24} md={8}>
                <div className="footer-brand">
                  <div className="footer-logo">
                    <ThunderboltOutlined />
                    <span>Molyte</span>
                  </div>
                  <p className="footer-desc">
                    新一代电解液研发模拟平台，以科技创新驱动新能源行业发展
                  </p>
                </div>
              </Col>
              <Col xs={24} md={8}>
                <div className="footer-links">
                  <h4 className="footer-title">快速链接</h4>
                  <a onClick={() => scrollToSection('features')}>功能介绍</a>
                  <a onClick={() => navigate('/research')}>数据库</a>
                  <a onClick={() => scrollToSection('news')}>新闻动态</a>
                  <a onClick={() => scrollToSection('about')}>关于我们</a>
                </div>
              </Col>
              <Col xs={24} md={8}>
                <div className="footer-contact">
                  <h4 className="footer-title">联系我们</h4>
                  <p><MailOutlined /> contact@molyte.com</p>
                  <p><PhoneOutlined /> +86 400-XXX-XXXX</p>
                  <p><EnvironmentOutlined /> 中国·北京</p>
                </div>
              </Col>
            </Row>
          </div>
          <div className="footer-bottom">
            <p>© 2025 Molyte. All rights reserved. | 隐私政策 | 服务条款</p>
          </div>
        </div>
      </footer>
    </div>
  );
}

