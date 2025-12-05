/**
 * 路由配置
 */
import { createBrowserRouter, Navigate } from 'react-router-dom';
import Home from './pages/Home';
import Login from './pages/Login';
import Dashboard from './pages/Dashboard';
import Projects from './pages/Projects';
import ProjectDetail from './pages/ProjectDetail';
import Electrolytes from './pages/Electrolytes';
import Jobs from './pages/Jobs';
import JobCreate from './pages/JobCreate';
import JobSubmit from './pages/JobSubmit';
import JobDetail from './pages/JobDetail';
import QCJobs from './pages/QCJobs';
import QCJobDetail from './pages/QCJobDetail';
import PostProcessAnalysis from './pages/PostProcessAnalysis';
import PostProcessDetail from './pages/PostProcessDetail';
import ChangePassword from './pages/ChangePassword';
import Profile from './pages/Profile';
import Recharge from './pages/Recharge';
import Layout from './components/Layout';
import PrivateRoute from './components/PrivateRoute';
import AdminDashboard from './pages/admin/AdminDashboard';
import UserManagement from './pages/admin/UserManagement';
import UserDetail from './pages/admin/UserDetail';
import AuditLogs from './pages/admin/AuditLogs';
import DataVisibilityAdmin from './pages/admin/DataVisibilityAdmin';
import Research from './pages/Research';
import DataVisibilityManager from './components/DataVisibilityManager';
import PublicResearch from './pages/PublicResearch';
import PublicResultDetail from './pages/PublicResultDetail';
import Guide from './pages/Guide';


const router = createBrowserRouter([
  {
    path: '/',
    element: <Home />,
  },
  {
    path: '/guide',
    element: <Guide />,
  },
  {
    path: '/research',
    element: <PublicResearch />,
  },
  {
    path: '/research/result/:jobId',
    element: <PublicResultDetail />,
  },
  {
    path: '/login',
    element: <Login />,
  },
  {
    path: '/workspace',
    element: (
      <PrivateRoute>
        <Layout />
      </PrivateRoute>
    ),
    children: [
      {
        path: 'dashboard',
        element: <Dashboard />,
      },
      {
        path: 'projects',
        element: <Projects />,
      },
      {
        path: 'projects/:id',
        element: <ProjectDetail />,
      },
      {
        path: 'electrolytes',
        element: <Electrolytes />,
      },
      // ===== 溶液电解质模块 =====
      {
        path: 'liquid-electrolyte',
        children: [
          // 默认跳转到 MD 模拟
          {
            index: true,
            element: <Navigate to="md" replace />,
          },
          // MD 模拟
          {
            path: 'md',
            element: <Jobs />,
          },
          {
            path: 'md/create/:systemId',
            element: <JobCreate />,
          },
          {
            path: 'md/:jobId/submit',
            element: <JobSubmit />,
          },
          {
            path: 'md/:id',
            element: <JobDetail />,
          },
          // 后处理分析
          {
            path: 'analysis',
            element: <PostProcessAnalysis />,
          },
          {
            path: 'analysis/create',
            element: <PostProcessDetail />,
          },
          {
            path: 'analysis/:id',
            element: <PostProcessDetail />,
          },
          // QC 任务
          {
            path: 'qc',
            element: <QCJobs />,
          },
          {
            path: 'qc/:id',
            element: <QCJobDetail />,
          },
        ],
      },
      // ===== 旧路由重定向（兼容性）=====
      {
        path: 'jobs',
        element: <Navigate to="/workspace/liquid-electrolyte/md" replace />,
      },
      {
        path: 'jobs/create/:systemId',
        element: <JobCreate />,  // 保持兼容
      },
      {
        path: 'jobs/:jobId/submit',
        element: <JobSubmit />,  // 保持兼容
      },
      {
        path: 'jobs/:id/detail',
        element: <Navigate to={`/workspace/liquid-electrolyte/md/${window.location.pathname.split('/')[3]}`} replace />,
      },
      {
        path: 'qc-jobs',
        element: <Navigate to="/workspace/liquid-electrolyte/qc" replace />,
      },
      {
        path: 'qc-jobs/:id',
        element: <QCJobDetail />,  // 保持兼容
      },
      {
        path: 'change-password',
        element: <ChangePassword />,
      },
      {
        path: 'profile',
        element: <Profile />,
      },
      {
        path: 'recharge',
        element: <Recharge />,
      },
      {
        path: 'research',
        element: <Research />,
      },
      // Admin routes
      {
        path: 'admin',
        element: <AdminDashboard />,
      },
      {
        path: 'admin/users',
        element: <UserManagement />,
      },
      {
        path: 'admin/users/:id',
        element: <UserDetail />,
      },
      {
        path: 'admin/logs',
        element: <AuditLogs />,
      },
      {
        path: 'admin/visibility',
        element: <DataVisibilityAdmin />,
      },
      {
        path: 'data-visibility',
        element: <DataVisibilityManager />,
      },
    ],
  },
  {
    path: '*',
    element: <div>404 - 页面不存在</div>,
  },
]);

export default router;

