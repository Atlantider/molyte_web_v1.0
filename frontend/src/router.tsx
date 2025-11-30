/**
 * 路由配置
 */
import { createBrowserRouter } from 'react-router-dom';
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
      {
        path: 'jobs',
        element: <Jobs />,
      },
      {
        path: 'jobs/create/:systemId',
        element: <JobCreate />,
      },
      {
        path: 'jobs/:jobId/submit',
        element: <JobSubmit />,
      },
      {
        path: 'jobs/:id/detail',
        element: <JobDetail />,
      },
      // QC量子化学计算路由
      {
        path: 'qc-jobs',
        element: <QCJobs />,
      },
      {
        path: 'qc-jobs/:id',
        element: <QCJobDetail />,
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

