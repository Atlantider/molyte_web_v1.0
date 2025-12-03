/**
 * 主应用组件
 */
import { useEffect } from 'react';
import { RouterProvider } from 'react-router-dom';
import { ConfigProvider } from 'antd';
import zhCN from 'antd/locale/zh_CN';
import router from './router';
import { useThemeStore, getAntdTheme } from './stores/themeStore';

function App() {
  const { mode } = useThemeStore();
  const themeConfig = getAntdTheme(mode);

  // 更新 body 背景色和主题属性
  useEffect(() => {
    document.body.style.backgroundColor = mode === 'dark' ? '#141414' : '#F5F7FA';
    document.body.style.transition = 'background-color 0.3s';
    document.documentElement.setAttribute('data-theme', mode);
  }, [mode]);

  return (
    <ConfigProvider locale={zhCN} theme={themeConfig}>
      <RouterProvider router={router} />
    </ConfigProvider>
  );
}

export default App;

