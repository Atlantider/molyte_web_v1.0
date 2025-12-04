/**
 * 主题状态管理
 * 支持 light/dark 模式切换，并持久化到 localStorage
 */
import { create } from 'zustand';
import { persist } from 'zustand/middleware';
import { theme } from 'antd';

export type ThemeMode = 'light' | 'dark';

interface ThemeState {
  mode: ThemeMode;
  isDark: boolean;
  toggleTheme: () => void;
  setTheme: (mode: ThemeMode) => void;
}

// 柔和的主色调 - 降低饱和度
export const themeTokens = {
  light: {
    colorPrimary: '#5B8DEF',  // 柔和蓝色
    colorSuccess: '#52C41A',
    colorWarning: '#FAAD14',
    colorError: '#F5222D',
    colorBgContainer: '#FFFFFF',
    colorBgLayout: '#F5F7FA',
    colorText: '#333333',
    colorTextSecondary: '#666666',
    colorBorder: '#E8E8E8',
    borderRadius: 8,
  },
  dark: {
    colorPrimary: '#6B9AFF',  // 暗模式主色
    colorSuccess: '#73D13D',
    colorWarning: '#FFC53D',
    colorError: '#FF4D4F',
    colorBgContainer: '#1F1F1F',
    colorBgLayout: '#141414',
    colorText: '#E8E8E8',
    colorTextSecondary: '#A0A0A0',
    colorBorder: '#303030',
    borderRadius: 8,
  },
};

// 获取 Ant Design 主题配置
export const getAntdTheme = (mode: ThemeMode) => {
  const tokens = themeTokens[mode];
  return {
    algorithm: mode === 'dark' ? theme.darkAlgorithm : theme.defaultAlgorithm,
    token: {
      colorPrimary: tokens.colorPrimary,
      colorSuccess: tokens.colorSuccess,
      colorWarning: tokens.colorWarning,
      colorError: tokens.colorError,
      colorBgContainer: tokens.colorBgContainer,
      colorBgLayout: tokens.colorBgLayout,
      colorText: tokens.colorText,
      colorTextSecondary: tokens.colorTextSecondary,
      colorBorder: tokens.colorBorder,
      borderRadius: tokens.borderRadius,
    },
    components: {
      Table: {
        headerBg: mode === 'dark' ? '#1a1a1a' : '#fafafa',
        headerColor: mode === 'dark' ? '#E8E8E8' : '#333333',
        rowHoverBg: mode === 'dark' ? 'rgba(255, 255, 255, 0.08)' : 'rgba(0, 0, 0, 0.04)',
        bodySortBg: mode === 'dark' ? '#1a1a1a' : '#fafafa',
      },
      Card: {
        headerBg: mode === 'dark' ? 'transparent' : '#fafafa',
      },
    },
  };
};

export const useThemeStore = create<ThemeState>()(
  persist(
    (set) => ({
      mode: 'light',
      isDark: false,
      toggleTheme: () => set((state) => ({
        mode: state.mode === 'light' ? 'dark' : 'light',
        isDark: state.mode === 'light',
      })),
      setTheme: (mode) => set({ mode, isDark: mode === 'dark' }),
    }),
    {
      name: 'molyte-theme',
    }
  )
);

