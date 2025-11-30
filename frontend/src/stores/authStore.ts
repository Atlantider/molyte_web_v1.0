/**
 * 认证状态管理
 */
import { create } from 'zustand';
import type { User } from '../types/api';
import * as authApi from '../api/auth';

interface AuthState {
  user: User | null;
  token: string | null;
  isAuthenticated: boolean;
  isLoading: boolean;
  error: string | null;

  // Actions
  login: (username: string, password: string) => Promise<void>;
  register: (
    email: string,
    username: string,
    password: string,
    user_type?: string,
    organization?: string,
    department?: string
  ) => Promise<void>;
  logout: () => void;
  loadUser: () => Promise<void>;
  clearError: () => void;
}

export const useAuthStore = create<AuthState>((set) => ({
  user: null,
  token: localStorage.getItem('access_token'),
  isAuthenticated: !!localStorage.getItem('access_token'),
  isLoading: false,
  error: null,

  login: async (username: string, password: string) => {
    set({ isLoading: true, error: null });
    try {
      const response = await authApi.login({ username, password });
      const token = response.access_token;
      
      // 保存 token
      localStorage.setItem('access_token', token);
      
      // 获取用户信息
      const user = await authApi.getCurrentUser();
      
      set({
        user,
        token,
        isAuthenticated: true,
        isLoading: false,
      });
    } catch (error: any) {
      const errorMessage = error.response?.data?.detail || '登录失败';
      set({
        error: errorMessage,
        isLoading: false,
        isAuthenticated: false,
      });
      throw error;
    }
  },

  register: async (
    email: string,
    username: string,
    password: string,
    user_type: string = 'STUDENT',
    organization: string = '',
    department?: string
  ) => {
    set({ isLoading: true, error: null });
    try {
      await authApi.register({
        email,
        username,
        password,
        user_type,
        organization,
        department
      });

      // 注册成功后自动登录（使用 username 而不是 email）
      await useAuthStore.getState().login(username, password);
    } catch (error: any) {
      const errorMessage = error.response?.data?.detail || '注册失败';
      set({
        error: errorMessage,
        isLoading: false,
      });
      throw error;
    }
  },

  logout: () => {
    authApi.logout();
    set({
      user: null,
      token: null,
      isAuthenticated: false,
    });
  },

  loadUser: async () => {
    const token = localStorage.getItem('access_token');
    if (!token) {
      set({ isAuthenticated: false });
      return;
    }

    set({ isLoading: true });
    try {
      const user = await authApi.getCurrentUser();
      set({
        user,
        token,
        isAuthenticated: true,
        isLoading: false,
      });
    } catch (error) {
      // Token 无效，清除
      localStorage.removeItem('access_token');
      set({
        user: null,
        token: null,
        isAuthenticated: false,
        isLoading: false,
      });
    }
  },

  clearError: () => set({ error: null }),
}));

