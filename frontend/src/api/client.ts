/**
 * Axios HTTP 客户端配置
 */
import axios from 'axios';
import type { AxiosError, AxiosResponse } from 'axios';

// 创建 axios 实例
const client = axios.create({
  baseURL: '/api/v1',
  timeout: 30000,
  headers: {
    'Content-Type': 'application/json',
  },
});

// 请求拦截器 - 自动添加 token
client.interceptors.request.use(
  (config) => {
    const token = localStorage.getItem('access_token');
    if (token) {
      config.headers.Authorization = `Bearer ${token}`;
    }
    return config;
  },
  (error) => {
    return Promise.reject(error);
  }
);

// 响应拦截器 - 统一错误处理
client.interceptors.response.use(
  (response: AxiosResponse) => {
    return response;
  },
  (error: AxiosError) => {
    if (error.response) {
      // 服务器返回错误状态码
      const status = error.response.status;
      
      if (status === 401) {
        // 未授权，清除 token 并跳转到登录页
        localStorage.removeItem('access_token');
        window.location.href = '/login';
      } else if (status === 403) {
        console.error('权限不足');
      } else if (status === 404) {
        console.error('资源不存在');
      } else if (status >= 500) {
        console.error('服务器错误');
      }
    } else if (error.request) {
      // 请求已发送但没有收到响应
      console.error('网络错误，请检查网络连接');
    } else {
      // 其他错误
      console.error('请求失败:', error.message);
    }
    
    return Promise.reject(error);
  }
);

export default client;

