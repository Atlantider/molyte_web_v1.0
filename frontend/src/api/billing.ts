/**
 * 计费和充值 API
 */
import client from './client';

// ============ 类型定义 ============

export interface BalanceInfo {
  balance: number;
  frozen: number;
  debt: number;
  available: number;
  price_per_hour: number;
}

export interface RechargeOrder {
  id: number;
  order_no: string;
  amount: number;
  cpu_hours: number;
  price_per_hour: number;
  payment_method: string;
  payment_status: string;
  created_at: string;
  paid_at?: string;
}

export interface Transaction {
  id: number;
  type: string;
  amount: number;
  balance_before: number;
  balance_after: number;
  description?: string;
  created_at: string;
}

export interface PricingConfig {
  cpu_hour_price: number;
  min_recharge_amount: number;
  max_debt_cpu_hours: number;
}

export interface CreateOrderRequest {
  amount: number;
  payment_method?: string;
}

// ============ 用户端 API ============

/**
 * 获取当前用户余额信息
 */
export async function getBalance(): Promise<BalanceInfo> {
  const response = await client.get('/billing/balance');
  return response.data;
}

/**
 * 检查是否可以提交任务
 */
export async function checkCanSubmit(): Promise<{ can_submit: boolean; reason: string }> {
  const response = await client.get('/billing/can-submit');
  return response.data;
}

/**
 * 创建充值订单
 */
export async function createOrder(data: CreateOrderRequest): Promise<RechargeOrder> {
  const response = await client.post('/billing/orders', data);
  return response.data;
}

/**
 * 获取充值订单列表
 */
export async function getOrders(skip = 0, limit = 20): Promise<RechargeOrder[]> {
  const response = await client.get('/billing/orders', { params: { skip, limit } });
  return response.data;
}

/**
 * 获取配额变更记录
 */
export async function getTransactions(skip = 0, limit = 20): Promise<Transaction[]> {
  const response = await client.get('/billing/transactions', { params: { skip, limit } });
  return response.data;
}

/**
 * 模拟支付
 */
export async function simulatePayment(orderId: number): Promise<{ success: boolean; message: string; balance: BalanceInfo }> {
  const response = await client.post('/billing/simulate-pay', { order_id: orderId });
  return response.data;
}

// ============ 管理端 API ============

/**
 * 获取定价配置
 */
export async function getPricingConfig(): Promise<PricingConfig> {
  const response = await client.get('/billing/admin/pricing');
  return response.data;
}

/**
 * 更新定价配置
 */
export async function updatePricingConfig(config: PricingConfig): Promise<{ success: boolean; message: string }> {
  const response = await client.put('/billing/admin/pricing', config);
  return response.data;
}

/**
 * 管理员调整用户余额
 */
export async function adminAdjustBalance(userId: number, amount: number, reason: string): Promise<{ success: boolean; message: string }> {
  const response = await client.post('/billing/admin/adjust', { user_id: userId, amount, reason });
  return response.data;
}

/**
 * 获取所有充值订单（管理员）
 */
export async function getAllOrders(skip = 0, limit = 50, status?: string): Promise<RechargeOrder[]> {
  const response = await client.get('/billing/admin/orders', { params: { skip, limit, status } });
  return response.data;
}

