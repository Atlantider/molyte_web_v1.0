"""
计费和充值相关模型
"""
from sqlalchemy import Column, Integer, String, Float, DateTime, Boolean, ForeignKey, Text, Enum as SQLEnum
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
from app.database import Base
import enum


class PaymentMethod(str, enum.Enum):
    """支付方式"""
    WECHAT = "wechat"
    ALIPAY = "alipay"
    ADMIN = "admin"  # 管理员手动充值
    SIMULATED = "simulated"  # 模拟支付（测试用）


class PaymentStatus(str, enum.Enum):
    """支付状态"""
    PENDING = "pending"      # 待支付
    PAID = "paid"            # 已支付
    FAILED = "failed"        # 支付失败
    CANCELLED = "cancelled"  # 已取消
    REFUNDED = "refunded"    # 已退款


class TransactionType(str, enum.Enum):
    """配额变更类型"""
    RECHARGE = "recharge"        # 充值
    CONSUME = "consume"          # 消费（任务扣费）
    REFUND = "refund"            # 退款
    ADMIN_ADJUST = "admin_adjust"  # 管理员调整
    DEBT_REPAY = "debt_repay"    # 偿还欠费


class SystemConfig(Base):
    """系统配置表"""
    __tablename__ = "system_configs"

    id = Column(Integer, primary_key=True, index=True)
    key = Column(String(100), unique=True, nullable=False, index=True)
    value = Column(Text, nullable=False)
    description = Column(String(500))
    
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())
    updated_by = Column(Integer, ForeignKey("users.id", ondelete="SET NULL"), nullable=True)

    def __repr__(self):
        return f"<SystemConfig(key={self.key}, value={self.value})>"


class RechargeOrder(Base):
    """充值订单表"""
    __tablename__ = "recharge_orders"

    id = Column(Integer, primary_key=True, index=True)
    order_no = Column(String(64), unique=True, nullable=False, index=True)  # 订单号
    user_id = Column(Integer, ForeignKey("users.id", ondelete="CASCADE"), nullable=False)
    
    # 金额信息
    amount = Column(Float, nullable=False)  # 充值金额（元）
    cpu_hours = Column(Float, nullable=False)  # 对应机时
    price_per_hour = Column(Float, nullable=False)  # 下单时的单价（元/核时）
    
    # 支付信息 - 使用 String 存储枚举值以避免 PostgreSQL 枚举大小写问题
    payment_method = Column(String(20), nullable=False)
    payment_status = Column(String(20), default="pending", nullable=False)
    transaction_id = Column(String(128))  # 第三方支付交易号
    
    # 时间戳
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    paid_at = Column(DateTime(timezone=True))
    expired_at = Column(DateTime(timezone=True))  # 订单过期时间
    
    # 备注
    remark = Column(String(500))
    
    # 关系
    user = relationship("User", back_populates="recharge_orders")

    def __repr__(self):
        return f"<RechargeOrder(order_no={self.order_no}, amount={self.amount}, status={self.payment_status})>"


class QuotaTransaction(Base):
    """配额变更流水表"""
    __tablename__ = "quota_transactions"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id", ondelete="CASCADE"), nullable=False)
    
    # 变更信息 - 使用 String 存储枚举值
    type = Column(String(20), nullable=False)
    amount = Column(Float, nullable=False)  # 变更的机时数（正为增加，负为消耗）
    balance_before = Column(Float, nullable=False)  # 变更前余额
    balance_after = Column(Float, nullable=False)  # 变更后余额
    
    # 关联信息
    reference_id = Column(Integer)  # 关联的订单ID或任务ID
    reference_type = Column(String(50))  # order / job
    
    description = Column(String(500))
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    
    # 关系
    user = relationship("User", back_populates="quota_transactions")

    def __repr__(self):
        return f"<QuotaTransaction(user_id={self.user_id}, type={self.type}, amount={self.amount})>"

