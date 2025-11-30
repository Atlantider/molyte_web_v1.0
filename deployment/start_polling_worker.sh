#!/bin/bash
#
# 启动轮询 Worker
#
# 用法: bash deployment/start_polling_worker.sh
#

set -e

echo "=========================================="
echo "  启动混合云轮询 Worker"
echo "=========================================="
echo ""

# 项目根目录
PROJECT_ROOT="/public/home/xiaoji/molyte_web"
cd "$PROJECT_ROOT"

# 激活 Conda 环境
echo "1️⃣  激活 Conda 环境..."
source /public/software/anaconda3/bin/activate molyte
echo "✅ Conda 环境已激活"
echo ""

# 检查配置文件
echo "2️⃣  检查配置文件..."
if [ ! -f "deployment/polling_worker_config.yaml" ]; then
    echo "❌ 配置文件不存在: deployment/polling_worker_config.yaml"
    echo "请先复制并编辑配置文件:"
    echo "  cp deployment/polling_worker_config.yaml.example deployment/polling_worker_config.yaml"
    exit 1
fi
echo "✅ 配置文件存在"
echo ""

# 检查依赖
echo "3️⃣  检查 Python 依赖..."
python -c "import oss2" 2>/dev/null || {
    echo "⚠️  OSS SDK 未安装，正在安装..."
    pip install oss2
}
python -c "import yaml" 2>/dev/null || {
    echo "⚠️  PyYAML 未安装，正在安装..."
    pip install pyyaml
}
echo "✅ 依赖检查完成"
echo ""

# 停止旧的 Worker
echo "4️⃣  停止旧的 Worker 进程..."
pkill -f "polling_worker.py" || echo "   (没有运行中的 Worker)"
sleep 2
echo "✅ 旧进程已停止"
echo ""

# 启动 Worker
echo "5️⃣  启动 Worker..."
nohup python deployment/polling_worker.py \
    --config deployment/polling_worker_config.yaml \
    > /tmp/polling_worker_stdout.log 2>&1 &

WORKER_PID=$!
echo "✅ Worker 已启动 (PID: $WORKER_PID)"
echo ""

# 等待启动
echo "6️⃣  等待 Worker 启动..."
sleep 3

# 检查进程
if ps -p $WORKER_PID > /dev/null; then
    echo "✅ Worker 运行正常"
    echo ""
    echo "=========================================="
    echo "  Worker 启动成功！"
    echo "=========================================="
    echo ""
    echo "📝 日志文件:"
    echo "   - Worker 日志: /tmp/polling_worker.log"
    echo "   - 标准输出: /tmp/polling_worker_stdout.log"
    echo ""
    echo "📊 查看日志:"
    echo "   tail -f /tmp/polling_worker.log"
    echo ""
    echo "🛑 停止 Worker:"
    echo "   bash deployment/stop_polling_worker.sh"
    echo ""
else
    echo "❌ Worker 启动失败"
    echo ""
    echo "请检查日志:"
    echo "   tail -50 /tmp/polling_worker_stdout.log"
    exit 1
fi

