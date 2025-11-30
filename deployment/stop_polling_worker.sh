#!/bin/bash
#
# 停止轮询 Worker
#
# 用法: bash deployment/stop_polling_worker.sh
#

echo "=========================================="
echo "  停止混合云轮询 Worker"
echo "=========================================="
echo ""

echo "正在停止 Worker 进程..."
pkill -f "polling_worker.py"

sleep 2

# 检查是否还有进程
if pgrep -f "polling_worker.py" > /dev/null; then
    echo "⚠️  Worker 进程仍在运行，强制终止..."
    pkill -9 -f "polling_worker.py"
    sleep 1
fi

if pgrep -f "polling_worker.py" > /dev/null; then
    echo "❌ 无法停止 Worker 进程"
    exit 1
else
    echo "✅ Worker 已停止"
fi

echo ""
echo "=========================================="
echo "  Worker 已成功停止"
echo "=========================================="

