#!/bin/bash
#
# 启动轮询 Worker
#
# 用法:
#   bash deployment/start_polling_worker.sh                    # 使用默认配置（阿里云）
#   bash deployment/start_polling_worker.sh tencent            # 使用腾讯云配置
#   bash deployment/start_polling_worker.sh /path/to/config    # 使用自定义配置
#

set -e

echo "=========================================="
echo "  启动混合云轮询 Worker"
echo "=========================================="
echo ""

# 项目根目录
PROJECT_ROOT="/public/home/xiaoji/molyte_web"
cd "$PROJECT_ROOT"

# 确定配置文件
if [ -z "$1" ]; then
    # 默认配置（阿里云）
    CONFIG_FILE="deployment/polling_worker_config.yaml"
    echo "使用默认配置（阿里云）: $CONFIG_FILE"
elif [ "$1" == "tencent" ]; then
    # 腾讯云配置
    CONFIG_FILE="deployment/polling_worker_config_tencent.yaml"
    echo "使用腾讯云配置: $CONFIG_FILE"
else
    # 自定义配置文件
    CONFIG_FILE="$1"
    echo "使用自定义配置: $CONFIG_FILE"
fi
echo ""

# 激活 Conda 环境
echo "1️⃣  激活 Conda 环境..."
source /public/software/anaconda3/bin/activate molyte
echo "✅ Conda 环境已激活"
echo ""

# 检查配置文件
echo "2️⃣  检查配置文件..."
if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ 配置文件不存在: $CONFIG_FILE"
    echo ""
    echo "可用的配置文件模板："
    echo "  - deployment/polling_worker_config.yaml (阿里云)"
    echo "  - deployment/polling_worker_config_tencent.yaml (腾讯云)"
    echo ""
    echo "请先编辑配置文件，填入正确的 API Token 和对象存储配置"
    exit 1
fi
echo "✅ 配置文件存在: $CONFIG_FILE"
echo ""

# 检查依赖
echo "3️⃣  检查 Python 依赖..."

# 检查 PyYAML
python -c "import yaml" 2>/dev/null || {
    echo "⚠️  PyYAML 未安装，正在安装..."
    pip install pyyaml
}

# 检查对象存储 SDK（根据配置文件判断）
if grep -q "^cos:" "$CONFIG_FILE"; then
    echo "检测到腾讯云 COS 配置..."
    python -c "from qcloud_cos import CosConfig, CosS3Client" 2>/dev/null || {
        echo "⚠️  腾讯云 COS SDK 未安装，正在安装..."
        pip install cos-python-sdk-v5
    }
elif grep -q "^oss:" "$CONFIG_FILE"; then
    echo "检测到阿里云 OSS 配置..."
    python -c "import oss2" 2>/dev/null || {
        echo "⚠️  阿里云 OSS SDK 未安装，正在安装..."
        pip install oss2
    }
fi

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
    --config "$CONFIG_FILE" \
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

