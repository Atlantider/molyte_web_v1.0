#!/bin/bash

# 完整启动脚本 - 启动所有 Molyte Web 服务

set -e

PROJECT_ROOT="/opt/molyte_web_v1.0"
CONDA_ENV="molyte"

echo ""
echo "╔════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                            ║"
echo "║              🚀 Molyte Web 完整启动 - 启动所有服务                         ║"
echo "║                                                                            ║"
echo "╚════════════════════════════════════════════════════════════════════════════╝"
echo ""

# 检测环境并激活
echo "📦 检测 Python 环境..."
if [ -d "$PROJECT_ROOT/backend/venv" ]; then
    # 腾讯云环境：使用 venv
    echo "使用虚拟环境: $PROJECT_ROOT/backend/venv"
    PYTHON_BIN="$PROJECT_ROOT/backend/venv/bin/python"
    UVICORN_BIN="$PROJECT_ROOT/backend/venv/bin/uvicorn"
    CELERY_BIN="$PROJECT_ROOT/backend/venv/bin/celery"
elif [ -f "/public/software/anaconda3/bin/activate" ]; then
    # 校园网环境：使用 conda
    echo "使用 Conda 环境: $CONDA_ENV"
    source /public/software/anaconda3/bin/activate $CONDA_ENV
    PYTHON_BIN="python"
    UVICORN_BIN="uvicorn"
    CELERY_BIN="celery"
else
    echo "❌ 未找到 Python 环境"
    exit 1
fi
echo "✅ Python 环境已准备"
echo ""

# 1. 启动 Redis
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "1️⃣  启动 Redis..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if ps aux | grep -v grep | grep redis-server > /dev/null; then
    echo "✅ Redis 已经在运行"
else
    /public/home/xiaoji/software/redis/src/redis-server /public/home/xiaoji/software/redis/redis.conf --daemonize yes
    sleep 1
    if /public/home/xiaoji/software/redis/src/redis-cli ping > /dev/null 2>&1; then
        echo "✅ Redis 启动成功"
    else
        echo "❌ Redis 启动失败"
        exit 1
    fi
fi
echo ""

# 2. 启动 Celery Worker
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "2️⃣  启动 Celery Worker..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if ps aux | grep -v grep | grep "celery.*worker" > /dev/null; then
    echo "⚠️  Celery Worker 已经在运行，跳过"
else
    cd $PROJECT_ROOT/backend
    nohup celery -A app.celery_app worker \
        --loglevel=info \
        --concurrency=4 \
        --logfile=/tmp/celery_worker.log \
        --pidfile=/tmp/celery_worker.pid \
        > /tmp/celery_worker_stdout.log 2>&1 &
    sleep 3
    if ps aux | grep -v grep | grep "celery.*worker" > /dev/null; then
        echo "✅ Celery Worker 启动成功"
    else
        echo "❌ Celery Worker 启动失败"
        exit 1
    fi
fi
echo ""

# 3. 启动 Celery Beat
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "3️⃣  启动 Celery Beat..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if ps aux | grep -v grep | grep "celery.*beat" > /dev/null; then
    echo "⚠️  Celery Beat 已经在运行，跳过"
else
    cd $PROJECT_ROOT/backend
    nohup celery -A app.celery_app beat \
        --loglevel=info \
        --logfile=/tmp/celery_beat.log \
        --pidfile=/tmp/celery_beat.pid \
        > /tmp/celery_beat_stdout.log 2>&1 &
    sleep 2
    if ps aux | grep -v grep | grep "celery.*beat" > /dev/null; then
        echo "✅ Celery Beat 启动成功"
    else
        echo "❌ Celery Beat 启动失败"
        exit 1
    fi
fi
echo ""

# 4. 启动 FastAPI 后端
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "4️⃣  启动 FastAPI 后端..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if ps aux | grep -v grep | grep "uvicorn.*app.main:app" > /dev/null; then
    echo "⚠️  FastAPI 后端已经在运行，跳过"
else
    cd $PROJECT_ROOT/backend
    nohup /public/software/anaconda3/envs/molyte/bin/python -m uvicorn app.main:app \
        --host 0.0.0.0 \
        --port 8000 \
        --reload \
        > /tmp/fastapi.log 2>&1 &
    sleep 3
    if ps aux | grep -v grep | grep "uvicorn.*app.main:app" > /dev/null; then
        echo "✅ FastAPI 后端启动成功"
    else
        echo "❌ FastAPI 后端启动失败"
        exit 1
    fi
fi
echo ""

# 5. 启动 React 前端
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "5️⃣  启动 React 前端..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if ps aux | grep -v grep | grep "vite" > /dev/null; then
    echo "⚠️  前端已经在运行，跳过"
else
    cd $PROJECT_ROOT/frontend
    
    # 检查 node_modules 是否存在
    if [ ! -d "node_modules" ]; then
        echo "📦 安装前端依赖..."
        npm install > /tmp/npm_install.log 2>&1
        if [ $? -ne 0 ]; then
            echo "❌ 前端依赖安装失败"
            tail -20 /tmp/npm_install.log
            exit 1
        fi
        echo "✅ 前端依赖安装完成"
    fi
    
    nohup npm run dev > /tmp/vite.log 2>&1 &
    sleep 3
    if ps aux | grep -v grep | grep "vite" > /dev/null; then
        echo "✅ React 前端启动成功"
    else
        echo "❌ React 前端启动失败"
        tail -20 /tmp/vite.log
        exit 1
    fi
fi
echo ""

# 显示启动完成信息
echo "╔════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                            ║"
echo "║                   ✅ 所有服务启动完成！                                    ║"
echo "║                                                                            ║"
echo "╚════════════════════════════════════════════════════════════════════════════╝"
echo ""

echo "📊 服务状态："
echo "   Redis:         $(ps aux | grep -v grep | grep redis-server | wc -l) 个进程"
echo "   Celery Worker: $(ps aux | grep -v grep | grep 'celery.*worker' | wc -l) 个进程"
echo "   Celery Beat:   $(ps aux | grep -v grep | grep 'celery.*beat' | wc -l) 个进程"
echo "   FastAPI:       $(ps aux | grep -v grep | grep 'uvicorn.*app.main:app' | wc -l) 个进程"
echo "   Vite:          $(ps aux | grep -v grep | grep 'vite' | wc -l) 个进程"
echo ""

echo "🌐 访问地址："
echo "   前端:          http://localhost:5173"
echo "   后端 API:      http://localhost:8000"
echo "   API 文档:      http://localhost:8000/docs"
echo ""

echo "📝 日志文件："
echo "   FastAPI:       /tmp/fastapi.log"
echo "   Celery Worker: /tmp/celery_worker.log"
echo "   Celery Beat:   /tmp/celery_beat.log"
echo "   Vite:          /tmp/vite.log"
echo ""

echo "🔍 查看日志："
echo "   tail -f /tmp/fastapi.log"
echo "   tail -f /tmp/celery_worker.log"
echo "   tail -f /tmp/celery_beat.log"
echo "   tail -f /tmp/vite.log"
echo ""

echo "⏹️  停止所有服务："
echo "   bash $PROJECT_ROOT/stop_all_services.sh"
echo ""

