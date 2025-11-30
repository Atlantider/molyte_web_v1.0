#!/bin/bash

# 完整停止脚本 - 停止所有 Molyte Web 服务

echo ""
echo "╔════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                            ║"
echo "║              ⏹️  Molyte Web 完整停止 - 停止所有服务                        ║"
echo "║                                                                            ║"
echo "╚════════════════════════════════════════════════════════════════════════════╝"
echo ""

# 1. 停止 Vite 前端
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "1️⃣  停止 Vite 前端..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if ps aux | grep -v grep | grep "vite" > /dev/null; then
    pkill -f "vite"
    sleep 1
    if ps aux | grep -v grep | grep "vite" > /dev/null; then
        echo "⚠️  Vite 进程仍在运行，强制杀死..."
        pkill -9 -f "vite"
    fi
    echo "✅ Vite 前端已停止"
else
    echo "ℹ️  Vite 前端未运行"
fi
echo ""

# 2. 停止 FastAPI 后端
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "2️⃣  停止 FastAPI 后端..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if ps aux | grep -v grep | grep "uvicorn.*app.main:app" > /dev/null; then
    pkill -f "uvicorn.*app.main:app"
    sleep 1
    if ps aux | grep -v grep | grep "uvicorn.*app.main:app" > /dev/null; then
        echo "⚠️  FastAPI 进程仍在运行，强制杀死..."
        pkill -9 -f "uvicorn.*app.main:app"
    fi
    echo "✅ FastAPI 后端已停止"
else
    echo "ℹ️  FastAPI 后端未运行"
fi
echo ""

# 3. 停止 Celery Worker
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "3️⃣  停止 Celery Worker..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if ps aux | grep -v grep | grep "celery.*worker" > /dev/null; then
    pkill -f "celery.*worker"
    sleep 1
    if ps aux | grep -v grep | grep "celery.*worker" > /dev/null; then
        echo "⚠️  Celery Worker 进程仍在运行，强制杀死..."
        pkill -9 -f "celery.*worker"
    fi
    echo "✅ Celery Worker 已停止"
else
    echo "ℹ️  Celery Worker 未运行"
fi
echo ""

# 4. 停止 Celery Beat
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "4️⃣  停止 Celery Beat..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if ps aux | grep -v grep | grep "celery.*beat" > /dev/null; then
    pkill -f "celery.*beat"
    sleep 1
    if ps aux | grep -v grep | grep "celery.*beat" > /dev/null; then
        echo "⚠️  Celery Beat 进程仍在运行，强制杀死..."
        pkill -9 -f "celery.*beat"
    fi
    echo "✅ Celery Beat 已停止"
else
    echo "ℹ️  Celery Beat 未运行"
fi
echo ""

# 5. 停止 Redis
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "5️⃣  停止 Redis..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if ps aux | grep -v grep | grep redis-server > /dev/null; then
    /public/home/xiaoji/software/redis/src/redis-cli shutdown
    sleep 1
    if ps aux | grep -v grep | grep redis-server > /dev/null; then
        echo "⚠️  Redis 进程仍在运行，强制杀死..."
        pkill -9 -f redis-server
    fi
    echo "✅ Redis 已停止"
else
    echo "ℹ️  Redis 未运行"
fi
echo ""

# 显示停止完成信息
echo "╔════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                            ║"
echo "║                   ✅ 所有服务已停止！                                      ║"
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

echo "💡 提示："
echo "   要重新启动所有服务，运行:"
echo "   bash /public/home/xiaoji/molyte_web/start_complete.sh"
echo ""

