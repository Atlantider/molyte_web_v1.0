#!/bin/bash
# 代码同步脚本
# 用途：在校园网集群和腾讯云之间同步代码

set -e  # 遇到错误立即退出

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 配置
TENCENT_SERVER="root@42.193.107.130"
TENCENT_PATH="/root/molyte_web"  # 修改为实际路径
LOCAL_PATH="/public/home/xiaoji/molyte_web"

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  Molyte Web 代码同步工具${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# 显示使用说明
show_usage() {
    echo "用法: $0 [选项]"
    echo ""
    echo "选项:"
    echo "  push      - 推送本地代码到 GitHub（校园网 → GitHub）"
    echo "  pull      - 从 GitHub 拉取代码到本地（GitHub → 校园网）"
    echo "  deploy    - 部署到腾讯云（GitHub → 腾讯云）"
    echo "  status    - 检查本地和远程状态"
    echo "  sync      - 完整同步（本地 → GitHub → 腾讯云）"
    echo ""
    echo "示例:"
    echo "  $0 push     # 提交本地修改到 GitHub"
    echo "  $0 deploy   # 部署最新代码到腾讯云"
    echo "  $0 sync     # 完整同步流程"
}

# 检查 Git 状态
check_git_status() {
    echo -e "${YELLOW}📊 检查 Git 状态...${NC}"
    
    cd "$LOCAL_PATH"
    
    if [[ -n $(git status -s) ]]; then
        echo -e "${RED}⚠️  有未提交的修改:${NC}"
        git status -s
        return 1
    else
        echo -e "${GREEN}✅ 工作目录干净${NC}"
        return 0
    fi
}

# 推送到 GitHub
push_to_github() {
    echo -e "${YELLOW}📤 推送代码到 GitHub...${NC}"
    
    cd "$LOCAL_PATH"
    
    # 检查是否有修改
    if [[ -z $(git status -s) ]]; then
        echo -e "${GREEN}✅ 没有需要提交的修改${NC}"
        return 0
    fi
    
    # 显示修改
    echo -e "${BLUE}修改的文件:${NC}"
    git status -s
    echo ""
    
    # 确认
    read -p "是否提交这些修改? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo -e "${YELLOW}❌ 取消提交${NC}"
        return 1
    fi
    
    # 输入提交信息
    read -p "请输入提交信息: " commit_msg
    
    # 提交
    git add .
    git commit -m "$commit_msg"
    git push origin main
    
    echo -e "${GREEN}✅ 代码已推送到 GitHub${NC}"
}

# 从 GitHub 拉取
pull_from_github() {
    echo -e "${YELLOW}📥 从 GitHub 拉取代码...${NC}"
    
    cd "$LOCAL_PATH"
    
    # 检查是否有未提交的修改
    if [[ -n $(git status -s) ]]; then
        echo -e "${RED}⚠️  有未提交的修改，请先提交或暂存${NC}"
        git status -s
        return 1
    fi
    
    # 拉取
    git pull origin main
    
    echo -e "${GREEN}✅ 代码已更新${NC}"
}

# 部署到腾讯云
deploy_to_tencent() {
    echo -e "${YELLOW}🚀 部署到腾讯云...${NC}"
    
    # 确认
    read -p "是否部署到腾讯云服务器? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo -e "${YELLOW}❌ 取消部署${NC}"
        return 1
    fi
    
    # SSH 到腾讯云执行部署
    ssh "$TENCENT_SERVER" << 'EOF'
        set -e
        
        echo "📂 进入项目目录..."
        cd /root/molyte_web  # 修改为实际路径
        
        echo "📥 拉取最新代码..."
        git pull origin main
        
        echo "🔄 重启后端服务..."
        # 根据实际情况选择重启方式
        if systemctl is-active --quiet molyte-backend; then
            sudo systemctl restart molyte-backend
            echo "✅ 后端服务已重启 (systemd)"
        elif command -v pm2 &> /dev/null; then
            pm2 restart molyte-backend
            echo "✅ 后端服务已重启 (PM2)"
        else
            echo "⚠️  请手动重启后端服务"
        fi
        
        echo "✅ 部署完成"
EOF
    
    echo -e "${GREEN}✅ 腾讯云部署完成${NC}"
}

# 检查状态
check_status() {
    echo -e "${YELLOW}📊 检查代码状态...${NC}"
    echo ""
    
    cd "$LOCAL_PATH"
    
    echo -e "${BLUE}=== 本地状态 ===${NC}"
    git status -s
    echo ""
    
    echo -e "${BLUE}=== 本地分支 ===${NC}"
    git branch -vv
    echo ""
    
    echo -e "${BLUE}=== 远程状态 ===${NC}"
    git fetch origin
    LOCAL=$(git rev-parse @)
    REMOTE=$(git rev-parse @{u})
    BASE=$(git merge-base @ @{u})
    
    if [ $LOCAL = $REMOTE ]; then
        echo -e "${GREEN}✅ 本地和远程同步${NC}"
    elif [ $LOCAL = $BASE ]; then
        echo -e "${YELLOW}⚠️  远程有新提交，需要 pull${NC}"
    elif [ $REMOTE = $BASE ]; then
        echo -e "${YELLOW}⚠️  本地有新提交，需要 push${NC}"
    else
        echo -e "${RED}⚠️  本地和远程分叉，需要合并${NC}"
    fi
    echo ""
    
    echo -e "${BLUE}=== 最近 5 次提交 ===${NC}"
    git log --oneline -5
}

# 完整同步
full_sync() {
    echo -e "${YELLOW}🔄 开始完整同步...${NC}"
    echo ""
    
    # 1. 检查状态
    check_status
    echo ""
    
    # 2. 推送到 GitHub
    if [[ -n $(git status -s) ]]; then
        push_to_github || return 1
    else
        echo -e "${GREEN}✅ 本地无修改，跳过推送${NC}"
    fi
    echo ""
    
    # 3. 部署到腾讯云
    deploy_to_tencent
    echo ""
    
    echo -e "${GREEN}✅ 完整同步完成！${NC}"
}

# 主逻辑
case "${1:-}" in
    push)
        push_to_github
        ;;
    pull)
        pull_from_github
        ;;
    deploy)
        deploy_to_tencent
        ;;
    status)
        check_status
        ;;
    sync)
        full_sync
        ;;
    *)
        show_usage
        exit 1
        ;;
esac

