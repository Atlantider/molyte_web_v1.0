/**
 * 项目卡片组件
 */
import { Card, Space, Button, Popconfirm, Typography } from 'antd';
import { EditOutlined, DeleteOutlined, FolderOutlined, ProjectOutlined, CalendarOutlined } from '@ant-design/icons';
import type { Project } from '../types';
import dayjs from 'dayjs';

const { Text, Paragraph } = Typography;

interface ProjectCardProps {
  project: Project;
  onEdit: (project: Project) => void;
  onDelete: (id: number) => void;
  onView: (project: Project) => void;
}

export default function ProjectCard({ project, onEdit, onDelete, onView }: ProjectCardProps) {
  // 处理卡片点击 - 整个卡片可点击跳转
  const handleCardClick = (e: React.MouseEvent) => {
    // 如果点击的是按钮或其子元素，不触发卡片跳转
    const target = e.target as HTMLElement;
    if (target.closest('button') || target.closest('.ant-popconfirm') || target.closest('.ant-card-actions')) {
      return;
    }
    onView(project);
  };

  return (
    <Card
      hoverable
      onClick={handleCardClick}
      style={{
        height: '100%',
        cursor: 'pointer',
        transition: 'all 0.3s ease',
        border: 'none',
        borderRadius: 12,
        boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
      }}
      styles={{
        body: { padding: '20px' },
      }}
      actions={[
        <Button type="link" icon={<FolderOutlined />} onClick={(e) => { e.stopPropagation(); onView(project); }}>
          查看详情
        </Button>,
        <Button type="link" icon={<EditOutlined />} onClick={(e) => { e.stopPropagation(); onEdit(project); }}>
          编辑
        </Button>,
        <Popconfirm
          title="确定要删除这个项目吗？"
          description="删除后将无法恢复，相关的电解质体系和任务也会被删除。"
          onConfirm={() => onDelete(project.id)}
          okText="确定"
          cancelText="取消"
        >
          <Button type="link" danger icon={<DeleteOutlined />} onClick={(e) => e.stopPropagation()}>
            删除
          </Button>
        </Popconfirm>,
      ]}
    >
      <div style={{ display: 'flex', gap: 16 }}>
        {/* 左侧图标 */}
        <div style={{
          width: 48,
          height: 48,
          borderRadius: 12,
          background: 'linear-gradient(135deg, #4facfe 0%, #00f2fe 100%)',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          flexShrink: 0,
        }}>
          <ProjectOutlined style={{ fontSize: 24, color: '#fff' }} />
        </div>

        {/* 右侧内容 */}
        <div style={{ flex: 1, minWidth: 0 }}>
          <Text strong style={{ fontSize: 16, display: 'block', marginBottom: 8 }}>
            {project.name}
          </Text>
          {project.description ? (
            <Paragraph
              ellipsis={{ rows: 2 }}
              style={{ marginBottom: 12, color: '#666', fontSize: 13 }}
            >
              {project.description}
            </Paragraph>
          ) : (
            <Text type="secondary" style={{ fontSize: 13, display: 'block', marginBottom: 12 }}>
              暂无描述
            </Text>
          )}
          <Space size={4}>
            <CalendarOutlined style={{ color: '#bfbfbf', fontSize: 12 }} />
            <Text type="secondary" style={{ fontSize: 12 }}>
              {dayjs(project.created_at).format('YYYY-MM-DD HH:mm')}
            </Text>
          </Space>
        </div>
      </div>
    </Card>
  );
}

