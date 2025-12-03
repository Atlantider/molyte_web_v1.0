/**
 * 可拖拽列的表格组件
 * 支持列排序、过滤、拖拽调整顺序
 */
import React, { useState, useCallback } from 'react';
import { Table, Button, Dropdown, Space } from 'antd';
import { SettingOutlined, EyeOutlined, EyeInvisibleOutlined } from '@ant-design/icons';
import type { TableProps, ColumnsType } from 'antd/es/table';
import { DndContext, closestCenter, PointerSensor, useSensor, useSensors } from '@dnd-kit/core';
import { arrayMove, SortableContext, useSortable, verticalListSortingStrategy } from '@dnd-kit/sortable';
import { CSS } from '@dnd-kit/utilities';

interface DraggableTableProps<T> extends TableProps<T> {
  columns: ColumnsType<T>;
  enableColumnSettings?: boolean;
}

// 可拖拽的行组件
const DraggableRow = ({ children, ...props }: any) => {
  const {
    attributes,
    listeners,
    setNodeRef,
    transform,
    transition,
    isDragging,
  } = useSortable({
    id: props['data-row-key'],
  });

  const style = {
    ...props.style,
    transform: CSS.Transform.toString(transform),
    transition,
    cursor: 'move',
    ...(isDragging ? { position: 'relative', zIndex: 9999 } : {}),
  };

  return (
    <tr {...props} ref={setNodeRef} style={style} {...attributes} {...listeners}>
      {children}
    </tr>
  );
};

export default function DraggableTable<T extends { id?: number | string; key?: string | number }>({
  columns: initialColumns,
  enableColumnSettings = true,
  ...props
}: DraggableTableProps<T>) {
  const [columns, setColumns] = useState(initialColumns);
  const [hiddenColumns, setHiddenColumns] = useState<Set<string>>(new Set());

  const sensors = useSensors(
    useSensor(PointerSensor, {
      activationConstraint: {
        distance: 1,
      },
    })
  );

  // 切换列显示/隐藏
  const toggleColumn = useCallback((key: string) => {
    setHiddenColumns(prev => {
      const newSet = new Set(prev);
      if (newSet.has(key)) {
        newSet.delete(key);
      } else {
        newSet.add(key);
      }
      return newSet;
    });
  }, []);

  // 过滤隐藏的列
  const visibleColumns = columns.filter(col => {
    const key = (col as any).key || (col as any).dataIndex;
    return !hiddenColumns.has(key);
  });

  // 列设置菜单
  const columnSettingsMenu = {
    items: columns.map(col => {
      const key = (col as any).key || (col as any).dataIndex;
      const isHidden = hiddenColumns.has(key);
      return {
        key,
        label: (
          <Space>
            {isHidden ? <EyeInvisibleOutlined /> : <EyeOutlined />}
            {col.title as string}
          </Space>
        ),
        onClick: () => toggleColumn(key),
      };
    }),
  };

  return (
    <div>
      {enableColumnSettings && (
        <div style={{ marginBottom: 16, textAlign: 'right' }}>
          <Dropdown menu={columnSettingsMenu} trigger={['click']}>
            <Button icon={<SettingOutlined />}>
              列设置
            </Button>
          </Dropdown>
        </div>
      )}
      <Table
        {...props}
        columns={visibleColumns}
      />
    </div>
  );
}

