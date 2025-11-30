"""
RDF 绘图工具

生成精美的 RDF 和配位数图表
"""

import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from typing import List, Dict, Any
import logging

logger = logging.getLogger(__name__)

# 设置中文字体和样式
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'

# 使用 seaborn 样式
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)


class RDFPlotter:
    """RDF 绘图器"""
    
    def __init__(self):
        # 定义颜色方案
        self.colors = sns.color_palette("husl", 10)
        
    def plot_rdf_combined(
        self,
        rdf_data: List[Dict[str, Any]],
        output_file: Path,
        title: str = "Radial Distribution Functions"
    ):
        """
        绘制组合的 RDF 图（所有 RDF 在一张图上）
        
        Args:
            rdf_data: RDF 数据列表
            output_file: 输出文件路径
            title: 图表标题
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        
        # 绘制 RDF
        for i, data in enumerate(rdf_data):
            label = f"{data['center_label']} → {data['target_label']}"
            color = self.colors[i % len(self.colors)]
            ax1.plot(data['r'], data['g_r'], label=label, color=color, linewidth=2, alpha=0.8)
        
        ax1.set_ylabel('g(r)', fontsize=14, fontweight='bold')
        ax1.set_ylim(bottom=0)
        ax1.legend(loc='upper right', frameon=True, shadow=True, fontsize=10)
        ax1.grid(True, alpha=0.3)
        ax1.set_title(title, fontsize=16, fontweight='bold', pad=20)
        
        # 绘制配位数
        for i, data in enumerate(rdf_data):
            if data.get('coordination_number'):
                label = f"{data['center_label']} → {data['target_label']}"
                color = self.colors[i % len(self.colors)]
                ax2.plot(data['r'], data['coordination_number'], label=label, color=color, linewidth=2, alpha=0.8)
        
        ax2.set_xlabel('r (Å)', fontsize=14, fontweight='bold')
        ax2.set_ylabel('Coordination Number', fontsize=14, fontweight='bold')
        ax2.set_ylim(bottom=0)
        ax2.legend(loc='lower right', frameon=True, shadow=True, fontsize=10)
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Saved combined RDF plot to {output_file}")
    
    def plot_rdf_individual(
        self,
        rdf_data: List[Dict[str, Any]],
        output_dir: Path,
        prefix: str = "rdf"
    ):
        """
        为每个 RDF 对绘制单独的图
        
        Args:
            rdf_data: RDF 数据列表
            output_dir: 输出目录
            prefix: 文件名前缀
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        
        for i, data in enumerate(rdf_data):
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
            
            label = f"{data['center_label']} → {data['target_label']}"
            color = self.colors[i % len(self.colors)]
            
            # 绘制 RDF
            ax1.plot(data['r'], data['g_r'], color=color, linewidth=2.5)
            ax1.fill_between(data['r'], data['g_r'], alpha=0.3, color=color)
            ax1.set_ylabel('g(r)', fontsize=14, fontweight='bold')
            ax1.set_ylim(bottom=0)
            ax1.grid(True, alpha=0.3)
            ax1.set_title(label, fontsize=14, fontweight='bold', pad=15)
            
            # 绘制配位数
            if data.get('coordination_number'):
                ax2.plot(data['r'], data['coordination_number'], color=color, linewidth=2.5)
                ax2.fill_between(data['r'], data['coordination_number'], alpha=0.3, color=color)
                
                # 标注最终配位数
                final_cn = data['coordination_number'][-1]
                ax2.axhline(y=final_cn, color='red', linestyle='--', alpha=0.5, linewidth=1)
                ax2.text(
                    data['r'][-1] * 0.7, final_cn * 1.1,
                    f'CN = {final_cn:.2f}',
                    fontsize=12, fontweight='bold', color='red'
                )
            
            ax2.set_xlabel('r (Å)', fontsize=14, fontweight='bold')
            ax2.set_ylabel('Coordination Number', fontsize=14, fontweight='bold')
            ax2.set_ylim(bottom=0)
            ax2.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            # 生成文件名
            safe_label = label.replace(' ', '_').replace('→', 'to').replace('(', '').replace(')', '')
            output_file = output_dir / f"{prefix}_{i+1}_{safe_label}.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Saved individual RDF plot to {output_file}")

    def plot_rdf_categorized(
        self,
        rdf_data: List[Dict[str, Any]],
        output_file: Path,
        title: str = "Radial Distribution Functions"
    ):
        """
        绘制分类的 RDF 图（典型 RDF 和其他 RDF 分开显示）

        Args:
            rdf_data: RDF 数据列表
            output_file: 输出文件路径
            title: 图表标题
        """
        # 分类 RDF
        typical_rdfs = []
        other_rdfs = []

        for data in rdf_data:
            center = data['center_label']
            target = data['target_label']

            # 判断是否是典型 RDF（阳离子到配位原子）
            is_typical = any(center.startswith(cat) for cat in ['Li_', 'Na_', 'K_', 'Mg_', 'Ca_'])

            if is_typical:
                typical_rdfs.append(data)
            else:
                other_rdfs.append(data)

        # 创建子图
        n_plots = (1 if typical_rdfs else 0) + (1 if other_rdfs else 0)
        if n_plots == 0:
            logger.warning("No RDF data to plot")
            return

        fig = plt.figure(figsize=(12, 6 * n_plots))
        plot_idx = 1

        # 绘制典型 RDF
        if typical_rdfs:
            ax1 = plt.subplot(n_plots, 1, plot_idx)
            ax2 = ax1.twinx()

            for i, data in enumerate(typical_rdfs):
                label = f"{data['center_label']} → {data['target_label']}"
                color = self.colors[i % len(self.colors)]

                # RDF 在左轴
                ax1.plot(data['r'], data['g_r'], label=label, color=color, linewidth=2, alpha=0.8)

                # 配位数在右轴
                if data.get('coordination_number'):
                    ax2.plot(data['r'], data['coordination_number'],
                            linestyle='--', color=color, linewidth=1.5, alpha=0.6)

            ax1.set_xlabel('r (Å)', fontsize=14, fontweight='bold')
            ax1.set_ylabel('g(r)', fontsize=14, fontweight='bold', color='black')
            ax1.set_ylim(bottom=0)
            ax1.legend(loc='upper left', frameon=True, shadow=True, fontsize=10)
            ax1.grid(True, alpha=0.3)
            ax1.set_title('Typical RDFs (Cation-Coordination)', fontsize=14, fontweight='bold', pad=15)

            ax2.set_ylabel('Coordination Number', fontsize=14, fontweight='bold', color='gray')
            ax2.set_ylim(bottom=0)
            ax2.tick_params(axis='y', labelcolor='gray')

            plot_idx += 1

        # 绘制其他 RDF
        if other_rdfs:
            ax3 = plt.subplot(n_plots, 1, plot_idx)

            for i, data in enumerate(other_rdfs):
                label = f"{data['center_label']} → {data['target_label']}"
                color = self.colors[(i + len(typical_rdfs)) % len(self.colors)]
                ax3.plot(data['r'], data['g_r'], label=label, color=color, linewidth=2, alpha=0.8)

            ax3.set_xlabel('r (Å)', fontsize=14, fontweight='bold')
            ax3.set_ylabel('g(r)', fontsize=14, fontweight='bold')
            ax3.set_ylim(bottom=0)
            ax3.legend(loc='upper right', frameon=True, shadow=True, fontsize=10)
            ax3.grid(True, alpha=0.3)
            ax3.set_title('Other RDFs', fontsize=14, fontweight='bold', pad=15)

        plt.suptitle(title, fontsize=16, fontweight='bold', y=0.995)
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        logger.info(f"Saved categorized RDF plot to {output_file}")

