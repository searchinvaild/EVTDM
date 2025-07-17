#!/usr/bin/env python3
"""
增强版：读取OpenFOAM transportProperties文件并分析kEffective的演化
支持参数设置：最大时间值、子图显示开关、NuMax显示开关
"""

import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

# 参数配置 - 用户可以在这里修改设置
CONFIG = {
    # 计算的最大时间值（秒）
    'max_time': 1000,
    
    # 四个子图的显示开关（True/False）
    'show_log_plot': True,       # 对数坐标图
    'show_linear_plot': True,    # 线性坐标图
    'show_rate_plot': True,      # 变化率图
    'show_normalized_plot': True, # 归一化图
    
    # 是否显示NuMax参考线
    'show_nuMax': True
}

def read_transport_properties(filepath):
    """读取transportProperties文件并提取参数"""
    with open(filepath, 'r') as f:
        content = f.read()
    
    pattern = r'timeVaryingGroutCoeffs\s*\{([^}]+)\}'
    match = re.search(pattern, content, re.DOTALL)
    
    if not match:
        raise ValueError("未找到timeVaryingGroutCoeffs块")
    
    coeffs_block = match.group(1)
    params = {}
    param_pattern = r'(\w+)\s+([0-9.e+-]+);'
    matches = re.findall(param_pattern, coeffs_block)
    
    for param_name, param_value in matches:
        params[param_name] = float(param_value)
    
    return params

def analyze_kEffective_evolution(params, max_time):
    """
    详细分析kEffective的演化特征
    """
    k = params['k']
    nuMax = params['nuMax']
    timeCoeff = params['timeCoeff']
    
    # 关键时间点
    t_50_percent = np.log(0.5 * nuMax / k) / timeCoeff
    t_90_percent = np.log(0.9 * nuMax / k) / timeCoeff
    t_99_percent = np.log(0.99 * nuMax / k) / timeCoeff
    
    # 创建时间数组
    time_array = np.arange(0, max_time + 1, 1)
    
    # 计算kEffective
    kEffective = np.minimum(nuMax, k * np.exp(timeCoeff * time_array))
    
    # 计算变化率
    dkEff_dt = np.gradient(kEffective, time_array)
    
    return {
        'time': time_array,
        'kEffective': kEffective,
        'dkEff_dt': dkEff_dt,
        't_50': t_50_percent,
        't_90': t_90_percent,
        't_99': t_99_percent
    }

def create_comprehensive_plot(results, params, config, save_prefix='kEffective'):
    """
    创建综合分析图，根据配置显示子图
    """
    # 计算需要显示的子图数量
    num_plots = sum([
        config['show_log_plot'],
        config['show_linear_plot'],
        config['show_rate_plot'],
        config['show_normalized_plot']
    ])
    
    if num_plots == 0:
        print("警告: 所有子图显示开关均为关闭状态，跳过绘图")
        return
    
    # 创建子图布局
    if num_plots == 1:
        fig, axes = plt.subplots(1, 1, figsize=(10, 7))
        axes = [[axes]]  # 统一为二维数组格式
    elif num_plots == 2:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        axes = [axes]  # 统一为二维数组格式
    elif num_plots == 3:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        # 隐藏最后一个子图
        axes[1, 1].axis('off')
    else:  # 4个子图
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 当前子图索引
    plot_index = 0
    
    # 子图1: kEffective vs 时间（对数刻度）
    if config['show_log_plot']:
        ax = axes[plot_index // 2][plot_index % 2] if num_plots > 1 else axes[0][0]
        plot_index += 1
        
        # 跳过0点避免对数刻度问题
        valid_idx = results['time'] > 0
        ax.semilogy(results['time'][valid_idx], results['kEffective'][valid_idx], 'b-', linewidth=2)
        
        # 显示NuMax参考线
        if config['show_nuMax']:
            ax.axhline(y=params['nuMax'], color='r', linestyle='--', label='nuMax')
        
        # 只显示在0-max_time秒范围内的关键时间点
        max_time = config['max_time']
        if 0 < results['t_50'] <= max_time:
            ax.axvline(x=results['t_50'], color='g', linestyle=':', alpha=0.7, label='50% nuMax')
        if 0 < results['t_90'] <= max_time:
            ax.axvline(x=results['t_90'], color='orange', linestyle=':', alpha=0.7, label='90% nuMax')
        if 0 < results['t_99'] <= max_time:
            ax.axvline(x=results['t_99'], color='purple', linestyle=':', alpha=0.7, label='99% nuMax')
        
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('kEffective (m²/s)')
        ax.set_title('kEffective Evolution (Log scale)')
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.set_xlim(0, max_time)  # 设置固定时间范围
    
    # 子图2: kEffective vs 时间（线性刻度）
    if config['show_linear_plot']:
        ax = axes[plot_index // 2][plot_index % 2] if num_plots > 1 else axes[0][0]
        plot_index += 1
        
        ax.plot(results['time'], results['kEffective'], 'b-', linewidth=2)
        
        # 显示NuMax参考线
        if config['show_nuMax']:
            ax.axhline(y=params['nuMax'], color='r', linestyle='--', label='nuMax')
        
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('kEffective (m²/s)')
        ax.set_title('kEffective Evolution (Linear scale)')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, config['max_time'])  # 设置固定时间范围
        if config['show_nuMax']:
            ax.legend()
    
    # 子图3: 变化率
    if config['show_rate_plot']:
        ax = axes[plot_index // 2][plot_index % 2] if num_plots > 1 else axes[0][0]
        plot_index += 1
        
        ax.plot(results['time'], np.abs(results['dkEff_dt']), 'g-', linewidth=2)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('|dkEffective/dt| (m²/s²)')
        ax.set_title('Rate of Change')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, config['max_time'])  # 设置固定时间范围
    
    # 子图4: 归一化图
    if config['show_normalized_plot']:
        ax = axes[plot_index // 2][plot_index % 2] if num_plots > 1 else axes[0][0]
        plot_index += 1
        
        kEff_normalized = results['kEffective'] / params['nuMax']
        ax.plot(results['time'], kEff_normalized, 'b-', linewidth=2)
        ax.axhline(y=0.5, color='g', linestyle=':', alpha=0.7, label='50%')
        ax.axhline(y=0.9, color='orange', linestyle=':', alpha=0.7, label='90%')
        ax.axhline(y=0.99, color='purple', linestyle=':', alpha=0.7, label='99%')
        ax.axhline(y=1.0, color='r', linestyle='--', label='100%')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('kEffective / nuMax')
        ax.set_title('Normalized kEffective')
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.set_xlim(0, config['max_time'])  # 设置固定时间范围
        ax.set_ylim(0, 1.1)
    
    # 添加总标题和参数信息
    fig.suptitle(f'kEffective Evolution Analysis (0-{config["max_time"]}s)\n' + 
                 f'k={params["k"]:.3e}, timeCoeff={params["timeCoeff"]:.3e}, nuMax={params["nuMax"]:.3e}',
                 fontsize=14)
    
    plt.tight_layout()
    plt.savefig(f'{save_prefix}_analysis.png', dpi=300, bbox_inches='tight')
    print(f"综合分析图已保存至: {save_prefix}_analysis.png")
    
    # 保存数据到CSV
    df = pd.DataFrame({
        'time_s': results['time'],
        'kEffective_m2s': results['kEffective'],
        'dkEff_dt_m2s2': results['dkEff_dt'],
        'normalized_kEff': results['kEffective'] / params['nuMax']
    })
    csv_path = f'{save_prefix}_data.csv'
    df.to_csv(csv_path, index=False)
    print(f"数据已保存至: {csv_path}")

def main():
    """
    主函数
    """
    # 打印当前配置
    print("当前配置:")
    for key, value in CONFIG.items():
        print(f"  {key}: {value}")
    
    transport_file = Path('constant/transportProperties')
    
    if not transport_file.exists():
        print(f"错误: 找不到文件 {transport_file}")
        return
    
    try:
        # 读取参数
        params = read_transport_properties(transport_file)
        print("\n成功读取参数:")
        for key, value in params.items():
            print(f"  {key}: {value:.3e}")
        
        # 分析kEffective演化
        results = analyze_kEffective_evolution(params, CONFIG['max_time'])
        
        # 打印关键信息
        print(f"\n关键时间点:")
        print(f"  达到50% nuMax: {results['t_50']:.1f} s")
        print(f"  达到90% nuMax: {results['t_90']:.1f} s")
        print(f"  达到99% nuMax: {results['t_99']:.1f} s")
        
        # 创建图表
        create_comprehensive_plot(results, params, CONFIG)
        
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()