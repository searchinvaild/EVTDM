#!/usr/bin/env python3
"""
增强版：读取OpenFOAM transportProperties文件并分析kEffective的演化
"""

import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

def read_transport_properties(filepath):
    """同上"""
    # [保持相同的实现]
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

def analyze_kEffective_evolution(params):
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
    
    # 创建详细的时间数组
    t_max = max(t_99_percent * 1.2, 10000)
    
    # 对数间隔的时间点，在早期有更多采样点
    t_log = np.logspace(-2, np.log10(t_max), 500)
    t_linear = np.linspace(0, t_max, 500)
    time_array = np.unique(np.concatenate([t_log, t_linear]))
    time_array = np.sort(time_array)
    
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

def create_comprehensive_plot(results, params, save_prefix='kEffective'):
    """
    创建综合分析图
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 子图1: kEffective vs 时间（对数刻度）
    ax1 = axes[0, 0]
    ax1.semilogy(results['time'], results['kEffective'], 'b-', linewidth=2)
    ax1.axhline(y=params['nuMax'], color='r', linestyle='--', label='nuMax')
    ax1.axvline(x=results['t_50'], color='g', linestyle=':', alpha=0.7, label='50% nuMax')
    ax1.axvline(x=results['t_90'], color='orange', linestyle=':', alpha=0.7, label='90% nuMax')
    ax1.axvline(x=results['t_99'], color='purple', linestyle=':', alpha=0.7, label='99% nuMax')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('kEffective (m²/s)')
    ax1.set_title('kEffective Evolution (Log scale)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # 子图2: kEffective vs 时间（线性刻度）
    ax2 = axes[0, 1]
    ax2.plot(results['time'], results['kEffective'], 'b-', linewidth=2)
    ax2.axhline(y=params['nuMax'], color='r', linestyle='--', label='nuMax')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('kEffective (m²/s)')
    ax2.set_title('kEffective Evolution (Linear scale)')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, results['t_99'] * 1.2)
    
    # 子图3: 变化率
    ax3 = axes[1, 0]
    ax3.semilogy(results['time'], np.abs(results['dkEff_dt']), 'g-', linewidth=2)
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('|dkEffective/dt| (m²/s²)')
    ax3.set_title('Rate of Change')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, results['t_99'] * 1.2)
    
    # 子图4: 归一化图
    ax4 = axes[1, 1]
    kEff_normalized = results['kEffective'] / params['nuMax']
    ax4.plot(results['time'], kEff_normalized, 'b-', linewidth=2)
    ax4.axhline(y=0.5, color='g', linestyle=':', alpha=0.7, label='50%')
    ax4.axhline(y=0.9, color='orange', linestyle=':', alpha=0.7, label='90%')
    ax4.axhline(y=0.99, color='purple', linestyle=':', alpha=0.7, label='99%')
    ax4.axhline(y=1.0, color='r', linestyle='--', label='100%')
    ax4.set_xlabel('Time (s)')
    ax4.set_ylabel('kEffective / nuMax')
    ax4.set_title('Normalized kEffective')
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    ax4.set_xlim(0, results['t_99'] * 1.2)
    ax4.set_ylim(0, 1.1)
    
    # 添加总标题和参数信息
    fig.suptitle(f'kEffective Evolution Analysis\n' + 
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
    transport_file = Path('constant/transportProperties')
    
    if not transport_file.exists():
        print(f"错误: 找不到文件 {transport_file}")
        return
    
    try:
        # 读取参数
        params = read_transport_properties(transport_file)
        print("成功读取参数:")
        for key, value in params.items():
            print(f"  {key}: {value:.3e}")
        
        # 分析kEffective演化
        results = analyze_kEffective_evolution(params)
        
        # 打印关键信息
        print(f"\n关键时间点:")
        print(f"  达到50% nuMax: {results['t_50']:.1f} s")
        print(f"  达到90% nuMax: {results['t_90']:.1f} s")
        print(f"  达到99% nuMax: {results['t_99']:.1f} s")
        
        # 创建图表
        create_comprehensive_plot(results, params)
        
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
