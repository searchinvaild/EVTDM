#!/usr/bin/env python3
"""
计算timeVaryingGrout模型的nuCalc值
根据公式: nuCalc = min(nuMax, (tau0 + k*exp(timeCoeff*t)*pow(sr, n)) / srLimited)
"""

import os
import re
import numpy as np
from pathlib import Path

class TimeVaryingGroutCalculator:
    def __init__(self, case_dir="."):
        self.case_dir = Path(case_dir)
        self.coeffs = {}
        self.read_transport_properties()
        
    def read_transport_properties(self):
        """读取transportProperties文件中的系数"""
        transport_file = self.case_dir / "constant" / "transportProperties"
        
        if not transport_file.exists():
            raise FileNotFoundError(f"找不到文件: {transport_file}")
            
        with open(transport_file, 'r') as f:
            content = f.read()
            
        # 查找timeVaryingGroutCoeffs块
        coeffs_match = re.search(r'timeVaryingGroutCoeffs\s*\{([^}]+)\}', content, re.DOTALL)
        if not coeffs_match:
            raise ValueError("找不到timeVaryingGroutCoeffs")
            
        coeffs_text = coeffs_match.group(1)
        
        # 提取各个系数
        patterns = {
            'k': r'k\s+([0-9.e+-]+)',
            'n': r'n\s+([0-9.e+-]+)',
            'tau0': r'tau0\s+([0-9.e+-]+)',
            'nuMax': r'nuMax\s+([0-9.e+-]+)',
            'timeCoeff': r'timeCoeff\s+([0-9.e+-]+)'
        }
        
        for key, pattern in patterns.items():
            match = re.search(pattern, coeffs_text)
            if match:
                self.coeffs[key] = float(match.group(1))
            else:
                raise ValueError(f"找不到系数: {key}")
                
        print("读取的系数:")
        for key, value in self.coeffs.items():
            print(f"  {key} = {value}")
            
    def read_debug1_file(self, file_path):
        """读取Debug1文件中的srLimited值"""
        with open(file_path, 'r') as f:
            content = f.read()
            
        # 查找internalField部分
        internal_match = re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(([\s\S]+?)\)\s*;', content)
        if not internal_match:
            raise ValueError(f"无法解析文件: {file_path}")
            
        n_cells = int(internal_match.group(1))
        values_text = internal_match.group(2)
        
        # 提取数值
        values = []
        for line in values_text.strip().split('\n'):
            if line.strip():
                try:
                    values.append(float(line.strip()))
                except ValueError:
                    continue
                    
        if len(values) != n_cells:
            print(f"警告: 期望{n_cells}个值，实际读取{len(values)}个")
            
        return np.array(values)
        
    def calculate_nu(self, sr_limited, time):
        """根据公式计算nu值"""
        k = self.coeffs['k']
        n = self.coeffs['n']
        tau0 = self.coeffs['tau0']
        nuMax = self.coeffs['nuMax']
        timeCoeff = self.coeffs['timeCoeff']
        
        # 计算nu
        # nuCalc = min(nuMax, (tau0 + k*exp(timeCoeff*t)*pow(sr, n)) / srLimited)
        numerator = tau0 + k * np.exp(timeCoeff * time) * np.power(sr_limited, n)
        nu_calc = numerator / sr_limited
        
        # 应用最大值限制
        nu_calc = np.minimum(nu_calc, nuMax)
        
        return nu_calc
        
    def write_nu_field(self, nu_values, time_dir, original_file):
        """将计算的nu值写入OpenFOAM格式的文件"""
        output_file = time_dir / "timeVaryingGrout_nuCalc"
        
        # 读取原始文件的头部信息
        with open(original_file, 'r') as f:
            lines = f.readlines()
            
        # 找到dimensions行
        dim_line_idx = None
        for i, line in enumerate(lines):
            if line.strip().startswith('dimensions'):
                dim_line_idx = i
                break
                
        # 构建输出内容
        header = lines[:dim_line_idx]
        
        # 修改object名称
        for i, line in enumerate(header):
            if 'object' in line:
                header[i] = '    object      timeVaryingGrout_nuCalc;\n'
                
        output_lines = header
        output_lines.append('dimensions      [0 2 -1 0 0 0 0];\n')  # 动力粘度的量纲
        output_lines.append('\n')
        output_lines.append(f'internalField   nonuniform List<scalar>\n')
        output_lines.append(f'{len(nu_values)}\n')
        output_lines.append('(\n')
        
        # 写入数值
        for value in nu_values:
            output_lines.append(f'{value:.6e}\n')
            
        output_lines.append(');\n')
        output_lines.append('\n')
        
        # 添加边界条件（从原文件复制）
        in_boundary = False
        for line in lines[dim_line_idx+1:]:
            if 'boundaryField' in line:
                in_boundary = True
            if in_boundary:
                output_lines.append(line)
                
        # 写入文件
        with open(output_file, 'w') as f:
            f.writelines(output_lines)
            
        return output_file
        
    def process_all_timesteps(self):
        """处理所有时间步"""
        # 获取所有时间目录
        time_dirs = []
        for item in self.case_dir.iterdir():
            if item.is_dir():
                try:
                    time = float(item.name)
                    time_dirs.append((time, item))
                except ValueError:
                    continue
                    
        time_dirs.sort(key=lambda x: x[0])
        
        print(f"\n找到{len(time_dirs)}个时间目录")
        
        for time, time_dir in time_dirs:
            # 查找Debug1文件
            debug1_file = None
            for file in time_dir.glob("timeVaryingGrout_Debug1*"):
                debug1_file = file
                break
                
            if debug1_file:
                print(f"\n处理时间 t={time}:")
                print(f"  读取文件: {debug1_file}")
                
                try:
                    # 读取srLimited值
                    sr_limited = self.read_debug1_file(debug1_file)
                    print(f"  读取{len(sr_limited)}个网格单元")
                    
                    # 计算nu值
                    nu_calc = self.calculate_nu(sr_limited, time)
                    
                    # 统计信息
                    mask = sr_limited > 1e-10  # 只统计非零值
                    if np.any(mask):
                        print(f"  srLimited范围: {np.min(sr_limited[mask]):.6e} - {np.max(sr_limited[mask]):.6e}")
                        print(f"  nuCalc范围: {np.min(nu_calc[mask]):.6e} - {np.max(nu_calc[mask]):.6e}")
                        print(f"  非零单元数: {np.sum(mask)}")
                    
                    # 写入结果
                    output_file = self.write_nu_field(nu_calc, time_dir, debug1_file)
                    print(f"  输出文件: {output_file}")
                    
                except Exception as e:
                    print(f"  处理失败: {e}")
                    
    def plot_statistics(self):
        """绘制统计图表（可选）"""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("\n未安装matplotlib，跳过绘图")
            return
            
        times = []
        nu_means = []
        nu_maxs = []
        sr_means = []
        
        # 收集数据
        for item in self.case_dir.iterdir():
            if item.is_dir():
                try:
                    time = float(item.name)
                    debug1_file = None
                    for file in item.glob("timeVaryingGrout_Debug1*"):
                        debug1_file = file
                        break
                        
                    if debug1_file:
                        sr_limited = self.read_debug1_file(debug1_file)
                        nu_calc = self.calculate_nu(sr_limited, time)
                        
                        mask = sr_limited > 1e-10
                        if np.any(mask):
                            times.append(time)
                            nu_means.append(np.mean(nu_calc[mask]))
                            nu_maxs.append(np.max(nu_calc[mask]))
                            sr_means.append(np.mean(sr_limited[mask]))
                            
                except:
                    continue
                    
        if times:
            # 排序
            sorted_indices = np.argsort(times)
            times = np.array(times)[sorted_indices]
            nu_means = np.array(nu_means)[sorted_indices]
            nu_maxs = np.array(nu_maxs)[sorted_indices]
            sr_means = np.array(sr_means)[sorted_indices]
            
            # 绘图
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
            
            ax1.plot(times, nu_means, 'b-', label='Mean nu')
            ax1.plot(times, nu_maxs, 'r--', label='Max nu')
            ax1.set_xlabel('Time (s)')
            ax1.set_ylabel('Viscosity (Pa.s)')
            ax1.set_title('Viscosity Evolution')
            ax1.legend()
            ax1.grid(True)
            ax1.set_yscale('log')
            
            ax2.plot(times, sr_means, 'g-')
            ax2.set_xlabel('Time (s)')
            ax2.set_ylabel('Mean Strain Rate (1/s)')
            ax2.set_title('Strain Rate Evolution')
            ax2.grid(True)
            ax2.set_yscale('log')
            
            plt.tight_layout()
            plt.savefig('viscosity_evolution.png', dpi=150)
            print("\n统计图表已保存为: viscosity_evolution.png")
            

def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='计算timeVaryingGrout粘度场')
    parser.add_argument('--case', '-c', default='.', help='OpenFOAM案例目录')
    parser.add_argument('--plot', '-p', action='store_true', help='生成统计图表')
    
    args = parser.parse_args()
    
    # 创建计算器
    calculator = TimeVaryingGroutCalculator(args.case)
    
    # 处理所有时间步
    calculator.process_all_timesteps()
    
    # 生成图表
    if args.plot:
        calculator.plot_statistics()
        
    print("\n处理完成!")
    

if __name__ == "__main__":
    main()
