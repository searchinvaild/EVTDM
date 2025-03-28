import numpy as np
import matplotlib.pyplot as plt

class TimeSlurryViscosity:
    def __init__(self, params):
        """
        初始化黏度模型参数
        Parameters:
            params : dict
                k: 幂律系数 [m^2/s^n]
                n: 幂律指数 [-]
                tau0: 屈服应力系数 [Pa]
                nuMax: 最大黏度 [m^2/s]
                timeCoeff: 时间影响指数 [-]
                gamma_dot: 剪切速率 [1/s] (可设置为数组进行多工况分析)
        """
        self.k = params['k']
        self.n = params['n']
        self.tau0 = params['tau0']
        self.nuMax = params['nuMax']
        self.timeCoeff = params['timeCoeff']
        self.gamma_dot = params['gamma_dot']
        
        # 单位转换常数（根据实际需要调整）
        self.t_scale = 1.0  # 时间单位换算系数
        self.stress_scale = 1.0  # 应力单位换算系数

    def calculate_viscosity(self, t):
        """
        计算指定时间的黏度
        Parameters:
            t : float/array
                时间值（单位：秒）
        Returns:
            nu : float/array
                计算得到的运动黏度 [m^2/s]
        """
        t = np.array(t) * self.t_scale
        
        # 计算时间相关项
        time_term = self.k * (t**self.timeCoeff)
        
        # 计算剪切速率项（处理零剪切速率情况）
        gamma = np.where(self.gamma_dot > 1e-6, self.gamma_dot, 1e-6)
        
        # 核心计算公式
        nu = (self.tau0 + time_term * (gamma**self.n)) / gamma
        
        # 应用黏度上限
        return np.minimum(nu, self.nuMax)

    def plot_viscosity_evolution(self, t_range=(0, 100), dt=1):
        """
        绘制黏度随时间变化曲线
        Parameters:
            t_range : tuple
                时间范围 (start, end) [秒]
            dt : float
                时间步长 [秒]
        """
        t = np.arange(t_range[0], t_range[1]+dt, dt)
        nu = self.calculate_viscosity(t)
        
        plt.figure(figsize=(10, 6))
        plt.plot(t, nu, 'b-', linewidth=2, 
                label=f'γ̇={self.gamma_dot} 1/s\nk={self.k}, n={self.n}')
        
        plt.axhline(y=self.nuMax, color='r', linestyle='--', 
                   label=f'nuMax={self.nuMax}')
        
        plt.title("Time-dependent Viscosity Evolution")
        plt.xlabel("Time (s)", fontsize=12)
        plt.ylabel("Kinematic Viscosity (m²/s)", fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        
        # 自动保存图片
        plt.savefig('viscosity_evolution.png', dpi=300)
        plt.show()

if __name__ == "__main__":
    # ========== 参数设置 ==========
    # 示例参数（需要根据实际情况修改）
    params = {
        'k': 0.01136,     # 示例值
        'n': 0.8751,       # 示例值
        'tau0': 0.00000178571, # 示例值
        'nuMax': 1e-1,     # 示例值
        'timeCoeff': 1.23, # 示例值
        'gamma_dot': 1.0  # 剪切速率 [1/s]
    }
    
    # ========== 执行计算 ==========
    model = TimeSlurryViscosity(params)
    
    # 计算单个时间点
    t_test = 10.0
    print(f"At t={t_test}s, nu={model.calculate_viscosity(t_test):.2e} m²/s")
    
    # 绘制演化曲线（时间范围可调）
    model.plot_viscosity_evolution(t_range=(0, 20), dt=0.1)
