import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, integrate, nsolve

P = symbols('P')
T = symbols('T')

class Minecraft:

    def __init__(self, abcd: tuple, h: float, s: float, v: float):
        a, b, c, d = abcd
        self.abcd = (a, b*1e-5, c, d)
        self.h = h
        self.v = v
        self.s = s*1e-3
        self.transtorm = {
            'T_values': [],
            'P_values': []
        }

    def get_delta_abce_transform_to(self, after: 'Minecraft') -> tuple:
        return (
            after.abcd[0] - self.abcd[0],
            after.abcd[1] - self.abcd[1],
            after.abcd[2] - self.abcd[2],
            after.abcd[3] - self.abcd[3]
        )

    def get_delta_h_transform_to(self, after: 'Minecraft') -> float:
        return after.h - self.h
    
    def get_delta_s_transform_to(self, after: 'Minecraft') -> float:
        return after.s - self.s
    
    def get_delta_v_transform_to(self, after: 'Minecraft') -> float:
        return after.v - self.v
    
    @classmethod
    def get_expr_of_delta_g(cls, delta_abcd: tuple, delta_h: float, delta_s: float, delta_v: float) -> any:
        a, b, c, d = delta_abcd
        exp_delta_Cp = a + b*T + c*(T**-2) + d*(T**-0.5)
        expr = (delta_h + integrate(exp_delta_Cp , (T,298,T)) + integrate(delta_v ,(P,1,P))) - T*(delta_s + integrate((exp_delta_Cp)/T ,(T,298,T)))
        return expr
    
    @classmethod
    def sub_solve(cls, expr, eqauls, start, stop):
        T_values = np.linspace(start, stop, stop-start)
        P_values = []

        # 逐个 T 解出 P
        for t in T_values:
            try:
                P_sol = nsolve(expr.subs(T, t)-eqauls, P, 1.0)  # 初始猜测 P = 1.0
                P_values.append(float(P_sol))
                print(f't={t},p={P_sol}; ')
            except:
                P_values.append(np.nan)  # 解失败时插入 NaN
        T_values = [x - 273.15 for x in T_values]
        return (T_values, P_values)


andalusite = Minecraft(abcd=(0.2773, -0.6588, -1914.1, -2.2656), h=-2588.77, s=92.70, v=5.153)
kyanite = Minecraft(abcd=(0.2794, -0.7124, -2055.6, -2.2894), h=-2593.13, s=83.50, v=4.414)
sillimanite = Minecraft(abcd=(0.2802, -0.6900, -1375.7, -2.3994), h=-2585.89, s=95.50, v=4.986)

# andalusite to kyanite
expr = Minecraft.get_expr_of_delta_g(
    delta_abcd = andalusite.get_delta_abce_transform_to(kyanite),
    delta_h = andalusite.get_delta_h_transform_to(kyanite),
    delta_s = andalusite.get_delta_s_transform_to(kyanite),
    delta_v = andalusite.get_delta_v_transform_to(kyanite)
)
andalusite.transtorm['T_values'], andalusite.transtorm['P_values'] = Minecraft.sub_solve(expr, eqauls=0, start=200, stop=1200)

# kyanite to sillimanite
expr = Minecraft.get_expr_of_delta_g(
    delta_abcd = kyanite.get_delta_abce_transform_to(sillimanite),
    delta_h = kyanite.get_delta_h_transform_to(sillimanite),
    delta_s = kyanite.get_delta_s_transform_to(sillimanite),
    delta_v = kyanite.get_delta_v_transform_to(sillimanite)
)
kyanite.transtorm['T_values'], kyanite.transtorm['P_values'] = Minecraft.sub_solve(expr, eqauls=0, start=200, stop=1200)

# andalusite to sillimanite
expr = Minecraft.get_expr_of_delta_g(
    delta_abcd = andalusite.get_delta_abce_transform_to(sillimanite),
    delta_h = andalusite.get_delta_h_transform_to(sillimanite),
    delta_s = andalusite.get_delta_s_transform_to(sillimanite),
    delta_v = andalusite.get_delta_v_transform_to(sillimanite)
)
sillimanite.transtorm['T_values'], sillimanite.transtorm['P_values'] = Minecraft.sub_solve(expr, eqauls=0, start=200, stop=1200)


# 画图
plt.figure(figsize=(9, 16))
plt.plot(andalusite.transtorm['T_values'], andalusite.transtorm['P_values'], label='andalusite to kyanite', color='green', linestyle=':')
plt.plot(kyanite.transtorm['T_values'], kyanite.transtorm['P_values'], label='kyanite to sillimanite', color='blue', linestyle=':')
plt.plot(sillimanite.transtorm['T_values'], sillimanite.transtorm['P_values'], label='andalusite to sillimanite', color='red', linestyle=':')
plt.xlim(0, 1000) 
plt.ylim(0, 12) 
plt.yticks(np.arange(0, 13, 2))
plt.xlabel('Temperature (°C)')
# plt.xlabel('T (K)')
plt.ylabel('Pressure (kbar)')
plt.title('P-T Relationship')
plt.grid(True)
plt.legend()
plt.show()

