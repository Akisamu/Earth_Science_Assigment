import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, integrate, Eq, nsolve

P = symbols('P')
T = symbols('T')

# Cp 参数
aq, bq, cq, dq = 0.0965, -0.0577e-5, -444.8, -0.7982
ac, bc, cc, dc = 0.1107, -0.5189-5, 0.0, -1.1283
a, b, c, d = aq-ac, (bq-bc)*1e-5, cq-cc, dq-dc

# ΔH(298,1) 和 ΔS(298,1)，单位可根据体系调整
delta_H_298_1 =  (-905.52) - (-910.88) # J/mol
delta_S_298_1 = -1.1e-3   # J/mol·K
delta_v_298_1 = 2.064 - 2.269   # kJ/kbar

exp_delta_Cp = a*T + b*T + c*T**(-2) + d*T**(-0.5)
expr = (delta_H_298_1 + integrate(exp_delta_Cp , (T,298,T)) + integrate(delta_v_298_1 ,(P,1,P))) - T*((delta_S_298_1) + integrate(exp_delta_Cp/T ,(T,298,T)))

eq = Eq(expr, 0)

# 构造 T 值范围（从 200 到 1000，按需修改）
T_values = np.linspace(300, 2000, 1700)
P_values = []

# 逐个 T 解出 P
for t in T_values:
    try:
        P_sol = nsolve(expr.subs(T, t), P, 1.0)  # 初始猜测 P = 1.0
        P_values.append(float(P_sol)*1e-3)
    except:
        P_values.append(np.nan)  # 解失败时插入 NaN

# 画图
plt.figure(figsize=(16, 9))
plt.plot(T_values, P_values, label='P vs T', color='green')
plt.xlim(500, 2000)  # 横轴：300 到 2000
plt.ylim(0, 80)      # 纵轴：1 到 10
plt.yticks(np.arange(0, 81, 5))
plt.xlabel('T (k)')
plt.ylabel('P (kbar)')
plt.title('P-T Relationship from Given Equation')
plt.grid(True)
plt.legend()
plt.show()

