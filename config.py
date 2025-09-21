T = 1        # 模擬一年
N = 252      # 模擬252個交易日
M = 5        # 模擬5條路徑
seed = 42    # 亂數種子
r = 0.03     # Risk-free rate
K = 100

S0 = 100     # 初始價格
mu = 0.05    # 年化報酬率
sigma = 0.2  # 年化波動率

v0 = 0.04      # 初始方差 (即 20% 年化波動率^2)
kappa = 2.0    # mean reversion speed
theta = 0.04   # long-term variance (20% 年化波動率^2)
sigma_v = 0.3  # vol-of-vol
rho = -0.7     # 負相關 (常見於市場)