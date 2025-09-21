import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from numpy.random import default_rng
import config
import utils

def col_name (num):
    List_temp = ['t'] + list(range(1, num+1))
    return List_temp

def plot_paths(t, S, title="Geometric Brownian Motion Simulation"):
    """
    plotting the path of GBM sim
    parameters:
        t : ndarray (N+1,)   - 時間序列
        S : ndarray (N+1, M) - 模擬的價格路徑
        title : str          - 圖片標題
    return: None
    """
    plt.figure(figsize=(10, 6))
    for m in range(S.shape[1]):
        plt.plot(t, S[:, m], lw=1)
    plt.title(title)
    plt.xlabel("Time (Years)")
    plt.ylabel("Price")
    plt.grid(True)
    plt.show()
    return None

def simulate_gbm(S0, mu, sigma, T, N, M=1, seed=None):
    """
    GBM simulation
    Parameters:
        S0 : float   - 初始價格
        mu : float   - 漲幅率 (drift)
        sigma : float - 波動率 (volatility)
        T : float    - 總時間 (通常以年計)
        N : int      - 時間步數
        M : int      - 模擬路徑數 (default=1)
        seed : int   - 隨機種子 (optional)
    return:
        S : ndarray (N+1, M) - 模擬的價格路徑
    """
    if seed is not None: np.random.seed(seed)

    dt = T / N
    t = np.linspace(0, T, N+1) #0到T，共N+1步

    # 初始化價格矩陣
    S = np.zeros((N+1, M))
    S[0] = S0

    # 模擬
    for m in range(M):
        for i in range(1, N+1):
            z = np.random.normal()
            S[i, m] = S[i-1, m] * np.exp((mu - 0.5 * sigma**2) * dt + sigma * np.sqrt(dt) * z)

    return t, S

def simulate_heston(S0, v0, mu, kappa, theta, sigma_v, rho, T, N, M=1, seed=None):
    """
    模擬 Heston Model

    參數:
        S0 : float    - 初始價格
        v0 : float    - 初始方差
        mu : float    - 資產報酬率
        kappa : float - 均值回歸速度
        theta : float - 長期均值 (長期方差)
        sigma_v : float - 方差過程波動率
        rho : float   - 資產價格與方差之間的相關係數
        T : float     - 總時間 (通常以年計)
        N : int       - 時間步數
        M : int       - 模擬路徑數
        seed : int    - 隨機種子 (optional)

    回傳:
        t : ndarray (N+1,)
        S : ndarray (N+1, M) - 模擬的價格路徑
    """
    if seed is not None: np.random.seed(seed)

    dt = T / N
    t = np.linspace(0, T, N+1)

    # 初始化矩陣
    S = np.zeros((N+1, M))
    v = np.zeros((N+1, M))
    S[0, :] = S0
    v[0, :] = v0

    # 模擬路徑
    for i in range(1, N+1):
        # 產生相關常態亂數
        z1 = np.random.normal(size=M)
        z2 = np.random.normal(size=M)
        dW_v = np.sqrt(dt) * z1
        dW_s = np.sqrt(dt) * (rho * z1 + np.sqrt(1 - rho**2) * z2)

        # 更新 v_t (CIR process, 保證非負: max(...,0))
        v[i, :] = np.maximum(
            v[i-1, :] + kappa * (theta - v[i-1, :]) * dt + sigma_v * np.sqrt(np.maximum(v[i-1, :], 0)) * dW_v,
            0
        )

        # 更新 S_t
        S[i, :] = S[i-1, :] * np.exp((mu - 0.5 * v[i-1, :]) * dt + np.sqrt(v[i-1, :]) * dW_s)

    return t, S
