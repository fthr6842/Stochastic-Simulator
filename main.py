import pandas as pd
import numpy as np
from config import S0, mu, sigma, T, N, seed, v0, kappa, theta, sigma_v, rho, M
from utils import col_name,\
                    plot_paths,\
                    simulate_gbm,\
                    simulate_heston

if __name__ == "__main__":
    
    # GBM 模擬數據
    t_GBM, array_GBM = simulate_gbm(S0, mu, sigma, T, N, M, seed)
    plot_paths(t_GBM, array_GBM, title="GBM Simulation")
    
    # Heston model 模擬數據
    t_HM, array_HM = simulate_heston(S0, v0, mu, kappa, theta, sigma_v, rho, T, N, M, seed)
    plot_paths(t_HM, array_HM, title="Heston model Simulation")
    
    # 資料整型
    combined_GBM = np.column_stack((t_GBM, array_GBM))
    df_GBM = pd.DataFrame(combined_GBM, columns=col_name(M))
    
    combined_HM = np.column_stack((t_HM, array_HM))
    df_HM = pd.DataFrame(combined_HM, columns=col_name(M))
    
    # 儲存資料
    df_GBM.to_csv(f"output_GBM_{N+1}by{M}.csv", index=True, encoding="utf-8-sig")
    df_HM.to_csv(f"output_HM_{N+1}by{M}.csv", index=True, encoding="utf-8-sig")
    pass

