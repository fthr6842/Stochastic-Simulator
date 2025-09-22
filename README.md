# Stochastic-Simulator
#### as title...
#### 目標是收錄GBM, Heston model, SVJ, SVCJ的模擬、繪圖、參數估計、MC訂價，並共用統一規格。
## 隨機過程模擬器
### 1. main.py: 主程式。可產生模擬數據、圖片，並將結果儲存。
### 2. config.py: 參數檔。
### 3. utils.py: 函數檔。目前已收錄GBM, Heston model的模擬、繪圖、MC訂價。
## CRR bin tree with Heston model
#### u_{i,t} = e^{\sqrt{v_t \delta t}}, d_{i,t}=u_{i,t}^{-1}
#### p_{i,t} = \frac{e^{r\delta t}-d_{i,t}}{u_{i,t}-d_{i,t}}
