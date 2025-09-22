# Stochastic-Simulator
#### as title...
#### 目標是收錄GBM, Heston model, SVJ, SVCJ的模擬、繪圖、參數估計、MC訂價，並共用統一規格。
## 隨機過程模擬器
### 1. main.py: 主程式。可產生模擬數據、圖片，並將結果儲存。
### 2. config.py: 參數檔。
### 3. utils.py: 函數檔。目前已收錄GBM, Heston model的模擬、繪圖、MC訂價。
## CRR bin tree with Heston model
原始 CRR binomial tree 假設標的資產價格服從GBM，現在考慮資產價格服從Heston model。  
先回顧 GBM 假設下的CRR推導:  
假設資產價格為:  
### $dS_t = \mu S_t dt + \sigma S_t dW_t$
基於 Ito's lemma，我們可得:  
### $\ln(\frac{S_{t+\Delta t}}{S_{t}})\sim N((\mu - \frac{1}{2} \sigma^2) \Delta t, \sigma^2 \Delta t)$
因此設定上升幅度與下降幅度為:  
### $u = e^{\sigma \Delta t}, \quad d = e^{- \sigma \Delta t} = u^{-1}$
以及上升機率為:  
### $p = \frac{e^{r\Delta t} - d}{u - d}$
先回顧Heston model:  
### $dS_t = \mu S_t dt + \sqrt{v_t} S_t dW_t^{(S)}$
### $dv_t = \kappa (\theta - v_t) dt + \xi \sqrt{v_t} dW_t^{(v)}$
其中  
### $dW_t^{(S)} \cdot dW_t^{(v)} = \rho dt$


### $u_{i,t} = e^{\sqrt{v_t \Delta t}}, d_{i,t}=u_{i,t}^{-1}$
### $p_{i,t} = \frac{e^{r\Delta t}-d_{i,t}}{u_{i,t}-d_{i,t}}$
