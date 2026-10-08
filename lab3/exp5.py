import math

# data 
L = 78.30 # cm
wl = 6.346e-5 # cm

# 2/2/2
dist_fine_222 = 0.1 # cm
th_fine_222 = 0.2 # cm
m_fine_222 = 8
th_env_222 = 4.4 # cm
m_env_222 = 1

# 2/1/2
dist_fine_212 = 0.1 # cm
th_fine_212 = 0.3 # cm
m_fine_212 = 7
th_env_212 = 5.1 # cm
m_env_212 = 1

# 2/2/2 analysis
d_222 = (wl * L * 2 * 1)  / th_fine_222
err_d_222 = d_222 * math.sqrt(((0.05/L) ** 2) + ((0.05/th_fine_222) ** 2))
a_222 = (wl * L * 2) / th_env_222
err_a_222 = a_222 * math.sqrt(((0.05/L) ** 2) + ((0.05/th_env_222) ** 2))
miss_222 = d_222 / a_222 
err_miss_222 = miss_222 * math.sqrt(((err_a_222/a_222) ** 2) + ((err_d_222/d_222) ** 2))
print(f"Slit thickness: {a_222:.6f} +/- {err_a_222}, center-to-center slit seperations: {d_222:.6f} +/- {err_d_222}, missing orders: {miss_222:.6f} +/- {err_miss_222}")

# 2/1/2 analysis
d_212 = (wl * L * 2 * 1)  / th_fine_212
err_d_212 = d_212 * math.sqrt(((0.05/L) ** 2) + ((0.05/th_fine_212) ** 2))
a_212 = (wl * L * 2) / th_env_212
err_a_212 = a_212 * math.sqrt(((0.05/L) ** 2) + ((0.05/th_env_212) ** 2))
miss_212 = d_212 / a_212 
err_miss_212 = miss_212 * math.sqrt(((err_a_212/a_212) ** 2) + ((err_d_212/d_212) ** 2))
print(f"Slit thickness: {a_212} +/- {err_a_212}, center-to-center slit seperations: {d_212} +/- {err_d_212}, missing orders: {miss_212} +/- {err_miss_212}")



