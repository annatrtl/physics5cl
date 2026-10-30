import numpy as np, pandas as pd
import matplotlib.pyplot as plt

c = 4.184  # J/g.K
df = pd.read_csv("exp3data.csv")       
trial = lambda n: df[[f"time{n}", f"temp{n}"]].dropna().values.T

L_wall= 0.2e-2                
L_base = 0.2e-4
r_bot = 5.60e-2/2 - L_wall

# cup trials 1-3: (trial, water height, inner radius at water level, bath temp C)
cup = [(1, 2.60e-2, 6.20e-2/2 - L_wall, 45.126),
       (2, 2.45e-2, 6.20e-2/2 - L_wall, 45.522),
       (3, 6.15e-2, 6.95e-2/2 - L_wall, 45.773)]

# metal can trials 4-6: (trial, water height, bath temp C)
r_metal = 3.5e-2 - 0.03e-2           
metal = [(4, 6.85e-2 - 0.03e-2, 45.989), (5, 6.90e-2 - 0.03e-2, 45.091), (6, 6.20e-2 - 0.03e-2, 44.715)]
skip = 10                                  

skip_metal = 40

mass = lambda r1, r2, h: 1e6*np.pi*h/3*(r1**2 + r1*r2 + r2**2)
def g(h, r2):                                    # G = k*g
    wall = np.pi*(r_bot + r2 + L_wall)*np.hypot(h, r2 - r_bot)/L_wall
    return wall + np.pi*r_bot**2/L_base

def decay(n, Tb, skip=skip):
    #Fit ln(Tb - T) = b - t/tau
    t, T = trial(n)
    s = t >= t[T.argmin()] + skip
    m = np.polyfit(t[s], np.log(abs(Tb - T[s])), 1)[0]
    r2 = np.corrcoef(t[s], np.log(abs(Tb - T[s])))[0, 1]**2
    return -m, r2                                # rate (1/s), R^2

# styrofoam cup. G = m*c*rate, k = G/g
k, Gc = [], []
for n, h, r2, Tb in cup:
    rate, R2 = decay(n, Tb)
    Gc.append(rate*mass(r_bot, r2, h)*c); k.append(Gc[-1]/g(h, r2))
    print(f"cup trial {n}: tau = {1/rate/60:.0f} min, G = {Gc[-1]:.3f} W/K, k = {k[-1]:.4f}, R2 = {R2:.3f}")
k = np.array(k)
print(f"k = {k.mean():.4f} +/- {k.std(ddof=1)/len(k)**.5:.4f} (SEM) W/m.K, {100*(k.mean()-0.033)/0.033:+.3f}% vs 0.033\n")

# leak = G * (T_inside - T_room) * duration
print("Leak for dT = 15 K over 300 s:", round(np.mean(Gc)*15*300), "J\n")

# metal can
Gm = []
for n, h, Tb in metal:
    rate, R2 = decay(n, Tb, skip_metal)
    Gm.append(rate*1e6*np.pi*r_metal**2*h*c)
    print(f"metal trial {n}: tau = {1/rate:.0f} s, G = {Gm[-1]:.2f} W/K, log-fit R2 = {R2:.3f}")
print(f"Metal G = {np.mean(Gm):.2f} W/K, about {np.mean(Gm)/np.mean(Gc):.0f}x the cup")

sigma_T = 0.02                                    
plt.figure(figsize=(6, 4.5))
for n, h, Tb in metal:
    t, T = trial(n)
    s = t >= t[T.argmin()] + skip_metal
    x, dT = t[s], abs(Tb - T[s])
    y, sig = np.log(dT), sigma_T/dT              
    m, b = np.polyfit(x, y, 1, w=1/sig)          
    chi2 = np.sum(((y - (m*x + b))/sig)**2)/(len(x) - 2)
    print(f"metal trial {n}: reduced chi2 = {chi2:.2f}")
    plt.plot(x, y, ".", ms=3, label=f"trial {n}")
    plt.plot(x, m*x + b, "k-", lw=0.8)
plt.xlabel("time (s)"); plt.ylabel("ln(T_bath - T)"); plt.title("Metal can: log of ΔT vs time")
plt.legend(); plt.tight_layout(); plt.savefig("metal_log_fit.png", dpi=150); plt.show()
