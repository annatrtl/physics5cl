import numpy as np, pandas as pd
from scipy.optimize import least_squares, brentq
df = pd.read_csv("exp4.csv")               

rods = [("Al",    1, 20.20e-2, 1.15e-2, 3.15e-2, 2700*897),
        ("Cu",    2, 20.45e-2, 2.25e-2, 3.45e-2, 8960*385),
        ("Brass", 3, 20.45e-2, 2.20e-2, 3.25e-2, 8470*380)]
dL_wet = 0.5e-2 

def T_model(t, alpha, L, x, N=40): 
    j = 2*np.arange(N) + 1
    return np.sum(4/(j*np.pi)*np.sin(j*np.pi*x/(2*L))
                  * np.exp(-(j*np.pi/(2*L))**2*alpha*np.clip(t, 0, None)[:, None]), axis=1)

def fit_alpha(t, T, L, x):        
    f = lambda p: p[1]*T_model(t - p[2], p[0]*1e-5, L, x) - T
    r = min((least_squares(f, [a, T[:5].mean(), d], bounds=([1e-3, 0, -200], [50, 60, 170]))
             for a in (0.1, 0.5, 2, 10) for d in (0, 20, 60)), key=lambda r: r.cost)
    return r.x[0]*1e-5, np.sqrt(2*r.cost/len(t))

k = {}
for name, n, Lfull, Lwet, pdep, rc in rods:
    t, T = df[[f"time{n}", f"temp{n}"]].dropna().values.T
    if name != "Brass":
        a, rms = fit_alpha(t, T, Lfull - Lwet, Lfull - Lwet - pdep)
        k[name] = a*rc
        # sensitivity to the length in water
        ks = [fit_alpha(t, T, Lfull - Lwet - d, Lfull - Lwet - d - pdep)[0]*rc for d in (-dL_wet, dL_wet)]
        print(f"{name}: alpha = {a:.2e} m2/s, k = {k[name]:.0f} W/m.K "
              f"(range {min(ks+[k[name]]):.0f}-{max(ks+[k[name]]):.0f} for length in water +/-{dL_wet*100:.1f} cm), rms = {rms:.2f} C")
    else:                                    # no decay in 180 s: only an upper limit
        L, x = Lfull - Lwet, Lfull - Lwet - pdep
        for drop in (0.01, 0.03, 0.05):      # largest fall allowed by the data at 180 s
            a = brentq(lambda a: T_model(np.array([180.0]), a, L, x)[0] - (1 - drop), 1e-9, 1e-3)
            print(f"Brass: fall < {drop:.0%} by 180 s -> alpha < {a:.1e}, k < {a*rc:.0f} W/m.K, "
                  f"k/k_Al < {a*rc/k['Al']:.2f}")

# ratio uncertainty: refit Al and Cu with all combinations of the length-in-water shift
Al, Cu = rods[0], rods[1]
ratios = []
for dA in (-dL_wet, 0, dL_wet):
    for dC in (-dL_wet, 0, dL_wet):
        out = []
        for (name, n, Lf, Lw, pd_, rc), d in ((Al, dA), (Cu, dC)):
            t, T = df[[f"time{n}", f"temp{n}"]].dropna().values.T
            out.append(fit_alpha(t, T, Lf - Lw - d, Lf - Lw - d - pd_)[0]*rc)
        ratios.append(out[1]/out[0])
print(f"k_Cu/k_Al = {k['Cu']/k['Al']:.2f}  (range {min(ratios):.2f}-{max(ratios):.2f})")
