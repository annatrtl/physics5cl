import numpy as np
import matplotlib.pyplot as plt

# data
wl_err = 0 # nm ; monocromatic lasers
c = 3e8 # m/s


vio_wl = 405 # nm

wls = np.array([405, 532, 635]) * 1e-9  # m ;  violet, green, red in order
f = c / wls
f_err = 0

vio_v = np.array([0.0, -200.0, -400.0, -600.0, -800.0, -1000.0, -1100.0, -1160.0, -1170.0, 100.0, 0.0, 200]) # mV
vio_v /= 1000 # V
vio_curr = np.array([(2221+2224)/2, (1457+1468)/2, (847+850)/2, (422+426)/2, (165+169)/2, (40)/2, (10+18)/2, (0+3)/2, (0+1)/2, 2623, (2210+2221)/2, (3024+3028)/2]) # nA

green_wl = 532 # nm
green_v = np.array([0.0, -100.0, -200.0, -300.0, -400.0, -500.0, -550.0, -545.0, -535.0]) # mV
green_v /= 1000 # V
green_curr = np.array([(423+431)/2, 282, (156+163)/2, (76+80)/2, (23+34)/2, (5+9)/2, (3+0)/2, (0+5)/2, (1+5)/2]) # nA

red_wl = 635 # nm
red_v = np.array([0.0, -50.0, -100.0, -120.0, -125.0, -132.0, -185.0, 0.0, 50.0, 100.0]) # mV
red_v /= 1000 # V
red_curr = np.array([34, (16+20)/2, (5+9)/2, (1+5)/2, (0+3)/2, (0+3)/2, (0+3)/2, (30+38)/2, (52+60)/2, (80+83)/2]) # nA

v = [vio_v, green_v, red_v]
i = [vio_curr, green_curr, red_curr]

# more careful analysis
windows = [(4, 450), (4, 200), (2.5, 40)]
ke, ke_err = [], []
for V, I, (lo, hi) in zip(v, i, windows):
    mask = (V <= 0) & (I > lo) & (I <= hi)
    (a, b0), cov = np.polyfit(V[mask], I[mask], 1, cov=True)
    vs = -b0 / a
    J = np.array([b0 / a**2, -1 / a])
    ke.append(-vs)
    ke_err.append(np.sqrt(J @ cov @ J + 0.001**2))
ke, ke_err = np.array(ke), np.array(ke_err)
print("Stopping voltages [V]:", ke, "+/-", ke_err)

m = 4e-15
for _ in range(10):
    sig = np.sqrt(ke_err**2 + (m * f_err)**2)
    (m, b), cov = np.polyfit(f, ke, 1, w=1/sig, cov="unscaled")

h, phi, f0 = m, -b, -b / m
dh, dphi = np.sqrt(np.diag(cov))
J = np.array([b / m**2, -1 / m])
df0 = np.sqrt(J @ cov @ J)

print(f"h: {h:.3e} +/- {dh:.3e} eV*s  (accepted 4.136e-15)")
print(f"phi: {phi:.5f} +/- {dphi:.5f} eV")
print(f"f0: {f0:.3e} +/- {df0:.3e} Hz")

fig = plt.figure()
plt.xlabel("Frequency [Hz]")
plt.ylabel("KEmax [eV]")
plt.title("Stopping energy vs frequency")

plt.errorbar(f, ke, xerr=f_err, yerr=ke_err, fmt="o", capsize=3)
plt.plot(f, m * f + b, label="Weighted fit")

plt.legend()
plt.show()
