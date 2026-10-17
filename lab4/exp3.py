import numpy as np
import matplotlib.pyplot as plt

# data
wl_err = 10 # nm
c = 3e8 # m/s


blue_wl = 462 # nm

wls = np.array([462, 533, 618]) * 1e-9  # m ;  blue, green, red in order
f = c / wls
f_err = c * 10e-9 / wls**2

blue_v = np.array([0.0, -500.0, -700.0, -900.0, -930.0, -951.0, -960.0, -955.0, 200.0, 400.0, 450.0]) # mV
blue_v /= 1000 # V
blue_curr = np.array([(8360+8350)/2, (1890+1870)/2, (580+575)/2, (51+47)/2, (25+18)/2, (12+5)/2, 0, (3+0)/2, 1000 * (11.73+11.72)/2, 1000 * (15.15+15.14)/2, 1000 * (16.00+15.99)/2]) # nA

green_wl = 533 # nm
green_v = np.array([0.0, -200.0, -300.0, -400.0, -500.0, -600.0, -700.0, -721.0, -724.0, -722.0, 0.0, 200.0, 400.0, 450.0]) # mV
green_v /= 1000 # V
green_curr = np.array([(3265+3287)/2, (1652+1663)/2, (1007+1011)/2, (521+524)/2, 200, (58+65)/2, (7+10)/2, (3+7)/2, (0+3)/2, (0+3)/2, (3320+3330)/2, (5229+5240)/2, (7172+7180)/2, (7655+7658)/2]) # nA

red_wl = 618 # nm
red_v = np.array([0.0, -100.0, -200.0, -300.0, -350.0, -355.0, -360.0, -361.0, 0.0, 100.0, 200.0, 250.0]) # mV
red_v /= 1000 # V
red_curr = np.array([(819+823)/2, (379+393)/2, (112+123)/2, (14+18)/2, (3+7)/2, (3+10)/2, (0+3)/2, (0+3)/2, (849+852)/2, (1370+1377)/2, (1953+1957)/2, (2241+2244)/2]) # nA

v = [blue_v, green_v, red_v]
i = [blue_curr, green_curr, red_curr]

# linear fit for each LED
blue_coef = np.polyfit(blue_curr, blue_v, 1) 

green_coef = np.polyfit(green_curr, green_v, 1) 

red_coef = np.polyfit(red_curr, red_v, 1) 

# KEmax = hc / wl - phi rough plot/fit
ke_rough = np.array([np.abs(V[(V < 0) & (I <= 3)]).min() for V, I in zip(v, i)])
freq = c / wls
photoelec_coef = np.polyfit(freq, ke_rough, 1)
photoelec_lin_func = np.poly1d(photoelec_coef)
photoelec_yfit = photoelec_lin_func(freq)

f_cutoff = -1 * photoelec_coef[1] / photoelec_coef[0]

fig = plt.figure()
plt.title("Stopping energy vs frequency")
plt.xlabel("Frequency [Hz]")
plt.ylabel("KEmax [eV]")

plt.errorbar(freq, ke_rough, fmt=".")
plt.plot(freq, photoelec_yfit, label="Linear fit")

print("Rough calc")
print(f"Stopping voltages: {ke_rough}")
print(f"h: {photoelec_coef[0]:.3e} eV*s")
print(f"phi: {-1 * photoelec_coef[1]:.5f} eV")
print(f"f0: {f_cutoff:.3e} Hz")

# more careful analysis
ke, ke_err = [], []
for V, I in zip(v, i):
    mask = (V < 0) & (I > 4) & (I <= 300)
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
