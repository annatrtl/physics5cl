import numpy as np
import matplotlib.pyplot as plt

# data
wl_err = 10 # nm
c = 3e8 # m/s

blue_wl = 462 # nm
blue_v = np.array([0, -500, -700, -900, -930, -951, -960, -955, 200, 400, 450]) # mV
blue_curr = np.array([(8360+8350)/2, (1890+1870)/2, (580+575)/2, (51+47)/2, (25+18)/2, (12+5)/2, 0, (3+0)/2, 1000 * (11.73+11.72)/2, 1000 * (15.15+15.14)/2, 1000 * (16.00+15.99)/2]) # nA

green_wl = 533 # nm
green_v = np.array([0, -200, -300, -400, -500, -600, -700, -721, -724, -722, 0, 200, 400, 450]) # mV
green_curr = np.array([(3265+3287)/2, (1652+1663)/2, (1007+1011)/2, (521+524)/2, 200, (58+65)/2, (7+10)/2, (3+7)/2, (0+3)/2, (0+3)/2, (3320+3330)/2, (5229+5240)/2, (7172+7180)/2, (7655+7658)/2]) # nA

red_wl = 618 # nm
red_v = np.array([0, -100, -200, -300, -350, -355, -360, -361, 0, 100, 200, 250]) # mV
red_curr = np.array([(819+823)/2, (379+393)/2, (112+123)/2, (14+18)/2, (3+7)/2, (3+10)/2, (0+3)/2, (0+3)/2, (849+852)/2, (1370+1377)/2, (1953+1957)/2, (2241+2244)/2]) # nA

# linear fit for each LED
blue_coef = np.polyfit(blue_curr, blue_v, 1) 

green_coef = np.polyfit(green_curr, green_v, 1) 

red_coef = np.polyfit(red_curr, red_v, 1) 

# KEmax = hc / wl - phi rough plot/fit
ke_max = np.array([blue_coef[0], green_coef[0], red_coef[0]])
wls = np.array([blue_wl, green_wl, red_wl]) * 1e-9
freq = c / wls
photoelec_coef = np.polyfit(freq, ke_max, 1)
photoelec_lin_func = np.poly1d(photoelec_coef)
photoelec_yfit = photoelec_lin_func(freq)

f_cutoff = -1 * photoelec_coef[1] / photoelec_coef[0]

fig = plt.figure()
plt.title("Photoelectric fit")
plt.xlabel("Frequency [1/s]")
plt.ylabel("KEmax [eV]")

plt.errorbar(freq, ke_max, fmt=".")
plt.plot(freq, photoelec_yfit, label="idk")

print(f"Plank's const: {-1 * photoelec_coef[0]} eV*s; work func: {photoelec_coef[1]} eV; cutoff frequency: {f_cutoff} Hz")

# more careful analysis
ke_err = np.array([0.009, 0.003, 0.001]) # found through looking at eV steps around I=0 value
freq_err = (c / (wls ** 2)) * wl_err
degree = 1
weights = 1 / ke_err
coef, cov = np.polyfit(freq, ke_max, degree, w=weights, cov=True)
h = coef[0]
phi = coef[1]
h_err = np.sqrt(cov[0, 0])
phi_err = np.sqrt(cov[1, 1])
freq_cutoff = phi / h
freq_cutoff_err = freq_cutoff * np.sqrt((phi_err / phi) ** 2 + (h_err / h) ** 2)
print(f"plank's const: {h} +/- {h_err} eV * s; work func: {phi} +/- {phi_err} eV; cutoff frequency: {freq_cutoff} +/- {freq_cutoff_err} Hz")  

plt.legend()
plt.show()
