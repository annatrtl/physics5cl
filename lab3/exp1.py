import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# 1.1 
L = 87.10 # cm
w = 0.65 # cm
wavelength = 6.35e-5 # cm

a = (2 * L * wavelength) / w
a_err = a * np.sqrt(((0.05 / L) ** 2) + ((0.05 / w) ** 2))

print(f"Slit width: {a:.5f} +/- {a_err:.5f} cm")

# 1.3
mth_order = np.array([1, 2, 3, 4])
d_min_to_min = np.array([0.45, 1.05, 1.60, 2.10]) # cm
calc_first_minima = d_min_to_min / (2 * mth_order)
print(calc_first_minima)

# 1.4
m = 6
d = 2.8 # cm
a = 4 * 0.04393 * 0.1 # cm
L = 84.50 # cm
wavelength = (a * d) / (2 * m * L)
w_err = wavelength * np.sqrt(((0.05 / d) ** 2) + ((0.05 / L) ** 2))
w_acc = 5.32e-5 # cm
print(f"Calculated wavelength: {wavelength:.9f} +/- {w_err:.9f} cm ; accepted wavelength: {w_acc:.9f} cm")

# 1.4 agreement test
print("Agreement test: |exp - acc| <= err")
print(f"Agreement test: |{wavelength:.9f} - {w_acc:.9f}| <= {w_err:.9f}; |{abs(wavelength - w_acc):.9f}| <= {w_err:.9f}")
print(f"Values agree? {abs(wavelength - w_acc) < w_err}")
