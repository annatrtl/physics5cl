import numpy as np
import math

# data
purple = np.array([12.60, 12.55]) # cm
turq = np.array([14.05, 14.05]) # cm
red = np.array([19.85, 19.45]) # cm

L = 52.10 # cm
thickness = 0.25 # cm
d = 1/5000

# rough wavelengths
rough_p_wl0 = d * (purple[0]/math.sqrt(purple[0] ** 2 + L ** 2))
rough_p_wl1 = d * (purple[1]/math.sqrt(purple[1] ** 2 + L ** 2))

rough_t_wl0 = d * (turq[0]/math.sqrt(turq[0] ** 2 + L ** 2))
rough_t_wl1 = d * (turq[1]/math.sqrt(turq[1] ** 2 + L ** 2))

rough_r_wl0 = d * (red[0]/math.sqrt(red[0] ** 2 + L ** 2))
rough_r_wl1 = d * (red[1]/math.sqrt(red[1] ** 2 + L ** 2))

print(f"Wavelengths without uncertainty --> purple: {rough_p_wl0 * 1e7:.3f} nm and {rough_p_wl1 * 1e7:.3f} nm, turquoise: {rough_t_wl0 * 1e7:.3f} nm and {rough_t_wl1 * 1e7:.3f} nm, red: {rough_r_wl0 * 1e7:.3f} nm and {rough_r_wl1 * 1e7:.3f} nm")

# predicted wavelengths w/ Rydberg formula
R = 1.097e7
pred_p_wl = ((R * (1/(2**2) - 1/(5**2))) ** -1) * 1e9 # nm

pred_t_wl = ((R * (1/(2**2) - 1/(4**2))) ** -1) * 1e9 # nm

pred_r_wl = ((R * (1/(2**2) - 1/(3**2))) ** -1) * 1e9 # nm

print(f"Predicted wavelengths w/ Rydberg --> purple: {pred_p_wl:.3f} nm, turquoise: {pred_t_wl:.3f} nm, red: {pred_r_wl:.3f} nm") 




