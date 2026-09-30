import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# data
intensity = np.array([4648.823, 4105.042, 3210.140, 2021.742, 783.725, 212.835, 117.059]) # lux ; 0-6,000 lux range
int_err = np.array([1.815, 2.508, 2.092, 1.648, 1.742, 0.977, 0.227]) # lux
angle = np.array([0, 15, 30, 45, 60, 75, 90]) # degrees
angle_err = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])

# degrees to radians
angle_rad = np.radians(angle)
angle_rad_err = np.radians(angle_err)

# malus's law
def malus(theta, i0, theta0, ibg):
    return i0 * (np.cos(theta-theta0)) ** 2 + ibg

# initial fit for curve steepness (di_dtheta)
guesses = [4600, 0, 100]
popt_i, pcov_i = curve_fit(malus, angle_rad, intensity, sigma=int_err, absolute_sigma=True, p0=guesses)
i0_i, theta0_i, ibg_i = popt_i

# propagating uncertainty
di_dtheta = -2 * i0_i * np.cos(angle_rad - theta0_i) * np.sin(angle_rad - theta0_i)
eff_err = np.sqrt(int_err**2 + (di_dtheta * angle_rad_err)**2)

# final fit w/ effective uncertainties
popt, pcov = curve_fit(malus, angle_rad, intensity, sigma=eff_err, absolute_sigma=True, p0=popt_i)
io_fit, theta0_fit, ibg_fit = popt
perr = np.sqrt(np.diag(pcov))

# radians to degrees
theta0_deg = np.degrees(theta0_fit)
theta_err_deg = np.degrees(perr[1])

# continuous lines
angle_fit_deg = np.linspace(-5, 95, 200)
angle_fit_rad = np.radians(angle_fit_deg)
intensity_fit = malus(angle_fit_rad, *popt)

# reduced chi squared
residuals = intensity - malus(angle_rad, *popt)
chi_sq = np.sum((residuals / eff_err) ** 2)
red_chi_sq = chi_sq / (len(intensity) - 3)
print(f"Reduced chi squared: {red_chi_sq}")

# plot
fig = plt.figure()
plt.xlabel("Angle [degrees]")
plt.ylabel("Intensity [lux]")
plt.title("Intensity vs. Angle w/ Malus Law fit")

plt.errorbar(angle, intensity, xerr=angle_err, yerr=int_err, fmt=".")
plt.plot(angle_fit_deg, intensity_fit, label=f"Malus's law fit; reduced chi squared: {red_chi_sq:.3f}")
text = f"X error: 0.5 degrees, Y error: {[1.815, 2.508, 2.092, 1.648, 1.742, 0.977, 0.227]}"
fig = fig.text(0.5, 0.02, text, wrap=True, horizontalalignment='center')


plt.legend()
plt.show()
