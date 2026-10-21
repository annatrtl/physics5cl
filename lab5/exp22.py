import numpy as np
import matplotlib.pyplot as plt
import math

# data - exp1
d3 = np.array([12.20, 12.55, 12.60, 11.90, 12.20*36.40/36.60, 12.25*36.40/36.60]) # cm

L = 36.40 # cm
thickness = 0.25 # cm
d = 1/5000

# data - exp2
purple = np.array([12.60, 12.55]) # cm
turq = np.array([14.05, 14.05]) # cm
red = np.array([19.85, 19.45]) # cm

L1 = 52.10 # cm

# standard deviation & standard error
mean = np.mean(d3)
stdev = np.sqrt(np.sum((d3 - mean) ** 2) / (len(d3) - 1))
stderr = stdev / math.sqrt(len(d3))
print(f"Standard deviation: {stdev:.3f}, standard error: {stderr:.3f}")

# scatter plot
ni = np.array([5, 4, 3]) # purple, turqouise, red
dist = np.array([np.mean(purple), np.mean(turq), np.mean(red)])
wls = d * dist / np.sqrt(dist ** 2 + L1 ** 2)
wlserr = d*L1**2/(dist**2 + L1**2)**1.5 * stderr

# linear regression
Racc = 1.097e5 # 1/cm
b = Racc / 4
x = 1/ni**2
y = 1/wls
yerror = wlserr / wls ** 2
w = 1 / yerror ** 2

m = np.sum(w*x*(y-b)) / np.sum(w*x**2)
merr=  1 / np.sqrt(np.sum(w*x**2))
R = -m
Rerr = merr

chi_sq = np.sum(w*(y-b-m*x)**2)
dof = len(x) - 1

print(f"R = {R:.4e} ± {Rerr:.2e} m^-1, chi_squared/dof = {chi_sq/dof:.2f}")
print(f"t = {abs(R-Racc)/Rerr:.2f}, % diff = {100*(R-Racc)/Racc:.2f}%")

# plot
fig = plt.figure()
plt.xlabel("1/ni**2")
plt.ylabel("1/wavelength [1/cm^2]")
plt.title("")
plt.errorbar(x, y, yerr=yerror, label="data", fmt=".")
plt.plot(x, b + m*x, label="weighted fit")

plt.legend()
plt.show()
