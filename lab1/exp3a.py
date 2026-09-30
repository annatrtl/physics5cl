import matplotlib.pyplot as plt
import numpy as np
import math

# data 
d_o_dist = np.array([5.10, 5.75, 7.40, 8.05, 8.90]) # lens D object dist (cm)
a_d_dist = np.array([26.15, 26.55, (32.25+27.10)/2, (32.20+27.40)/2, (30.10+26.50)/2]) # dist between lens D and A (cm)
a_i_dist = np.array([40.45, 39.50, (31.80+37.00)/2, (32.45+36.40)/2, (32.65+36.25)/2]) # lens A image dist (cm)
a_f = 18 # focal length of lens A (cm)

# analysis
a_o_dist = (a_f * a_i_dist / (a_i_dist - a_f)) # lens A object dist (cm)
d_i_dist = -1 * np.abs(a_d_dist - a_o_dist) # lens D image dist (cm)

# plot calcs
inv_d_o = 1 / d_o_dist
inv_d_i = 1 / d_i_dist

coef = np.polyfit(inv_d_o, inv_d_i, 1)
lin_func = np.poly1d(coef)
yfit = lin_func(inv_d_o)

# plot
fig = plt.figure()
plt.xlabel("Inverse d_o [cm^-1]")
plt.ylabel("Inverse d_i [cm^-1]")
plt.title("Inverse d_i vs. inverse d_o")

plt.errorbar(inv_d_o, inv_d_i, fmt=".")
plt.plot(inv_d_o, yfit, label=f"y={coef[0]:.3f}x + {coef[1]:.3f}", color="red") 


# agreement test
exp_f = 1 / coef[1]
acc_f = -15
error = 0.05

print("Agreement test: |exp - acc| <= 2 * sqrt(err^2 + err^2)")
print(f"Agreement test: |{exp_f:.3f} - {acc_f}| < = 2 * sqrt({error}^2 + {error}^2)")
print(f"Values agree? {abs(exp_f - acc_f):.3f} / {2*math.sqrt(2 * (error ** 2)):.3f} ; Values agree? {(abs(exp_f - acc_f) / (2* math.sqrt(2 * (error ** 2))) < 1)}")

# residuals
fig = plt.figure()
plt.xlabel("Inverse d_o [cm^-1]")
plt.ylabel("Residuals for inverse d_i [cm^-1]")
plt.title("Residuals")

plt.scatter(inv_d_o, yfit - inv_d_i)

plt.legend()
plt.show()






