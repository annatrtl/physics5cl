import matplotlib.pyplot as plt
import numpy as np
import math

# data 
o_dist = np.array([31.10, 38.90, 35.75, 54.20, 44.50]) # object dist (cm)
i_dist = np.array([44.65, 36.20, 31.30, 30.40, 38.10]) # image dist (cm)

# plot calcs
inv_d_o = 1 / o_dist
inv_d_i = 1 / i_dist

#coef = np.polyfit(inv_d_o, inv_d_i, 1)
#lin_func = np.poly1d(coef)
#yfit = lin_func(inv_d_o)

slope = -1
inter = np.mean(inv_d_i - (slope * inv_d_o))
yfit = (slope * inv_d_o) + inter

# plot
fig = plt.figure()
plt.xlabel("Inverse d_o [cm^-1]")
plt.ylabel("Inverse d_i [cm^-1]")
plt.title("Inverse d_i vs. inverse d_o")

plt.errorbar(inv_d_o, inv_d_i, fmt=".")
plt.plot(inv_d_o, yfit, label=f"y={slope:.3f}x + {inter:.3f}", color="red") 

"""
# agreement test
exp_f = 1 / coef[1]
acc_f = -15
error = 0.05

print("Agreement test: |exp - acc| <= 2 * sqrt(err^2 + err^2)")
print(f"Agreement test: |{exp_f:.3f} - {acc_f}| < = 2 * sqrt({error}^2 + {error}^2)")
print(f"Values agree? {abs(exp_f - acc_f):.3f} / {2*math.sqrt(2 * (error ** 2)):.3f} ; Values agree? {(abs(exp_f - acc_f) / (2* math.sqrt(2 * (error ** 2))) < 1)}")

"""


plt.legend()
plt.show()






