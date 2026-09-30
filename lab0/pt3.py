import matplotlib.pyplot as plt
import numpy as np

# data 
do = np.array([34.95, 33.85, 32.05, 30.40, 27.80, 26.15]) # cm
di = np.array([27.85, 28.65, 30.35, 33.10, 34.85, 36.95]) # cm
hi = np.array([1.20, 1.30, 1.40, 1.60, 1.85, 2.10]) # cm
# y = -x + b ; x = 1/do ; y = 1/di
y = 1 / do
x = 1 / di

# a) covariance and coef of linear correlation
mx = np.mean(x)
my = np.mean(y)
co_var = np.sum((x - mx) * (y - my)) / (len(x) - 1)

var_x = np.sqrt(np.sum((x - mx) ** 2) / (len(x) - 1))
var_y = np.sqrt(np.sum((y - my) ** 2) / (len(x) - 1))
lin_corr = co_var / (var_x * var_y)
print(f"cv: {co_var}, lc: {lin_corr}")

# b) ay_equiv
error = 0.05 #cm
ax = error / (do ** 2)
ay = error / (di ** 2)
ay_equiv = np.sqrt((ay ** 2) + (ax ** 2))
print(f"ay_equiv: {ay_equiv}")

# c) weighted least-squares approach
w = 1 / (ay_equiv ** 2)
m_x_weighed = np.sum(w * x) / np.sum(w)
m_y_weighed = np.sum(w * y) / np.sum(w)

b_fit = m_x_weighed + m_y_weighed
a_b = (np.sum(w) ** (-1 / 2))
print(f"b: {b_fit}, a_b: {a_b}") 

# d)
y_exp = -x + b_fit
chi_sq = np.sum(((y - y_exp) / ay_equiv) ** 2)
red_chi_sq = chi_sq / (len(x) - 1) #parameters for dof = 1 bc slope fixed
print(f"reduced chi sq: {red_chi_sq}")

# e) focal length
f = 1 / b_fit
a_f = a_b / (b_fit ** 2)
print(f"f is: {f} +/- {a_f}")

# f) magnification eq
ratio = di / do
ho_fit = np.sum(ratio * hi) / np.sum(ratio ** 2)
a_m = ratio * np.sqrt((error / do) ** 2 + (error / di) ** 2)
a_total = np.sqrt(error ** 2 + (ho_fit * a_m) ** 2)
hi_exp = ho_fit * ratio

chi_sq = np.sum(((hi - hi_exp) / error) ** 2)
red_chi_sq = chi_sq / (len(hi) - 1)
print(f"reduced chi squared: {red_chi_sq}")


# plot
"""
fig = plt.figure()
plt.xlabel("Incident angle [Degrees]")
plt.ylabel("Reflected angle [Degrees]")
#plt.title("Reflected vs incident angles for law of reflection experiment")
plt.title("Reflected vs incident angles with expected law of reflection line") # compared w/ law of reflection line

#coef = np.polyfit(inc, ref, 1)
coef = np.polyfit(inc, inc, 1) # b) compare w/ law of reflection line
lin_func = np.poly1d(coef)
yfit = lin_func(inc)

plt.errorbar(inc, ref, xerr=error, yerr=error, fmt=".")
#plt.plot(inc, yfit, label=f"y={coef[0]:.3f} + {coef[1]:.3f}")
plt.plot(inc, yfit, label=f"law of reflection; y={coef[0]:.3f} + {coef[1]:.3f}", color="red") # compared w/ law of reflection line
text = f"Error bars: 0.5 degrees"
fig = fig.text(0.5, 0.02, text, wrap=True, horizontalalignment='center')


# reduced chi squared calc
unc = 0.5
chi_sq = np.sum(((ref - yfit) / unc) ** 2)
red_chi_sq = chi_sq / (len(inc) - 2)
print(f"Reduced chi squared: {red_chi_sq:.9f}")

plt.legend()
plt.show()
"""
