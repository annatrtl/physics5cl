import matplotlib.pyplot as plt
import numpy as np

# data 
inc = np.array([45.5, 40.0, 19.5, 30.0, 49.5, 60.0, 70.0, 80.0, 10.0, 25.0]) # degrees
ref = np.array([45.0, 40.5, 20.0, 30.0, 50.0, 60.0, 69.5, 80.0, 11.0, 26.0]) # degrees
error = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]) # degrees

# plot
coef = np.polyfit(inc, inc, 1) # b) compare w/ law of reflection line
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
