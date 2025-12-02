import numpy as np
from scipy.stats import halfnorm
#
# Set the precision for printing
#
np.set_printoptions(precision=16)
#
# Half-Normal distribution
# https://docs.scipy.org/doc/scipy-1.16.1/reference/generated/scipy.stats.halfnorm.html
#
print("\n\n\n")
print("Default Half-Normal distribution")
pset = [0.001, 0.01, 0.025, 0.05, 0.1, 0.999,0.990, 0.975, 0.950, 0.90]
lb, ub = halfnorm.support()
vals   = halfnorm.ppf(pset) # same as quantiles
print("Lower bound for the support: ",lb)
print("Upper bound for the support: ",ub)
print("Probability set: ")
print(pset)
print("Related percentiles: ")
print(vals)
print("Cumulative distribution: ")
print(halfnorm.cdf(vals))
print("Probability density: ")
print(halfnorm.pdf(vals))
print("Check if the results are consistent: ",np.allclose(pset, halfnorm.cdf(vals)))

print("\n\n\n")
print("Default Half-Normal distribution")
vals = [0.5, 1.0, 2.0, 4.0]
print("Evaluation points: ")
print(vals)
print("Probability density: ")
print(halfnorm.pdf(vals))
mu  = halfnorm.mean()
sigma = halfnorm.std()
print("Mean: ",mu)
print("Standard deviation: ",sigma)
xset = [mu - sigma, mu, mu + sigma, mu + 2 * sigma,  mu + 3 * sigma, mu + 4 * sigma, mu + 5 * sigma];
vals = halfnorm.cdf(xset)
print("Evaluation points: ")
print(xset)
print("Cumulative distribution: ")
print(vals)

