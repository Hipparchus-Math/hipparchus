#
# Licensed to the Hipparchus project under one or more
# contributor license agreements.  See the NOTICE file distributed with
# this work for additional information regarding copyright ownership.
# The Hipparchus project licenses this file to You under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with
# the License.  You may obtain a copy of the License at
#
#      https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
 */
import numpy as np
from scipy.stats import halfcauchy
#
# Set the precision for printing
#
np.set_printoptions(precision=16)
#
# Half-Cauchy distribution
# https://docs.scipy.org/doc/scipy-1.16.1/reference/generated/scipy.stats.halfcauchy.html
#
print("\n\n\n")
print("Default Half-Cauchy distribution")
pset = [0.001, 0.01, 0.025, 0.05, 0.1, 0.999,0.990, 0.975, 0.950, 0.90]
lb, ub = halfcauchy.support()
vals   = halfcauchy.ppf(pset) # same as quantiles
print("Lower bound for the support: ",lb)
print("Upper bound for the support: ",ub)
print("Probability set: ")
print(pset)
print("Related percentiles: ")
print(vals)
print("Cumulative distribution: ")
print(halfcauchy.cdf(vals))
print("Probability density: ")
print(halfcauchy.pdf(vals))
print("Check if the results are consistent: ",np.allclose(pset, halfcauchy.cdf(vals)))

