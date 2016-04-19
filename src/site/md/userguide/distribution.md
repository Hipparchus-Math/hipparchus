# 8 Probability Distributions
## 8.1 Overview
The distributions package provides a framework and implementations for some commonly used
probability distributions. Continuous univariate distributions are represented by implementations of
the [RealDistribution](../apidocs/org/hipparchus/distribution/RealDistribution.html)
interface.  Discrete distributions implement
[IntegerDistribution](../apidocs/org/hipparchus/distribution/IntegerDistribution.html)
(values must be mapped to integers) and there is an
[EnumeratedDistribution](../apidocs/org/hipparchus/distribution/EnumeratedDistribution.html)
class representing discrete distributions with a finite, enumerated set of values.  Finally, multivariate
real-valued distributions can be represented via the
[MultivariateRealDistribution](../apidocs/org/hipparchus/distribution/MultiVariateRealDistribution.html)
interface.

An overview of available continuous distributions:<br/>
![Overview of continuous distributions](../images/userguide/real_distribution_examples.png)


## 8.2 Distribution Framework
The distribution framework provides the means to compute probability density
functions (`density(&middot;)`), probability mass functions
(`probability(&middot;)`) and distribution functions
(`cumulativeProbability(&middot;)`) for both
discrete (integer-valued) and continuous probability distributions.
The framework also allows for the computation of inverse cumulative probabilities
and sampling from distributions.

For an instance `f` of a distribution `F`,
and a domain value, `x`, `f.cumulativeProbability(x)`
computes `P(X &lt;= x)` where `X` is a random variable distributed
as `f`, i.e., `f.cumulativeProbability(&middot;)` represents
the distribution function of `f`. If `f` is continuous,
(implementing the `RealDistribution` interface) the probability density
function of `f` is represented by `f.density(&middot;)`.
For discrete `f` (implementing `IntegerDistribution`), the probability
mass function is represented by `f.probability(&middot;)`.  Continuous
distributions also implement `probability(&middot;)` with the same
definition (`f.probability(x)` represents `P(X = x)`
where `X` is distributed as `f`), though in the continuous
case, this will usually be identically 0.


    TDistribution t = new TDistribution(29);
    double lowerTail = t.cumulativeProbability(-2.656);     // P(T(29) &lt;= -2.656)
    double upperTail = 1.0 - t.cumulativeProbability(2.75); // P(T(29) &gt;= 2.75)
All distributions implement a `sample()` method to support random sampling from the
distribution. Implementation classes expose constructors allowing the default
[RandomGenerator](../apidocs/org/hipparchus/random/RandomGenerator.html)
used by the sampling algorithm to be overridden.  If sampling is not going to be used, providing
a null `RandomGenerator` constructor argument will avoid the overhead of initializing
the default generator.

Inverse distribution functions can be computed using the
`inverseCumulativeProbability` methods.  For continuous `f`
and `p` a probability, `f.inverseCumulativeProbability(p)` returns
* inf{x in R | P(X≤x) ≥ p} for 0 &lt; p &lt; 1},
* inf{x in R | P(X≤x) &gt; 0} for p = 0}.

For discrete `f`, the definition is the same, with `Z` (the integers)
in place of `R`.  Note that in the discrete case, the &ge; in the definition
can make a difference when `p` is an attained value of the distribution.


<!--
TODO: add section on multivariate distributions
-->
## 8.3 User Defined Distributions
User-defined distributions can be implemented using
[RealDistribution](../apidocs/org/hipparchus/distribution/RealDistribution.html),
[IntegerDistribution](../apidocs/org/hipparchus/distribution/IntegerDistribution.html) and
[MultivariateRealDistribution](../apidocs/org/hipparchus/distribution/MultivariateRealDistribution.html)
interfaces serve as base types.  These serve as the basis for all the distributions directly supported by
Apache Commons Math.  To aid in implementing distributions,
the [AbstractRealDistribution](../apidocs/org/hipparchus/distribution/AbstractRealDistribution.html),
[AbstractIntegerDistribution](../apidocs/org/hipparchus/distribution/AbstractIntegerDistribution.html) and
[AbstractMultivariateRealDistribution](../apidocs/org/hipparchus/distribution/AbstractMultivariateRealDistribution.html)
provide implementation building blocks and offer basic distribution functionality.
By extending these abstract classes directly, much of the repetitive distribution
implementation is already developed and should save time and effort in developing
user-defined distributions.


