<!--
 Licensed to the Hipparchus project under one or more
 contributor license agreements.  See the NOTICE file distributed with
 this work for additional information regarding copyright ownership.
 The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 (the "License"); you may not use this file except in compliance with
 the License.  You may obtain a copy of the License at

      https://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
-->
# Curve Fitting
## Overview
The fitting package deals with curve fitting for univariate real functions.
When a univariate real function y = f(x) does depend on some unknown parameters
p<sub>0</sub>, p<sub>1</sub> ... p<sub>n-1</sub>, curve fitting can be used to
find these parameters. It does this by <em>fitting</em> the curve so it remains
very close to a set of observed points (x<sub>0</sub>, y<sub>0</sub>),
(x<sub>1</sub>, y<sub>1</sub>) ... (x<sub>k-1</sub>, y<sub>k-1</sub>). This
fitting is done by finding the parameters values that minimizes the objective
function Σ(y<sub>i</sub> - f(x<sub>i</sub>))<sup>2</sup>. This is actually a
least-squares problem.

For all provided curve fitters, the operating principle is the same.
Users must

* create an instance of the fitter using the `create` factory method of the appropriate class,
* call the [fit](../apidocs/org/hipparchus/fitting/AbstractCurveFitter) with a `Collection` of [ observed data points](../apidocs/org/hipparchus/fitting/WeightedObservedPoint.html) as argument, which will return an array with the parameters that best fit the given data points.

The list of observed data points to be passed to `fit` can be built by incrementally
adding instances to an instance of [WeightedObservedPoints](../apidocs/org/hipparchus/fitting/WeightedObservedPoints.html),
and then retrieve the list of `WeightedObservedPoint` by calling the `toList`
method on that container.
A weight can be associated with each observed point; it allows to take into account uncertainty
on some points when they come from noisy measurements for example. If no such information exist and
all points should be treated the same, it is safe to put 1.0 as the weight for all points (and this
is the default when no weight is passed to the `add`.


Some fitters require that initial values for the parameters are provided by the user,
through the `withStartPoint` method, before attempting to perform the fit.
When that's the case the fitter class usually defines a guessing procedure within a
`ParameterGuesser` inner class, that attempts to provide appropriate initial
values based on the user-supplied data.
When initial values are required but are not provided, the `fit` method will
internally call the guessing procedure.

## Implemented Functions

Fitting of specific functions are provided through the following classes:

* create an instance of the fitter using the `create` factory method of the appropriate class,
* call the [fit](../apidocs/org/hipparchus/fitting/AbstractCurveFitter) with a `Collection` of [ observed data points](../apidocs/org/hipparchus/fitting/WeightedObservedPoint.html) as argument, which will return an array with the parameters that best fit the given data points.

The following example shows how to fit data with a polynomial function.

    // Collect data.
    final WeightedObservedPoints obs = new WeightedObservedPoints();
    obs.add(-1.00, 2.021170021833143);
    obs.add(-0.99, 2.221135431136975);
    obs.add(-0.98, 2.09985277659314);
    obs.add(-0.97, 2.0211192647627025);
    // ... Lots of lines omitted ...
    obs.addt(0.99, -2.4345814727089854);
    
    // Instantiate a third-degree polynomial fitter.
    final PolynomialCurveFitter fitter = PolynomialCurveFitter.create(3);
    
    // Retrieve fitted parameters (coefficients of the polynomial function).
    final double[] coeff = fitter.fit(obs.toList());


## General Case

The [AbstractCurveFitter](../apidocs/org/hipparchus/fitting/AbstractCurveFitter.html)
class provides the basic functionality for implementing other curve fitting classes.
Users must provide their own implementation of the curve template as a class that implements
the [ParametricUnivariateFunction](../apidocs/org/hipparchus/analysis/ParametricUnivariateFunction.html)
interface.
