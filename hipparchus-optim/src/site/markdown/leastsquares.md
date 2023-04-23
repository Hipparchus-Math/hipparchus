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
# Least squares

## Overview

The least squares package fits a non-linear parametric model to a set of observed
values by minimizing a cost function with a specific form.
The fitting consists in finding the values
for some parameters \\( p \\) that minimize a cost function
$$ J = \sum_i w_i(o_i - f_i)^2 $$ where \\( w_i \\) are the weights, \\( o_i \\) are the
observed values and \\( f_i = f_i( p ) \\) are the values computed from the model.
The various \\( r_i = o_i - f_i \\)
terms are the residuals which quantify the deviation between the set of
observed values and theoretical values computed from
measurement models depending on the free parameters.

Among other properties, the non-linear least squares is the
Maximum Likelihood Estimator (MLE) if the observations have normally
distributed errors and the weights are set to the reciprocal of their
variance. That is $$ o_i = N(f_i(p), \sigma) $$ and
$$ w_i = \frac{1}{\sigma^2} $$ where \\( N(\\mu, \\sigma) \\) denotes a normal distribution
with mean \\( \\mu \\) and standard deviation \\( \\sigma \\). For more on non-linear least
squares, its theory, and its properties the reader is encouraged to consult a text book
such as [1](#ref1).

Two engines devoted to least-squares problems are available. The first one is
based on the [Gauss-Newton](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/GaussNewtonOptimizer.html)
method. The second one is the
[Levenberg-Marquardt](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/LevenbergMarquardtOptimizer.html) method.


## LeastSquaresBuilder and LeastSquaresFactory

In order to solve a least-squares fitting problem, the user must provide the following elements:

* an implementation of the measurement model \\( f(p) \\) and its Jacobian, \\( \\frac{\\partial f}{\\partial p} \\). This is best done by implementing [MultivariateJacobianFunction](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/MultivariateJacobianFunction.html).
* the observed (or target) values: \\( o \\).
* the start values for all parameters: \\( s \\).
* optionally a validator for the parameters \\( p \\). This is an implementation of [ParameterValidator](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/ParameterValidator.html).
* optionally weights for sample point: \\( w \\), this defaults to 1.0 if not provided.
* a maximum number of iterations.
* a maximum number of model evaluations, which may be different than iterations for Levenberg-Marquardt.
* a convergence criterion, which is an implementation of [ConvergenceChecker&lt;Evaluation&gt;](../apidocs/org/hipparchus/optim/ConvergenceChecker.html)

The elements of the list above can be provided as an implementation of the
[LeastSquaresProblem](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/LeastSquaresProblem.html)
interface. However, this may be cumbersome to do directly, so some
helper classes are available. The first helper is a mutable builder:
[LeastSquaresBuilder](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/LeastSquaresBuilder.html).
The second helper is an utility factory:
[LeastSquaresFactory](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/LeastSquaresFactory.html).

The builder class is better suited when setting the various elements of the least squares
problem is done progressively in different places in the user code. In this case, the user
would create first an empty builder and configure it progressively by calling its methods
(`start`, `target`, `model`, ...). Once the configuration
is complete, calling the `build` method would create the least squares problem.

The factory utility is better suited when the various elements of the least squares
problem are all known at one place and the problem can be built in just one sweep, calling
to one of the static `LeastSquaresFactory.create` method.


## Model Function

The model function is used by the least squares engine to evaluate the model
\\( f(p) \\). It is therefore a multivariate
function (it depends on the various \\( p_k \\)) and it is vector-valued (it has several
components \\( f_i \\)). There must be exactly one model function \\( f_i \\) for
each observed (or target) value \\( o_i \\). In order for the problem to be well defined, the
number of parameters must be less than the number of observations.
Failing to ensure this may lead to the engine throwing an exception as the underlying linear
algebra operations may encounter singular matrices. It is not unusual to have a large number
of observations (several thousands) and only a dozen parameters. There are no limitations on these
numbers, though.

As the least squares engine uses the Jacobians matrices of the model function, both
its value and its derivatives <em>with respect to the parameters, \\( p \\)</em>,
must be available. There are two ways to provide this:

* an implementation of the measurement model \\( f(p) \\) and its Jacobian, \\( \\frac{\\partial f}{\\partial p} \\). This is best done by implementing [MultivariateJacobianFunction](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/MultivariateJacobianFunction.html).
* the observed (or target) values: \\( o \\).
* the start values for all parameters: \\( s \\).
* optionally a validator for the parameters \\( p \\). This is an implementation of [ParameterValidator](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/ParameterValidator.html).
* optionally weights for sample point: \\( w \\), this defaults to 1.0 if not provided.
* a maximum number of iterations.
* a maximum number of model evaluations, which may be different than iterations for Levenberg-Marquardt.
* a convergence criterion, which is an implementation of [ConvergenceChecker&lt;Evaluation&gt;](../apidocs/org/hipparchus/optim/ConvergenceChecker.html)

The first alternative is best suited for models which are not computationally intensive
as it allows more modularized code with one method for each type of computation. The second
alternative is best suited for models which are computationally intensive and evaluating
both the values and derivatives in one sweep saves a lot of work.

The `point` parameter of the `value` methods in the
[MultivariateVectorFunction](../apidocs/org/hipparchus/analysis/MultivariateVectorFunction.html),
[MultivariateMatrixFunction](../apidocs/org/hipparchus/analysis/MultivariateMatrixFunction.html),
or [MultivariateJacobianFunction](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/MultivariateJacobianFunction.html)
interfaces will contain the parameter vector \\( p \\). The values will be the
model values \\( f(p) \\) and the derivatives will be the derivatives of
the model values with respect to the parameters
\\( \\frac{\\partial f(p)}{\\partial p} \\).

There are no requirements on how to compute value and derivatives. The
[DerivativeStructure](../apidocs/org/hipparchus/analysis/differentiation/DerivativeStructure.html)
class may be useful to compute analytical derivatives without
pencil and paper, but this class is not mandated by the API which only expects the
derivatives as a Jacobian matrix containing primitive double entries.

One non-obvious feature provided by both the builder and the factory is lazy evaluation. This feature
allows to defer calls to the model functions until they are really needed by the engine. This
can save some calls for engines that evaluate the value and the Jacobians in different loops
(this is the case for Levenberg-Marquardt). However, lazy evaluation is possible <em>only</em>
if the model functions are themselves separated, i.e. it can be used only with the first
alternative above. Setting up the `lazyEvaluation` flag to `true` in the builder
or factory and setting up the model function as one
[MultivariateJacobianFunction](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/MultivariateJacobianFunction.html)
instance at the same time will trigger an illegal state exception telling that the model function
misses required functionality.


## Parameters Validation

In some cases, the model function requires parameters to lie within a specific domain. For example
a parameter may be used in a square root and needs to be positive, or another parameter represents
the sine of an angle and should be within -1 and +1, or several parameters may need to remain in
the unit circle and the sum of their squares must be smaller than 1. The least square solvers available
in Hipparchus currently don't allow to set up constraints on the parameters. This is a
known missing feature. There are two ways to circumvent this.

Both ways are achieved by setting up a
[ParameterValidator](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/ParameterValidator.html)
instance. The input of the value and jacobian model functions will always be the output of
the parameter validator if one exists.

One way to constrain parameters is to use a continuous mapping between the parameters that the
least squares solver will handle and the real parameters of the mathematical model. Using mapping
functions like `logit` and `sigmoid`, one can map a finite range to the
infinite real line. Using mapping functions based on `log` and `exp`, one
can map a semi-infinite range to the infinite real line. It is possible to use such a mapping so
that the engine will always see unbounded parameters, whereas on the other side of the mapping the
mathematical model will always see parameters mapped correctly to the expected range. Care must be
taken with derivatives as one must remember that the parameters have been mapped. Care must also
be taken with convergence status. This may be tricky.

Another way to constrain parameters is to simply truncate the parameters back to the domain when
one search point escapes from it and not care about derivatives. This works <em>only</em> if the
solution is expected to be inside the domain and not at the boundary, as points out of the domain
will only be temporary test points with a cost function higher than the real solution and will soon
be dropped by the underlying engine. As a rule of thumb, these conditions are met only when the
domain boundaries correspond to unrealistic values that will never be achieved (null distances,
negative masses, ...) but they will not be met when the domain boundaries are more operational
limits (a maximum weight that can be handled by a device, a minimum temperature that can be
sustained by an instrument, ...).


## Tuning

Among the elements to be provided to the least squares problem builder or factory
are some tuning parameters for the solver.

The maximum number of iterations refers to the engine algorithm main loop, whereas the
maximum number of evaluations refers to the number of calls to evaluate the model. Some
algorithms (like Levenberg-Marquardt) have two embedded loops, with iteration number
being incremented at outer loop level, but a new evaluation being done at each inner
loop. In this case, the number of evaluations will be greater than the number of iterations.
Other algorithms (like Gauss-Newton) have only one level of loops. In this case, the
number of evaluations will equal to the number of iterations. In any case, the maximum
numbers are really only intended as safeguard to prevent infinite loops, so the exact
value of the limit is not important so it is common to select some almost arbitrary number
much larger than the expected number of evaluations and use it for both
`maxIterations` and `maxEvaluations`. As an example, if the least
squares solver usually finds a solution in 50 iterations, setting a maximum value to 1000
is probably safe and prevents infinite loops. If the least squares solver needs several
hundreds of evaluations, it would probably be safer to set the maximum value to 10000 or
even 1000000 to avoid failures in slightly more demanding cases. Very fine tuning of
these maximum numbers is often worthless, they are only intended as safeguards. One can
think of these parameters as the divergence criteria.

Convergence checking is delegated to a dedicated interface from the `optim`
package: [ConvergenceChecker](../apidocs/org/hipparchus/optim/ConvergenceChecker.html),
parameterized with either the specific
[Evaluation](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/LeastSquaresProblem.Evaluation.html)
class used for least squares problems or the general
[PointVectorValuePair](../apidocs/org/hipparchus/optim/PointVectorValuePair.html).
Each time convergence is checked, both the previous
and the current evaluations of the least squares problem are provided, so the checker can
compare them and decide whether convergence has been reached or not. The predefined convergence
checker implementations that can be useful for least squares fitting are:

* an implementation of the measurement model \\( f(p) \\) and its Jacobian, \\( \\frac{\\partial f}{\\partial p} \\). This is best done by implementing [MultivariateJacobianFunction](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/MultivariateJacobianFunction.html).
* the observed (or target) values: \\( o \\).
* the start values for all parameters: \\( s \\).
* optionally a validator for the parameters \\( p \\). This is an implementation of [ParameterValidator](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/ParameterValidator.html).
* optionally weights for sample point: \\( w \\), this defaults to 1.0 if not provided.
* a maximum number of iterations.
* a maximum number of model evaluations, which may be different than iterations for Levenberg-Marquardt.
* a convergence criterion, which is an implementation of [ConvergenceChecker&lt;Evaluation&gt;](../apidocs/org/hipparchus/optim/ConvergenceChecker.html)

Of course, users can also provide their own implementation of the
[ConvergenceChecker](../apidocs/org/hipparchus/optim/ConvergenceChecker.html) interface.


## Optimization Engine

Once the least squares problem has been created, using either the builder or the factory,
it is passed to an optimization engine for solving. Two engines devoted to least-squares
problems are available. The first one is
based on the [Gauss-Newton](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/GaussNewtonOptimizer.html)
method. The second one is the [Levenberg-Marquardt](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/LevenbergMarquardtOptimizer.html)
method. For both increased readability and in order to leverage
possible future changes in the configuration, it is recommended to use the fluent-style API to
build and configure the optimizers. This means creating a first temporary version of the optimizer
with a default parameterless constructor, and then to set up the various configuration parameters
using the available `withXxx` methods that all return a new optimizer instance. Only the
final fully configured instance is used. As an example, setting up a Levenberg-Marquardt with
all configuration set to default except the cost relative tolerance and parameter relative tolerance
would be done as follows:

      LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer().
                                        withCostRelativeTolerance(1.0e-12).
                                        withParameterRelativeTolerance(1.0e-12);

As another example, setting up a Gauss-Newton optimizer and forcing the decomposition to SVD (the
default is QR decomposition) would be done as follows:

      LeastSquaresOptimizer optimizer = new GaussNewtonOptimizer().
                                        withDecomposition(GaussNewtonOptimizer.Decomposition.SVD);


## Solving

Solving the least squares problem is done by calling the `optimize` method of the
optimizer and passing the least squares problem as the single parameter:

      LeastSquaresOptimizer.Optimum optimum = optimizer.optimize(leastSquaresProblem);

The
[LeastSquaresOptimizer.Optimum](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/LeastSquaresOptimizer.Optimum.html)
class is a specialized
[Evaluation](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/LeastSquaresProblem.Evaluation.html)
with additional methods te retrieve the number of evaluations and number of iterations performed.
The most important methods are inherited from the
[Evaluation](../apidocs/org/hipparchus/optim/nonlinear/vector/leastsquares/LeastSquaresProblem.Evaluation.html)
class and correspond to the point (i.e. the parameters), cost, Jacobian, RMS, covariance ...


## Example

The following simple example shows how to find the center of a circle of known radius to
to best fit observed 2D points. It is a simplified version of one of the JUnit test cases.
In the complete test case, both the circle center and its radius are fitted, here the
radius is fixed. From the optimizer's perspective the parameters are the location (x, y)
of the center of the circle, the observed values are all the assumed radius of the circle:
70.0 and the computed values are the distance from each of the observedPoints to center of
the circle (as given by the current value of the parameters).

      final double radius = 70.0;
      final Vector2D[] observedPoints = new Vector2D[] {
          new Vector2D( 30.0,  68.0),
          new Vector2D( 50.0,  -6.0),
          new Vector2D(110.0, -20.0),
          new Vector2D( 35.0,  15.0),
          new Vector2D( 45.0,  97.0)
      };
    
      // the model function components are the distances to current estimated center,
      // they should be as close as possible to the specified radius
      MultivariateJacobianFunction distancesToCurrentCenter = new MultivariateJacobianFunction() {
          public Pair&lt;RealVector, RealMatrix&gt; value(final RealVector point) {
    
              Vector2D center = new Vector2D(point.getEntry(0), point.getEntry(1));
    
              RealVector value = new ArrayRealVector(observedPoints.length);
              RealMatrix jacobian = new Array2DRowRealMatrix(observedPoints.length, 2);
    
              for (int i = 0; i &lt; observedPoints.length; ++i) {
                  Vector2D o = observedPoints[i];
                  Vector2D diff = center.subtract(o);
                  double modelI = diff.getNorm();
                  value.setEntry(i, modelI);
                  // derivative with respect to p0 = x center
                  jacobian.setEntry(i, 0, diff.getX() / modelI);
                  // derivative with respect to p1 = y center
                  jacobian.setEntry(i, 1, diff.getY() / modelI);
              }
    
              return new Pair&lt;RealVector, RealMatrix&gt;(value, jacobian);
    
          }
      };
    
      // the target is to have all points at the specified radius from the center
      double[] prescribedDistances = new double[observedPoints.length];
      Arrays.fill(prescribedDistances, radius);
    
      // least squares problem to solve : modeled radius should be close to target radius
      LeastSquaresProblem problem = new LeastSquaresBuilder().
                                    start(new double[] { 100.0, 50.0 }).
                                    model(distancesToCurrentCenter).
                                    target(prescribedDistances).
                                    lazyEvaluation(false).
                                    maxEvaluations(1000).
                                    maxIterations(1000).
                                    build();
      LeastSquaresOptimizer.Optimum optimum = new LevenbergMarquardtOptimizer().optimize(problem);
      Vector2D fittedCenter = new Vector2D(optimum.getPoint().getEntry(0), optimum.getPoint().getEntry(1));
      System.out.println("fitted center: " + fittedCenter.getX() + " " + fittedCenter.getY());
      System.out.println("RMS: "           + optimum.getRMS());
      System.out.println("evaluations: "   + optimum.getEvaluations());
      System.out.println("iterations: "    + optimum.getIterations());


## References

<a name="ref1"/>[1] Crassidis, John L., and John L. Junkins. Optimal estimation of dynamic
systems. Boca Raton, FL: CRC Press, 2012.
