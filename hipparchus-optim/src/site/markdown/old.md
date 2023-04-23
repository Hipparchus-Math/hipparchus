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
# Optimization
<em>The contents of this section currently describes deprecated classes.</em>
Please refer to the new [API description](../apidocs/org/hipparchus/optim/package-summary.html).

Least squares optimizers are not in this package anymore, they have been moved
in a dedicated least-squares sub-package described in the [least squares](./leastsquares.html)
section.


## Overview
The optimization package provides algorithms to optimize (i.e. either minimize
or maximize) some objective or cost function. The package is split in several
sub-packages dedicated to different kind of functions or algorithms.

* the univariate package handles univariate scalar functions,
* the linear package handles multivariate vector linear functions with linear constraints,
* the direct package handles multivariate scalar functions using direct search methods (i.e. not using derivatives),
* the general package handles multivariate scalar or vector functions using derivatives.
* the fitting package handles curve fitting by univariate real functions


The top level optimization package provides common interfaces for the optimization
algorithms provided in sub-packages. The main interfaces defines defines optimizers
and convergence checkers. The functions that are optimized by the algorithms provided
by this package and its sub-packages are a subset of the one defined in the
`analysis` package, namely the real and vector valued functions. These
functions are called objective function here. When the goal is to minimize, the
functions are often called cost function, this name is not used in this package.

The type of goal, i.e. minimization or maximization, is defined by the enumerated
[        GoalType](../apidocs/org/hipparchus/optimization/GoalType.html)
which has only two values: `MAXIMIZE` and `MINIMIZE`.

Optimizers are the algorithms that will either minimize or maximize, the objective
function by changing its input variables set until an optimal set is found. There
are only four interfaces defining the common behavior of optimizers, one for each
supported type of objective function:

* the univariate package handles univariate scalar functions,
* the linear package handles multivariate vector linear functions with linear constraints,
* the direct package handles multivariate scalar functions using direct search methods (i.e. not using derivatives),
* the general package handles multivariate scalar or vector functions using derivatives.
* the fitting package handles curve fitting by univariate real functions



Despite there are only four types of supported optimizers, it is possible to optimize
a transform a <a
href="../apidocs/org.hipparchus/analysis/MultivariateVectorFunction.html">
non-differentiable multivariate vectorial function</a> by converting it to a <a
href="../apidocs/org.hipparchus/analysis/MultivariateFunction.html">
non-differentiable multivariate real function</a> thanks to the <a
href="../apidocs/org.hipparchus/optimization/LeastSquaresConverter.html">
LeastSquaresConverter</a> helper class. The transformed function can be optimized using
any implementation of the <a
href="../apidocs/org.hipparchus/optimization/MultivariateOptimizer.html">
MultivariateOptimizer</a> interface.


For each of the four types of supported optimizers, there is a special implementation
which wraps a classical optimizer in order to add it a multi-start feature. This feature
call the underlying optimizer several times in sequence with different starting points
and returns the best optimum found or all optima if desired. This is a classical way to
prevent being trapped into a local extremum when looking for a global one.


## Univariate Functions
A [          UnivariateOptimizer](../apidocs/org/hipparchus/optimization/univariate/UnivariateOptimizer.html)
is used to find the minimal values of a univariate real-valued
function `f`.

These algorithms usage is very similar to root-finding algorithms usage explained
in the analysis package. The main difference is that the `solve` methods in root
finding algorithms is replaced by `optimize` methods.


## Linear Programming
This package provides an implementation of George Dantzig's simplex algorithm
for solving linear optimization problems with linear equality and inequality
constraints.


## Direct Methods
Direct search methods only use cost function values, they don't
need derivatives and don't either try to compute approximation of
the derivatives. According to a 1996 paper by Margaret H. Wright
([Direct Search Methods: Once Scorned, Now Respectable](http://cm.bell-labs.com/cm/cs/doc/96/4-02.ps.gz)
), they are used
when either the computation of the derivative is impossible (noisy
functions, unpredictable discontinuities) or difficult (complexity,
computation cost). In the first cases, rather than an optimum, a
<em>not too bad</em> point is desired. In the latter cases, an
optimum is desired but cannot be reasonably found. In all cases
direct search methods can be useful.

Simplex-based direct search methods are based on comparison of
the cost function values at the vertices of a simplex (which is a
set of n+1 points in dimension n) that is updated by the algorithms
steps.

The instances can be built either in single-start or in
multi-start mode. Multi-start is a traditional way to try to avoid
being trapped in a local minimum and miss the global minimum of a
function. It can also be used to verify the convergence of an
algorithm. In multi-start mode, the `minimizes`method
returns the best minimum found after all starts, and the `etMinima`
method can be used to retrieve all minima from all starts (including the one
already provided by the `minimizes` method).

The `direct` package provides four solvers:

* the univariate package handles univariate scalar functions,
* the linear package handles multivariate vector linear functions with linear constraints,
* the direct package handles multivariate scalar functions using direct search methods (i.e. not using derivatives),
* the general package handles multivariate scalar or vector functions using derivatives.
* the fitting package handles curve fitting by univariate real functions


The first two simplex-based methods do not handle simple bounds constraints by themselves.
However there are two adapters(<a
href="../apidocs/org.hipparchus/optimization/direct/MultivariateFunctionMappingAdapter.html">
MultivariateFunctionMappingAdapter</a> and <a
href="../apidocs/org.hipparchus/optimization/direct/MultivariateFunctionPenaltyAdapter.html">
MultivariateFunctionPenaltyAdapter</a>) that can be used to wrap the user function in
such a way the wrapped function is unbounded and can be used with these optimizers, despite
the fact the underlying function is still bounded and will be called only with feasible
points that fulfill the constraints. Note however that using these adapters are only a
poor man solutions to simple bounds optimization constraints. Better solutions are to use an
optimizer that directly supports simple bounds. Some caveats of the mapping adapter
solution are that

* the univariate package handles univariate scalar functions,
* the linear package handles multivariate vector linear functions with linear constraints,
* the direct package handles multivariate scalar functions using direct search methods (i.e. not using derivatives),
* the general package handles multivariate scalar or vector functions using derivatives.
* the fitting package handles curve fitting by univariate real functions

One caveat of penalty adapter is that if start point or start simplex is outside of the allowed
range, only the penalty function is used, and the optimizer may converge without ever entering
the allowed range.

The last methods do handle simple bounds constraints directly, so the adapters are not needed
with them.


## General Case
The general package deals with non-linear vectorial optimization problems when
the partial derivatives of the objective function are available.

One important class of estimation problems is weighted least
squares problems. They basically consist in finding the values
for some parameters p<sub>k</sub> such that a cost function
J = sum(w<sub>i</sub>(mes<sub>i</sub> - mod<sub>i</sub>)<sup>2</sup>) is
minimized. The various (target<sub>i</sub> - model<sub>i</sub>(p<sub>k</sub>))
terms are called residuals. They represent the deviation between a set of
target values target<sub>i</sub> and theoretical values computed from
models model<sub>i</sub> depending on free parameters p<sub>k</sub>.
The w<sub>i</sub> factors are weights. One classical use case is when the
target values are experimental observations or measurements.

Solving a least-squares problem is finding the free parameters p<sub>k</sub>
of the theoretical models such that they are close to the target values, i.e.
when the residual are small.

Two optimizers are available in the general package, both devoted to least-squares
problems. The first one is based on the <a
href="../apidocs/org.hipparchus/optimization/general/GaussNewtonOptimizer.html">
Gauss-Newton</a> method. The second one is the <a
href="../apidocs/org.hipparchus/optimization/general/LevenbergMarquardtOptimizer.html">
Levenberg-Marquardt</a> method.

In order to solve a vectorial optimization problem, the user must provide it as
an object implementing the <a
href="../apidocs/org.hipparchus/analysis/DifferentiableMultivariateVectorFunction.html">
DifferentiableMultivariateVectorFunction</a> interface. The object will be provided to
the `estimate` method of the optimizer, along with the target and weight arrays,
thus allowing the optimizer to compute the residuals at will. The last parameter to the
`estimate` method is the point from which the optimizer will start its
search for the optimal point.


## Quadratic Problem Example


We are looking to find the best parameters [a, b, c] for the quadratic function
\\(f(x) = a x^2 + b x + c\\).
The data set below was generated using [a = 8, b = 10, c = 16].
A random number between zero and one was added to each y value calculated.

| X | Y |
| --- |
| 1 | 34.234064369 |
| 2 | 68.2681162306108 |
| 3 | 118.615899084602 |
| 4 | 184.138197238557 |
| 5 | 266.599877916276 |
| 6 | 364.147735251579 |
| 7 | 478.019226091914 |
| 8 | 608.140949270688 |
| 9 | 754.598868667148 |
| 10 | 916.128818085883 |

First we need to implement the interface [DifferentiableMultivariateVectorFunction](../apidocs/org/hipparchus/analysis/DifferentiableMultivariateVectorFunction.html).
This requires the implementation of the method signatures:

* the univariate package handles univariate scalar functions,
* the linear package handles multivariate vector linear functions with linear constraints,
* the direct package handles multivariate scalar functions using direct search methods (i.e. not using derivatives),
* the general package handles multivariate scalar or vector functions using derivatives.
* the fitting package handles curve fitting by univariate real functions


We'll tackle the implementation of the `MultivariateMatrixFunction jacobian()` method first.  You may wish to familiarize yourself with what a [ Jacobian Matrix](http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant) is.
In this case the Jacobian is the partial derivative of the function with respect
to the parameters a, b and c.  These derivatives are computed as follows:

* the univariate package handles univariate scalar functions,
* the linear package handles multivariate vector linear functions with linear constraints,
* the direct package handles multivariate scalar functions using direct search methods (i.e. not using derivatives),
* the general package handles multivariate scalar or vector functions using derivatives.
* the fitting package handles curve fitting by univariate real functions

For a quadratic which has three variables the Jacobian Matrix will have three columns, one for each variable, and the number
of rows will equal the number of rows in our data set, which in this case is ten.  So for example for <tt>[a = 1, b = 1, c = 1]</tt>, the Jacobian Matrix is (excluding the first column which shows the value of x):


| x | \\(\\frac{\\partial{(ax^2 + bx + c)}}{\\partial{a}}\\) | \\(\\frac{\\partial{(ax^2 + bx + c)}}{\\partial{b}}\\) | \\(\\frac{\\partial{(ax^2 + bx + c)}}{\\partial{c}}\\) |
| --- |
| 1 | 1 | 1 | 1 |
| 2 | 4 | 2 | 1 |
| 3 | 9 | 3 | 1 |
| 4 | 16 | 4 | 1 |
| 5 | 25 | 5 | 1 |
| 6 | 36 | 6 | 1 |
| 7 | 49 | 7 | 1 |
| 8 | 64 | 8 | 1 |
| 9 | 81 | 9 | 1 |
| 10 | 100 | 10 | 1 |

The implementation of the `MultivariateMatrixFunction jacobian()` for this problem looks like this (The `x`
parameter is an ArrayList containing the independent values of the data set):



     private double[][] jacobian(double[] variables) {
         double[][] jacobian = new double[x.size()][3];
         for (int i = 0; i &lt; jacobian.length; ++i) {
             jacobian[i][0] = x.get(i) * x.get(i);
             jacobian[i][1] = x.get(i);
             jacobian[i][2] = 1.0;
         }
         return jacobian;
     }
    
     public MultivariateMatrixFunction jacobian() {
         return new MultivariateMatrixFunction() {
             private static final long serialVersionUID = -8673650298627399464L;
             public double[][] value(double[] point) {
                 return jacobian(point);
             }
         };
     }

Note that if for some reason the derivative of the objective function with respect
to its variables is difficult to obtain,
[Numerical differentiation](http://en.wikipedia.org/wiki/Numerical_differentiation) can be used.



The implementation of the `double[] value(double[] point)` method, which returns
a `double` array containing the
values the objective function returns per given independent value
and the current set of variables or parameters,
can be seen below:

    public double[] value(double[] variables) {
        double[] values = new double[x.size()];
        for (int i = 0; i &lt; values.length; ++i) {
            values[i] = (variables[0] * x.get(i) + variables[1]) * x.get(i) + variables[2];
        }
        return values;
    }

Below is the the class containing all the implementation details
(Taken from the Hipparchus <b>org.hipparchus.optimization.general.LevenbergMarquardtOptimizerTest</b>):



    private static class QuadraticProblem
        implements DifferentiableMultivariateVectorFunction, Serializable {
    
        private static final long serialVersionUID = 7072187082052755854L;
        private List&lt;Double&gt; x;
        private List&lt;Double&gt; y;
    
        public QuadraticProblem() {
            x = new ArrayList&lt;Double&gt;();
            y = new ArrayList&lt;Double&gt;();
        }
    
        public void addPoint(double x, double y) {
            this.x.add(x);
            this.y.add(y);
        }
    
        public double[] calculateTarget() {
            double[] target = new double[y.size()];
            for (int i = 0; i &lt; y.size(); i++) {
                target[i] = y.get(i).doubleValue();
            }
            return target;
        }
    
        private double[][] jacobian(double[] variables) {
            double[][] jacobian = new double[x.size()][3];
            for (int i = 0; i &lt; jacobian.length; ++i) {
                jacobian[i][0] = x.get(i) * x.get(i);
                jacobian[i][1] = x.get(i);
                jacobian[i][2] = 1.0;
            }
            return jacobian;
        }
    
        public double[] value(double[] variables) {
            double[] values = new double[x.size()];
            for (int i = 0; i &lt; values.length; ++i) {
                values[i] = (variables[0] * x.get(i) + variables[1]) * x.get(i) + variables[2];
            }
            return values;
        }
    
        public MultivariateMatrixFunction jacobian() {
            return new MultivariateMatrixFunction() {
                private static final long serialVersionUID = -8673650298627399464L;
                public double[][] value(double[] point) {
                    return jacobian(point);
                }
            };
        }
    }

The below code shows how to go about using the above class
and a LevenbergMarquardtOptimizer instance to produce an
optimal set of quadratic curve fitting parameters:

    QuadraticProblem problem = new QuadraticProblem();
    
    problem.addPoint(1, 34.234064369);
    problem.addPoint(2, 68.2681162306);
    problem.addPoint(3, 118.6158990846);
    problem.addPoint(4, 184.1381972386);
    problem.addPoint(5, 266.5998779163);
    problem.addPoint(6, 364.1477352516);
    problem.addPoint(7, 478.0192260919);
    problem.addPoint(8, 608.1409492707);
    problem.addPoint(9, 754.5988686671);
    problem.addPoint(10, 916.1288180859);
    
    LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
    
    final double[] weights = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    
    final double[] initialSolution = {1, 1, 1};
    
    PointVectorValuePair optimum = optimizer.optimize(100,
                                                      problem,
                                                      problem.calculateTarget(),
                                                      weights,
                                                      initialSolution);
   
    final double[] optimalValues = optimum.getPoint();
   
    System.out.println(&quot;A: &quot; + optimalValues[0]);
    System.out.println(&quot;B: &quot; + optimalValues[1]);
    System.out.println(&quot;C: &quot; + optimalValues[2]);


If you run the above sample you will see the following printed by the console:

    A: 7.998832172372726
    B: 10.001841530162448
    C: 16.324008168386605


In addition to least squares solving, the <a
href="../apidocs/org.hipparchus/optimization/general/NonLinearConjugateGradientOptimizer.html">
NonLinearConjugateGradientOptimizer</a> class provides a non-linear conjugate gradient algorithm
to optimize <a
href="../apidocs/org.hipparchus/analysis/DifferentiableMultivariateFunction.html">
DifferentiableMultivariateFunction</a>. Both the Fletcher-Reeves and the Polak-Ribi&#232;re
search direction update methods are supported. It is also possible to set up a preconditioner
or to change the line-search algorithm of the inner loop if desired (the default one is a Brent
solver).

The [PowellOptimizer](../apidocs/org/hipparchus/optimization/direct/PowellOptimizer.html)
provides an optimization method for non-differentiable functions.


