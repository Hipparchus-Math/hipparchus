<!--
 Licensed to the Hipparchus project under one or more
 contributor license agreements.  See the NOTICE file distributed with
 this work for additional information regarding copyright ownership.
 The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 (the "License"); you may not use this file except in compliance with
 the License.  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
-->
# Utilities

## Overview

The [org.hipparchus.util](../apidocs/org/hipparchus/util/package-summary.html)
package collects a group of array utilities, value transformers, and numerical
routines used by implementation classes in Hipparchus.


## Double array utilities

To maintain statistics based on a "rolling" window of values, a resizable
array implementation was developed and is provided for reuse in the
`util` package.  The core functionality provided is described in
the documentation for the interface,
[DoubleArray](../apidocs/org/hipparchus/util/DoubleArray.html). 
This interface adds one method, `addElementRolling(double)` to basic list accessors.
The `addElementRolling` method adds an element (the actual parameter) to the end
of the list and removes the first element in the list.

The [ResizableDoubleArray](../apidocs/org/hipparchus/util/ResizableDoubleArray.html)
class provides a configurable, array-backed
implementation of the `DoubleArray` interface.
When `addElementRolling` is invoked, the underlying
array is expanded if necessary, the new element is added to the end of the
array and the "usable window" of the array is moved forward, so that
the first element is effectively discarded, what was the second becomes the
first, and so on.  To efficiently manage storage, two maintenance
operations need to be periodically performed -- orphaned elements at the
beginning of the array need to be reclaimed and space for new elements at
the end needs to be created.  Both of these operations are handled
automatically, with frequency / effect driven by the configuration
properties `expansionMode`, `expansionFactor` and
`contractionCriteria.`  See
[ResizableDoubleArray](../apidocs/org/hipparchus/util/ResizableDoubleArray.html)
for details.


## Primitive int/double hash map

The [OpenIntToDoubleHashMap](../apidocs/org/hipparchus/util/OpenIntToDoubleHashMap.html)
class provides a specialized hash map implementation for int/double. This implementation
has a much smaller memory overhead than standard `java.util.HashMap` class.
It uses open addressing and primitive arrays, which greatly reduces the number of
intermediate objects and improve data locality.


## Continued Fractions

The [ContinuedFraction](../apidocs/org/hipparchus/util/ContinuedFraction.html)
class provides a generic way to create and evaluate continued fractions.
The easiest way to create a continued fraction is to subclass `ContinuedFraction`
and override the `getA` and `getB` methods which return
the continued fraction terms.  The precise definition of these terms is
explained in [Continued Fraction, equation (1)](http://mathworld.wolfram.com/ContinuedFraction.html)
from MathWorld.

As an example, the constant Pi can be computed using a
[continued fraction](http://functions.wolfram.com/Constants/Pi/10/0002/).
The following anonymous class provides the implementation:

    ContinuedFraction c = new ContinuedFraction() {
        public double getA(int n, double x) {
            switch(n) {
                case 0: return 3.0;
                default: return 6.0;
            }
        }
        
        public double getB(int n, double x) {
            double y = (2.0 * n) - 1.0;
            return y * y;
        }
    }

Then, to evaluate Pi, simply call any of the `evaluate` methods
(Note, the point of evaluation in this example is meaningless since Pi is a constant).

For a more practical use of continued fractions, consider the
[exponential function](http://functions.wolfram.com/ElementaryFunctions/Exp/10/).
The following anonymous class provides its implementation:

    ContinuedFraction c = new ContinuedFraction() {
        public double getA(int n, double x) {
            if (n % 2 == 0) {
                switch(n) {
                    case 0: return 1.0;
                    default: return 2.0;
                }
            } else {
                return n;
            }
        }
        
        public double getB(int n, double x) {
            if (n % 2 == 0) {
                return -x;
            } else {
                return x;
            }
        }
    }

Then, to evaluate *e*<sup>x</sup> for any value x, simply call any of the
`evaluate` methods.


## Fast mathematical functions

Hipparchus provides a faster, more accurate, portable alternative
to the regular `Math` and `StrictMath`classes for large scale computation.

FastMath is a drop-in replacement for both Math and StrictMath. This
means that for any method in Math (say `Math.sin(x)` or
`Math.cbrt(y)`), user can directly change the class and use the
methods as is (using `FastMath.sin(x)` or `FastMath.cbrt(y)`
in the previous example).

FastMath speed is achieved by relying heavily on optimizing compilers to
native code present in many JVM todays and use of large tables. Precomputed
literal arrays are provided in this class to speed up load time. These
precomputed tables are used in the default configuration, to improve speed
even at first use of the class. If users prefer to compute the tables
automatically at load time, they can change a compile-time constant. This will
increase class load time at first use, but this overhead will occur only once
per run, regardless of the number of subsequent calls to computation methods.
Note that FastMath is extensively used inside Hipparchus, so by
calling some algorithms, the one-shot overhead when the constant is set to
false will occur regardless of the end-user calling FastMath methods directly
or not. Performance figures for a specific JVM and hardware can be evaluated by
running the FastMathTestPerformance tests in the test directory of the source
distribution.

FastMath accuracy should be mostly independent of the JVM as it relies only
on IEEE-754 basic operations and on embedded tables. Almost all operations
are accurate to about 0.5 ulp throughout the domain range. This statement, of
course is only a rough global observed behavior, it is <em>not</em> a guarantee
for <em>every</em> double numbers input (see William Kahan's <a
href="http://en.wikipedia.org/wiki/Rounding#The_table-maker.27s_dilemma">Table
Maker's Dilemma</a>).

FastMath additionally implements the following methods not found in Math/StrictMath:

* asinh(double)
* acosh(double)
* atanh(double)
* pow(double,int)

The following methods are found in Math/StrictMath since 1.6 only, they are provided by FastMath even in 1.5 Java virtual machines

* copySign(double, double)
* getExponent(double)
* nextAfter(double,double)
* nextUp(double)
* scalb(double, int)
* copySign(float, float)
* getExponent(float)
* nextAfter(float,double)
* nextUp(float)
* scalb(float, int)


## Miscellaneous

A collection of reusable math functions is provided in the
[ArithmeticUtils](../apidocs/org/hipparchus/util/ArithmeticUtils.html)

* Binomial coefficients -- "n choose k" available as an (exact) long value, `binomialCoefficient(int, int)` for small n, k; as a double, `binomialCoefficientDouble(int, int)` for larger values; and in a "super-sized" version, `binomialCoefficientLog(int, int)` that returns the natural logarithm of the value.
* Stirling numbers of the second kind -- S(n,k) as an exact long value `stirlingS2(int, int)` for small n, k.
* Factorials -- like binomial coefficients, these are available as exact long values, `factorial(int)`; doubles, `factorialDouble(int)`; or logs, `factorialLog(int)`.
* Least common multiple and greatest common denominator functions.

The [MultidimensionalCounter](../apidocs/org/hipparchus/util/MultidimensionalCounter.html)
is a utility class that converts a set of indices (identifying points in a multidimensional
space) to a single index (e.g. identifying a location in a one-dimensional array.

