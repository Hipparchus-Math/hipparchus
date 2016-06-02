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
# Overview

## About The User Guide
This guide is intended to help programmers quickly find what they need to develop
solutions using Hipparchus.  It also provides a supplement to the javadoc API documentation,
providing a little more explanation of the mathematical objects and functions included
in the package.



## What's in hipparchus
Hipparchus is made up of a small set of math/stat utilities addressing
programming problems like the ones in the list below.  This list is not exhaustive,
it's just meant to give a feel for the kinds of things that Hipparchus provides.

* Computing means, variances and other summary statistics for a list of numbers
* Fitting a line to a set of data points using linear regression
* Fitting a curve to a set of data points
* Finding a smooth curve that passes through a collection of points (interpolation)
* Fitting a parametric model to a set of measurements using least-squares methods
* Solving equations involving real-valued functions (i.e. root-finding)
* Solving systems of linear equations
* Solving Ordinary Differential Equations
* Minimizing multi-dimensional functions
* Generating random numbers with more restrictions (e.g distribution, range) than what is possible using the JDK
* Generating random samples and/or datasets that are "like" the data in an input file
* Performing statistical significance tests
* Miscellaneous mathematical functions such as factorials, binomial coefficients and "special functions" (e.g. gamma, beta functions)

We are actively seeking ideas for additional components that fit into the
[Hipparchus vision](index.html#summary) of a set of lightweight,
self-contained math/stat components useful for solving common programming problems.
Suggestions for new components or enhancements to existing functionality are always welcome!
All feedback/suggestions for improvement should be sent to the
[developers mailing list](mail-lists.html).



## How Hipparchus is organized
Hipparchus is divided into sixteen subpackages, based on functionality provided.

* Computing means, variances and other summary statistics for a list of numbers
* Fitting a line to a set of data points using linear regression
* Fitting a curve to a set of data points
* Finding a smooth curve that passes through a collection of points (interpolation)
* Fitting a parametric model to a set of measurements using least-squares methods
* Solving equations involving real-valued functions (i.e. root-finding)
* Solving systems of linear equations
* Solving Ordinary Differential Equations
* Minimizing multi-dimensional functions
* Generating random numbers with more restrictions (e.g distribution, range) than what is possible using the JDK
* Generating random samples and/or datasets that are "like" the data in an input file
* Performing statistical significance tests
* Miscellaneous mathematical functions such as factorials, binomial coefficients and "special functions" (e.g. gamma, beta functions)

Package javadocs are [here](apidocs/index.html)



## How interface contracts are specified in hipparchus javadoc
You should always read the javadoc class and method comments carefully when using
Hipparchus components in your programs.  The javadoc provides references to the algorithms
that are used, usage notes about limitations, performance, etc. as well as interface contracts.
Interface contracts are specified in terms of preconditions (what has to be true in order
for the method to return valid results), special values returned (e.g. Double.NaN)
or exceptions that may be thrown if the preconditions are not met, and definitions for returned
values/objects or state changes.

When the actual parameters provided to a method or the internal state of an object
make a computation meaningless, a
[MathIllegalArgumentException](apidocs/org/hipparchus/exception/MathIllegalArgumentException.html)
or
[MathIllegalStateException](apidocs/org/hipparchus/exception/MathIllegalStateException.html)
may be thrown. Exact conditions under which runtime
exceptions (and any other exceptions) are thrown are specified in the javadoc method
comments.
In some cases, to be consistent with the [IEEE 754 standard](http://grouper.ieee.org/groups/754/)
for floating point arithmetic and with java.lang.Math, Hipparchus
methods return `Double.NaN` values. Conditions under which `Double.NaN`
or other special values are returned are fully specified in the javadoc method comments.

The policy for dealing with null references is as
follows: When an argument is unexpectedly null, a
[NullArgumentException](apidocs/org/hipparchus/exception/NullArgumentException.html)
is raised to signal the illegal argument. Note that this
class does not inherit from the standard `NullPointerException` but is a subclass
of `MathIllegalArgumentException`.



## Dependencies
Hipparchus requires JDK 1.8+ and has no runtime dependencies.



## License
Hipparchus is distributed under the terms of the [Apache License, Version 2.0](license.html)

