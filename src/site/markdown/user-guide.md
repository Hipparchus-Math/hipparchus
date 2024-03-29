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
# User Guide

* [About the User Guide](overview.html#About_the_User_Guide)
* [What's in Hipparchus](overview.html#Whats_in_Hipparchus)
* [How Hipparchus is organized](overview.html#How_Hipparchus_is_organized)
* [How interface contracts are specified in Hipparchus javadoc](overview.html#How_interface_contracts_are_specified_in_Hipparchus_javadoc)
* [Dependencies](overview.html#Dependencies)
* [Statistics](hipparchus-stat/index.html#Overview)
 * [Descriptive statistics](hipparchus-stat/index.html#Descriptive_statistics)
 * [Frequency distributions](hipparchus-stat/index.html#Frequency_distributions)
 * [Simple regression](hipparchus-stat/index.html#Simple_regression)
 * [Multiple Regression](hipparchus-stat/index.html#Multiple_linear_regression)
 * [Rank transformations](hipparchus-stat/index.html#Rank_transformations)
 * [Covariance and correlation](hipparchus-stat/index.html#Covariance_and_correlation)
 * [Statistical Tests](hipparchus-stat/index.html#Statistical_tests)
* [Random Data Generation](hipparchus-core/random.html#Overview)
 * [Random numbers](hipparchus-core/random.html#Random_numbers)
 * [Random Vectors](hipparchus-core/random.html#Random_Vectors)
 * [Random Strings](hipparchus-core/random.html#Random_Strings)
 * [Random permutations, combinations, sampling](hipparchus-core/random.html#Random_permutations_combinations_sampling)
 * [Generating data 'like' an input file](hipparchus-core/random.html#Generating_data_like_an_input_file)
 * [PRNG Pluggability](hipparchus-core/random.html#PRNG_Pluggability)
* [Linear Algebra](hipparchus-core/linear.html#Overview)
 * [Real matrices](hipparchus-core/linear.html#Real_matrices)
 * [Real vectors](hipparchus-core/linear.html#Real_vectors)
 * [Solving linear systems](hipparchus-core/linear.html#Solving_linear_systems)
 * [Eigenvalues/eigenvectors and singular values/singular vectors](hipparchus-core/linear.html#Eigenvalueseigenvectors_and_singular_valuessingular_vectors)
 * [Non-real fields (complex, fractions ...)](hipparchus-core/linear.html#Non-real_fields_complex_fractions_...)
* [Numerical Analysis](hipparchus-core/analysis.html#Overview)
 * [Error handling](hipparchus-core/analysis.html#Error-handling)
 * [Root-finding](hipparchus-core/analysis.html#Root-finding)
 * [Interpolation](hipparchus-core/analysis.html#Interpolation)
 * [Integration](hipparchus-core/analysis.html#Integration)
 * [Polynomials](hipparchus-core/analysis.html#Polynomials)
 * [Differentiation](hipparchus-core/analysis.html#Differentiation)
* [Special functions](hipparchus-core/special.html#Overview)
 * [Erf functions](hipparchus-core/special.html#Erf_functions)
 * [Gamma functions](hipparchus-core/special.html#Gamma_functions)
 * [Beta funtions](hipparchus-core/special.html#Beta_funtions)
* [Utilities](hipparchus-core/utilities.html#Overview)
 * [Double array utilities](hipparchus-core/utilities.html#Double_array_utilities)
 * [int/double hash map](hipparchus-core/utilities.html#intdouble_hash_map)
 * [Continued Fractions](hipparchus-core/utilities.html#Continued_Fractions)
 * [Binomial coefficients, factorials, Stirling numbers and other common math functions](hipparchus-core/utilities.html#binomial_coefficients_factorials_Stirling_numbers_and_other_common_math_functions)
 * [Fast mathematical functions](hipparchus-core/utilities.html#fast_math)
 * [Miscellaneous](hipparchus-core/utilities.html#miscellaneous)
* [Complex Analysis](hipparchus-core/complex.html#Overview)
 * [Complex Numbers](hipparchus-core/complex.html#Complex_Numbers)
 * [Complex Transcendental Functions](hipparchus-core/complex.html#Complex_Transcendental_Functions)
 * [Complex Formatting and Parsing](hipparchus-core/complex.html#Complex_Formatting_and_Parsing)
* [Probability Distributions](hipparchus-core/distribution.html#Overview)
 * [Distribution Framework](hipparchus-core/distribution.html#Distribution_Framework)
 * [User Defined Distributions](hipparchus-core/distribution.html#User_Defined_Distributions)
* [Exact Fractions](hipparchus-core/fraction.html#Overview)
 * [Fraction Numbers](hipparchus-core/fraction.html#Fraction_Numbers)
 * [Fraction Formatting and Parsing](hipparchus-core/fraction.html#Fraction_Formatting_and_Parsing)
* [Transform methods](hipparchus-fft/index.html)
* [Geometry](hipparchus-geometry/index.html#Overview)
 * [Euclidean spaces](hipparchus-geometry/index.html#Euclidean_spaces)
 * [n-Sphere](hipparchus-geometry/index.html#n-Sphere)
 * [Binary Space Partitioning](hipparchus-geometry/index.html#Binary_Space_Partitioning)
 * [Regions](hipparchus-geometry/index.html#Regions)
* [Optimization](hipparchus-optim/index.html#Overview)
 * [Univariate Functions](hipparchus-optim/index.html#Univariate_Functions)
 * [Linear Programming](hipparchus-optim/index.html#Linear_Programming)
 * [Direct Methods](hipparchus-optim/index.html#Direct_Methods)
 * [General Case](hipparchus-optim/index.html#General_Case)
* [Fitting](hipparchus-fitting/fitting.html#Overview)
 * [Implemented Functions](hipparchus-fitting/fitting.html#Implemented_Functions)
 * [General Case](hipparchus-fitting/fitting.html#General_Case)
 * [Least Squares](hipparchus-optim/leastsquares.html#Overview)
 * [LeastSquaresBuilder and LeastSquaresFactory](hipparchus-optim/leastsquares.html#LeastSquaresBuilder_and_LeastSquaresFactory)
 * [Model Function](hipparchus-optim/leastsquares.html#Model_Function)
 * [Parameters Validation](hipparchus-optim/leastsquares.html#Parameters_Validation)
 * [Tuning](hipparchus-optim/leastsquares.html#Tuning)
 * [Optimization Engine](hipparchus-optim/leastsquares.html#Optimization_Engine)
 * [Solving](hipparchus-optim/leastsquares.html#Solving)
 * [Example](hipparchus-optim/leastsquares.html#Example)
* [Ordinary Differential Equations](hipparchus-ode/index.html#Overview)
 * [Continuous Output](hipparchus-ode/index.html#Continuous_Output)
 * [Discrete Events Handling](hipparchus-ode/index.html#Discrete_Events_Handling)
* [Integration](hipparchus-ode/index.html#Available_Integrators)
* [Derivatives](hipparchus-ode/index.html#Derivatives)
* [Clustering](hipparchus-clustering/index.html#overview)
 * [Clustering algorithms and distance measures](hipparchus-clustering/index.html#clustering)
 * [Implementation](hipparchus-clustering/index.html#implementation)
* [Exceptions](hipparchus-core/exceptions.html#Overview)
 * [Unchecked Exceptions](hipparchus-core/exceptions.html#Unchecked_Exceptions)
 * [Hierarchies](hipparchus-core/exceptions.html#Hierarchies)
 * [Features](hipparchus-core/exceptions.html#Features)
