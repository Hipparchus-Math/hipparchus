/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
/**
 *
 * <p>
 *   This package holds the main interfaces and basic building block classes
 *   dealing with differentiation.
 *   The core class is {@link org.hipparchus.analysis.differentiation.DerivativeStructure
 *   DerivativeStructure} which holds the value and the differentials of a function. This class
 *   handles some arbitrary number of free parameters and arbitrary differentiation order. It is used
 *   both as the input and the output type for the {@link
 *   org.hipparchus.analysis.differentiation.UnivariateDifferentiableFunction
 *   UnivariateDifferentiableFunction} interface. Any differentiable function should implement this
 *   interface.
 * </p>
 * <p>
 *   The {@link org.hipparchus.analysis.differentiation.UnivariateDerivative1 UnivariateDerivative1},
 *   {@link org.hipparchus.analysis.differentiation.UnivariateDerivative2 UnivariateDerivative2} and
 *   {@link org.hipparchus.analysis.differentiation.Gradient Gradient} classes are more restricted
 *   implementation of classes holding the value and the differentials of a function. These classes
 *   handle only either one free parameter (i.e. univariate functions) with derivation orders 1 or 2,
 *   or several free parameters with derivation order 1. As they arefar less general than {@link
 *   org.hipparchus.analysis.differentiation.DerivativeStructure DerivativeStructure}, they have less
 *   overhead and are more efficient in their respective domains.
 * </p>
 * <p>
 *   The {@link org.hipparchus.analysis.differentiation.UnivariateFunctionDifferentiator
 *   UnivariateFunctionDifferentiator} interface defines a way to differentiate a simple {@link
 *   org.hipparchus.analysis.UnivariateFunction UnivariateFunction} and get a {@link
 *   org.hipparchus.analysis.differentiation.UnivariateDifferentiableFunction
 *   UnivariateDifferentiableFunction}.
 * </p>
 * <p>
 *   Similar interfaces also exist for multivariate functions and for vector or matrix valued functions.
 * </p>
 *
 */
package org.hipparchus.analysis.differentiation;
