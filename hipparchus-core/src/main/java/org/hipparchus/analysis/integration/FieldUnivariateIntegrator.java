/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.analysis.integration;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;

/**
 * Interface for univariate real integration algorithms.
 * @param <T> Type of the field elements.
 * @since 2.0
 */
public interface FieldUnivariateIntegrator<T extends CalculusFieldElement<T>> {

    /**
     * Get the relative accuracy.
     *
     * @return the accuracy
     */
    double getRelativeAccuracy();

    /**
     * Get the absolute accuracy.
     *
     * @return the accuracy
     */
    double getAbsoluteAccuracy();

    /**
     * Get the min limit for the number of iterations.
     *
     * @return the actual min limit
     */
    int getMinimalIterationCount();

    /**
     * Get the upper limit for the number of iterations.
     *
     * @return the actual upper limit
     */
    int getMaximalIterationCount();

    /**
     * Integrate the function in the given interval.
     *
     * @param maxEval Maximum number of evaluations.
     * @param f the integrand function
     * @param min the lower bound for the interval
     * @param max the upper bound for the interval
     * @return the value of integral
     * @throws MathIllegalStateException if the maximum number of function
     * evaluations is exceeded
     * @throws MathIllegalStateException if the maximum iteration count is exceeded
     * or the integrator detects convergence problems otherwise
     * @throws MathIllegalArgumentException if {@code min > max} or the endpoints do not
     * satisfy the requirements specified by the integrator
     * @throws NullArgumentException if {@code f} is {@code null}.
     */
    T integrate(int maxEval, CalculusFieldUnivariateFunction<T> f, T min, T max)
        throws MathIllegalArgumentException, MathIllegalStateException, NullArgumentException;

    /**
     * Get the number of function evaluations of the last run of the integrator.
     *
     * @return number of function evaluations
     */
    int getEvaluations();

    /**
     * Get the number of iterations of the last run of the integrator.
     *
     * @return number of iterations
     */
    int getIterations();

}
