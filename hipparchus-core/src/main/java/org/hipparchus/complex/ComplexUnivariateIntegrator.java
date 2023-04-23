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
package org.hipparchus.complex;

import java.util.function.DoubleFunction;

import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.analysis.integration.UnivariateIntegrator;

/**
 * Wrapper to perform univariate complex integration using an underlying real integration algorithms.
 * @since 2.0
 */
public class ComplexUnivariateIntegrator  {

    /** Underlying real integrator. */
    private UnivariateIntegrator integrator;

    /** Crate a complex integrator from a real integrator.
     * @param integrator underlying real integrator to use
     */
    public ComplexUnivariateIntegrator(final UnivariateIntegrator integrator) {
        this.integrator = integrator;
    }

    /**
     * Integrate a function along a straight path between points.
     *
     * @param maxEval maximum number of evaluations (real and imaginary
     * parts are evaluated separately, so up to twice this number may be used)
     * @param f the integrand function
     * @param start start point of the integration path
     * @param end end point of the integration path
     * @return the value of integral along the straight path
     */
    public Complex integrate(final int maxEval, final CalculusFieldUnivariateFunction<Complex> f,
                             final Complex start, final Complex end) {

        // linear mapping from real interval [0; 1] to function value along complex straight path from start to end
        final Complex                 rate   = end.subtract(start);
        final DoubleFunction<Complex> mapped = t -> f.value(start.add(rate.multiply(t)));

        // integrate real and imaginary parts separately
        final double real      = integrator.integrate(maxEval, t -> mapped.apply(t).getRealPart(),      0.0, 1.0);
        final double imaginary = integrator.integrate(maxEval, t -> mapped.apply(t).getImaginaryPart(), 0.0, 1.0);

        // combine integrals
        return new Complex(real, imaginary).multiply(rate);

    }

    /**
     * Integrate a function along a polyline path between any number of points.
     *
     * @param maxEval maximum number of evaluations (real and imaginary
     * parts are evaluated separately and each path segments are also evaluated
     * separately, so up to 2n times this number may be used for n segments)
     * @param f the integrand function
     * @param start start point of the integration path
     * @param path successive points defining the path vertices
     * @return the value of integral along the polyline path
     */
    public Complex integrate(final int maxEval, final CalculusFieldUnivariateFunction<Complex> f,
                             final Complex start, final Complex...path) {
        Complex sum      = Complex.ZERO;
        Complex previous = start;
        for (final Complex current : path) {
            sum = sum.add(integrate(maxEval, f, previous, current));
            previous = current;
        }
        return sum;
    }

}
