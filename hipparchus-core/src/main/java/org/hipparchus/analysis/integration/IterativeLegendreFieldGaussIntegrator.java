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
package org.hipparchus.analysis.integration;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.analysis.integration.gauss.FieldGaussIntegrator;
import org.hipparchus.analysis.integration.gauss.FieldGaussIntegratorFactory;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;

/**
 * This algorithm divides the integration interval into equally-sized
 * sub-interval and on each of them performs a
 * <a href="http://mathworld.wolfram.com/Legendre-GaussQuadrature.html">
 * Legendre-Gauss</a> quadrature.
 * Because of its <em>non-adaptive</em> nature, this algorithm can
 * converge to a wrong value for the integral (for example, if the
 * function is significantly different from zero toward the ends of the
 * integration interval).
 * In particular, a change of variables aimed at estimating integrals
 * over infinite intervals as proposed
 * <a href="http://en.wikipedia.org/w/index.php?title=Numerical_integration#Integrals_over_infinite_intervals">
 *  here</a> should be avoided when using this class.
 *
 * @param <T> Type of the field elements.
 * @since 2.0
 */
public class IterativeLegendreFieldGaussIntegrator<T extends CalculusFieldElement<T>>
    extends BaseAbstractFieldUnivariateIntegrator<T> {

    /** Factory that computes the points and weights. */
    private final FieldGaussIntegratorFactory<T> factory;

    /** Number of integration points (per interval). */
    private final int numberOfPoints;

    /**
     * Builds an integrator with given accuracies and iterations counts.
     *
     * @param field field to which function argument and value belong
     * @param n Number of integration points.
     * @param relativeAccuracy Relative accuracy of the result.
     * @param absoluteAccuracy Absolute accuracy of the result.
     * @param minimalIterationCount Minimum number of iterations.
     * @param maximalIterationCount Maximum number of iterations.
     * @throws MathIllegalArgumentException if minimal number of iterations
     * or number of points are not strictly positive.
     * @throws MathIllegalArgumentException if maximal number of iterations
     * is smaller than or equal to the minimal number of iterations.
     */
    public IterativeLegendreFieldGaussIntegrator(final Field<T> field, final int n,
                                                 final double relativeAccuracy,
                                                 final double absoluteAccuracy,
                                                 final int minimalIterationCount,
                                                 final int maximalIterationCount)
        throws MathIllegalArgumentException {
        super(field, relativeAccuracy, absoluteAccuracy, minimalIterationCount, maximalIterationCount);
        if (n <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_OF_POINTS, n);
        }
        factory = new FieldGaussIntegratorFactory<>(field);
        numberOfPoints = n;
    }

    /**
     * Builds an integrator with given accuracies.
     *
     * @param field field to which function argument and value belong
     * @param n Number of integration points.
     * @param relativeAccuracy Relative accuracy of the result.
     * @param absoluteAccuracy Absolute accuracy of the result.
     * @throws MathIllegalArgumentException if {@code n < 1}.
     */
    public IterativeLegendreFieldGaussIntegrator(final Field<T> field, final int n,
                                                 final double relativeAccuracy,
                                                 final double absoluteAccuracy)
        throws MathIllegalArgumentException {
        this(field, n, relativeAccuracy, absoluteAccuracy,
             DEFAULT_MIN_ITERATIONS_COUNT, DEFAULT_MAX_ITERATIONS_COUNT);
    }

    /**
     * Builds an integrator with given iteration counts.
     *
     * @param field field to which function argument and value belong
     * @param n Number of integration points.
     * @param minimalIterationCount Minimum number of iterations.
     * @param maximalIterationCount Maximum number of iterations.
     * @throws MathIllegalArgumentException if minimal number of iterations
     * is not strictly positive.
     * @throws MathIllegalArgumentException if maximal number of iterations
     * is smaller than or equal to the minimal number of iterations.
     * @throws MathIllegalArgumentException if {@code n < 1}.
     */
    public IterativeLegendreFieldGaussIntegrator(final Field<T> field, final int n,
                                                 final int minimalIterationCount,
                                                 final int maximalIterationCount)
                                                                 throws MathIllegalArgumentException {
        this(field, n, DEFAULT_RELATIVE_ACCURACY, DEFAULT_ABSOLUTE_ACCURACY,
             minimalIterationCount, maximalIterationCount);
    }

    /** {@inheritDoc} */
    @Override
    protected T doIntegrate()
        throws MathIllegalArgumentException, MathIllegalStateException {
        // Compute first estimate with a single step.
        T oldt = stage(1);

        int n = 2;
        while (true) {
            // Improve integral with a larger number of steps.
            final T t = stage(n);

            // Estimate the error.
            final double delta = FastMath.abs(t.subtract(oldt)).getReal();
            final double limit =
                FastMath.max(getAbsoluteAccuracy(),
                             FastMath.abs(oldt).add(FastMath.abs(t)).multiply(0.5 * getRelativeAccuracy()).getReal());

            // check convergence
            if (iterations.getCount() + 1 >= getMinimalIterationCount() &&
                delta <= limit) {
                return t;
            }

            // Prepare next iteration.
            final double ratio = FastMath.min(4, FastMath.pow(delta / limit, 0.5 / numberOfPoints));
            n = FastMath.max((int) (ratio * n), n + 1);
            oldt = t;
            iterations.increment();
        }
    }

    /**
     * Compute the n-th stage integral.
     *
     * @param n Number of steps.
     * @return the value of n-th stage integral.
     * @throws MathIllegalStateException if the maximum number of evaluations
     * is exceeded.
     */
    private T stage(final int n)
        throws MathIllegalStateException {

        final T min = getMin();
        final T max = getMax();
        final T step = max.subtract(min).divide(n);

        T sum = getField().getZero();
        for (int i = 0; i < n; i++) {
            // Integrate over each sub-interval [a, b].
            final T a = min.add(step.multiply(i));
            final T b = a.add(step);
            final FieldGaussIntegrator<T> g = factory.legendre(numberOfPoints, a, b);
            sum = sum.add(g.integrate(super::computeObjectiveValue));
        }

        return sum;
    }

}
