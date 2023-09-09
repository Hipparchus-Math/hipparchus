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
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;

/**
 * Implements <a href="http://mathworld.wolfram.com/SimpsonsRule.html">
 * Simpson's Rule</a> for integration of real univariate functions. For
 * reference, see <b>Introduction to Numerical Analysis</b>, ISBN 038795452X,
 * chapter 3.
 * <p>
 * This implementation employs the basic trapezoid rule to calculate Simpson's
 * rule.</p>
 * @param <T> Type of the field elements.
 * @since 2.0
 */
public class FieldSimpsonIntegrator<T extends CalculusFieldElement<T>> extends BaseAbstractFieldUnivariateIntegrator<T> {

    /** Maximal number of iterations for Simpson. */
    public static final int SIMPSON_MAX_ITERATIONS_COUNT = 64;

    /**
     * Build a Simpson integrator with given accuracies and iterations counts.
     * @param field field to which function argument and value belong
     * @param relativeAccuracy relative accuracy of the result
     * @param absoluteAccuracy absolute accuracy of the result
     * @param minimalIterationCount minimum number of iterations
     * @param maximalIterationCount maximum number of iterations
     * (must be less than or equal to {@link #SIMPSON_MAX_ITERATIONS_COUNT})
     * @exception MathIllegalArgumentException if minimal number of iterations
     * is not strictly positive
     * @exception MathIllegalArgumentException if maximal number of iterations
     * is lesser than or equal to the minimal number of iterations
     * @exception MathIllegalArgumentException if maximal number of iterations
     * is greater than {@link #SIMPSON_MAX_ITERATIONS_COUNT}
     */
    public FieldSimpsonIntegrator(final Field<T> field,
                                  final double relativeAccuracy,
                                  final double absoluteAccuracy,
                                  final int minimalIterationCount,
                                  final int maximalIterationCount)
        throws MathIllegalArgumentException {
        super(field, relativeAccuracy, absoluteAccuracy, minimalIterationCount, maximalIterationCount);
        if (maximalIterationCount > SIMPSON_MAX_ITERATIONS_COUNT) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_LARGE_BOUND_EXCLUDED,
                                                   maximalIterationCount, SIMPSON_MAX_ITERATIONS_COUNT);
        }
    }

    /**
     * Build a Simpson integrator with given iteration counts.
     * @param field field to which function argument and value belong
     * @param minimalIterationCount minimum number of iterations
     * @param maximalIterationCount maximum number of iterations
     * (must be less than or equal to {@link #SIMPSON_MAX_ITERATIONS_COUNT})
     * @exception MathIllegalArgumentException if minimal number of iterations
     * is not strictly positive
     * @exception MathIllegalArgumentException if maximal number of iterations
     * is lesser than or equal to the minimal number of iterations
     * @exception MathIllegalArgumentException if maximal number of iterations
     * is greater than {@link #SIMPSON_MAX_ITERATIONS_COUNT}
     */
    public FieldSimpsonIntegrator(final Field<T> field,
                                  final int minimalIterationCount,
                                  final int maximalIterationCount)
        throws MathIllegalArgumentException {
        super(field, minimalIterationCount, maximalIterationCount);
        if (maximalIterationCount > SIMPSON_MAX_ITERATIONS_COUNT) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_LARGE_BOUND_EXCLUDED,
                                                   maximalIterationCount, SIMPSON_MAX_ITERATIONS_COUNT);
        }
    }

    /**
     * Construct an integrator with default settings.
     * @param field field to which function argument and value belong
     * (max iteration count set to {@link #SIMPSON_MAX_ITERATIONS_COUNT})
     */
    public FieldSimpsonIntegrator(final Field<T> field) {
        super(field, DEFAULT_MIN_ITERATIONS_COUNT, SIMPSON_MAX_ITERATIONS_COUNT);
    }

    /** {@inheritDoc} */
    @Override
    protected T doIntegrate()
        throws MathIllegalStateException {

        FieldTrapezoidIntegrator<T> qtrap = new FieldTrapezoidIntegrator<>(getField());
        if (getMinimalIterationCount() == 1) {
            return qtrap.stage(this, 1).multiply(4).subtract(qtrap.stage(this, 0)).divide(3.0);
        }

        // Simpson's rule requires at least two trapezoid stages.
        T olds = getField().getZero();
        T oldt = qtrap.stage(this, 0);
        while (true) {
            final T t = qtrap.stage(this, iterations.getCount());
            iterations.increment();
            final T s = t.multiply(4).subtract(oldt).divide(3.0);
            if (iterations.getCount() >= getMinimalIterationCount()) {
                final double delta = FastMath.abs(s.subtract(olds)).getReal();
                final double rLimit = FastMath.abs(olds).add(FastMath.abs(s)).multiply(0.5 * getRelativeAccuracy()).getReal();
                if ((delta <= rLimit) || (delta <= getAbsoluteAccuracy())) {
                    return s;
                }
            }
            olds = s;
            oldt = t;
        }

    }

}
