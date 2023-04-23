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
package org.hipparchus.analysis.integration.gauss;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Pair;

/**
 * This class's implements {@link #integrate(CalculusFieldUnivariateFunction) integrate}
 * method assuming that the integral is symmetric about 0.
 * This allows to reduce numerical errors.
 *
 * @param <T> Type of the field elements.
 * @since 2.0
 */
public class SymmetricFieldGaussIntegrator<T extends CalculusFieldElement<T>> extends FieldGaussIntegrator<T> {
    /**
     * Creates an integrator from the given {@code points} and {@code weights}.
     * The integration interval is defined by the first and last value of
     * {@code points} which must be sorted in increasing order.
     *
     * @param points Integration points.
     * @param weights Weights of the corresponding integration nodes.
     * @throws MathIllegalArgumentException if the {@code points} are not
     * sorted in increasing order.
     * @throws MathIllegalArgumentException if points and weights don't have the same length
     */
    public SymmetricFieldGaussIntegrator(T[] points, T[] weights)
        throws MathIllegalArgumentException {
        super(points, weights);
    }

    /**
     * Creates an integrator from the given pair of points (first element of
     * the pair) and weights (second element of the pair.
     *
     * @param pointsAndWeights Integration points and corresponding weights.
     * @throws MathIllegalArgumentException if the {@code points} are not
     * sorted in increasing order.
     *
     * @see #SymmetricFieldGaussIntegrator(CalculusFieldElement[], CalculusFieldElement[])
     */
    public SymmetricFieldGaussIntegrator(Pair<T[], T[]> pointsAndWeights)
        throws MathIllegalArgumentException {
        this(pointsAndWeights.getFirst(), pointsAndWeights.getSecond());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public T integrate(CalculusFieldUnivariateFunction<T> f) {
        final int ruleLength = getNumberOfPoints();

        final T zero = getPoint(0).getField().getZero();
        if (ruleLength == 1) {
            return getWeight(0).multiply(f.value(zero));
        }

        final int iMax = ruleLength / 2;
        T s = zero;
        T c = zero;
        for (int i = 0; i < iMax; i++) {
            final T p = getPoint(i);
            final T w = getWeight(i);

            final T f1 = f.value(p);
            final T f2 = f.value(p.negate());

            final T y = w.multiply(f1.add(f2)).subtract(c);
            final T t = s.add(y);

            c = t.subtract(s).subtract(y);
            s = t;
        }

        if (ruleLength % 2 != 0) {
            final T w = getWeight(iMax);

            final T y = w.multiply(f.value(zero)).subtract(c);
            final T t = s.add(y);

            s = t;
        }

        return s;
    }
}
