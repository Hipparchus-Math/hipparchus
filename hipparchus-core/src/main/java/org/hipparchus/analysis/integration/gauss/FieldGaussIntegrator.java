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
package org.hipparchus.analysis.integration.gauss;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.Pair;

/**
 * Class that implements the Gaussian rule for
 * {@link #integrate(CalculusFieldUnivariateFunction) integrating} a weighted
 * function.
 *
 * @param <T> Type of the field elements.
 * @since 2.0
 */
public class FieldGaussIntegrator<T extends CalculusFieldElement<T>> {
    /** Nodes. */
    private final T[] points;
    /** Nodes weights. */
    private final T[] weights;

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
    public FieldGaussIntegrator(T[] points, T[] weights)
        throws MathIllegalArgumentException {

        MathArrays.checkEqualLength(points, weights);
        MathArrays.checkOrder(points, MathArrays.OrderDirection.INCREASING, true, true);

        this.points = points.clone();
        this.weights = weights.clone();
    }

    /**
     * Creates an integrator from the given pair of points (first element of
     * the pair) and weights (second element of the pair.
     *
     * @param pointsAndWeights Integration points and corresponding weights.
     * @throws MathIllegalArgumentException if the {@code points} are not
     * sorted in increasing order.
     *
     * @see #FieldGaussIntegrator(CalculusFieldElement[], CalculusFieldElement[])
     */
    public FieldGaussIntegrator(Pair<T[], T[]> pointsAndWeights)
        throws MathIllegalArgumentException {
        this(pointsAndWeights.getFirst(), pointsAndWeights.getSecond());
    }

    /**
     * Returns an estimate of the integral of {@code f(x) * w(x)},
     * where {@code w} is a weight function that depends on the actual
     * flavor of the Gauss integration scheme.
     * The algorithm uses the points and associated weights, as passed
     * to the {@link #FieldGaussIntegrator(CalculusFieldElement[], CalculusFieldElement[]) constructor}.
     *
     * @param f Function to integrate.
     * @return the integral of the weighted function.
     */
    public T integrate(CalculusFieldUnivariateFunction<T> f) {
        T s = points[0].getField().getZero();
        T c = s;
        for (int i = 0; i < points.length; i++) {
            final T x = points[i];
            final T w = weights[i];
            final T y = w.multiply(f.value(x)).subtract(c);
            final T t = s.add(y);
            c = t.subtract(s).subtract(y);
            s = t;
        }
        return s;
    }

    /** Get order of the integration rule.
     * @return the order of the integration rule (the number of integration
     * points).
     */
    public int getNumberOfPoints() {
        return points.length;
    }

    /**
     * Gets the integration point at the given index.
     * The index must be in the valid range but no check is performed.
     * @param index index of the integration point
     * @return the integration point.
     */
    public T getPoint(int index) {
        return points[index];
    }

    /**
     * Gets the weight of the integration point at the given index.
     * The index must be in the valid range but no check is performed.
     * @param index index of the integration point
     * @return the weight.
     */
    public T getWeight(int index) {
        return weights[index];
    }
}
