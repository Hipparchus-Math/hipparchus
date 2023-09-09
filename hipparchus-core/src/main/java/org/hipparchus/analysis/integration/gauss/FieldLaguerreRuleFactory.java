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
import org.hipparchus.Field;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.Pair;

/**
 * Factory that creates Gauss-type quadrature rule using Laguerre polynomials.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature">Gauss-Laguerre quadrature (Wikipedia)</a>
 * @param <T> Type of the number used to represent the points and weights of
 * the quadrature rules.
 * @since 2.0
 */
public class FieldLaguerreRuleFactory<T extends CalculusFieldElement<T>> extends FieldAbstractRuleFactory<T> {

    /** Simple constructor
     * @param field field to which rule coefficients belong
     */
    public FieldLaguerreRuleFactory(final Field<T> field) {
        super(field);
    }

    /** {@inheritDoc} */
    @Override
    public Pair<T[], T[]> computeRule(int numberOfPoints)
        throws MathIllegalArgumentException {

        final Field<T> field = getField();

        // find nodes as roots of Laguerre polynomial
        final Laguerre<T> p      =  new Laguerre<>(numberOfPoints);
        final T[]      points = findRoots(numberOfPoints, p::ratio);

        // compute weights
        final T[] weights = MathArrays.buildArray(field, numberOfPoints);
        final int      n1         = numberOfPoints + 1;
        final long     n1Squared  = n1 * (long) n1;
        final Laguerre<T> laguerreN1 = new Laguerre<>(n1);
        for (int i = 0; i < numberOfPoints; i++) {
            final T y = laguerreN1.value(points[i]);
            weights[i] = points[i].divide(y.multiply(y).multiply(n1Squared));
        }

        return new Pair<>(points, weights);

    }

    /** Laguerre polynomial.
     * @param <T> Type of the field elements.
     */
    private static class Laguerre<T extends CalculusFieldElement<T>> {

        /** Degree. */
        private int degree;

        /** Simple constructor.
         * @param degree polynomial degree
         */
        Laguerre(int degree) {
            this.degree = degree;
        }

        /** Evaluate polynomial.
         * @param x point at which polynomial must be evaluated
         * @return value of the polynomial
         */
        public T value(final T x) {
            return lNlNm1(x)[0];
        }

        /** Compute ratio L(x)/L'(x).
         * @param x point at which ratio must be computed
         * @return ratio L(x)/L'(x)
         */
        public T ratio(T x) {
            T[] l = lNlNm1(x);
            return x.multiply(l[0]).divide(l[0].subtract(l[1]).multiply(degree));
        }

        /** Compute Lₙ(x) and Lₙ₋₁(x).
         * @param x point at which polynomials are evaluated
         * @return array containing Lₙ(x) at index 0 and Lₙ₋₁(x) at index 1
         */
        private T[] lNlNm1(final T x) {
            T[] l = MathArrays.buildArray(x.getField(), 2);
            l[0] = x.subtract(1).negate();
            l[1] = x.getField().getOne();
            for (int n = 1; n < degree; n++) {
                // apply recurrence relation (n+1) Lₙ₊₁(x) = (2n + 1 - x) Lₙ(x) - n Lₙ₋₁(x)
                final T lp = l[0].multiply(x.negate().add(2 * n + 1)).subtract(l[1].multiply(n)).divide(n + 1);
                l[1] = l[0];
                l[0] = lp;
            }
            return l;
        }

    }

}
