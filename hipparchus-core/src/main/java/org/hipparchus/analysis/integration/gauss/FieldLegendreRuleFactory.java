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
 * Factory that creates Gauss-type quadrature rule using Legendre polynomials.
 * In this implementation, the lower and upper bounds of the natural interval
 * of integration are -1 and 1, respectively.
 * The Legendre polynomials are evaluated using the recurrence relation
 * presented in <a href="http://en.wikipedia.org/wiki/Abramowitz_and_Stegun">
 * Abramowitz and Stegun, 1964</a>.
 *
 * @param <T> Type of the number used to represent the points and weights of
 * the quadrature rules.
 * @since 2.0
 */
public class FieldLegendreRuleFactory<T extends CalculusFieldElement<T>> extends FieldAbstractRuleFactory<T> {

    /** Simple constructor
     * @param field field to which rule coefficients belong
     */
    public FieldLegendreRuleFactory(final Field<T> field) {
        super(field);
    }

    /** {@inheritDoc} */
    @Override
    public Pair<T[], T[]> computeRule(int numberOfPoints)
        throws MathIllegalArgumentException {

        final Field<T> field = getField();

        if (numberOfPoints == 1) {
            // Break recursion.
            final T[] points  = MathArrays.buildArray(field, numberOfPoints);
            final T[] weights = MathArrays.buildArray(field, numberOfPoints);
            points[0]  = field.getZero();
            weights[0] = field.getZero().newInstance(2);
            return new Pair<>(points, weights);
        }

        // find nodes as roots of Legendre polynomial
        final Legendre<T> p      =  new Legendre<>(numberOfPoints);
        final T[]         points = findRoots(numberOfPoints, p::ratio);
        enforceSymmetry(points);

        // compute weights
        final T[] weights = MathArrays.buildArray(field, numberOfPoints);
        for (int i = 0; i <= numberOfPoints / 2; i++) {
            final T c = points[i];
            final T[] pKpKm1 = p.pNpNm1(c);
            final T d = pKpKm1[1].subtract(c.multiply(pKpKm1[0])).multiply(numberOfPoints);
            weights[i] = c.multiply(c).subtract(1).multiply(-2).divide(d.multiply(d));

            // symmetrical point
            final int idx = numberOfPoints - i - 1;
            weights[idx]  = weights[i];

        }

        return new Pair<>(points, weights);

    }

    /** Legendre polynomial.
     * @param <T> Type of the field elements.
     */
    private static class Legendre<T extends CalculusFieldElement<T>> {

        /** Degree. */
        private int degree;

        /** Simple constructor.
         * @param degree polynomial degree
         */
        Legendre(int degree) {
            this.degree = degree;
        }

        /** Compute ratio P(x)/P'(x).
         * @param x point at which ratio must be computed
         * @return ratio P(x)/P'(x)
         */
        public T ratio(T x) {
            T pm = x.getField().getOne();
            T p  = x;
            T d  = x.getField().getOne();
            for (int n = 1; n < degree; n++) {
                // apply recurrence relations (n+1) Pₙ₊₁(x)  = (2n+1) x Pₙ(x) - n Pₙ₋₁(x)
                // and                              P'ₙ₊₁(x) = (n+1) Pₙ(x) + x P'ₙ(x)
                final T pp = p.multiply(x.multiply(2 * n + 1)).subtract(pm.multiply(n)).divide(n + 1);
                d  = p.multiply(n + 1).add(d.multiply(x));
                pm = p;
                p  = pp;
            }
            return p.divide(d);
        }

        /** Compute Pₙ(x) and Pₙ₋₁(x).
         * @param x point at which polynomials are evaluated
         * @return array containing Pₙ(x) at index 0 and Pₙ₋₁(x) at index 1
         */
        private T[] pNpNm1(final T x) {
            T[] p = MathArrays.buildArray(x.getField(), 2);
            p[0] = x;
            p[1] = x.getField().getOne();
            for (int n = 1; n < degree; n++) {
                // apply recurrence relation (n+1) Pₙ₊₁(x) = (2n+1) x Pₙ(x) - n Pₙ₋₁(x)
                final T pp = p[0].multiply(x.multiply(2 * n + 1)).subtract(p[1].multiply(n)).divide(n + 1);
                p[1] = p[0];
                p[0] = pp;
            }
            return p;
        }

    }

}
