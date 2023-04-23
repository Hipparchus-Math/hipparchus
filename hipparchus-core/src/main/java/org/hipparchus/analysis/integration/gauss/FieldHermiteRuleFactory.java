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
 * Factory that creates a
 * <a href="http://en.wikipedia.org/wiki/Gauss-Hermite_quadrature">
 * Gauss-type quadrature rule using Hermite polynomials</a>
 * of the first kind.
 * Such a quadrature rule allows the calculation of improper integrals
 * of a function
 * <p>
 *  \(f(x) e^{-x^2}\)
 * </p><p>
 * Recurrence relation and weights computation follow
 * <a href="http://en.wikipedia.org/wiki/Abramowitz_and_Stegun">
 * Abramowitz and Stegun, 1964</a>.
 * </p><p>
 * The coefficients of the standard Hermite polynomials grow very rapidly.
 * In order to avoid overflows, each Hermite polynomial is normalized with
 * respect to the underlying scalar product.
 * @param <T> Type of the number used to represent the points and weights of
 * the quadrature rules.
 * @since 2.0
 */
public class FieldHermiteRuleFactory<T extends CalculusFieldElement<T>> extends FieldAbstractRuleFactory<T> {

    /** Simple constructor
     * @param field field to which rule coefficients belong
     */
    public FieldHermiteRuleFactory(final Field<T> field) {
        super(field);
    }

    /** {@inheritDoc} */
    @Override
    protected Pair<T[], T[]> computeRule(int numberOfPoints)
        throws MathIllegalArgumentException {

        final Field<T> field  = getField();
        final T        sqrtPi = field.getZero().getPi().sqrt();

        if (numberOfPoints == 1) {
            // Break recursion.
            final T[] points  = MathArrays.buildArray(field, numberOfPoints);
            final T[] weights = MathArrays.buildArray(field, numberOfPoints);
            points[0]  = field.getZero();
            weights[0] = sqrtPi;
            return new Pair<>(points, weights);
        }

        // find nodes as roots of Hermite polynomial
        final T[] points = findRoots(numberOfPoints, new Hermite<>(field, numberOfPoints)::ratio);
        enforceSymmetry(points);

        // compute weights
        final T[] weights = MathArrays.buildArray(field, numberOfPoints);
        final Hermite<T> hm1 = new Hermite<>(field, numberOfPoints - 1);
        for (int i = 0; i < numberOfPoints; i++) {
            final T y = hm1.hNhNm1(points[i])[0];
            weights[i] = sqrtPi.divide(y.multiply(y).multiply(numberOfPoints));
        }

        return new Pair<>(points, weights);

    }

    /** Hermite polynomial, normalized to avoid overflow.
     * <p>
     * The regular Hermite polynomials and associated weights are given by:
     *   <pre>
     *     H₀(x)   = 1
     *     H₁(x)   = 2 x
     *     Hₙ₊₁(x) = 2x Hₙ(x) - 2n Hₙ₋₁(x), and H'ₙ(x) = 2n Hₙ₋₁(x)
     *     wₙ(xᵢ) = [2ⁿ⁻¹ n! √π]/[n Hₙ₋₁(xᵢ)]²
     *   </pre>
     * </p>
     * <p>
     * In order to avoid overflow with normalize the polynomials hₙ(x) = Hₙ(x) / √[2ⁿ n!]
     * so the recurrence relations and weights become:
     *   <pre>
     *     h₀(x)   = 1
     *     h₁(x)   = √2 x
     *     hₙ₊₁(x) = [√2 x hₙ(x) - √n hₙ₋₁(x)]/√(n+1), and h'ₙ(x) = 2n hₙ₋₁(x)
     *     uₙ(xᵢ) = √π/[n Nₙ₋₁(xᵢ)²]
     *   </pre>
     * </p>
     * @param <T> Type of the field elements.
     */
    private static class Hermite<T extends CalculusFieldElement<T>> {

        /** √2. */
        private final T sqrt2;

        /** Degree. */
        private final int degree;

        /** Simple constructor.
         * @param field field to which rule coefficients belong
         * @param degree polynomial degree
         */
        Hermite(Field<T> field, int degree) {
            this.sqrt2  = field.getZero().newInstance(2).sqrt();
            this.degree = degree;
        }

        /** Compute ratio H(x)/H'(x).
         * @param x point at which ratio must be computed
         * @return ratio H(x)/H'(x)
         */
        public T ratio(T x) {
            T[] h = hNhNm1(x);
            return h[0].divide(h[1].multiply(2 * degree));
        }

        /** Compute Nₙ(x) and Nₙ₋₁(x).
         * @param x point at which polynomials are evaluated
         * @return array containing Nₙ(x) at index 0 and Nₙ₋₁(x) at index 1
         */
        private T[] hNhNm1(final T x) {
            T[] h = MathArrays.buildArray(x.getField(), 2);
            h[0] = sqrt2.multiply(x);
            h[1] = x.getField().getOne();
            T sqrtN = x.getField().getOne();
            for (int n = 1; n < degree; n++) {
                // apply recurrence relation hₙ₊₁(x) = [√2 x hₙ(x) - √n hₙ₋₁(x)]/√(n+1)
                final T sqrtNp = x.getField().getZero().newInstance(n + 1).sqrt();
                final T hp = (h[0].multiply(x).multiply(sqrt2).subtract(h[1].multiply(sqrtN))).divide(sqrtNp);
                h[1]  = h[0];
                h[0]  = hp;
                sqrtN = sqrtNp;
            }
            return h;
        }

    }

}
