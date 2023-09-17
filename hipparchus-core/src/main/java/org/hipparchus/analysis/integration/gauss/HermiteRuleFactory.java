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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
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
 * </p>
 * <p>
 * Recurrence relation and weights computation follow
 * <a href="http://en.wikipedia.org/wiki/Abramowitz_and_Stegun">
 * Abramowitz and Stegun, 1964</a>.
 * </p>
 *
 */
public class HermiteRuleFactory extends AbstractRuleFactory {

    /** √π. */
    private static final double SQRT_PI = 1.77245385090551602729;

    /** Empty constructor.
     * <p>
     * This constructor is not strictly necessary, but it prevents spurious
     * javadoc warnings with JDK 18 and later.
     * </p>
     * @since 3.0
     */
    public HermiteRuleFactory() { // NOPMD - unnecessary constructor added intentionally to make javadoc happy
        // nothing to do
    }

    /** {@inheritDoc} */
    @Override
    protected Pair<double[], double[]> computeRule(int numberOfPoints)
        throws MathIllegalArgumentException {

        if (numberOfPoints == 1) {
            // Break recursion.
            return new Pair<>(new double[] { 0 } , new double[] { SQRT_PI });
        }

        // find nodes as roots of Hermite polynomial
        final double[] points = findRoots(numberOfPoints, new Hermite(numberOfPoints)::ratio);
        enforceSymmetry(points);

        // compute weights
        final double[] weights = new double[numberOfPoints];
        final Hermite hm1 = new Hermite(numberOfPoints - 1);
        for (int i = 0; i < numberOfPoints; i++) {
            final double y = hm1.hNhNm1(points[i])[0];
            weights[i] = SQRT_PI / (numberOfPoints * y * y);
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
     */
    private static class Hermite {

        /** √2. */
        private static final double SQRT2 = FastMath.sqrt(2);

        /** Degree. */
        private final int degree;

        /** Simple constructor.
         * @param degree polynomial degree
         */
        Hermite(int degree) {
            this.degree = degree;
        }

        /** Compute ratio H(x)/H'(x).
         * @param x point at which ratio must be computed
         * @return ratio H(x)/H'(x)
         */
        public double ratio(double x) {
            double[] h = hNhNm1(x);
            return h[0] / (h[1] * 2 * degree);
        }

        /** Compute Nₙ(x) and Nₙ₋₁(x).
         * @param x point at which polynomials are evaluated
         * @return array containing Nₙ(x) at index 0 and Nₙ₋₁(x) at index 1
         */
        private double[] hNhNm1(final double x) {
            double[] h = { SQRT2 * x, 1 };
            double sqrtN = 1;
            for (int n = 1; n < degree; n++) {
                // apply recurrence relation hₙ₊₁(x) = [√2 x hₙ(x) - √n hₙ₋₁(x)]/√(n+1)
                final double sqrtNp = FastMath.sqrt(n + 1);
                final double hp = (h[0] * x * SQRT2 - h[1] * sqrtN) / sqrtNp;
                h[1]  = h[0];
                h[0]  = hp;
                sqrtN = sqrtNp;
            }
            return h;
        }

    }

}
