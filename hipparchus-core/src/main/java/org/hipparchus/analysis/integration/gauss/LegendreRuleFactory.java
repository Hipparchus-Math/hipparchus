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
import org.hipparchus.util.Pair;

/**
 * Factory that creates Gauss-type quadrature rule using Legendre polynomials.
 * In this implementation, the lower and upper bounds of the natural interval
 * of integration are -1 and 1, respectively.
 * The Legendre polynomials are evaluated using the recurrence relation
 * presented in <a href="http://en.wikipedia.org/wiki/Abramowitz_and_Stegun">
 * Abramowitz and Stegun, 1964</a>.
 *
 */
public class LegendreRuleFactory extends AbstractRuleFactory {

    /** Empty constructor.
     * <p>
     * This constructor is not strictly necessary, but it prevents spurious
     * javadoc warnings with JDK 18 and later.
     * </p>
     * @since 3.0
     */
    public LegendreRuleFactory() { // NOPMD - unnecessary constructor added intentionally to make javadoc happy
        // nothing to do
    }

    /** {@inheritDoc} */
    @Override
    protected Pair<double[], double[]> computeRule(int numberOfPoints)
        throws MathIllegalArgumentException {

        if (numberOfPoints == 1) {
            // Break recursion.
           return new Pair<>(new double[] { 0 } , new double[] { 2 });
        }

        // find nodes as roots of Legendre polynomial
        final Legendre p      =  new Legendre(numberOfPoints);
        final double[] points = findRoots(numberOfPoints, p::ratio);
        enforceSymmetry(points);

        // compute weights
        final double[] weights = new double[numberOfPoints];
        for (int i = 0; i <= numberOfPoints / 2; i++) {
            final double c = points[i];
            final double[] pKpKm1 = p.pNpNm1(c);
            final double d = numberOfPoints * (pKpKm1[1] - c * pKpKm1[0]);
            weights[i] = 2 * (1 - c * c) / (d * d);

            // symmetrical point
            final int idx = numberOfPoints - i - 1;
            weights[idx]  = weights[i];

        }

        return new Pair<>(points, weights);

    }

    /** Legendre polynomial. */
    private static class Legendre {

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
        public double ratio(double x) {
            double pm = 1;
            double p  = x;
            double d  = 1;
            for (int n = 1; n < degree; n++) {
                // apply recurrence relations (n+1) Pₙ₊₁(x)  = (2n+1) x Pₙ(x) - n Pₙ₋₁(x)
                // and                              P'ₙ₊₁(x) = (n+1) Pₙ(x) + x P'ₙ(x)
                final double pp = (p * (x * (2 * n + 1)) - pm * n) / (n + 1);
                d  = p * (n + 1) + d * x;
                pm = p;
                p  = pp;
            }
            return p / d;
        }

        /** Compute Pₙ(x) and Pₙ₋₁(x).
         * @param x point at which polynomials are evaluated
         * @return array containing Pₙ(x) at index 0 and Pₙ₋₁(x) at index 1
         */
        private double[] pNpNm1(final double x) {
            double[] p = { x, 1 };
            for (int n = 1; n < degree; n++) {
                // apply recurrence relation (n+1) Pₙ₊₁(x) = (2n+1) x Pₙ(x) - n Pₙ₋₁(x)
                final double pp = (p[0] * (x * (2 * n + 1)) - p[1] * n) / (n + 1);
                p[1] = p[0];
                p[0] = pp;
            }
            return p;
        }

    }

}
