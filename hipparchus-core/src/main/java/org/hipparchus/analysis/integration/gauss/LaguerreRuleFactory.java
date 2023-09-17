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

import org.hipparchus.util.Pair;

/**
 * Factory that creates Gauss-type quadrature rule using Laguerre polynomials.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature">Gauss-Laguerre quadrature (Wikipedia)</a>
 */
public class LaguerreRuleFactory extends AbstractRuleFactory {

    /** Empty constructor.
     * <p>
     * This constructor is not strictly necessary, but it prevents spurious
     * javadoc warnings with JDK 18 and later.
     * </p>
     * @since 3.0
     */
    public LaguerreRuleFactory() { // NOPMD - unnecessary constructor added intentionally to make javadoc happy
        // nothing to do
    }

    /** {@inheritDoc} */
    @Override
    protected Pair<double[], double[]> computeRule(int numberOfPoints) {

        // find nodes as roots of Laguerre polynomial
        final double[] points  = findRoots(numberOfPoints, new Laguerre(numberOfPoints)::ratio);

        // compute weights
        final double[] weights    = new double[numberOfPoints];
        final int      n1         = numberOfPoints + 1;
        final long     n1Squared  = n1 * (long) n1;
        final Laguerre laguerreN1 = new Laguerre(n1);
        for (int i = 0; i < numberOfPoints; i++) {
            final double val = laguerreN1.value(points[i]);
            weights[i] = points[i] / (n1Squared * val * val);
        }

        return new Pair<>(points, weights);

    }

    /** Laguerre polynomial. */
    private static class Laguerre {

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
        public double value(final double x) {
            return lNlNm1(x)[0];
        }

        /** Compute ratio L(x)/L'(x).
         * @param x point at which ratio must be computed
         * @return ratio L(x)/L'(x)
         */
        public double ratio(double x) {
            double[] l = lNlNm1(x);
            return x * l[0] / (degree * (l[0] - l[1]));
        }

        /** Compute Lₙ(x) and Lₙ₋₁(x).
         * @param x point at which polynomials are evaluated
         * @return array containing Lₙ(x) at index 0 and Lₙ₋₁(x) at index 1
         */
        private double[] lNlNm1(final double x) {
            double[] l = { 1 - x, 1 };
            for (int n = 1; n < degree; n++) {
                // apply recurrence relation (n+1) Lₙ₊₁(x) = (2n + 1 - x) Lₙ(x) - n Lₙ₋₁(x)
                final double lp = (l[0] * (2 * n + 1 - x) - l[1] * n) / (n + 1);
                l[1] = l[0];
                l[0] = lp;
            }
            return l;
        }

    }

}
