/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
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
import org.hipparchus.Field;
import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.analysis.solvers.AllowedSolution;
import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver;
import org.hipparchus.analysis.solvers.FieldBracketingNthOrderBrentSolver;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
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
 * The initial interval for the application of the bisection method is
 * based on the roots of the previous Hermite polynomial (interlacing).
 * Upper and lower bounds of these roots are provided by </p>
 * <blockquote>
 *  I. Krasikov,
 *  <em>Nonnegative quadratic forms and bounds on orthogonal polynomials</em>,
 *  Journal of Approximation theory <b>111</b>, 31-49
 * </blockquote>
 * @param <T> Type of the number used to represent the points and weights of
 * the quadrature rules.
 * @since 2.0
 */
public class FieldHermiteRuleFactory<T extends CalculusFieldElement<T>> extends FieldAbstractRuleFactory<T> {
    /** &pi;<sup>1/2</sup> */
    private final T sqrtPi;
    /** &pi;<sup>-1/4</sup> */
    private final T h0;
    /** &pi;<sup>-1/4</sup> &radic;2 */
    private final T h1;

    /** Simple constructor
     * @param field field to which rule coefficients belong
     */
    public FieldHermiteRuleFactory(final Field<T> field) {
        sqrtPi = field.getZero().getPi().sqrt();
        h0     = sqrtPi.sqrt().reciprocal();
        h1     = h0.multiply(field.getOne().newInstance(2).sqrt());
    }

    /** {@inheritDoc} */
    @Override
    protected Pair<T[], T[]> computeRule(int numberOfPoints)
        throws MathIllegalArgumentException {

        final Field<T> field   = sqrtPi.getField();
        final T        zero    = field.getZero();
        final T[]      points  = MathArrays.buildArray(field, numberOfPoints);
        final T[]      weights = MathArrays.buildArray(field, numberOfPoints);

        if (numberOfPoints == 1) {
            // Break recursion.
            points[0]  = zero;
            weights[0] = sqrtPi;
            return new Pair<>(points, weights);
        }

        // Get previous rule.
        // If it has not been computed yet it will trigger a recursive call
        // to this method.
        final int lastNumPoints = numberOfPoints - 1;
        final T[] previousPoints = getRule(lastNumPoints).getFirst();
        final NormalizedHermite hm = new NormalizedHermite(numberOfPoints - 1);
        final NormalizedHermite h  = new NormalizedHermite(numberOfPoints);
        final T tol = field.getOne().ulp().multiply(10);
        final BracketedRealFieldUnivariateSolver<T> solver = new FieldBracketingNthOrderBrentSolver<>(tol, tol, tol, 5);

        final T sqrtTwoTimesLastNumPoints = FastMath.sqrt(zero.newInstance(2 * lastNumPoints));
        final T sqrtTwoTimesNumPoints = FastMath.sqrt(zero.newInstance(2 * numberOfPoints));

        // Find i-th root of H[n+1] by bracketing.
        final int iMax = numberOfPoints / 2;
        for (int i = 0; i < iMax; i++) {
            // Lower-bound of the interval.
            T a = (i == 0) ? sqrtTwoTimesLastNumPoints.negate() : previousPoints[i - 1];
            // Upper-bound of the interval.
            T b = (iMax == 1) ? zero.newInstance(-0.5) : previousPoints[i];
            // find root
            final T c = solver.solve(1000, h, a, b, AllowedSolution.ANY_SIDE);
            points[i] = c;

            final T d = sqrtTwoTimesNumPoints.multiply(hm.value(c));
            weights[i] = d.multiply(d).reciprocal().multiply(2);

            // symmetrical point
            final int idx = lastNumPoints - i;
            points[idx]   = c.negate();
            weights[idx]  = weights[i];

        }

        // If "numberOfPoints" is odd, 0 is a root.
        // Note: as written, the test for oddness will work for negative
        // integers too (although it is not necessary here), preventing
        // a FindBugs warning.
        if (numberOfPoints % 2 != 0) {
            T hmz = h0;
            for (int j = 1; j < numberOfPoints; j += 2) {
                final double jp1 = j + 1;
                hmz = hmz.multiply(zero.newInstance(j).divide(jp1).sqrt()).negate();
            }
            final T d = sqrtTwoTimesNumPoints.multiply(hmz);
            points[iMax] = zero;
            weights[iMax] = d.multiply(d).divide(2).reciprocal();

        }

        return new Pair<>(points, weights);

    }

    /** Hermite polynomial, normalized to avoid overflow. */
    private class NormalizedHermite implements CalculusFieldUnivariateFunction<T> {

        /** Degree. */
        private int degree;

        /** Simple constructor.
         * @param degree polynomial degree
         */
        NormalizedHermite(int degree) {
            this.degree = degree;
        }

        /** {@inheritDoc} */
        @Override
        public T value(T x) {
            T hm = h0;
            T h  = h1.multiply(x);
            for (int j = 1; j < degree; j++) {
                // Compute H[j](x)
                final double jp1 = j + 1;
                final T hp = x.multiply(h).multiply(x.newInstance(2).divide(jp1).sqrt()).
                             subtract(hm.multiply(x.newInstance(j).divide(jp1).sqrt()));
                hm = h;
                h  = hp;
            }
            return h;
        }

    }

}
