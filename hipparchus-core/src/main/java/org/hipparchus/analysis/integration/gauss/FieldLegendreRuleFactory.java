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

    /** Field to which rule coefficients belong. */
    private final Field<T> field;

    /** Simple constructor
     * @param field field to which rule coefficients belong
     */
    public FieldLegendreRuleFactory(final Field<T> field) {
        this.field = field;
    }

    /** {@inheritDoc} */
    @Override
    public Pair<T[], T[]> computeRule(int numberOfPoints)
        throws MathIllegalArgumentException {

        final T[]      points  = MathArrays.buildArray(field, numberOfPoints);
        final T[]      weights = MathArrays.buildArray(field, numberOfPoints);

        if (numberOfPoints == 1) {
            // Break recursion.
            points[0]  = field.getZero();
            weights[0] = field.getZero().newInstance(2);
            return new Pair<>(points, weights);
        }

        // Get previous rule.
        // If it has not been computed yet it will trigger a recursive call
        // to this method.
        final T[] previousPoints = getRule(numberOfPoints - 1).getFirst();
        final Legendre pm = new Legendre(numberOfPoints - 1);
        final Legendre p  = new Legendre(numberOfPoints);
        final T tol = field.getOne().ulp().multiply(10);
        final BracketedRealFieldUnivariateSolver<T> solver = new FieldBracketingNthOrderBrentSolver<>(tol, tol, tol, 5);

        // Find i-th root of P[n+1]
        final int iMax = numberOfPoints / 2;
        for (int i = 0; i < iMax; i++) {
            // Lower-bound of the interval.
            final T a = (i == 0) ? field.getOne().negate() : previousPoints[i - 1];
            // Upper-bound of the interval.
            final T b = previousPoints[i];
            // find root
            final T c = solver.solve(1000, p, a, b, AllowedSolution.ANY_SIDE);
            points[i] = c;

            final T d = pm.value(c).subtract(c.multiply(p.value(c))).multiply(numberOfPoints);
            weights[i] = field.getOne().subtract(c.multiply(c)).multiply(2).divide(d.multiply(d));

            // symmetrical point
            final int idx = numberOfPoints - i - 1;
            points[idx]   = c.negate();
            weights[idx]  = weights[i];

        }

        // If "numberOfPoints" is odd, 0 is a root.
        // Note: as written, the test for oddness will work for negative
        // integers too (although it is not necessary here), preventing
        // a FindBugs warning.
        if (numberOfPoints % 2 != 0) {
            T pmz = field.getOne();
            for (int j = 1; j < numberOfPoints; j += 2) {
                // pmc = -j * pmc / (j + 1);
                pmz = pmz.multiply(-j).divide(j + 1);
            }

            final T d   = pmz.multiply(numberOfPoints);
            points[iMax]  = field.getZero();
            weights[iMax] = d.multiply(d).reciprocal().multiply(2);
        }

        return new Pair<>(points, weights);

    }

    /** Legendre polynomial. */
    private class Legendre implements CalculusFieldUnivariateFunction<T> {

        /** Degree. */
        private int degree;

        /** Simple constructor.
         * @param degree polynomial degree
         */
        Legendre(int degree) {
            this.degree = degree;
        }

        /** {@inheritDoc} */
        @Override
        public T value(T x) {
            T pm = x.getField().getOne();
            T p  = x;
            for (int k = 1; k < degree; k++) {
                // apply recurrence relation (k+1) P_{k+1}(x) = (2k+1) x P_k(x) - k P_{k-1}(x)
                final T pp = p.multiply(x.multiply(2 * k + 1)).subtract(pm.multiply(k)).divide(k + 1);
                pm = p;
                p  = pp;
            }
            return p;
        }

    }

}
