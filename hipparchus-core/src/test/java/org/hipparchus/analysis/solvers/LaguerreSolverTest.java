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
package org.hipparchus.analysis.solvers;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.analysis.polynomials.PolynomialFunction;
import org.hipparchus.complex.Complex;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Test case for Laguerre solver.
 * <p>
 * Laguerre's method is very efficient in solving polynomials. Test runs
 * show that for a default absolute accuracy of 1E-6, it generally takes
 * less than 5 iterations to find one root, provided solveAll() is not
 * invoked, and 15 to 20 iterations to find all roots for quintic function.
 *
 */
public final class LaguerreSolverTest {
    /**
     * Test of solver for the linear function.
     */
    @Test
    public void testLinearFunction() {
        double min, max, expected, result, tolerance;

        // p(x) = 4x - 1
        double[] coefficients = { -1.0, 4.0 };
        PolynomialFunction f = new PolynomialFunction(coefficients);
        LaguerreSolver solver = new LaguerreSolver();

        min = 0.0; max = 1.0; expected = 0.25;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assertions.assertEquals(expected, result, tolerance);
    }

    /**
     * Test of solver for the quadratic function.
     */
    @Test
    public void testQuadraticFunction() {
        double min, max, expected, result, tolerance;

        // p(x) = 2x^2 + 5x - 3 = (x+3)(2x-1)
        double[] coefficients = { -3.0, 5.0, 2.0 };
        PolynomialFunction f = new PolynomialFunction(coefficients);
        LaguerreSolver solver = new LaguerreSolver();

        min = 0.0; max = 2.0; expected = 0.5;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assertions.assertEquals(expected, result, tolerance);

        min = -4.0; max = -1.0; expected = -3.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assertions.assertEquals(expected, result, tolerance);
    }

    /**
     * Test of solver for the quintic function.
     */
    @Test
    public void testQuinticFunction() {
        double min, max, expected, result, tolerance;

        // p(x) = x^5 - x^4 - 12x^3 + x^2 - x - 12 = (x+1)(x+3)(x-4)(x^2-x+1)
        double[] coefficients = { -12.0, -1.0, 1.0, -12.0, -1.0, 1.0 };
        PolynomialFunction f = new PolynomialFunction(coefficients);
        LaguerreSolver solver = new LaguerreSolver();

        min = -2.0; max = 2.0; expected = -1.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assertions.assertEquals(expected, result, tolerance);

        min = -5.0; max = -2.5; expected = -3.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assertions.assertEquals(expected, result, tolerance);

        min = 3.0; max = 6.0; expected = 4.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assertions.assertEquals(expected, result, tolerance);
    }

    /**
     * Test of solver for the quintic function using
     * {@link LaguerreSolver#solveAllComplex(double[],double) solveAllComplex}.
     */
    @Test
    public void testQuinticFunction2() {
        // p(x) = x^5 + 4x^3 + x^2 + 4 = (x+1)(x^2-x+1)(x^2+4)
        final double[] coefficients = { 4.0, 0.0, 1.0, 4.0, 0.0, 1.0 };
        final LaguerreSolver solver = new LaguerreSolver();
        final Complex[] result = solver.solveAllComplex(coefficients, 0);

        for (Complex expected : new Complex[] { new Complex(0, -2),
                                                new Complex(0, 2),
                                                new Complex(0.5, 0.5 * FastMath.sqrt(3)),
                                                new Complex(-1, 0),
                                                new Complex(0.5, -0.5 * FastMath.sqrt(3.0)) }) {
            final double tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                                                  FastMath.abs(expected.norm() * solver.getRelativeAccuracy()));
            UnitTestUtils.assertContains(result, expected, tolerance);
        }
    }

    /**
     * Test of parameters for the solver.
     */
    @Test
    public void testParameters() {
        double[] coefficients = { -3.0, 5.0, 2.0 };
        PolynomialFunction f = new PolynomialFunction(coefficients);
        LaguerreSolver solver = new LaguerreSolver();

        try {
            // bad interval
            solver.solve(100, f, 1, -1);
            Assertions.fail("Expecting MathIllegalArgumentException - bad interval");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            // no bracketing
            solver.solve(100, f, 2, 3);
            Assertions.fail("Expecting MathIllegalArgumentException - no bracketing");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testIssue177() {
        doTestIssue177(new double[] {-100.0,    0.0, 0.0, 0.0, 1.0}, FastMath.sqrt(10.0));
        doTestIssue177(new double[] {        -100.0, 0.0, 0.0, 1.0}, FastMath.cbrt(100.0));
        doTestIssue177(new double[] { -16.0,    0.0, 0.0, 0.0, 1.0}, 2.0);
    }

    private void doTestIssue177(final double[] coefficients, final double expected) {
        Complex[] roots = new LaguerreSolver(1.0e-5).solveAllComplex(coefficients, 0);
        Assertions.assertEquals(coefficients.length - 1, roots.length);
        for (final Complex root : roots) {
            Assertions.assertEquals(expected, root.norm(), 1.0e-15);
            Assertions.assertEquals(0.0, MathUtils.normalizeAngle(roots.length * root.getArgument(), 0.0), 1.0e-15);
        }
    }

}
