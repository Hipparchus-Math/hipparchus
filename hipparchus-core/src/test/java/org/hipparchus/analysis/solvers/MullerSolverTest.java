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

import org.hipparchus.analysis.QuinticFunction;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.function.Expm1;
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test case for {@link MullerSolver Muller} solver.
 * <p>
 * Muller's method converges almost quadratically near roots, but it can
 * be very slow in regions far away from zeros. Test runs show that for
 * reasonably good initial values, for a default absolute accuracy of 1E-6,
 * it generally takes 5 to 10 iterations for the solver to converge.
 * <p>
 * Tests for the exponential function illustrate the situations where
 * Muller solver performs poorly.
 *
 */
final class MullerSolverTest {
    /**
     * Test of solver for the sine function.
     */
    @Test
    void testSinFunction() {
        UnivariateFunction f = new Sin();
        UnivariateSolver solver = new MullerSolver();
        double min, max, expected, result, tolerance;

        min = 3.0; max = 4.0; expected = FastMath.PI;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        assertEquals(expected, result, tolerance);

        min = -1.0; max = 1.5; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        assertEquals(expected, result, tolerance);
    }

    /**
     * Test of solver for the quintic function.
     */
    @Test
    void testQuinticFunction() {
        UnivariateFunction f = new QuinticFunction();
        UnivariateSolver solver = new MullerSolver();
        double min, max, expected, result, tolerance;

        min = -0.4; max = 0.2; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        assertEquals(expected, result, tolerance);

        min = 0.75; max = 1.5; expected = 1.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        assertEquals(expected, result, tolerance);

        min = -0.9; max = -0.2; expected = -0.5;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        assertEquals(expected, result, tolerance);
    }

    /**
     * Test of solver for the exponential function.
     * <p>
     * It takes 10 to 15 iterations for the last two tests to converge.
     * In fact, if not for the bisection alternative, the solver would
     * exceed the default maximal iteration of 100.
     */
    @Test
    void testExpm1Function() {
        UnivariateFunction f = new Expm1();
        UnivariateSolver solver = new MullerSolver();
        double min, max, expected, result, tolerance;

        min = -1.0; max = 2.0; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        assertEquals(expected, result, tolerance);

        min = -20.0; max = 10.0; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        assertEquals(expected, result, tolerance);

        min = -50.0; max = 100.0; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        assertEquals(expected, result, tolerance);
    }

    /**
     * Test of parameters for the solver.
     */
    @Test
    void testParameters() {
        UnivariateFunction f = new Sin();
        UnivariateSolver solver = new MullerSolver();

        try {
            // bad interval
            double root = solver.solve(100, f, 1, -1);
            System.out.println("root=" + root);
            fail("Expecting MathIllegalArgumentException - bad interval");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            // no bracketing
            solver.solve(100, f, 2, 3);
            fail("Expecting MathIllegalArgumentException - no bracketing");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }
}
