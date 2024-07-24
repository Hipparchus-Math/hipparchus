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
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Test case for {@link RiddersSolver Ridders} solver.
 * <p>
 * Ridders' method converges superlinearly, more specific, its rate of
 * convergence is sqrt(2). Test runs show that for a default absolute
 * accuracy of 1E-6, it generally takes less than 5 iterations for close
 * initial bracket and 5 to 10 iterations for distant initial bracket
 * to converge.
 *
 */
public final class RiddersSolverTest {
    /**
     * Test of solver for the sine function.
     */
    @Test
    public void testSinFunction() {
        UnivariateFunction f = new Sin();
        UnivariateSolver solver = new RiddersSolver();
        double min, max, expected, result, tolerance;

        min = 3.0; max = 4.0; expected = FastMath.PI;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assertions.assertEquals(expected, result, tolerance);

        min = -1.0; max = 1.5; expected = 0.0;
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
        UnivariateFunction f = new QuinticFunction();
        UnivariateSolver solver = new RiddersSolver();
        double min, max, expected, result, tolerance;

        min = -0.4; max = 0.2; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assertions.assertEquals(expected, result, tolerance);

        min = 0.75; max = 1.5; expected = 1.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assertions.assertEquals(expected, result, tolerance);

        min = -0.9; max = -0.2; expected = -0.5;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assertions.assertEquals(expected, result, tolerance);
    }

    /**
     * Test of solver for the exponential function.
     */
    @Test
    public void testExpm1Function() {
        UnivariateFunction f = new Expm1();
        UnivariateSolver solver = new RiddersSolver();
        double min, max, expected, result, tolerance;

        min = -1.0; max = 2.0; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assertions.assertEquals(expected, result, tolerance);

        min = -20.0; max = 10.0; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assertions.assertEquals(expected, result, tolerance);

        min = -50.0; max = 100.0; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assertions.assertEquals(expected, result, tolerance);
    }

    /**
     * Test of parameters for the solver.
     */
    @Test
    public void testParameters() {
        UnivariateFunction f = new Sin();
        UnivariateSolver solver = new RiddersSolver();

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
}
