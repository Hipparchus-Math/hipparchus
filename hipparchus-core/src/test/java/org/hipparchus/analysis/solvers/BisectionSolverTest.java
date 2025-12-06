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
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 */
final class BisectionSolverTest {
    @Test
    void testSinZero() {
        UnivariateFunction f = new Sin();
        double result;

        BisectionSolver solver = new BisectionSolver();
        result = solver.solve(100, f, 3, 4);
        assertEquals(FastMath.PI, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 1, 4);
        assertEquals(FastMath.PI, result, solver.getAbsoluteAccuracy());
    }

    @Test
    void testQuinticZero() {
        UnivariateFunction f = new QuinticFunction();
        double result;

        BisectionSolver solver = new BisectionSolver();
        result = solver.solve(100, f, -0.2, 0.2);
        assertEquals(0, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, -0.1, 0.3);
        assertEquals(0, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, -0.3, 0.45);
        assertEquals(0, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.3, 0.7);
        assertEquals(0.5, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.2, 0.6);
        assertEquals(0.5, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.05, 0.95);
        assertEquals(0.5, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.85, 1.25);
        assertEquals(1.0, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.8, 1.2);
        assertEquals(1.0, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.85, 1.75);
        assertEquals(1.0, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.55, 1.45);
        assertEquals(1.0, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.85, 5);
        assertEquals(1.0, result, solver.getAbsoluteAccuracy());

        assertTrue(solver.getEvaluations() > 0);
    }

    @Test
    void testMath369() {
        UnivariateFunction f = new Sin();
        BisectionSolver solver = new BisectionSolver();
        assertEquals(FastMath.PI, solver.solve(100, f, 3.0, 3.2, 3.1), solver.getAbsoluteAccuracy());
    }

    @Test
    void testHipparchusGithub40() {
        assertThrows(MathIllegalArgumentException.class, () -> new BisectionSolver().solve(100, x -> Math.cos(x) + 2, 0.0, 5.0));
    }

}

