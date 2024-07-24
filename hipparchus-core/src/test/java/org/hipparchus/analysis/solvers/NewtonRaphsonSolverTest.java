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
import org.hipparchus.analysis.differentiation.UnivariateDifferentiableFunction;
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;


/**
 */
public final class NewtonRaphsonSolverTest {
    /**
     *
     */
    @Test
    public void testSinZero() {
        UnivariateDifferentiableFunction f = new Sin();
        double result;

        NewtonRaphsonSolver solver = new NewtonRaphsonSolver();
        result = solver.solve(100, f, 3, 4);
        Assertions.assertEquals(FastMath.PI, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 1, 4);
        Assertions.assertEquals(FastMath.PI, result, solver.getAbsoluteAccuracy());

        Assertions.assertTrue(solver.getEvaluations() > 0);
    }

    /**
     *
     */
    @Test
    public void testQuinticZero() {
        final UnivariateDifferentiableFunction f = new QuinticFunction();
        double result;

        NewtonRaphsonSolver solver = new NewtonRaphsonSolver();
        result = solver.solve(100, f, -0.2, 0.2);
        Assertions.assertEquals(0, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, -0.1, 0.3);
        Assertions.assertEquals(0, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, -0.3, 0.45);
        Assertions.assertEquals(0, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.3, 0.7);
        Assertions.assertEquals(0.5, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.2, 0.6);
        Assertions.assertEquals(0.5, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.05, 0.95);
        Assertions.assertEquals(0.5, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.85, 1.25);
        Assertions.assertEquals(1.0, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.8, 1.2);
        Assertions.assertEquals(1.0, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.85, 1.75);
        Assertions.assertEquals(1.0, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.55, 1.45);
        Assertions.assertEquals(1.0, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 0.85, 5);
        Assertions.assertEquals(1.0, result, solver.getAbsoluteAccuracy());
    }
}
