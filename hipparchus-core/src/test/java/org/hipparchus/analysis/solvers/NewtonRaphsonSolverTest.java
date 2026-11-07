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
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.analysis.differentiation.UnivariateDerivative1;
import org.hipparchus.analysis.differentiation.UnivariateDifferentiableFunction;
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.analysis.polynomials.PolynomialFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;


/**
 */
final class NewtonRaphsonSolverTest {

    @Test
    void testConstructor() {
        // GIVEN
        final double expectedRelativeAccuracy = 1.;
        final double expectedAbsoluteAccuracy = 2.;
        // WHEN
        final NewtonRaphsonSolver solver = new NewtonRaphsonSolver(expectedRelativeAccuracy, expectedAbsoluteAccuracy, 3.);
        // THEN
        Assertions.assertEquals(expectedAbsoluteAccuracy, solver.getAbsoluteAccuracy());
        Assertions.assertEquals(expectedRelativeAccuracy, solver.getRelativeAccuracy());
    }

    @Test
    @Deprecated
    void testDeprecated() {
        // GIVEN
        final NewtonRaphsonSolver solver = new NewtonRaphsonSolver();
        final Sin sin = new Sin();
        solver.setup(10, sin, -1., 1., 0.);
        // WHEN
        final DerivativeStructure derivativeStructure = solver.computeObjectiveValueAndDerivative(0.);
        // THEN
        final UnivariateDerivative1 derivative1 = solver.computeObjectiveValueAndDerivative1(0.);
        Assertions.assertEquals(derivative1.getDerivative(1), derivativeStructure.getPartialDerivative(1));
    }

    /**
     *
     */
    @Test
    void testSinZero() {
        UnivariateDifferentiableFunction f = new Sin();
        double result;

        NewtonRaphsonSolver solver = new NewtonRaphsonSolver();
        result = solver.solve(100, f, 3, 4);
        assertEquals(FastMath.PI, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 1, 4);
        assertEquals(FastMath.PI, result, solver.getAbsoluteAccuracy());

        assertTrue(solver.getEvaluations() > 0);
    }

    /**
     *
     */
    @Test
    void testQuinticZero() {
        final UnivariateDifferentiableFunction f = new QuinticFunction();
        double result;

        NewtonRaphsonSolver solver = new NewtonRaphsonSolver();
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
    }

    @Test
    @SuppressWarnings("deprecation")
    void testDeprecatedOverrideStillCalled() {
        // a 4.0.3 subclass overriding the deprecated hook must still be used by the solver
        final int[] calls = new int[1];
        final NewtonRaphsonSolver solver = new NewtonRaphsonSolver() {
            @Override
            protected DerivativeStructure computeObjectiveValueAndDerivative(final double point) {
                ++calls[0];
                return super.computeObjectiveValueAndDerivative(point);
            }
        };
        checkSqrt2(solver, calls);
    }

    @Test
    void testNewOverrideCalled() {
        final int[] calls = new int[1];
        final NewtonRaphsonSolver solver = new NewtonRaphsonSolver() {
            @Override
            protected UnivariateDerivative1 computeObjectiveValueAndDerivative1(final double point) {
                ++calls[0];
                return super.computeObjectiveValueAndDerivative1(point);
            }
        };
        checkSqrt2(solver, calls);
    }

    /** Solve x^2 - 2 = 0 and check the hook was called once per evaluation, with the 4.0.3 count. */
    private static void checkSqrt2(final NewtonRaphsonSolver solver, final int[] calls) {
        final double root = solver.solve(100, new PolynomialFunction(new double[] { -2, 0, 1 }), 0, 3);
        Assertions.assertEquals(FastMath.sqrt(2), root, 1.0e-6);
        Assertions.assertEquals(4, calls[0]);
        Assertions.assertEquals(4, solver.getEvaluations());
    }

}
