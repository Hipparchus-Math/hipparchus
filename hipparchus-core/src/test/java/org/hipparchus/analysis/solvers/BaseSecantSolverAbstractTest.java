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
import org.hipparchus.analysis.XMinus5Function;
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import java.lang.reflect.Field;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Base class for root-finding algorithms tests derived from
 * {@link BaseSecantSolver}.
 *
 */
public abstract class BaseSecantSolverAbstractTest {
    /** Returns the solver to use to perform the tests.
     * @return the solver to use to perform the tests
     */
    protected abstract UnivariateSolver getSolver();

    /** Returns the expected number of evaluations for the
     * {@link #testQuinticZero} unit test. A value of {@code -1} indicates that
     * the test should be skipped for that solver.
     * @return the expected number of evaluations for the
     * {@link #testQuinticZero} unit test
     */
    protected abstract int[] getQuinticEvalCounts();

    @Test
    public void testSinZero() {
        // The sinus function is behaved well around the root at pi. The second
        // order derivative is zero, which means linear approximating methods
        // still converge quadratically.
        UnivariateFunction f = new Sin();
        double result;
        UnivariateSolver solver = getSolver();

        result = solver.solve(100, f, 3, 4);
        //System.out.println(
        //    "Root: " + result + " Evaluations: " + solver.getEvaluations());
        assertEquals(FastMath.PI, result, solver.getAbsoluteAccuracy());
        assertTrue(solver.getEvaluations() <= 6);
        result = solver.solve(100, f, 1, 4);
        //System.out.println(
        //    "Root: " + result + " Evaluations: " + solver.getEvaluations());
        assertEquals(FastMath.PI, result, solver.getAbsoluteAccuracy());
        assertTrue(solver.getEvaluations() <= 7);
    }

    @Test
    public void testQuinticZero() {
        // The quintic function has zeros at 0, +-0.5 and +-1.
        // Around the root of 0 the function is well behaved, with a second
        // derivative of zero a 0.
        // The other roots are less well to find, in particular the root at 1,
        // because the function grows fast for x>1.
        // The function has extrema (first derivative is zero) at 0.27195613
        // and 0.82221643, intervals containing these values are harder for
        // the solvers.
        UnivariateFunction f = new QuinticFunction();
        double result;
        UnivariateSolver solver = getSolver();
        double atol = solver.getAbsoluteAccuracy();
        int[] counts = getQuinticEvalCounts();

        // Tests data: initial bounds, and expected solution, per test case.
        double[][] testsData = {{-0.2,  0.2,  0.0},
                                {-0.1,  0.3,  0.0},
                                {-0.3,  0.45, 0.0},
                                { 0.3,  0.7,  0.5},
                                { 0.2,  0.6,  0.5},
                                { 0.05, 0.95, 0.5},
                                { 0.85, 1.25, 1.0},
                                { 0.8,  1.2,  1.0},
                                { 0.85, 1.75, 1.0},
                                { 0.55, 1.45, 1.0},
                                { 0.85, 5.0,  1.0},
                               };
        int maxIter = 500;

        for(int i = 0; i < testsData.length; i++) {
            // Skip test, if needed.
            if (counts[i] == -1) continue;

            // Compute solution.
            double[] testData = testsData[i];
            result = solver.solve(maxIter, f, testData[0], testData[1]);
            //System.out.println(
            //    "Root: " + result + " Evaluations: " + solver.getEvaluations());

            // Check solution.
            assertEquals(result, testData[2], atol);
            assertTrue(solver.getEvaluations() <= counts[i] + 1,
                    "" + solver.getEvaluations() + " <= " + (counts[i] + 1));
        }
    }

    @Test
    public void testRootEndpoints() {
        UnivariateFunction f = new XMinus5Function();
        UnivariateSolver solver = getSolver();

        // End-point is root. This should be a special case in the solver, and
        // the initial end-point should be returned exactly.
        double result = solver.solve(100, f, 5.0, 6.0);
        assertEquals(5.0, result, 0.0);

        result = solver.solve(100, f, 4.0, 5.0);
        assertEquals(5.0, result, 0.0);

        result = solver.solve(100, f, 5.0, 6.0, 5.5);
        assertEquals(5.0, result, 0.0);

        result = solver.solve(100, f, 4.0, 5.0, 4.5);
        assertEquals(5.0, result, 0.0);
    }

    @Test
    public void testCloseEndpoints() {
        UnivariateFunction f = new XMinus5Function();
        UnivariateSolver solver = getSolver();

        double result = solver.solve(100, f, 5.0, FastMath.nextUp(5.0));
        assertEquals(5.0, result, 0.0);

        result = solver.solve(100, f, FastMath.nextDown(5.0), 5.0);
        assertEquals(5.0, result, 0.0);
    }

    @Test
    public void testBadEndpoints() {
        UnivariateFunction f = new Sin();
        UnivariateSolver solver = getSolver();
        try {  // bad interval
            solver.solve(100, f, 1, -1);
            fail("Expecting MathIllegalArgumentException - bad interval");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {  // no bracket
            solver.solve(100, f, 1, 1.5);
            fail("Expecting MathIllegalArgumentException - non-bracketing");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {  // no bracket
            solver.solve(100, f, 1, 1.5, 1.2);
            fail("Expecting MathIllegalArgumentException - non-bracketing");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testSolutionLeftSide() {
        UnivariateFunction f = new Sin();
        UnivariateSolver solver = getSolver();
        double left = -1.5;
        double right = 0.05;
        for(int i = 0; i < 10; i++) {
            // Test whether the allowed solutions are taken into account.
            double solution = getSolution(solver, 100, f, left, right, AllowedSolution.LEFT_SIDE);
            if (!Double.isNaN(solution)) {
                assertTrue(solution <= 0.0);
            }

            // Prepare for next test.
            left -= 0.1;
            right += 0.3;
        }
    }

    @Test
    public void testSolutionRightSide() {
        UnivariateFunction f = new Sin();
        UnivariateSolver solver = getSolver();
        double left = -1.5;
        double right = 0.05;
        for(int i = 0; i < 10; i++) {
            // Test whether the allowed solutions are taken into account.
            double solution = getSolution(solver, 100, f, left, right, AllowedSolution.RIGHT_SIDE);
            if (!Double.isNaN(solution)) {
                assertTrue(solution >= 0.0);
            }

            // Prepare for next test.
            left -= 0.1;
            right += 0.3;
        }
    }
    @Test
    public void testSolutionBelowSide() {
        UnivariateFunction f = new Sin();
        UnivariateSolver solver = getSolver();
        double left = -1.5;
        double right = 0.05;
        for(int i = 0; i < 10; i++) {
            // Test whether the allowed solutions are taken into account.
            double solution = getSolution(solver, 100, f, left, right, AllowedSolution.BELOW_SIDE);
            if (!Double.isNaN(solution)) {
                assertTrue(f.value(solution) <= 0.0);
            }

            // Prepare for next test.
            left -= 0.1;
            right += 0.3;
        }
    }

    @Test
    public void testSolutionAboveSide() {
        UnivariateFunction f = new Sin();
        UnivariateSolver solver = getSolver();
        double left = -1.5;
        double right = 0.05;
        for(int i = 0; i < 10; i++) {
            // Test whether the allowed solutions are taken into account.
            double solution = getSolution(solver, 100, f, left, right, AllowedSolution.ABOVE_SIDE);
            if (!Double.isNaN(solution)) {
                assertTrue(f.value(solution) >= 0.0);
            }

            // Prepare for next test.
            left -= 0.1;
            right += 0.3;
        }
    }

    private double getSolution(UnivariateSolver solver, int maxEval, UnivariateFunction f,
                               double left, double right, AllowedSolution allowedSolution) {
        try {
            @SuppressWarnings("unchecked")
            BracketedUnivariateSolver<UnivariateFunction> bracketing =
            (BracketedUnivariateSolver<UnivariateFunction>) solver;
            return bracketing.solve(100, f, left, right, allowedSolution);
        } catch (ClassCastException cce) {
            double baseRoot = solver.solve(maxEval, f, left, right);
            if ((baseRoot <= left) || (baseRoot >= right)) {
                // the solution slipped out of interval
                return Double.NaN;
            }
            PegasusSolver bracketing =
                    new PegasusSolver(solver.getRelativeAccuracy(), solver.getAbsoluteAccuracy(),
                                      solver.getFunctionValueAccuracy());
            return UnivariateSolverUtils.forceSide(maxEval - solver.getEvaluations(),
                                                       f, bracketing, baseRoot, left, right,
                                                       allowedSolution);
        }
    }

    protected void checktype(UnivariateSolver solver, BaseSecantSolver.Method expected) {
        try {
            Field methodField = BaseSecantSolver.class.getDeclaredField("method");
            methodField.setAccessible(true);
            BaseSecantSolver.Method method = (BaseSecantSolver.Method) methodField.get(solver);
            assertEquals(expected, method);
        } catch (IllegalAccessException | NoSuchFieldException | SecurityException e) {
            fail(e.getLocalizedMessage());
        }
    }

}
