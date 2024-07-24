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

import org.hipparchus.analysis.FunctionUtils;
import org.hipparchus.analysis.MonitoredFunction;
import org.hipparchus.analysis.QuinticFunction;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.UnivariateDifferentiableFunction;
import org.hipparchus.analysis.function.Constant;
import org.hipparchus.analysis.function.Inverse;
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.analysis.function.Sqrt;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test case for {@link BrentSolver Brent} solver.
 * Because Brent-Dekker is guaranteed to converge in less than the default
 * maximum iteration count due to bisection fallback, it is quite hard to
 * debug. I include measured iteration counts plus one in order to detect
 * regressions. On average Brent-Dekker should use 4..5 iterations for the
 * default absolute accuracy of 10E-8 for sinus and the quintic function around
 * zero, and 5..10 iterations for the other zeros.
 *
 */
final class BrentSolverTest {
    @Test
    void testSinZero() {
        // The sinus function is behaved well around the root at pi. The second
        // order derivative is zero, which means linar approximating methods will
        // still converge quadratically.
        UnivariateFunction f = new Sin();
        double result;
        UnivariateSolver solver = new BrentSolver();
        // Somewhat benign interval. The function is monotone.
        result = solver.solve(100, f, 3, 4);
        // System.out.println(
        //    "Root: " + result + " Evaluations: " + solver.getEvaluations());
        assertEquals(FastMath.PI, result, solver.getAbsoluteAccuracy());
        assertTrue(solver.getEvaluations() <= 7);
        // Larger and somewhat less benign interval. The function is grows first.
        result = solver.solve(100, f, 1, 4);
        // System.out.println(
        //    "Root: " + result + " Evaluations: " + solver.getEvaluations());
        assertEquals(FastMath.PI, result, solver.getAbsoluteAccuracy());
        assertTrue(solver.getEvaluations() <= 8);
    }

    @Test
    void testQuinticZero() {
        // The quintic function has zeros at 0, +-0.5 and +-1.
        // Around the root of 0 the function is well behaved, with a second derivative
        // of zero a 0.
        // The other roots are less well to find, in particular the root at 1, because
        // the function grows fast for x>1.
        // The function has extrema (first derivative is zero) at 0.27195613 and 0.82221643,
        // intervals containing these values are harder for the solvers.
        UnivariateFunction f = new QuinticFunction();
        double result;
        // Brent-Dekker solver.
        UnivariateSolver solver = new BrentSolver();
        // Symmetric bracket around 0. Test whether solvers can handle hitting
        // the root in the first iteration.
        result = solver.solve(100, f, -0.2, 0.2);
        //System.out.println(
        //    "Root: " + result + " Evaluations: " + solver.getEvaluations());
        assertEquals(0, result, solver.getAbsoluteAccuracy());
        assertTrue(solver.getEvaluations() <= 3);
        // 1 iterations on i586 JDK 1.4.1.
        // Asymmetric bracket around 0, just for fun. Contains extremum.
        result = solver.solve(100, f, -0.1, 0.3);
        //System.out.println(
        //    "Root: " + result + " Evaluations: " + solver.getEvaluations());
        assertEquals(0, result, solver.getAbsoluteAccuracy());
        // 5 iterations on i586 JDK 1.4.1.
        assertTrue(solver.getEvaluations() <= 7);
        // Large bracket around 0. Contains two extrema.
        result = solver.solve(100, f, -0.3, 0.45);
        //System.out.println(
        //    "Root: " + result + " Evaluations: " + solver.getEvaluations());
        assertEquals(0, result, solver.getAbsoluteAccuracy());
        // 6 iterations on i586 JDK 1.4.1.
        assertTrue(solver.getEvaluations() <= 8);
        // Benign bracket around 0.5, function is monotonous.
        result = solver.solve(100, f, 0.3, 0.7);
        //System.out.println(
        //    "Root: " + result + " Evaluations: " + solver.getEvaluations());
        assertEquals(0.5, result, solver.getAbsoluteAccuracy());
        // 6 iterations on i586 JDK 1.4.1.
        assertTrue(solver.getEvaluations() <= 9);
        // Less benign bracket around 0.5, contains one extremum.
        result = solver.solve(100, f, 0.2, 0.6);
        // System.out.println(
        //    "Root: " + result + " Evaluations: " + solver.getEvaluations());
        assertEquals(0.5, result, solver.getAbsoluteAccuracy());
        assertTrue(solver.getEvaluations() <= 10);
        // Large, less benign bracket around 0.5, contains both extrema.
        result = solver.solve(100, f, 0.05, 0.95);
        //System.out.println(
        //    "Root: " + result + " Evaluations: " + solver.getEvaluations());
        assertEquals(0.5, result, solver.getAbsoluteAccuracy());
        assertTrue(solver.getEvaluations() <= 11);
        // Relatively benign bracket around 1, function is monotonous. Fast growth for x>1
        // is still a problem.
        result = solver.solve(100, f, 0.85, 1.25);
        //System.out.println(
        //    "Root: " + result + " Evaluations: " + solver.getEvaluations());
        assertEquals(1.0, result, solver.getAbsoluteAccuracy());
        assertTrue(solver.getEvaluations() <= 11);
        // Less benign bracket around 1 with extremum.
        result = solver.solve(100, f, 0.8, 1.2);
        //System.out.println(
        //    "Root: " + result + " Evaluations: " + solver.getEvaluations());
        assertEquals(1.0, result, solver.getAbsoluteAccuracy());
        assertTrue(solver.getEvaluations() <= 11);
        // Large bracket around 1. Monotonous.
        result = solver.solve(100, f, 0.85, 1.75);
        //System.out.println(
        //    "Root: " + result + " Evaluations: " + solver.getEvaluations());
        assertEquals(1.0, result, solver.getAbsoluteAccuracy());
        assertTrue(solver.getEvaluations() <= 13);
        // Large bracket around 1. Interval contains extremum.
        result = solver.solve(100, f, 0.55, 1.45);
        //System.out.println(
        //    "Root: " + result + " Evaluations: " + solver.getEvaluations());
        assertEquals(1.0, result, solver.getAbsoluteAccuracy());
        assertTrue(solver.getEvaluations() <= 10);
        // Very large bracket around 1 for testing fast growth behaviour.
        result = solver.solve(100, f, 0.85, 5);
        //System.out.println(
       //     "Root: " + result + " Evaluations: " + solver.getEvaluations());
        assertEquals(1.0, result, solver.getAbsoluteAccuracy());
        assertTrue(solver.getEvaluations() <= 15);

        try {
            result = solver.solve(5, f, 0.85, 5);
            fail("Expected MathIllegalStateException");
        } catch (MathIllegalStateException e) {
            // Expected.
        }
    }

    @Test
    void testRootEndpoints() {
        UnivariateFunction f = new Sin();
        BrentSolver solver = new BrentSolver();

        // endpoint is root
        double result = solver.solve(100, f, FastMath.PI, 4);
        assertEquals(FastMath.PI, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 3, FastMath.PI);
        assertEquals(FastMath.PI, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, FastMath.PI, 4, 3.5);
        assertEquals(FastMath.PI, result, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 3, FastMath.PI, 3.07);
        assertEquals(FastMath.PI, result, solver.getAbsoluteAccuracy());
    }

    @Test
    void testBadEndpoints() {
        UnivariateFunction f = new Sin();
        BrentSolver solver = new BrentSolver();
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
    void testInitialGuess() {
        MonitoredFunction f = new MonitoredFunction(new QuinticFunction());
        BrentSolver solver = new BrentSolver();
        double result;

        // no guess
        result = solver.solve(100, f, 0.6, 7.0);
        assertEquals(1.0, result, solver.getAbsoluteAccuracy());
        int referenceCallsCount = f.getCallsCount();
        assertTrue(referenceCallsCount >= 13);

        // invalid guess (it *is* a root, but outside of the range)
        try {
          result = solver.solve(100, f, 0.6, 7.0, 0.0);
          fail("a MathIllegalArgumentException was expected");
        } catch (MathIllegalArgumentException iae) {
            // expected behaviour
        }

        // bad guess
        f.setCallsCount(0);
        result = solver.solve(100, f, 0.6, 7.0, 0.61);
        assertEquals(1.0, result, solver.getAbsoluteAccuracy());
        assertTrue(f.getCallsCount() > referenceCallsCount);

        // good guess
        f.setCallsCount(0);
        result = solver.solve(100, f, 0.6, 7.0, 0.999999);
        assertEquals(1.0, result, solver.getAbsoluteAccuracy());
        assertTrue(f.getCallsCount() < referenceCallsCount);

        // perfect guess
        f.setCallsCount(0);
        result = solver.solve(100, f, 0.6, 7.0, 1.0);
        assertEquals(1.0, result, solver.getAbsoluteAccuracy());
        assertEquals(1, solver.getEvaluations());
        assertEquals(1, f.getCallsCount());
    }

    @Test
    void testMath832() {
        final UnivariateFunction f = new UnivariateFunction() {
                private final UnivariateDifferentiableFunction sqrt = new Sqrt();
                private final UnivariateDifferentiableFunction inv = new Inverse();
                private final UnivariateDifferentiableFunction func
                    = FunctionUtils.add(FunctionUtils.multiply(new Constant(1e2), sqrt),
                                        FunctionUtils.multiply(new Constant(1e6), inv),
                                        FunctionUtils.multiply(new Constant(1e4),
                                                               FunctionUtils.compose(inv, sqrt)));
                private final DSFactory factory = new DSFactory(1, 1);

                public double value(double x) {
                    return func.value(factory.variable(0, x)).getPartialDerivative(1);
                }

            };

        BrentSolver solver = new BrentSolver();
        final double result = solver.solve(99, f, 1, 1e30, 1 + 1e-10);
        assertEquals(804.93558250, result, 1e-8);
    }
}
