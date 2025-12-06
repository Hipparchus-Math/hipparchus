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

import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.analysis.QuinticFunction;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 */
class UnivariateSolverUtilsTest {

    private UnivariateFunction sin = new Sin();
    private CalculusFieldUnivariateFunction<Binary64> fieldSin = Binary64::sin;

    @Test
    void testSolveNull() {
        assertThrows(NullArgumentException.class, () -> UnivariateSolverUtils.solve(null, 0.0, 4.0));
    }

    @Test
    void testSolveBadEndpoints() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            double root = UnivariateSolverUtils.solve(sin, 4.0, -0.1, 1e-6);
            System.out.println("root=" + root);
        });
    }

    @Test
    void testSolveBadAccuracy() {
        try { // bad accuracy
            UnivariateSolverUtils.solve(sin, 0.0, 4.0, 0.0);
//             Assertions.fail("Expecting MathIllegalArgumentException"); // TODO needs rework since convergence behaviour was changed
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testSolveSin() {
        double x = UnivariateSolverUtils.solve(sin, 1.0, 4.0);
        assertEquals(FastMath.PI, x, 1.0e-4);
    }

    @Test
    void testSolveAccuracyNull()  {
        assertThrows(NullArgumentException.class, () -> {
            double accuracy = 1.0e-6;
            UnivariateSolverUtils.solve(null, 0.0, 4.0, accuracy);
        });
    }

    @Test
    void testSolveAccuracySin() {
        double accuracy = 1.0e-6;
        double x = UnivariateSolverUtils.solve(sin, 1.0,
                4.0, accuracy);
        assertEquals(FastMath.PI, x, accuracy);
    }

    @Test
    void testSolveNoRoot() {
        assertThrows(MathIllegalArgumentException.class, () -> UnivariateSolverUtils.solve(sin, 1.0, 1.5));
    }

    @Test
    void testBracketSin() {
        double[] result = UnivariateSolverUtils.bracket(sin,
                0.0, -2.0, 2.0);
        assertTrue(sin.value(result[0]) < 0);
        assertTrue(sin.value(result[1]) > 0);
    }

    @Test
    void testBracketCentered() {
        double initial = 0.1;
        double[] result = UnivariateSolverUtils.bracket(sin, initial, -2.0, 2.0, 0.2, 1.0, 100);
        assertTrue(result[0] < initial);
        assertTrue(result[1] > initial);
        assertTrue(sin.value(result[0]) < 0);
        assertTrue(sin.value(result[1]) > 0);
    }

    @Test
    void testBracketLow() {
        double initial = 0.5;
        double[] result = UnivariateSolverUtils.bracket(sin, initial, -2.0, 2.0, 0.2, 1.0, 100);
        assertTrue(result[0] < initial);
        assertTrue(result[1] < initial);
        assertTrue(sin.value(result[0]) < 0);
        assertTrue(sin.value(result[1]) > 0);
    }

    @Test
    void testBracketHigh(){
        double initial = -0.5;
        double[] result = UnivariateSolverUtils.bracket(sin, initial, -2.0, 2.0, 0.2, 1.0, 100);
        assertTrue(result[0] > initial);
        assertTrue(result[1] > initial);
        assertTrue(sin.value(result[0]) < 0);
        assertTrue(sin.value(result[1]) > 0);
    }

    @Test
    void testBracketLinear(){
        assertThrows(MathIllegalArgumentException.class, () -> UnivariateSolverUtils.bracket(x -> 1 - x, 1000, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 100));
    }

    @Test
    void testBracketExponential(){
        double[] result = UnivariateSolverUtils.bracket(x -> 1 - x, 1000, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0, 2.0, 10);
        assertTrue(result[0] <= 1);
        assertTrue(result[1] >= 1);
    }

    @Test
    void testBracketEndpointRoot() {
        double[] result = UnivariateSolverUtils.bracket(sin, 1.5, 0, 2.0, 100);
        assertEquals(0.0, sin.value(result[0]), 1.0e-15);
        assertTrue(sin.value(result[1]) > 0);
    }

    @Test
    void testNullFunction() {
        assertThrows(NullArgumentException.class, () -> UnivariateSolverUtils.bracket(null, 1.5, 0, 2.0));
    }

    @Test
    void testBadInitial() {
        assertThrows(MathIllegalArgumentException.class, () -> UnivariateSolverUtils.bracket(sin, 2.5, 0, 2.0));
    }

    @Test
    void testBadAdditive() {
        assertThrows(MathIllegalArgumentException.class, () -> UnivariateSolverUtils.bracket(sin, 1.0, -2.0, 3.0, -1.0, 1.0, 100));
    }

    @Test
    void testIterationExceeded() {
        assertThrows(MathIllegalArgumentException.class, () -> UnivariateSolverUtils.bracket(sin, 1.0, -2.0, 3.0, 1.0e-5, 1.0, 100));
    }

    @Test
    void testBadEndpoints() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            // endpoints not valid
            UnivariateSolverUtils.bracket(sin, 1.5, 2.0, 1.0);
        });
    }

    @Test
    void testBadMaximumIterations() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            // bad maximum iterations
            UnivariateSolverUtils.bracket(sin, 1.5, 0, 2.0, 0);
        });
    }

    @Test
    void testFieldBracketSin() {
        Binary64[] result = UnivariateSolverUtils.bracket(fieldSin,new Binary64(0.0),
                                                           new Binary64(-2.0),new Binary64(2.0));
        assertTrue(fieldSin.value(result[0]).getReal() < 0);
        assertTrue(fieldSin.value(result[1]).getReal() > 0);
    }

    @Test
    void testFieldBracketCentered() {
        Binary64 initial = new Binary64(0.1);
        Binary64[] result = UnivariateSolverUtils.bracket(fieldSin, initial,
                                                           new Binary64(-2.0), new Binary64(2.0),
                                                           new Binary64(0.2), new Binary64(1.0),
                                                           100);
        assertTrue(result[0].getReal() < initial.getReal());
        assertTrue(result[1].getReal() > initial.getReal());
        assertTrue(fieldSin.value(result[0]).getReal() < 0);
        assertTrue(fieldSin.value(result[1]).getReal() > 0);
    }

    @Test
    void testFieldBracketLow() {
        Binary64 initial = new Binary64(0.5);
        Binary64[] result = UnivariateSolverUtils.bracket(fieldSin, initial,
                                                           new Binary64(-2.0), new Binary64(2.0),
                                                           new Binary64(0.2), new Binary64(1.0),
                                                           100);
        assertTrue(result[0].getReal() < initial.getReal());
        assertTrue(result[1].getReal() < initial.getReal());
        assertTrue(fieldSin.value(result[0]).getReal() < 0);
        assertTrue(fieldSin.value(result[1]).getReal() > 0);
    }

    @Test
    void testFieldBracketHigh(){
        Binary64 initial = new Binary64(-0.5);
        Binary64[] result = UnivariateSolverUtils.bracket(fieldSin, initial,
                                                           new Binary64(-2.0), new Binary64(2.0),
                                                           new Binary64(0.2), new Binary64(1.0),
                                                           100);
        assertTrue(result[0].getReal() > initial.getReal());
        assertTrue(result[1].getReal() > initial.getReal());
        assertTrue(fieldSin.value(result[0]).getReal() < 0);
        assertTrue(fieldSin.value(result[1]).getReal() > 0);
    }

    @Test
    void testFieldBracketLinear(){
        assertThrows(MathIllegalArgumentException.class, () -> UnivariateSolverUtils.bracket(x -> x.negate().add(1),
            new Binary64(1000),
            new Binary64(Double.NEGATIVE_INFINITY), new Binary64(Double.POSITIVE_INFINITY),
            new Binary64(1.0), new Binary64(1.0), 100));
    }

    @Test
    void testFieldBracketExponential(){
        Binary64[] result = UnivariateSolverUtils.bracket(x -> x.negate().add(1),
        new Binary64(1000),
        new Binary64(Double.NEGATIVE_INFINITY), new Binary64(Double.POSITIVE_INFINITY),
        new Binary64(1.0), new Binary64(2.0), 10);
        assertTrue(result[0].getReal() <= 1);
        assertTrue(result[1].getReal() >= 1);
    }

    @Test
    void testFieldBracketEndpointRoot() {
        Binary64[] result = UnivariateSolverUtils.bracket(fieldSin,
                                                           new Binary64(1.5), new Binary64(0),
                                                           new Binary64(2.0), 100);
        assertEquals(0.0, fieldSin.value(result[0]).getReal(), 1.0e-15);
        assertTrue(fieldSin.value(result[1]).getReal() > 0);
    }

    @Test
    void testFieldNullFunction() {
        assertThrows(NullArgumentException.class, () -> UnivariateSolverUtils.bracket(null, new Binary64(1.5), new Binary64(0), new Binary64(2.0)));
    }

    @Test
    void testFieldBadInitial() {
        assertThrows(MathIllegalArgumentException.class, () -> UnivariateSolverUtils.bracket(fieldSin, new Binary64(2.5), new Binary64(0), new Binary64(2.0)));
    }

    @Test
    void testFieldBadAdditive() {
        assertThrows(MathIllegalArgumentException.class, () -> UnivariateSolverUtils.bracket(fieldSin, new Binary64(1.0), new Binary64(-2.0), new Binary64(3.0),
            new Binary64(-1.0), new Binary64(1.0), 100));
    }

    @Test
    void testFieldIterationExceeded() {
        assertThrows(MathIllegalArgumentException.class, () -> UnivariateSolverUtils.bracket(fieldSin, new Binary64(1.0), new Binary64(-2.0), new Binary64(3.0),
            new Binary64(1.0e-5), new Binary64(1.0), 100));
    }

    @Test
    void testFieldBadEndpoints() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            // endpoints not valid
            UnivariateSolverUtils.bracket(fieldSin, new Binary64(1.5), new Binary64(2.0), new Binary64(1.0));
        });
    }

    @Test
    void testFieldBadMaximumIterations() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            // bad maximum iterations
            UnivariateSolverUtils.bracket(fieldSin, new Binary64(1.5), new Binary64(0), new Binary64(2.0), 0);
        });
    }

    /** check the search continues when a = lowerBound and b &lt; upperBound. */
    @Test
    void testBracketLoopConditionForB() {
        double[] result = UnivariateSolverUtils.bracket(sin, -0.9, -1, 1, 0.1, 1, 100);
        assertTrue(result[0] <= 0);
        assertTrue(result[1] >= 0);
    }

    @Test
    void testMisc() {
        UnivariateFunction f = new QuinticFunction();
        double result;
        // Static solve method
        result = UnivariateSolverUtils.solve(f, -0.2, 0.2);
        assertEquals(0, result, 1E-8);
        result = UnivariateSolverUtils.solve(f, -0.1, 0.3);
        assertEquals(0, result, 1E-8);
        result = UnivariateSolverUtils.solve(f, -0.3, 0.45);
        assertEquals(0, result, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.3, 0.7);
        assertEquals(0.5, result, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.2, 0.6);
        assertEquals(0.5, result, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.05, 0.95);
        assertEquals(0.5, result, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.85, 1.25);
        assertEquals(1.0, result, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.8, 1.2);
        assertEquals(1.0, result, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.85, 1.75);
        assertEquals(1.0, result, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.55, 1.45);
        assertEquals(1.0, result, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.85, 5);
        assertEquals(1.0, result, 1E-6);
    }
}
