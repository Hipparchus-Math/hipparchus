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

package org.hipparchus.analysis.solvers;

import org.hipparchus.analysis.QuinticFunction;
import org.hipparchus.analysis.RealFieldUnivariateFunction;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.Decimal64;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 */
public class UnivariateSolverUtilsTest {

    private UnivariateFunction sin = new Sin();
    private RealFieldUnivariateFunction<Decimal64> fieldSin = x -> x.sin();

    @Test(expected=NullArgumentException.class)
    public void testSolveNull() {
        UnivariateSolverUtils.solve(null, 0.0, 4.0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testSolveBadEndpoints() {
        double root = UnivariateSolverUtils.solve(sin, 4.0, -0.1, 1e-6);
        System.out.println("root=" + root);
    }

    @Test
    public void testSolveBadAccuracy() {
        try { // bad accuracy
            UnivariateSolverUtils.solve(sin, 0.0, 4.0, 0.0);
//             Assert.fail("Expecting MathIllegalArgumentException"); // TODO needs rework since convergence behaviour was changed
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testSolveSin() {
        double x = UnivariateSolverUtils.solve(sin, 1.0, 4.0);
        Assert.assertEquals(FastMath.PI, x, 1.0e-4);
    }

    @Test(expected=NullArgumentException.class)
    public void testSolveAccuracyNull()  {
        double accuracy = 1.0e-6;
        UnivariateSolverUtils.solve(null, 0.0, 4.0, accuracy);
    }

    @Test
    public void testSolveAccuracySin() {
        double accuracy = 1.0e-6;
        double x = UnivariateSolverUtils.solve(sin, 1.0,
                4.0, accuracy);
        Assert.assertEquals(FastMath.PI, x, accuracy);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testSolveNoRoot() {
        UnivariateSolverUtils.solve(sin, 1.0, 1.5);
    }

    @Test
    public void testBracketSin() {
        double[] result = UnivariateSolverUtils.bracket(sin,
                0.0, -2.0, 2.0);
        Assert.assertTrue(sin.value(result[0]) < 0);
        Assert.assertTrue(sin.value(result[1]) > 0);
    }

    @Test
    public void testBracketCentered() {
        double initial = 0.1;
        double[] result = UnivariateSolverUtils.bracket(sin, initial, -2.0, 2.0, 0.2, 1.0, 100);
        Assert.assertTrue(result[0] < initial);
        Assert.assertTrue(result[1] > initial);
        Assert.assertTrue(sin.value(result[0]) < 0);
        Assert.assertTrue(sin.value(result[1]) > 0);
    }

    @Test
    public void testBracketLow() {
        double initial = 0.5;
        double[] result = UnivariateSolverUtils.bracket(sin, initial, -2.0, 2.0, 0.2, 1.0, 100);
        Assert.assertTrue(result[0] < initial);
        Assert.assertTrue(result[1] < initial);
        Assert.assertTrue(sin.value(result[0]) < 0);
        Assert.assertTrue(sin.value(result[1]) > 0);
    }

    @Test
    public void testBracketHigh(){
        double initial = -0.5;
        double[] result = UnivariateSolverUtils.bracket(sin, initial, -2.0, 2.0, 0.2, 1.0, 100);
        Assert.assertTrue(result[0] > initial);
        Assert.assertTrue(result[1] > initial);
        Assert.assertTrue(sin.value(result[0]) < 0);
        Assert.assertTrue(sin.value(result[1]) > 0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testBracketLinear(){
        UnivariateSolverUtils.bracket(new UnivariateFunction() {
            public double value(double x) {
                return 1 - x;
            }
        }, 1000, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 100);
    }

    @Test
    public void testBracketExponential(){
        double[] result = UnivariateSolverUtils.bracket(new UnivariateFunction() {
            public double value(double x) {
                return 1 - x;
            }
        }, 1000, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0, 2.0, 10);
        Assert.assertTrue(result[0] <= 1);
        Assert.assertTrue(result[1] >= 1);
    }

    @Test
    public void testBracketEndpointRoot() {
        double[] result = UnivariateSolverUtils.bracket(sin, 1.5, 0, 2.0, 100);
        Assert.assertEquals(0.0, sin.value(result[0]), 1.0e-15);
        Assert.assertTrue(sin.value(result[1]) > 0);
    }

    @Test(expected=NullArgumentException.class)
    public void testNullFunction() {
        UnivariateSolverUtils.bracket(null, 1.5, 0, 2.0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testBadInitial() {
        UnivariateSolverUtils.bracket(sin, 2.5, 0, 2.0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testBadAdditive() {
        UnivariateSolverUtils.bracket(sin, 1.0, -2.0, 3.0, -1.0, 1.0, 100);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testIterationExceeded() {
        UnivariateSolverUtils.bracket(sin, 1.0, -2.0, 3.0, 1.0e-5, 1.0, 100);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testBadEndpoints() {
        // endpoints not valid
        UnivariateSolverUtils.bracket(sin, 1.5, 2.0, 1.0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testBadMaximumIterations() {
        // bad maximum iterations
        UnivariateSolverUtils.bracket(sin, 1.5, 0, 2.0, 0);
    }

    @Test
    public void testFieldBracketSin() {
        Decimal64[] result = UnivariateSolverUtils.bracket(fieldSin,new Decimal64(0.0),
                                                           new Decimal64(-2.0),new Decimal64(2.0));
        Assert.assertTrue(fieldSin.value(result[0]).getReal() < 0);
        Assert.assertTrue(fieldSin.value(result[1]).getReal() > 0);
    }

    @Test
    public void testFieldBracketCentered() {
        Decimal64 initial = new Decimal64(0.1);
        Decimal64[] result = UnivariateSolverUtils.bracket(fieldSin, initial,
                                                           new Decimal64(-2.0), new Decimal64(2.0),
                                                           new Decimal64(0.2), new Decimal64(1.0),
                                                           100);
        Assert.assertTrue(result[0].getReal() < initial.getReal());
        Assert.assertTrue(result[1].getReal() > initial.getReal());
        Assert.assertTrue(fieldSin.value(result[0]).getReal() < 0);
        Assert.assertTrue(fieldSin.value(result[1]).getReal() > 0);
    }

    @Test
    public void testFieldBracketLow() {
        Decimal64 initial = new Decimal64(0.5);
        Decimal64[] result = UnivariateSolverUtils.bracket(fieldSin, initial,
                                                           new Decimal64(-2.0), new Decimal64(2.0),
                                                           new Decimal64(0.2), new Decimal64(1.0),
                                                           100);
        Assert.assertTrue(result[0].getReal() < initial.getReal());
        Assert.assertTrue(result[1].getReal() < initial.getReal());
        Assert.assertTrue(fieldSin.value(result[0]).getReal() < 0);
        Assert.assertTrue(fieldSin.value(result[1]).getReal() > 0);
    }

    @Test
    public void testFieldBracketHigh(){
        Decimal64 initial = new Decimal64(-0.5);
        Decimal64[] result = UnivariateSolverUtils.bracket(fieldSin, initial,
                                                           new Decimal64(-2.0), new Decimal64(2.0),
                                                           new Decimal64(0.2), new Decimal64(1.0),
                                                           100);
        Assert.assertTrue(result[0].getReal() > initial.getReal());
        Assert.assertTrue(result[1].getReal() > initial.getReal());
        Assert.assertTrue(fieldSin.value(result[0]).getReal() < 0);
        Assert.assertTrue(fieldSin.value(result[1]).getReal() > 0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testFieldBracketLinear(){
        UnivariateSolverUtils.bracket(new RealFieldUnivariateFunction<Decimal64>() {
            public Decimal64 value(Decimal64 x) {
                return x.negate().add(1);
            }
        },
        new Decimal64(1000),
        new Decimal64(Double.NEGATIVE_INFINITY), new Decimal64(Double.POSITIVE_INFINITY),
        new Decimal64(1.0), new Decimal64(1.0), 100);
    }

    @Test
    public void testFieldBracketExponential(){
        Decimal64[] result = UnivariateSolverUtils.bracket(new RealFieldUnivariateFunction<Decimal64>() {
            public Decimal64 value(Decimal64 x) {
                return x.negate().add(1);
            }
        },
        new Decimal64(1000),
        new Decimal64(Double.NEGATIVE_INFINITY), new Decimal64(Double.POSITIVE_INFINITY),
        new Decimal64(1.0), new Decimal64(2.0), 10);
        Assert.assertTrue(result[0].getReal() <= 1);
        Assert.assertTrue(result[1].getReal() >= 1);
    }

    @Test
    public void testFieldBracketEndpointRoot() {
        Decimal64[] result = UnivariateSolverUtils.bracket(fieldSin,
                                                           new Decimal64(1.5), new Decimal64(0),
                                                           new Decimal64(2.0), 100);
        Assert.assertEquals(0.0, fieldSin.value(result[0]).getReal(), 1.0e-15);
        Assert.assertTrue(fieldSin.value(result[1]).getReal() > 0);
    }

    @Test(expected=NullArgumentException.class)
    public void testFieldNullFunction() {
        UnivariateSolverUtils.bracket(null, new Decimal64(1.5), new Decimal64(0), new Decimal64(2.0));
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testFieldBadInitial() {
        UnivariateSolverUtils.bracket(fieldSin, new Decimal64(2.5), new Decimal64(0), new Decimal64(2.0));
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testFieldBadAdditive() {
        UnivariateSolverUtils.bracket(fieldSin, new Decimal64(1.0), new Decimal64(-2.0), new Decimal64(3.0),
                                      new Decimal64(-1.0), new Decimal64(1.0), 100);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testFieldIterationExceeded() {
        UnivariateSolverUtils.bracket(fieldSin, new Decimal64(1.0), new Decimal64(-2.0), new Decimal64(3.0),
                                      new Decimal64(1.0e-5), new Decimal64(1.0), 100);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testFieldBadEndpoints() {
        // endpoints not valid
        UnivariateSolverUtils.bracket(fieldSin, new Decimal64(1.5), new Decimal64(2.0), new Decimal64(1.0));
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testFieldBadMaximumIterations() {
        // bad maximum iterations
        UnivariateSolverUtils.bracket(fieldSin, new Decimal64(1.5), new Decimal64(0), new Decimal64(2.0), 0);
    }

    /** check the search continues when a = lowerBound and b &lt; upperBound. */
    @Test
    public void testBracketLoopConditionForB() {
        double[] result = UnivariateSolverUtils.bracket(sin, -0.9, -1, 1, 0.1, 1, 100);
        Assert.assertTrue(result[0] <= 0);
        Assert.assertTrue(result[1] >= 0);
    }

    @Test
    public void testMisc() {
        UnivariateFunction f = new QuinticFunction();
        double result;
        // Static solve method
        result = UnivariateSolverUtils.solve(f, -0.2, 0.2);
        Assert.assertEquals(result, 0, 1E-8);
        result = UnivariateSolverUtils.solve(f, -0.1, 0.3);
        Assert.assertEquals(result, 0, 1E-8);
        result = UnivariateSolverUtils.solve(f, -0.3, 0.45);
        Assert.assertEquals(result, 0, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.3, 0.7);
        Assert.assertEquals(result, 0.5, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.2, 0.6);
        Assert.assertEquals(result, 0.5, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.05, 0.95);
        Assert.assertEquals(result, 0.5, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.85, 1.25);
        Assert.assertEquals(result, 1.0, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.8, 1.2);
        Assert.assertEquals(result, 1.0, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.85, 1.75);
        Assert.assertEquals(result, 1.0, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.55, 1.45);
        Assert.assertEquals(result, 1.0, 1E-6);
        result = UnivariateSolverUtils.solve(f, 0.85, 5);
        Assert.assertEquals(result, 1.0, 1E-6);
    }
}
