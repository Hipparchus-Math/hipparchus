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
package org.hipparchus.optim.univariate;


import org.hipparchus.UnitTestUtils;
import org.hipparchus.analysis.FunctionUtils;
import org.hipparchus.analysis.QuinticFunction;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.analysis.function.StepFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.optim.ConvergenceChecker;
import org.hipparchus.optim.MaxEval;
import org.hipparchus.optim.nonlinear.scalar.GoalType;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 */
final class BrentOptimizerTest {

    @Test
    void testSinMin() {
        UnivariateFunction f = new Sin();
        UnivariateOptimizer optimizer = new BrentOptimizer(1e-10, 1e-14);
        assertEquals(3 * FastMath.PI / 2, optimizer.optimize(new MaxEval(200),
                                                                    new UnivariateObjectiveFunction(f),
                                                                    GoalType.MINIMIZE,
                                                                    new SearchInterval(4, 5)).getPoint(), 1e-8);
        assertTrue(optimizer.getEvaluations() <= 50);
        assertEquals(200, optimizer.getMaxEvaluations());
        assertEquals(3 * FastMath.PI / 2, optimizer.optimize(new MaxEval(200),
                                                                    new UnivariateObjectiveFunction(f),
                                                                    GoalType.MINIMIZE,
                                                                    new SearchInterval(1, 5)).getPoint(), 1e-8);
        assertTrue(optimizer.getEvaluations() <= 100);
        assertTrue(optimizer.getEvaluations() >= 15);
        try {
            optimizer.optimize(new MaxEval(10),
                               new UnivariateObjectiveFunction(f),
                               GoalType.MINIMIZE,
                               new SearchInterval(4, 5));
            fail("an exception should have been thrown");
        } catch (MathIllegalStateException fee) {
            // expected
        }
    }

    @Test
    void testSinMinWithValueChecker() {
        final UnivariateFunction f = new Sin();
        final ConvergenceChecker<UnivariatePointValuePair> checker = new SimpleUnivariateValueChecker(1e-5, 1e-14);
        // The default stopping criterion of Brent's algorithm should not
        // pass, but the search will stop at the given relative tolerance
        // for the function value.
        final UnivariateOptimizer optimizer = new BrentOptimizer(1e-10, 1e-14, checker);
        final UnivariatePointValuePair result = optimizer.optimize(new MaxEval(200),
                                                                   new UnivariateObjectiveFunction(f),
                                                                   GoalType.MINIMIZE,
                                                                   new SearchInterval(4, 5));
        assertEquals(3 * FastMath.PI / 2, result.getPoint(), 1e-3);
    }

    @Test
    void testBoundaries() {
        final double lower = -1.0;
        final double upper = +1.0;
        UnivariateFunction f = new UnivariateFunction() {
            @Override
            public double value(double x) {
                if (x < lower) {
                    throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL,
                                                           x, lower);
                } else if (x > upper) {
                    throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_LARGE,
                                                           x, upper);
                } else {
                    return x;
                }
            }
        };
        UnivariateOptimizer optimizer = new BrentOptimizer(1e-10, 1e-14);
        assertEquals(lower,
                            optimizer.optimize(new MaxEval(100),
                                               new UnivariateObjectiveFunction(f),
                                               GoalType.MINIMIZE,
                                               new SearchInterval(lower, upper)).getPoint(),
                            1.0e-8);
        assertEquals(upper,
                            optimizer.optimize(new MaxEval(100),
                                               new UnivariateObjectiveFunction(f),
                                               GoalType.MAXIMIZE,
                                               new SearchInterval(lower, upper)).getPoint(),
                            1.0e-8);
    }

    @Test
    void testQuinticMin() {
        // The function has local minima at -0.27195613 and 0.82221643.
        UnivariateFunction f = new QuinticFunction();
        UnivariateOptimizer optimizer = new BrentOptimizer(1e-10, 1e-14);
        assertEquals(-0.27195613, optimizer.optimize(new MaxEval(200),
                                                            new UnivariateObjectiveFunction(f),
                                                            GoalType.MINIMIZE,
                                                            new SearchInterval(-0.3, -0.2)).getPoint(), 1.0e-8);
        assertEquals( 0.82221643, optimizer.optimize(new MaxEval(200),
                                                            new UnivariateObjectiveFunction(f),
                                                            GoalType.MINIMIZE,
                                                            new SearchInterval(0.3,  0.9)).getPoint(), 1.0e-8);
        assertTrue(optimizer.getEvaluations() <= 50);

        // search in a large interval
        assertEquals(-0.27195613, optimizer.optimize(new MaxEval(200),
                                                            new UnivariateObjectiveFunction(f),
                                                            GoalType.MINIMIZE,
                                                            new SearchInterval(-1.0, 0.2)).getPoint(), 1.0e-8);
        assertTrue(optimizer.getEvaluations() <= 50);
    }

    @Test
    void testQuinticMinStatistics() {
        // The function has local minima at -0.27195613 and 0.82221643.
        UnivariateFunction f = new QuinticFunction();
        UnivariateOptimizer optimizer = new BrentOptimizer(1e-11, 1e-14);

        final UnitTestUtils.SimpleStatistics[] stat = new UnitTestUtils.SimpleStatistics[2];
        for (int i = 0; i < stat.length; i++) {
            stat[i] = new UnitTestUtils.SimpleStatistics();
        }

        final double min = -0.75;
        final double max = 0.25;
        final int nSamples = 200;
        final double delta = (max - min) / nSamples;
        for (int i = 0; i < nSamples; i++) {
            final double start = min + i * delta;
            stat[0].addValue(optimizer.optimize(new MaxEval(40),
                                                new UnivariateObjectiveFunction(f),
                                                GoalType.MINIMIZE,
                                                new SearchInterval(min, max, start)).getPoint());
            stat[1].addValue(optimizer.getEvaluations());
        }

        final double meanOptValue = stat[0].getMean();
        final double medianEval = stat[1].getMedian();
        assertTrue(meanOptValue > -0.2719561281);
        assertTrue(meanOptValue < -0.2719561280);
        assertEquals(23, (int) medianEval);

        // MATH-1121: Ensure that the iteration counter is incremented.
        assertTrue(optimizer.getIterations() > 0);
    }

    @Test
    void testQuinticMax() {
        // The quintic function has zeros at 0, +-0.5 and +-1.
        // The function has a local maximum at 0.27195613.
        UnivariateFunction f = new QuinticFunction();
        UnivariateOptimizer optimizer = new BrentOptimizer(1e-12, 1e-14);
        assertEquals(0.27195613, optimizer.optimize(new MaxEval(100),
                                                           new UnivariateObjectiveFunction(f),
                                                           GoalType.MAXIMIZE,
                                                           new SearchInterval(0.2, 0.3)).getPoint(), 1e-8);
        try {
            optimizer.optimize(new MaxEval(5),
                               new UnivariateObjectiveFunction(f),
                               GoalType.MAXIMIZE,
                               new SearchInterval(0.2, 0.3));
            fail("an exception should have been thrown");
        } catch (MathIllegalStateException miee) {
            // expected
        }
    }

    @Test
    void testMinEndpoints() {
        UnivariateFunction f = new Sin();
        UnivariateOptimizer optimizer = new BrentOptimizer(1e-8, 1e-14);

        // endpoint is minimum
        double result = optimizer.optimize(new MaxEval(50),
                                           new UnivariateObjectiveFunction(f),
                                           GoalType.MINIMIZE,
                                           new SearchInterval(3 * FastMath.PI / 2, 5)).getPoint();
        assertEquals(3 * FastMath.PI / 2, result, 1e-6);

        result = optimizer.optimize(new MaxEval(50),
                                    new UnivariateObjectiveFunction(f),
                                    GoalType.MINIMIZE,
                                    new SearchInterval(4, 3 * FastMath.PI / 2)).getPoint();
        assertEquals(3 * FastMath.PI / 2, result, 1e-6);
    }

    @Test
    void testMath832() {
        final UnivariateFunction f = new UnivariateFunction() {
                @Override
                public double value(double x) {
                    final double sqrtX = FastMath.sqrt(x);
                    final double a = 1e2 * sqrtX;
                    final double b = 1e6 / x;
                    final double c = 1e4 / sqrtX;

                    return a + b + c;
                }
            };

        UnivariateOptimizer optimizer = new BrentOptimizer(1e-10, 1e-8);
        final double result = optimizer.optimize(new MaxEval(1483),
                                                 new UnivariateObjectiveFunction(f),
                                                 GoalType.MINIMIZE,
                                                 new SearchInterval(Double.MIN_VALUE,
                                                                    Double.MAX_VALUE)).getPoint();

        assertEquals(804.9355825, result, 1e-6);
    }

    /**
     * Contrived example showing that prior to the resolution of MATH-855
     * (second revision), the algorithm would not return the best point if
     * it happened to be the initial guess.
     */
    @Test
    void testKeepInitIfBest() {
        final double minSin = 3 * FastMath.PI / 2;
        final double offset = 1e-8;
        final double delta = 1e-7;
        final UnivariateFunction f1 = new Sin();
        final UnivariateFunction f2 = new StepFunction(new double[] { minSin, minSin + offset, minSin + 2 * offset},
                                                       new double[] { 0, -1, 0 });
        final UnivariateFunction f = FunctionUtils.add(f1, f2);
        // A slightly less stringent tolerance would make the test pass
        // even with the previous implementation.
        final double relTol = 1e-8;
        final UnivariateOptimizer optimizer = new BrentOptimizer(relTol, 1e-100);
        final double init = minSin + 1.5 * offset;
        final UnivariatePointValuePair result
            = optimizer.optimize(new MaxEval(200),
                                 new UnivariateObjectiveFunction(f),
                                 GoalType.MINIMIZE,
                                 new SearchInterval(minSin - 6.789 * delta,
                                                    minSin + 9.876 * delta,
                                                    init));

        final double sol = result.getPoint();
        final double expected = init;

//         System.out.println("numEval=" + numEval);
//         System.out.println("min=" + init + " f=" + f.value(init));
//         System.out.println("sol=" + sol + " f=" + f.value(sol));
//         System.out.println("exp=" + expected + " f=" + f.value(expected));

        assertTrue(f.value(sol) <= f.value(expected), "Best point not reported");
    }

    /**
     * Contrived example showing that prior to the resolution of MATH-855,
     * the algorithm, by always returning the last evaluated point, would
     * sometimes not report the best point it had found.
     */
    @Test
    void testMath855() {
        final double minSin = 3 * FastMath.PI / 2;
        final double offset = 1e-8;
        final double delta = 1e-7;
        final UnivariateFunction f1 = new Sin();
        final UnivariateFunction f2 = new StepFunction(new double[] { minSin, minSin + offset, minSin + 5 * offset },
                                                       new double[] { 0, -1, 0 });
        final UnivariateFunction f = FunctionUtils.add(f1, f2);
        final UnivariateOptimizer optimizer = new BrentOptimizer(1e-8, 1e-100);
        final UnivariatePointValuePair result
            = optimizer.optimize(new MaxEval(200),
                                 new UnivariateObjectiveFunction(f),
                                 GoalType.MINIMIZE,
                                 new SearchInterval(minSin - 6.789 * delta,
                                                    minSin + 9.876 * delta));

        final double sol = result.getPoint();
        final double expected = 4.712389027602411;

        // System.out.println("min=" + (minSin + offset) + " f=" + f.value(minSin + offset));
        // System.out.println("sol=" + sol + " f=" + f.value(sol));
        // System.out.println("exp=" + expected + " f=" + f.value(expected));

        assertTrue(f.value(sol) <= f.value(expected), "Best point not reported");
    }
}
