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

package org.hipparchus.optim.nonlinear.scalar.noderiv;

import org.hipparchus.analysis.MultivariateFunction;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.MaxEval;
import org.hipparchus.optim.PointValuePair;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.SimpleValueChecker;
import org.hipparchus.optim.nonlinear.scalar.GoalType;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class SimplexOptimizerMultiDirectionalTest {
    @Test
    void testBoundsUnsupported() {
        assertThrows(MathRuntimeException.class, () -> {
            SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
            final FourExtrema fourExtrema = new FourExtrema();

            optimizer.optimize(new MaxEval(100),
                new ObjectiveFunction(fourExtrema),
                GoalType.MINIMIZE,
                new InitialGuess(new double[]{-3, 0}),
                new NelderMeadSimplex(new double[]{0.2, 0.2}),
                new SimpleBounds(new double[]{-5, -1},
                    new double[]{5, 1}));
        });
    }

    @Test
    void testMinimize1() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-11, 1e-30);
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(200),
                                 new ObjectiveFunction(fourExtrema),
                                 GoalType.MINIMIZE,
                                 new InitialGuess(new double[] { -3, 0 }),
                                 new MultiDirectionalSimplex(new double[] { 0.2, 0.2 }));
        assertEquals(fourExtrema.xM, optimum.getPoint()[0], 4e-6);
        assertEquals(fourExtrema.yP, optimum.getPoint()[1], 3e-6);
        assertEquals(fourExtrema.valueXmYp, optimum.getValue(), 8e-13);
        assertTrue(optimizer.getEvaluations() > 120);
        assertTrue(optimizer.getEvaluations() < 150);

        // Check that the number of iterations is updated (MATH-949).
        assertTrue(optimizer.getIterations() > 0);
    }

    @Test
    void testMinimize2() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-11, 1e-30);
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(200),
                                 new ObjectiveFunction(fourExtrema),
                                 GoalType.MINIMIZE,
                                 new InitialGuess(new double[] { 1, 0 }),
                                 new MultiDirectionalSimplex(new double[] { 0.2, 0.2 }));
        assertEquals(fourExtrema.xP, optimum.getPoint()[0], 2e-8);
        assertEquals(fourExtrema.yM, optimum.getPoint()[1], 3e-6);
        assertEquals(fourExtrema.valueXpYm, optimum.getValue(), 2e-12);
        assertTrue(optimizer.getEvaluations() > 120);
        assertTrue(optimizer.getEvaluations() < 150);

        // Check that the number of iterations is updated (MATH-949).
        assertTrue(optimizer.getIterations() > 0);
    }

    @Test
    void testMaximize1() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-11, 1e-30);
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(200),
                                 new ObjectiveFunction(fourExtrema),
                                 GoalType.MAXIMIZE,
                                 new InitialGuess(new double[] { -3.0, 0.0 }),
                                 new MultiDirectionalSimplex(new double[] { 0.2, 0.2 }));
        assertEquals(fourExtrema.xM, optimum.getPoint()[0], 7e-7);
        assertEquals(fourExtrema.yM, optimum.getPoint()[1], 3e-7);
        assertEquals(fourExtrema.valueXmYm, optimum.getValue(), 2e-14);
        assertTrue(optimizer.getEvaluations() > 120);
        assertTrue(optimizer.getEvaluations() < 150);

        // Check that the number of iterations is updated (MATH-949).
        assertTrue(optimizer.getIterations() > 0);
    }

    @Test
    void testMaximize2() {
        SimplexOptimizer optimizer = new SimplexOptimizer(new SimpleValueChecker(1e-15, 1e-30));
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(200),
                                 new ObjectiveFunction(fourExtrema),
                                 GoalType.MAXIMIZE,
                                 new InitialGuess(new double[] { 1, 0 }),
                                 new MultiDirectionalSimplex(new double[] { 0.2, 0.2 }));
        assertEquals(fourExtrema.xP, optimum.getPoint()[0], 2e-8);
        assertEquals(fourExtrema.yP, optimum.getPoint()[1], 3e-6);
        assertEquals(fourExtrema.valueXpYp, optimum.getValue(), 2e-12);
        assertTrue(optimizer.getEvaluations() > 180);
        assertTrue(optimizer.getEvaluations() < 220);

        // Check that the number of iterations is updated (MATH-949).
        assertTrue(optimizer.getIterations() > 0);
    }

    @Test
    void testRosenbrock() {
        MultivariateFunction rosenbrock
            = new MultivariateFunction() {
                    public double value(double[] x) {
                        ++count;
                        double a = x[1] - x[0] * x[0];
                        double b = 1.0 - x[0];
                        return 100 * a * a + b * b;
                    }
                };

        count = 0;
        SimplexOptimizer optimizer = new SimplexOptimizer(-1, 1e-3);
        PointValuePair optimum
           = optimizer.optimize(new MaxEval(100),
                                new ObjectiveFunction(rosenbrock),
                                GoalType.MINIMIZE,
                                new InitialGuess(new double[] { -1.2, 1 }),
                                new MultiDirectionalSimplex(new double[][] {
                                        { -1.2,  1.0 },
                                        { 0.9, 1.2 },
                                        {  3.5, -2.3 } }));

        assertEquals(count, optimizer.getEvaluations());
        assertTrue(optimizer.getEvaluations() > 50);
        assertTrue(optimizer.getEvaluations() < 100);
        assertTrue(optimum.getValue() > 1e-2);
    }

    @Test
    void testPowell() {
        MultivariateFunction powell
            = new MultivariateFunction() {
                    public double value(double[] x) {
                        ++count;
                        double a = x[0] + 10 * x[1];
                        double b = x[2] - x[3];
                        double c = x[1] - 2 * x[2];
                        double d = x[0] - x[3];
                        return a * a + 5 * b * b + c * c * c * c + 10 * d * d * d * d;
                    }
                };

        count = 0;
        SimplexOptimizer optimizer = new SimplexOptimizer(-1, 1e-3);
        PointValuePair optimum
            = optimizer.optimize(new MaxEval(1000),
                                 new ObjectiveFunction(powell),
                                 GoalType.MINIMIZE,
                                 new InitialGuess(new double[] { 3, -1, 0, 1 }),
                                 new MultiDirectionalSimplex(4));
        assertEquals(count, optimizer.getEvaluations());
        assertTrue(optimizer.getEvaluations() > 800);
        assertTrue(optimizer.getEvaluations() < 900);
        assertTrue(optimum.getValue() > 1e-2);
    }

    @Test
    void testMath283() {
        // fails because MultiDirectional.iterateSimplex is looping forever
        // the while(true) should be replaced with a convergence check
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-14, 1e-14);
        final Gaussian2D function = new Gaussian2D(0, 0, 1);
        PointValuePair estimate = optimizer.optimize(new MaxEval(1000),
                                                     new ObjectiveFunction(function),
                                                     GoalType.MAXIMIZE,
                                                     new InitialGuess(function.getMaximumPosition()),
                                                     new MultiDirectionalSimplex(2));
        final double EPSILON = 1e-5;
        final double expectedMaximum = function.getMaximum();
        final double actualMaximum = estimate.getValue();
        assertEquals(expectedMaximum, actualMaximum, EPSILON);

        final double[] expectedPosition = function.getMaximumPosition();
        final double[] actualPosition = estimate.getPoint();
        assertEquals(expectedPosition[0], actualPosition[0], EPSILON );
        assertEquals(expectedPosition[1], actualPosition[1], EPSILON );
    }

    private static class FourExtrema implements MultivariateFunction {
        // The following function has 4 local extrema.
        final double xM = -3.841947088256863675365;
        final double yM = -1.391745200270734924416;
        final double xP =  0.2286682237349059125691;
        final double yP = -yM;
        final double valueXmYm = 0.2373295333134216789769; // Local maximum.
        final double valueXmYp = -valueXmYm; // Local minimum.
        final double valueXpYm = -0.7290400707055187115322; // Global minimum.
        final double valueXpYp = -valueXpYm; // Global maximum.

        public double value(double[] variables) {
            final double x = variables[0];
            final double y = variables[1];
            return (x == 0 || y == 0) ? 0 :
                FastMath.atan(x) * FastMath.atan(x + 2) * FastMath.atan(y) * FastMath.atan(y) / (x * y);
        }
    }

    private static class Gaussian2D implements MultivariateFunction {
        private final double[] maximumPosition;
        private final double std;

        public Gaussian2D(double xOpt, double yOpt, double std) {
            maximumPosition = new double[] { xOpt, yOpt };
            this.std = std;
        }

        public double getMaximum() {
            return value(maximumPosition);
        }

        public double[] getMaximumPosition() {
            return maximumPosition.clone();
        }

        public double value(double[] point) {
            final double x = point[0], y = point[1];
            final double twoS2 = 2.0 * std * std;
            return 1.0 / (twoS2 * FastMath.PI) * FastMath.exp(-(x * x + y * y) / twoS2);
        }
    }

    private int count;
}
