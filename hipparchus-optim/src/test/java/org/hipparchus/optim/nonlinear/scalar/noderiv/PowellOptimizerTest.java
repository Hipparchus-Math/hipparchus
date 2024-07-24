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
import org.hipparchus.analysis.SumSincFunction;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.MaxEval;
import org.hipparchus.optim.PointValuePair;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.GoalType;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Test for {@link PowellOptimizer}.
 */
public class PowellOptimizerTest {
    @Test
    public void testBoundsUnsupported() {
        assertThrows(MathRuntimeException.class, () -> {
            final MultivariateFunction func = new SumSincFunction(-1);
            final PowellOptimizer optim = new PowellOptimizer(1e-8, 1e-5,
                1e-4, 1e-4);

            optim.optimize(new MaxEval(100),
                new ObjectiveFunction(func),
                GoalType.MINIMIZE,
                new InitialGuess(new double[]{-3, 0}),
                new SimpleBounds(new double[]{-5, -1},
                    new double[]{5, 1}));
        });
    }

    @Test
    public void testSumSinc() {
        final MultivariateFunction func = new SumSincFunction(-1);

        int dim = 2;
        final double[] minPoint = new double[dim];
        for (int i = 0; i < dim; i++) {
            minPoint[i] = 0;
        }

        double[] init = new double[dim];

        // Initial is minimum.
        for (int i = 0; i < dim; i++) {
            init[i] = minPoint[i];
        }
        doTest(func, minPoint, init, GoalType.MINIMIZE, 1e-9, 1e-9);

        // Initial is far from minimum.
        for (int i = 0; i < dim; i++) {
            init[i] = minPoint[i] + 3;
        }
        doTest(func, minPoint, init, GoalType.MINIMIZE, 1e-9, 1e-5);
        // More stringent line search tolerance enhances the precision
        // of the result.
        doTest(func, minPoint, init, GoalType.MINIMIZE, 1e-9, 1e-9, 1e-7);
    }

    @Test
    public void testQuadratic() {
        final MultivariateFunction func = new MultivariateFunction() {
                public double value(double[] x) {
                    final double a = x[0] - 1;
                    final double b = x[1] - 1;
                    return a * a + b * b + 1;
                }
            };

        int dim = 2;
        final double[] minPoint = new double[dim];
        for (int i = 0; i < dim; i++) {
            minPoint[i] = 1;
        }

        double[] init = new double[dim];

        // Initial is minimum.
        for (int i = 0; i < dim; i++) {
            init[i] = minPoint[i];
        }
        doTest(func, minPoint, init, GoalType.MINIMIZE, 1e-9, 1e-8);

        // Initial is far from minimum.
        for (int i = 0; i < dim; i++) {
            init[i] = minPoint[i] - 20;
        }
        doTest(func, minPoint, init, GoalType.MINIMIZE, 1e-9, 1e-8);
    }

    @Test
    public void testMaximizeQuadratic() {
        final MultivariateFunction func = new MultivariateFunction() {
                public double value(double[] x) {
                    final double a = x[0] - 1;
                    final double b = x[1] - 1;
                    return -a * a - b * b + 1;
                }
            };

        int dim = 2;
        final double[] maxPoint = new double[dim];
        for (int i = 0; i < dim; i++) {
            maxPoint[i] = 1;
        }

        double[] init = new double[dim];

        // Initial is minimum.
        for (int i = 0; i < dim; i++) {
            init[i] = maxPoint[i];
        }
        doTest(func, maxPoint, init,  GoalType.MAXIMIZE, 1e-9, 1e-8);

        // Initial is far from minimum.
        for (int i = 0; i < dim; i++) {
            init[i] = maxPoint[i] - 20;
        }
        doTest(func, maxPoint, init, GoalType.MAXIMIZE, 1e-9, 1e-8);
    }

    /**
     * Ensure that we do not increase the number of function evaluations when
     * the function values are scaled up.
     * Note that the tolerances parameters passed to the constructor must
     * still hold sensible values because they are used to set the line search
     * tolerances.
     */
    @Test
    public void testRelativeToleranceOnScaledValues() {
        final MultivariateFunction func = new MultivariateFunction() {
                public double value(double[] x) {
                    final double a = x[0] - 1;
                    final double b = x[1] - 1;
                    return a * a * FastMath.sqrt(FastMath.abs(a)) + b * b + 1;
                }
            };

        int dim = 2;
        final double[] minPoint = new double[dim];
        for (int i = 0; i < dim; i++) {
            minPoint[i] = 1;
        }

        double[] init = new double[dim];
        // Initial is far from minimum.
        for (int i = 0; i < dim; i++) {
            init[i] = minPoint[i] - 20;
        }

        final double relTol = 1e-10;

        final int maxEval = 1000;
        // Very small absolute tolerance to rely solely on the relative
        // tolerance as a stopping criterion
        final PowellOptimizer optim = new PowellOptimizer(relTol, 1e-100);

        final PointValuePair funcResult = optim.optimize(new MaxEval(maxEval),
                                                         new ObjectiveFunction(func),
                                                         GoalType.MINIMIZE,
                                                         new InitialGuess(init));
        final double funcValue = func.value(funcResult.getPoint());
        final int funcEvaluations = optim.getEvaluations();

        final double scale = 1e10;
        final MultivariateFunction funcScaled = new MultivariateFunction() {
                public double value(double[] x) {
                    return scale * func.value(x);
                }
            };

        final PointValuePair funcScaledResult = optim.optimize(new MaxEval(maxEval),
                                                               new ObjectiveFunction(funcScaled),
                                                               GoalType.MINIMIZE,
                                                               new InitialGuess(init));
        final double funcScaledValue = funcScaled.value(funcScaledResult.getPoint());
        final int funcScaledEvaluations = optim.getEvaluations();

        // Check that both minima provide the same objective funciton values,
        // within the relative function tolerance.
        Assertions.assertEquals(1, funcScaledValue / (scale * funcValue), relTol);

        // Check that the numbers of evaluations are the same.
        Assertions.assertEquals(funcEvaluations, funcScaledEvaluations);
    }

    /**
     * @param func Function to optimize.
     * @param optimum Expected optimum.
     * @param init Starting point.
     * @param goal Minimization or maximization.
     * @param fTol Tolerance (relative error on the objective function) for
     * "Powell" algorithm.
     * @param pointTol Tolerance for checking that the optimum is correct.
     */
    private void doTest(MultivariateFunction func,
                        double[] optimum,
                        double[] init,
                        GoalType goal,
                        double fTol,
                        double pointTol) {
        final PowellOptimizer optim = new PowellOptimizer(fTol, FastMath.ulp(1d));

        final PointValuePair result = optim.optimize(new MaxEval(1000),
                                                     new ObjectiveFunction(func),
                                                     goal,
                                                     new InitialGuess(init));
        final double[] point = result.getPoint();

        for (int i = 0, dim = optimum.length; i < dim; i++) {
            Assertions.assertEquals(optimum[i], point[i], pointTol, "found[" + i + "]=" + point[i] + " value=" + result.getValue());
        }
    }

    /**
     * @param func Function to optimize.
     * @param optimum Expected optimum.
     * @param init Starting point.
     * @param goal Minimization or maximization.
     * @param fTol Tolerance (relative error on the objective function) for
     * "Powell" algorithm.
     * @param fLineTol Tolerance (relative error on the objective function)
     * for the internal line search algorithm.
     * @param pointTol Tolerance for checking that the optimum is correct.
     */
    private void doTest(MultivariateFunction func,
                        double[] optimum,
                        double[] init,
                        GoalType goal,
                        double fTol,
                        double fLineTol,
                        double pointTol) {
        final PowellOptimizer optim = new PowellOptimizer(fTol, FastMath.ulp(1d),
                                                          fLineTol, FastMath.ulp(1d));

        final PointValuePair result = optim.optimize(new MaxEval(1000),
                                                     new ObjectiveFunction(func),
                                                     goal,
                                                     new InitialGuess(init));
        final double[] point = result.getPoint();

        for (int i = 0, dim = optimum.length; i < dim; i++) {
            Assertions.assertEquals(optimum[i], point[i], pointTol, "found[" + i + "]=" + point[i] + " value=" + result.getValue());
        }

        Assertions.assertTrue(optim.getIterations() > 0);
    }
}
