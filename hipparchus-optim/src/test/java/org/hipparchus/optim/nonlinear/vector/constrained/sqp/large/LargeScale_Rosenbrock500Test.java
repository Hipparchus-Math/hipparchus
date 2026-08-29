/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements. See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.hipparchus.optim.nonlinear.vector.constrained.sqp.large;


import java.util.Locale;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.MaxIter;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

/**
 * Large-scale chained Rosenbrock benchmark with 500 variables.
 *
 * <p>This test reproduces the Rosenbrock 500 problem used by the ripopt
 * large-scale benchmark:</p>
 *
 * <pre>
 * minimize
 *
 *     f(x) = sum(i = 0, ..., n - 2) [
 *                 (1 - x[i])^2
 *               + 100 (x[i + 1] - x[i]^2)^2
 *            ]
 *
 * with
 *
 *     n = 500
 *     x0 = (-1.2, ..., -1.2)
 *
 * and no constraints or variable bounds.
 * </pre>
 *
 * <p>The unique reference minimizer is:</p>
 *
 * <pre>
 *     x* = (1, ..., 1)
 *     f* = 0
 * </pre>
 */
public class LargeScale_Rosenbrock500Test {

    /** Number of optimization variables. */
    private static final int DIMENSION = 500;

    /** Maximum number of SQP iterations used by the reference benchmark. */
    private static final int MAX_ITERATIONS = 3000;

    /** Initial value assigned to every variable. */
    private static final double INITIAL_VALUE = -1.2;

    /** Reference optimal objective value is 0 global  and 3.9866238543 local . */
    private static final double EXPECTED_OBJECTIVE = 3.9866238543;

    /**
     * Objective tolerance corresponding to the tolerance used by the
     * ripopt large-scale comparison.
     */
    private static final double OBJECTIVE_TOLERANCE = 1.0e-7;

    /**
     * Chained Rosenbrock objective function.
     */
    private static final class ChainedRosenbrock
        extends TwiceDifferentiableFunction {

        /** {@inheritDoc} */
        @Override
        public int dim() {
            return DIMENSION;
        }

        /**
         * Compute
         *
         * <pre>
         * f(x) = sum [
         *            (1 - x[i])^2
         *          + 100 (x[i + 1] - x[i]^2)^2
         *        ].
         * </pre>
         *
         * @param x current point
         * @return objective value
         */
        @Override
        public double value(final RealVector x) {

            double objective = 0.0;

            for (int i = 0; i < DIMENSION - 1; ++i) {

                final double xi = x.getEntry(i);
                final double next = x.getEntry(i + 1);

                final double firstResidual = 1.0 - xi;
                final double secondResidual = next - xi * xi;

                objective += firstResidual * firstResidual +
                             100.0 * secondResidual * secondResidual;
            }

            return objective;
        }

        /**
         * Compute the exact analytical gradient.
         *
         * <p>The gradient is accumulated term by term to reproduce the
         * implementation used by the reference benchmark.</p>
         *
         * @param x current point
         * @return exact gradient
         */
        @Override
        public RealVector gradient(final RealVector x) {

            final RealVector gradient =
                    new ArrayRealVector(DIMENSION);

            for (int i = 0; i < DIMENSION - 1; ++i) {

                final double xi = x.getEntry(i);
                final double next = x.getEntry(i + 1);
                final double residual = next - xi * xi;

                /*
                 * Derivative with respect to x[i]:
                 *
                 * -2 (1 - x[i])
                 * -400 x[i] (x[i + 1] - x[i]^2)
                 */
                gradient.addToEntry(
                        i,
                        -2.0 * (1.0 - xi) -
                        400.0 * xi * residual);

                /*
                 * Derivative with respect to x[i + 1]:
                 *
                 * 200 (x[i + 1] - x[i]^2)
                 */
                gradient.addToEntry(
                        i + 1,
                        200.0 * residual);
            }

            return gradient;
        }

        /**
         * Compute the exact Rosenbrock Hessian.
         *
         * <p>The Hessian is tridiagonal. It is represented here through a
         * Hipparchus {@link RealMatrix}; no primitive {@code double[][]}
         * algebra is used.</p>
         *
         * @param x current point
         * @return exact Hessian
         */
        @Override
        public RealMatrix hessian(final RealVector x) {

            final RealMatrix hessian =
                    MatrixUtils.createRealMatrix(
                            DIMENSION,
                            DIMENSION);

            for (int i = 0; i < DIMENSION - 1; ++i) {

                final double xi = x.getEntry(i);
                final double next = x.getEntry(i + 1);

                /*
                 * Contribution to H[i, i]:
                 *
                 * 2 + 1200 x[i]^2 - 400 x[i + 1]
                 */
                hessian.addToEntry(
                        i,
                        i,
                        2.0 +
                        1200.0 * xi * xi -
                        400.0 * next);

                /*
                 * Symmetric off-diagonal contribution:
                 *
                 * H[i, i + 1] = H[i + 1, i] = -400 x[i]
                 */
                final double offDiagonal =
                        -400.0 * xi;

                hessian.addToEntry(
                        i,
                        i + 1,
                        offDiagonal);

                hessian.addToEntry(
                        i + 1,
                        i,
                        offDiagonal);

                /*
                 * Contribution to H[i + 1, i + 1]:
                 *
                 * 200
                 */
                hessian.addToEntry(
                        i + 1,
                        i + 1,
                        200.0);
            }

            return hessian;
        }
    }

    /**
     * Create the original benchmark initial point:
     *
     * <pre>
     * x[i] = -1.2, i = 0, ..., 499.
     * </pre>
     *
     * @return initial point
     */
    private static double[] initialPoint() {

        final double[] initial = new double[DIMENSION];

        for (int i = 0; i < DIMENSION; ++i) {
            initial[i] = INITIAL_VALUE;
        }

        return initial;
    }

    /**
     * Execute the Rosenbrock 500 benchmark.
     */
    @Test
    public void testRosenbrock500() {

        final SQPOptimizerS2 optimizer =
                LargeScaleProblemTestUtils.newOptimizer();

        /*
         * Request the exact external analytical gradient supplied by the
         * objective function.
         */
        final SQPOption option = new SQPOption();
        option.setGradientMode(GradientMode.EXTERNAL);
        option.setMaxIteration(1000);
        final long startTime = System.nanoTime();

        final LagrangeSolution solution = optimizer.optimize(
                new MaxIter(MAX_ITERATIONS),
                new InitialGuess(initialPoint()),
                new ObjectiveFunction(new ChainedRosenbrock()),
                option);

        final double elapsedSeconds =
                (System.nanoTime() - startTime) * 1.0e-9;

        final double objective = solution.getValue();

        System.out.printf(
                Locale.US,
                "%nRosenbrock 500%n" +
                "----------------%n" +
                "variables : %d%n" +
                "objective : %.16e%n" +
                "target    : %.16e%n" +
                "time      : %.6f s%n",
                DIMENSION,
                objective,
                EXPECTED_OBJECTIVE,
                elapsedSeconds);

        LargeScaleProblemTestUtils.assertExpectedObjective(EXPECTED_OBJECTIVE,solution);

    }
}
