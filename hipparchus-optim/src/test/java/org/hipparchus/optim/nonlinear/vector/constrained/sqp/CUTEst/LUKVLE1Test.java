/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements. See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.CUTEst;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.OpenMapRealMatrix;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.MaxIter;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Disabled;

/**
 * CUTEst problem LUKVLE1.
 *
 * <p>Chained Rosenbrock objective with nonlinear equality constraints
 * containing cubic, trigonometric and exponential terms.</p>
 *
 * <p>The original SIF model is scalable. Change {@link #N} to 500 or 1000
 * after validating the conversion at the default size.</p>
 */
@Disabled
public class LUKVLE1Test {

    /** Number of variables. Must be at least 3. */
    private static final int N = 10000;

    /** Number of nonlinear equality constraints. */
    private static final int M = N - 2;

    /** Known global optimum. */
    private static final double EXPECTED_OBJECTIVE = 6.232458632;

    /** Equality feasibility tolerance used by this test. */
    private static final double FEASIBILITY_TOLERANCE = 1.0e-6;

    /**
     * Chained Rosenbrock objective:
     *
     * <pre>
     * f(x) = sum(i=0..N-2)
     *        [100 (x_i^2 - x_{i+1})^2 + (x_i - 1)^2].
     * </pre>
     */
    private static final class Objective extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector x) {
            double value = 0.0;

            for (int i = 0; i < N - 1; ++i) {
                final double xi = x.getEntry(i);
                final double residual1 = xi * xi - x.getEntry(i + 1);
                final double residual2 = xi - 1.0;

                value += 100.0 * residual1 * residual1 +
                         residual2 * residual2;
            }

            return value;
        }

        @Override
        public RealVector gradient(final RealVector x) {
            final RealVector gradient = new ArrayRealVector(N);

            for (int i = 0; i < N - 1; ++i) {
                final double xi = x.getEntry(i);
                final double residual = xi * xi - x.getEntry(i + 1);

                gradient.addToEntry(
                        i,
                        400.0 * xi * residual +
                        2.0 * (xi - 1.0));

                gradient.addToEntry(
                        i + 1,
                        -200.0 * residual);
            }

            return gradient;
        }

        @Override
        public RealMatrix hessian(final RealVector x) {
            final RealMatrix hessian =
                    new OpenMapRealMatrix(N, N);

            for (int i = 0; i < N - 1; ++i) {
                final double xi = x.getEntry(i);

                hessian.addToEntry(
                        i,
                        i,
                        1200.0 * xi * xi -
                        400.0 * x.getEntry(i + 1) +
                        2.0);

                hessian.addToEntry(
                        i + 1,
                        i + 1,
                        200.0);

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
            }

            return hessian;
        }
    }

    /**
     * Nonlinear equality constraints:
     *
     * <pre>
     * c_k(x) =
     *     3 x_{k+1}^3
     *   + sin(x_{k+1} - x_{k+2}) sin(x_{k+1} + x_{k+2})
     *   - x_k exp(x_k - x_{k+1})
     *   + 2 x_{k+2}
     *   + 4 x_{k+1}
     *   - 8
     *   = 0,
     *
     * k = 0,...,N-3.
     * </pre>
     */
    private static final class Equality extends EqualityConstraint {

        Equality() {
            super(new ArrayRealVector(M));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector x) {
            final RealVector constraints =
                    new ArrayRealVector(M);

            for (int k = 0; k < M; ++k) {
                final double x0 = x.getEntry(k);
                final double x1 = x.getEntry(k + 1);
                final double x2 = x.getEntry(k + 2);

                final double trigonometricTerm =
                        FastMath.sin(x1 - x2) *
                        FastMath.sin(x1 + x2);

                final double exponentialTerm =
                        x0 * FastMath.exp(x0 - x1);

                constraints.setEntry(
                        k,
                        3.0 * x1 * x1 * x1 +
                        trigonometricTerm -
                        exponentialTerm +
                        2.0 * x2 +
                        4.0 * x1 -
                        8.0);
            }

            return constraints;
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {
            final RealMatrix jacobian =
                    new OpenMapRealMatrix(M, N);

            for (int k = 0; k < M; ++k) {
                final double x0 = x.getEntry(k);
                final double x1 = x.getEntry(k + 1);
                final double x2 = x.getEntry(k + 2);

                final double exp =
                        FastMath.exp(x0 - x1);

                final double difference =
                        x1 - x2;

                final double sum =
                        x1 + x2;

                jacobian.setEntry(
                        k,
                        k,
                        -(1.0 + x0) * exp);

                jacobian.setEntry(
                        k,
                        k + 1,
                        9.0 * x1 * x1 +
                        FastMath.cos(difference) * FastMath.sin(sum) +
                        FastMath.sin(difference) * FastMath.cos(sum) +
                        x0 * exp +
                        4.0);

                jacobian.setEntry(
                        k,
                        k + 2,
                        2.0 -
                        FastMath.cos(difference) * FastMath.sin(sum) +
                        FastMath.sin(difference) * FastMath.cos(sum));
            }

            return jacobian;
        }
    }

    /**
     * Original alternating starting point:
     * x_1, x_3, ... = -1.2 and x_2, x_4, ... = 1.
     */
    private static double[] initialPoint() {
        final double[] initial = new double[N];

        for (int i = 0; i < N; ++i) {
            initial[i] = (i & 1) == 0 ? -1.2 : 1.0;
        }

        return initial;
    }
    
    @Test
    public void testLUKVLE1() {
        final Objective objective =
                new Objective();

        final Equality equality =
                new Equality();

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(5000),
                        new InitialGuess(initialPoint()),
                        new ObjectiveFunction(objective),
                        equality,
                        SimpleBounds.unbounded(N),
                        option);

       

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}