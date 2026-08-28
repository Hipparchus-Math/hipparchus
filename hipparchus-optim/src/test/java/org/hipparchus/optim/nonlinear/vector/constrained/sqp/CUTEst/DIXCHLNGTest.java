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
import org.junit.jupiter.api.Test;

/**
 * CUTEst problem {@code DIXCHLNG}.
 *
 * <p>This is the ten-variable equality-constrained challenge proposed by
 * L. C. W. Dixon for SQP methods. The standard CUTEst starting point lies
 * on the nonlinear constraint manifold and converges to the standard local
 * reference solution with objective approximately 2471.8978109.</p>
 */
public class DIXCHLNGTest {

    /** Number of variables. */
    private static final int N = 10;

    /** Number of equality constraints. */
    private static final int M = 5;

    /** Standard benchmark reference from the CUTEst starting point. */
    private static final double EXPECTED_OBJECTIVE =
            1.810282998041121e+03;

    /**
     * Dixon chained objective.
     */
    private static final class Objective
            extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {
            double value = 0.0;

            for (int i = 0; i < 7; ++i) {
                final double xi = point.getEntry(i);
                final double xi1 = point.getEntry(i + 1);
                final double xi2 = point.getEntry(i + 2);
                final double xi3 = point.getEntry(i + 3);

                final double r1 = xi1 - xi * xi;
                final double r2 = xi3 - xi2 * xi2;
                final double d0 = xi - 1.0;
                final double d1 = xi1 - 1.0;
                final double d2 = xi2 - 1.0;
                final double d3 = xi3 - 1.0;

                value += 100.0 * r1 * r1;
                value += d0 * d0;
                value += 90.0 * r2 * r2;
                value += d2 * d2;
                value += 10.1 * d1 * d1;
                value += 10.1 * d3 * d3;
                value += 19.8 * d1 * d3;
            }

            return value;
        }

        @Override
        public RealVector gradient(final RealVector point) {
            final RealVector gradient =
                    new ArrayRealVector(N);

            for (int i = 0; i < 7; ++i) {
                final double xi = point.getEntry(i);
                final double xi1 = point.getEntry(i + 1);
                final double xi2 = point.getEntry(i + 2);
                final double xi3 = point.getEntry(i + 3);

                addRosenbrockGradient(
                        gradient,
                        i,
                        i + 1,
                        xi,
                        xi1,
                        100.0);

                gradient.addToEntry(
                        i,
                        2.0 * (xi - 1.0));

                addRosenbrockGradient(
                        gradient,
                        i + 2,
                        i + 3,
                        xi2,
                        xi3,
                        90.0);

                gradient.addToEntry(
                        i + 2,
                        2.0 * (xi2 - 1.0));

                gradient.addToEntry(
                        i + 1,
                        20.2 * (xi1 - 1.0));

                gradient.addToEntry(
                        i + 3,
                        20.2 * (xi3 - 1.0));

                gradient.addToEntry(
                        i + 1,
                        19.8 * (xi3 - 1.0));

                gradient.addToEntry(
                        i + 3,
                        19.8 * (xi1 - 1.0));
            }

            return gradient;
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            final RealMatrix hessian =
                    new OpenMapRealMatrix(N, N);

            for (int i = 0; i < 7; ++i) {
                final double xi = point.getEntry(i);
                final double xi1 = point.getEntry(i + 1);
                final double xi2 = point.getEntry(i + 2);
                final double xi3 = point.getEntry(i + 3);

                addRosenbrockHessian(
                        hessian,
                        i,
                        i + 1,
                        xi,
                        xi1,
                        100.0);

                hessian.addToEntry(i, i, 2.0);

                addRosenbrockHessian(
                        hessian,
                        i + 2,
                        i + 3,
                        xi2,
                        xi3,
                        90.0);

                hessian.addToEntry(
                        i + 2,
                        i + 2,
                        2.0);

                hessian.addToEntry(
                        i + 1,
                        i + 1,
                        20.2);

                hessian.addToEntry(
                        i + 3,
                        i + 3,
                        20.2);

                hessian.addToEntry(
                        i + 1,
                        i + 3,
                        19.8);

                hessian.addToEntry(
                        i + 3,
                        i + 1,
                        19.8);
            }

            return hessian;
        }

        /**
         * Add derivatives of {@code scale * (v - u^2)^2}.
         */
        private static void addRosenbrockGradient(
                final RealVector gradient,
                final int uIndex,
                final int vIndex,
                final double u,
                final double v,
                final double scale) {

            final double residual = v - u * u;

            gradient.addToEntry(
                    uIndex,
                    -4.0 * scale * u * residual);

            gradient.addToEntry(
                    vIndex,
                    2.0 * scale * residual);
        }

        /**
         * Add the Hessian of {@code scale * (v - u^2)^2}.
         */
        private static void addRosenbrockHessian(
                final RealMatrix hessian,
                final int uIndex,
                final int vIndex,
                final double u,
                final double v,
                final double scale) {

            final double huu =
                    scale * (12.0 * u * u - 4.0 * v);

            final double huv =
                    -4.0 * scale * u;

            final double hvv =
                    2.0 * scale;

            hessian.addToEntry(
                    uIndex,
                    uIndex,
                    huu);

            hessian.addToEntry(
                    uIndex,
                    vIndex,
                    huv);

            hessian.addToEntry(
                    vIndex,
                    uIndex,
                    huv);

            hessian.addToEntry(
                    vIndex,
                    vIndex,
                    hvv);
        }
    }

    /**
     * Product equalities:
     *
     * <pre>
     * product(x_1,...,x_2)  = 1
     * product(x_1,...,x_4)  = 1
     * product(x_1,...,x_6)  = 1
     * product(x_1,...,x_8)  = 1
     * product(x_1,...,x_10) = 1
     * </pre>
     */
    private static final class ProductEqualities
            extends EqualityConstraint {

        ProductEqualities() {
            super(new ArrayRealVector(M));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector point) {
            final RealVector constraints =
                    new ArrayRealVector(M);

            double product = 1.0;
            int row = 0;

            for (int i = 0; i < N; ++i) {
                product *= point.getEntry(i);

                if ((i & 1) == 1) {
                    constraints.setEntry(
                            row++,
                            product - 1.0);
                }
            }

            return constraints;
        }

        @Override
        public RealMatrix jacobian(final RealVector point) {
            final RealMatrix jacobian =
                    new OpenMapRealMatrix(M, N);

            for (int row = 0; row < M; ++row) {
                final int endExclusive =
                        2 * (row + 1);

                for (int column = 0;
                     column < endExclusive;
                     ++column) {

                    double derivative = 1.0;

                    for (int k = 0;
                         k < endExclusive;
                         ++k) {

                        if (k != column) {
                            derivative *=
                                    point.getEntry(k);
                        }
                    }

                    jacobian.setEntry(
                            row,
                            column,
                            derivative);
                }
            }

            return jacobian;
        }
    }

    /**
     * Exact starting point from DIXCHLNG.SIF.
     */
    private static double[] initialPoint() {
        return new double[] {
            -2.0,
            -0.5,
             3.0,
             1.0 / 3.0,
            -4.0,
            -0.25,
             5.0,
             0.2,
            -6.0,
            -1.0 / 6.0
        };
    }

    @Test
    public void testDIXCHLNG() {
        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(5000),
                        new InitialGuess(initialPoint()),
                        new ObjectiveFunction(new Objective()),
                        new ProductEqualities(),
                        SimpleBounds.unbounded(N),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}