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
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

/**
 * CUTEst problem LUKVLE3.
 *
 * <p>Chained Powell singular objective with two nonlinear equality
 * constraints. The formulation follows the official CUTEst SIF model.</p>
 */
@Disabled
public class LUKVLE3Test {

    /** Number of variables. N must be even and at least 6. */
    private static final int N = 10000;

    /**
     * Standard CUTEst solution value from the standard starting point.
     * This value is independent of the inactive zero interior blocks.
     */
    private static final double EXPECTED_OBJECTIVE = 2.758658e+1;

    /**
     * Chained Powell singular objective.
     */
    private static final class Objective extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector x) {
            double value = 0.0;

            for (int i = 0; i < (N - 2) / 2; ++i) {
                final int j = 2 * i;

                final double a1 =
                        x.getEntry(j) +
                        10.0 * x.getEntry(j + 1);

                final double a2 =
                        x.getEntry(j + 2) -
                        x.getEntry(j + 3);

                final double a3 =
                        x.getEntry(j + 1) -
                        2.0 * x.getEntry(j + 2);

                final double a4 =
                        x.getEntry(j) -
                        x.getEntry(j + 3);

                value += a1 * a1 +
                         5.0 * a2 * a2 +
                         a3 * a3 * a3 * a3 +
                         10.0 * a4 * a4 * a4 * a4;
            }

            return value;
        }

        @Override
        public RealVector gradient(final RealVector x) {
            final RealVector gradient =
                    new ArrayRealVector(N);

            for (int i = 0; i < (N - 2) / 2; ++i) {
                final int j = 2 * i;

                final double a1 =
                        x.getEntry(j) +
                        10.0 * x.getEntry(j + 1);

                final double a2 =
                        x.getEntry(j + 2) -
                        x.getEntry(j + 3);

                final double a3 =
                        x.getEntry(j + 1) -
                        2.0 * x.getEntry(j + 2);

                final double a4 =
                        x.getEntry(j) -
                        x.getEntry(j + 3);

                final double a3Cube =
                        a3 * a3 * a3;

                final double a4Cube =
                        a4 * a4 * a4;

                gradient.addToEntry(
                        j,
                        2.0 * a1 +
                        40.0 * a4Cube);

                gradient.addToEntry(
                        j + 1,
                        20.0 * a1 +
                        4.0 * a3Cube);

                gradient.addToEntry(
                        j + 2,
                        10.0 * a2 -
                        8.0 * a3Cube);

                gradient.addToEntry(
                        j + 3,
                        -10.0 * a2 -
                        40.0 * a4Cube);
            }

            return gradient;
        }

        @Override
        public RealMatrix hessian(final RealVector x) {
            final RealMatrix hessian =
                    new OpenMapRealMatrix(N, N);

            for (int i = 0; i < (N - 2) / 2; ++i) {
                final int j = 2 * i;

                final double a3 =
                        x.getEntry(j + 1) -
                        2.0 * x.getEntry(j + 2);

                final double a4 =
                        x.getEntry(j) -
                        x.getEntry(j + 3);

                final double a3Squared =
                        a3 * a3;

                final double a4Squared =
                        a4 * a4;

                hessian.addToEntry(
                        j,
                        j,
                        2.0 + 120.0 * a4Squared);

                hessian.addToEntry(
                        j + 1,
                        j + 1,
                        200.0 + 12.0 * a3Squared);

                hessian.addToEntry(
                        j + 2,
                        j + 2,
                        10.0 + 48.0 * a3Squared);

                hessian.addToEntry(
                        j + 3,
                        j + 3,
                        10.0 + 120.0 * a4Squared);

                addSymmetric(
                        hessian,
                        j,
                        j + 1,
                        20.0);

                addSymmetric(
                        hessian,
                        j,
                        j + 3,
                        -120.0 * a4Squared);

                addSymmetric(
                        hessian,
                        j + 1,
                        j + 2,
                        -24.0 * a3Squared);

                addSymmetric(
                        hessian,
                        j + 2,
                        j + 3,
                        -10.0);
            }

            return hessian;
        }

        private static void addSymmetric(final RealMatrix matrix,
                                         final int row,
                                         final int column,
                                         final double value) {
            matrix.addToEntry(row, column, value);
            matrix.addToEntry(column, row, value);
        }
    }

    /**
     * Two nonlinear equality constraints:
     *
     * <pre>
     * c1(x) =
     *     3 x0^3
     *   + 2 x1
     *   - 5
     *   + sin(x0 - x1) sin(x0 + x1)
     *   = 0
     *
     * c2(x) =
     *     4 x(N-2)
     *   - x(N-3) exp(x(N-3) - x(N-2))
     *   - 3
     *   = 0
     * </pre>
     */
    private static final class Equality extends EqualityConstraint {

        Equality() {
            super(new ArrayRealVector(2));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector x) {
            final RealVector constraints =
                    new ArrayRealVector(2);

            final double x0 =
                    x.getEntry(0);

            final double x1 =
                    x.getEntry(1);

            constraints.setEntry(
                    0,
                    3.0 * x0 * x0 * x0 +
                    2.0 * x1 -
                    5.0 +
                    FastMath.sin(x0 - x1) *
                    FastMath.sin(x0 + x1));

            final double a =
                    x.getEntry(N - 2);

            final double b =
                    x.getEntry(N - 1);

            constraints.setEntry(
                    1,
                    4.0 * a -
                    a * FastMath.exp(a - b) -
                    3.0);

            return constraints;
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {
            final RealMatrix jacobian =
                    new OpenMapRealMatrix(2, N);

            final double x0 =
                    x.getEntry(0);

            final double x1 =
                    x.getEntry(1);

            final double difference =
                    x0 - x1;

            final double sum =
                    x0 + x1;

            jacobian.setEntry(
                    0,
                    0,
                    9.0 * x0 * x0 +
                    FastMath.cos(difference) *
                    FastMath.sin(sum) +
                    FastMath.sin(difference) *
                    FastMath.cos(sum));

            jacobian.setEntry(
                    0,
                    1,
                    2.0 -
                    FastMath.cos(difference) *
                    FastMath.sin(sum) +
                    FastMath.sin(difference) *
                    FastMath.cos(sum));

            final double a =
                    x.getEntry(N - 2);

            final double b =
                    x.getEntry(N - 1);

            final double exponential =
                    FastMath.exp(a - b);

            jacobian.setEntry(
                    1,
                    N - 2,
                    4.0 -
                    (1.0 + a) * exponential);

            jacobian.setEntry(
                    1,
                    N - 1,
                    a * exponential);

            return jacobian;
        }
    }

    /**
     * Standard CUTEst starting point:
     * (3, -1, 0, 1) repeated.
     */
    private static double[] initialPoint() {
        final double[] initial =
                new double[N];

        for (int i = 0; i < N; ++i) {
            switch (i & 3) {
                case 0:
                    initial[i] = 3.0;
                    break;
                case 1:
                    initial[i] = -1.0;
                    break;
                case 2:
                    initial[i] = 0.0;
                    break;
                default:
                    initial[i] = 1.0;
                    break;
            }
        }

        return initial;
    }

    @Test
    public void testLUKVLE3() {
        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newExternalOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(5000),
                        new InitialGuess(initialPoint()),
                        new ObjectiveFunction(new Objective()),
                        new Equality(),
                        SimpleBounds.unbounded(N),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}