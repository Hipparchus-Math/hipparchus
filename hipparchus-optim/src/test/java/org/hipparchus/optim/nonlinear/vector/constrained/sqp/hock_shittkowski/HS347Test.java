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

package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/**
 * Hock-Schittkowski test problem TP347.
 *
 * <p>The problem is:</p>
 *
 * <pre>
 * minimize
 *
 *     8204.37 log(H1 / H2)
 *   + 9008.72 log(H4 / H5)
 *   + 9330.46 log(H7 / H8)
 *
 * subject to
 *
 *     x1 + x2 + x3 = 1
 *
 *     0 <= x1 <= 1
 *     0 <= x2 <= 1
 *     0 <= x3 <= 1
 *
 * where
 *
 *     H1 = x1 + x2 + x3 + 0.03
 *     H2 = 0.09 x1 + x2 + x3 + 0.03
 *
 *     H4 = x2 + x3 + 0.03
 *     H5 = 0.07 x2 + x3 + 0.03
 *
 *     H7 = x3 + 0.03
 *     H8 = 0.13 x3 + 0.03
 * </pre>
 */
public class HS347Test {

    /** Problem dimension. */
    private static final int DIM = 3;

    /** Objective coefficients from the original Fortran problem. */
    private static final double A1 = 8204.37;

    /** Objective coefficients from the original Fortran problem. */
    private static final double A2 = 9008.72;

    /** Objective coefficients from the original Fortran problem. */
    private static final double A3 = 9330.46;

    /**
     * Lower protection used by the original Fortran {@code DMAX1}
     * expressions inside the logarithms.
     */
    private static final double LOG_ARGUMENT_MIN = 1.0e-4;

    /**
     * Container for the intermediate H values used by the objective
     * and its gradient.
     *
     * <p>The array deliberately retains the original one-based Fortran
     * indexing. Entry zero is unused.</p>
     */
    private static final class HValues {

        /** Intermediate values H(1), ..., H(8). */
        private final double[] h;

        /**
         * Compute the intermediate values.
         *
         * @param x current point
         */
        HValues(final RealVector x) {

            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);

            h = new double[9];

            h[1] = x1 + x2 + x3 + 0.03;
            h[2] = 0.09 * x1 + x2 + x3 + 0.03;
            h[3] = h[1] * h[2];

            h[4] = x2 + x3 + 0.03;
            h[5] = 0.07 * x2 + x3 + 0.03;
            h[6] = h[4] * h[5];

            h[7] = x3 + 0.03;
            h[8] = 0.13 * x3 + 0.03;
        }

        /**
         * Get an intermediate value using its original Fortran index.
         *
         * @param index one-based index
         * @return H(index)
         */
        double get(final int index) {
            return h[index];
        }
    }

    /** Objective function for TP347. */
    private static final class HS347Obj extends TwiceDifferentiableFunction {

        /** {@inheritDoc} */
        @Override
        public int dim() {
            return DIM;
        }

        /** {@inheritDoc} */
        @Override
        public double value(final RealVector x) {

            final HValues h = new HValues(x);

            /*
             * Exact equivalent of the original Fortran expressions:
             *
             * DLOG(DMAX1(H(1) / H(2), 1.D-4))
             * DLOG(DMAX1(H(4) / H(5), 1.D-4))
             * DLOG(DMAX1(H(7) / H(8), 1.D-4))
             */
            final double ratio1 =
                    FastMath.max(h.get(1) / h.get(2),
                                 LOG_ARGUMENT_MIN);

            final double ratio2 =
                    FastMath.max(h.get(4) / h.get(5),
                                 LOG_ARGUMENT_MIN);

            final double ratio3 =
                    FastMath.max(h.get(7) / h.get(8),
                                 LOG_ARGUMENT_MIN);

            return A1 * FastMath.log(ratio1) +
                   A2 * FastMath.log(ratio2) +
                   A3 * FastMath.log(ratio3);
        }

        /** {@inheritDoc} */
        @Override
        public RealVector gradient(final RealVector x) {

            final HValues h = new HValues(x);

            /*
             * Derivative of:
             *
             *     A log(u / v)
             *
             * is:
             *
             *     A (v du - u dv) / (u v).
             */

            final double gradient1 =
                    A1 *
                    (h.get(2) - 0.09 * h.get(1)) /
                    h.get(3);

            final double gradient2 =
                    A1 *
                    (h.get(2) - h.get(1)) /
                    h.get(3) +
                    A2 *
                    (h.get(5) - 0.07 * h.get(4)) /
                    h.get(6);

            final double gradient3 =
                    A1 *
                    (h.get(2) - h.get(1)) /
                    h.get(3) +
                    A2 *
                    (h.get(5) - h.get(4)) /
                    h.get(6) +
                    A3 *
                    (h.get(8) - 0.13 * h.get(7)) /
                    (h.get(7) * h.get(8));

            return new ArrayRealVector(
                    new double[] {
                        gradient1,
                        gradient2,
                        gradient3
                    },
                    false);
        }

        /** {@inheritDoc} */
        @Override
        public RealMatrix hessian(final RealVector x) {
            /*
             * The original Fortran problem does not provide an objective
             * Hessian. SQPOptimizerS2 therefore uses its quasi-Newton
             * Hessian approximation.
             */
            throw new UnsupportedOperationException(
                    "The exact objective Hessian is not provided for TP347.");
        }
    }

    /**
     * Generic equality constraint for:
     *
     * <pre>
     *     x1 + x2 + x3 - 1 = 0.
     * </pre>
     *
     * <p>The value returned by {@link #value(RealVector)} is already the
     * complete residual. Therefore the equality target passed to the
     * superclass must be zero.</p>
     */
    private static final class HS347Eq extends EqualityConstraint {

        /** Simple constructor. */
        HS347Eq() {
            super(new ArrayRealVector(
                    new double[] { 0.0 },
                    false));
        }

        /** {@inheritDoc} */
        @Override
        public int dim() {
            return DIM;
        }

        /** {@inheritDoc} */
        @Override
        public RealVector value(final RealVector x) {

            final double residual =
                    x.getEntry(0) +
                    x.getEntry(1) +
                    x.getEntry(2) -
                    1.0;

            return new ArrayRealVector(
                    new double[] { residual },
                    false);
        }

        /** {@inheritDoc} */
        @Override
        public RealMatrix jacobian(final RealVector x) {

            final RealMatrix jacobian =
                    MatrixUtils.createRealMatrix(1, DIM);

            jacobian.setEntry(0, 0, 1.0);
            jacobian.setEntry(0, 1, 1.0);
            jacobian.setEntry(0, 2, 1.0);

            return jacobian;
        }
    }

    /**
     * Get the original TP347 initial point.
     *
     * @return initial point
     */
    private static double[] start() {
        return new double[] {
            0.7,
            0.2,
            0.1
        };
    }

    /** Test TP347. */
    @Test
    public void testHS347() {

        final SQPOptimizerS2 optimizer =
                HSProblemTestUtils.newOptimizer();

        /*
         * Original variable bounds:
         *
         *     0 <= xi <= 1, i = 1, 2, 3.
         */
        final SimpleBounds bounds = new SimpleBounds(
                new double[] {
                    0.0,
                    0.0,
                    0.0
                },
                new double[] {
                    1.0,
                    1.0,
                    1.0
                });

        final LagrangeSolution solution = optimizer.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS347Obj()),
                new HS347Eq(),
                bounds);

        /*
         * Original validated objective value.
         *
         * The reference solution is approximately:
         *
         *     x = (0, 0, 1)
         *
         * with:
         *
         *     f(x) = 17374.625.
         */
        final double expectedObjective = 17374.625;

        HSProblemTestUtils.assertExpectedObjective(
                expectedObjective,
                solution);
    }
}