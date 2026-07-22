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

/**
 * CUTEst problem {@code CHAIN} with {@code NH = 400}.
 *
 * <p>The problem computes a uniform suspended chain of length four between
 * the fixed endpoints {@code X(0) = 1} and {@code X(400) = 3}, minimizing
 * its discretized potential energy. The two fixed endpoint variables are
 * eliminated exactly as in CUTEst, leaving 800 free variables and
 * 401 nonlinear equality constraints.</p>
 *
 * <p>All derivatives are evaluated by forward finite differences.</p>
 */
public class CHAINTest {

    /** Number of mesh intervals. */
    private static final int INTERVALS = 400;

    /** Number of mesh points. */
    private static final int POINTS = INTERVALS + 1;

    /** Number of free variables after eliminating the two fixed heights. */
    private static final int N = 2 * INTERVALS;

    /** 400 defect constraints plus one length constraint. */
    private static final int M = INTERVALS + 1;

    /** Uniform mesh spacing. */
    private static final double H =
            1.0 / INTERVALS;

    /** Prescribed total chain length. */
    private static final double LENGTH =
            4.0;

    /** Fixed left endpoint height. */
    private static final double LEFT_HEIGHT =
            1.0;

    /** Fixed right endpoint height. */
    private static final double RIGHT_HEIGHT =
            3.0;

    /** Standard CUTEst reference objective for NH = 400. */
    private static final double EXPECTED_OBJECTIVE =
            5.06862e+00;

    /**
     * Trapezoidal approximation of the chain potential energy:
     *
     * <pre>
     * H * [ 0.5 * X(0)   * sqrt(1 + U(0)^2)
     *       + sum(k=1..399) X(k) * sqrt(1 + U(k)^2)
     *       + 0.5 * X(400) * sqrt(1 + U(400)^2) ].
     * </pre>
     */
    private static final class Objective
            extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {
            double value =
                    0.5 *
                    LEFT_HEIGHT *
                    arcLengthFactor(u(point, 0));

            for (int k = 1; k < INTERVALS; ++k) {
                value +=
                        x(point, k) *
                        arcLengthFactor(u(point, k));
            }

            value +=
                    0.5 *
                    RIGHT_HEIGHT *
                    arcLengthFactor(u(point, INTERVALS));

            return H * value;
        }

        @Override
        public RealVector gradient(final RealVector point) {
            throw new UnsupportedOperationException(
                    "Gradient is evaluated by finite differences.");
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            throw new UnsupportedOperationException(
                    "Hessian is not required by the finite-difference option.");
        }
    }

    /**
     * Discretized derivative and total-length equalities.
     */
    private static final class ChainEqualities
            extends EqualityConstraint {

        ChainEqualities() {
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

            for (int k = 0; k < INTERVALS; ++k) {
                constraints.setEntry(
                        k,
                        x(point, k) -
                        x(point, k + 1) +
                        0.5 * H *
                        (u(point, k) +
                         u(point, k + 1)));
            }

            double discreteLength =
                    0.5 *
                    arcLengthFactor(u(point, 0));

            for (int k = 1; k < INTERVALS; ++k) {
                discreteLength +=
                        arcLengthFactor(u(point, k));
            }

            discreteLength +=
                    0.5 *
                    arcLengthFactor(u(point, INTERVALS));

            constraints.setEntry(
                    INTERVALS,
                    H * discreteLength - LENGTH);

            return constraints;
        }

        @Override
        public RealMatrix jacobian(final RealVector point) {
            throw new UnsupportedOperationException(
                    "Jacobian is evaluated by finite differences.");
        }
    }

    /**
     * Exact CUTEst starting point.
     *
     * <p>For {@code t = k / 400}:</p>
     *
     * <pre>
     * X(k) = 8 t (0.5 t - 0.25) + 1
     * U(k) = 8 (t - 0.25).
     * </pre>
     */
    private static double[] initialPoint() {
        final double[] initial =
                new double[N];

        initial[uIndex(0)] =
                initialU(0);

        for (int k = 1; k < INTERVALS; ++k) {
            initial[xIndex(k)] =
                    initialX(k);

            initial[uIndex(k)] =
                    initialU(k);
        }

        initial[uIndex(INTERVALS)] =
                initialU(INTERVALS);

        return initial;
    }

    /**
     * Return X(k), including the two eliminated fixed endpoints.
     */
    private static double x(final RealVector point,
                            final int k) {
        if (k == 0) {
            return LEFT_HEIGHT;
        } else if (k == INTERVALS) {
            return RIGHT_HEIGHT;
        }

        return point.getEntry(xIndex(k));
    }

    /**
     * Return U(k).
     */
    private static double u(final RealVector point,
                            final int k) {
        return point.getEntry(uIndex(k));
    }

    /**
     * CUTEst free-variable order after fixed-variable elimination:
     *
     * <pre>
     * U(0), X(1), U(1), ..., X(399), U(399), U(400).
     * </pre>
     */
    private static int xIndex(final int k) {
        return 2 * k - 1;
    }

    private static int uIndex(final int k) {
        if (k == INTERVALS) {
            return N - 1;
        }

        return 2 * k;
    }

    private static double initialX(final int k) {
        final double t =
                (double) k / INTERVALS;

        return 8.0 *
               t *
               (0.5 * t - 0.25) +
               LEFT_HEIGHT;
    }

    private static double initialU(final int k) {
        final double t =
                (double) k / INTERVALS;

        return 8.0 * (t - 0.25);
    }

    private static double arcLengthFactor(final double derivative) {
        return FastMath.sqrt(
                1.0 + derivative * derivative);
    }

    @Test
    public void testCHAIN() {
        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(5000),
                        new InitialGuess(initialPoint()),
                        new ObjectiveFunction(new Objective()),
                        new ChainEqualities(),
                        SimpleBounds.unbounded(N),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}