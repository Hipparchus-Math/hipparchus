/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
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
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/**
 * CUTEst problem {@code HALDMADS}.
 *
 * <p>This is the nonlinear minimax rational-approximation problem of
 * Hald and Madsen. Five coefficients define a rational approximation of
 * {@code exp(y)} at 21 equally spaced points in {@code [-1, 1]}. The sixth
 * variable is the minimax error bound.</p>
 *
 * <p>The original CUTEst model has six variables and 42 nonlinear
 * inequalities. Some solver tables report 48 variables because they include
 * one slack variable for each inequality.</p>
 */
public class HALDMADSTest {

    /** Number of rational-function coefficients plus the minimax variable. */
    private static final int N = 6;

    /** Number of sample points. */
    private static final int SAMPLE_COUNT = 21;

    /** Two minimax inequalities for every sample point. */
    private static final int M = 2 * SAMPLE_COUNT;

    /** Index of the minimax variable U. */
    private static final int U_INDEX = 5;

    /**
     * Objective reported for HALDMADS in the IPOPT CUTEr result table.
     *
     * <p>The historical SIF file also contains the lower-precision solution
     * annotation {@code 0.0001207}; the value below intentionally reproduces
     * the result in the solver table used for this test selection.</p>
     */
    private static final double EXPECTED_OBJECTIVE =
            1.223573743e-4;

    /** Official CUTEst starting point. */
    private static final double[] START = {
        0.5, 0.0, 0.0, 0.0, 0.0, 0.0
    };

    /** All variables are free in the original SIF model. */
    private static final double[] LOWER = {
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY
    };

    /** All variables are free in the original SIF model. */
    private static final double[] UPPER = {
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY
    };

    /** Sample abscissas y = -1.0, -0.9, ..., 1.0. */
    private static final double[] SAMPLE_Y =
            new double[SAMPLE_COUNT];

    /** Exact target values exp(y). */
    private static final double[] SAMPLE_EXP =
            new double[SAMPLE_COUNT];

    static {
        for (int i = 0; i < SAMPLE_COUNT; ++i) {
            final double y = -1.0 + 0.1 * i;
            SAMPLE_Y[i] = y;
            SAMPLE_EXP[i] = FastMath.exp(y);
        }
    }

    /**
     * Linear minimax objective f(x) = U.
     */
    private static final class Objective
            extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {
            return point.getEntry(U_INDEX);
        }

        @Override
        public RealVector gradient(final RealVector point) {
            final double[] gradient = new double[N];
            gradient[U_INDEX] = 1.0;
            return new ArrayRealVector(gradient, false);
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            return new OpenMapRealMatrix(N, N);
        }
    }

    /**
     * Canonical minimax inequalities.
     *
     * <p>For each sample point, CUTEst defines:</p>
     *
     * <pre>
     * rational(y) - exp(y) - U <= 0
     * exp(y) - rational(y) - U <= 0
     * </pre>
     *
     * <p>The Hipparchus inequality convention is {@code g(x) >= 0}, hence
     * the returned pair is:</p>
     *
     * <pre>
     * U - rational(y) + exp(y) >= 0
     * U + rational(y) - exp(y) >= 0
     * </pre>
     */
    private static final class MinimaxConstraints
            extends InequalityConstraint {

        MinimaxConstraints() {
            super(new ArrayRealVector(M));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector point) {

            final double[] constraints = new double[M];
            final double u = point.getEntry(U_INDEX);

            for (int i = 0; i < SAMPLE_COUNT; ++i) {
                final double approximation =
                        rationalValue(point, SAMPLE_Y[i]);
                final double target = SAMPLE_EXP[i];

                constraints[2 * i] =
                        u - approximation + target;

                constraints[2 * i + 1] =
                        u + approximation - target;
            }

            return new ArrayRealVector(constraints, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector point) {

            final RealMatrix jacobian =
                    new OpenMapRealMatrix(M, N);

            for (int i = 0; i < SAMPLE_COUNT; ++i) {

                final double y = SAMPLE_Y[i];
                final double y2 = y * y;
                final double y3 = y2 * y;

                final double numerator =
                        point.getEntry(0) +
                        y * point.getEntry(1);

                final double denominator =
                        1.0 +
                        y  * point.getEntry(2) +
                        y2 * point.getEntry(3) +
                        y3 * point.getEntry(4);

                final double denominatorSquared =
                        denominator * denominator;

                final double d0 = 1.0 / denominator;
                final double d1 = y / denominator;
                final double d2 =
                        -numerator * y /
                        denominatorSquared;
                final double d3 =
                        -numerator * y2 /
                        denominatorSquared;
                final double d4 =
                        -numerator * y3 /
                        denominatorSquared;

                final int firstRow = 2 * i;
                final int secondRow = firstRow + 1;

                jacobian.setEntry(firstRow, 0, -d0);
                jacobian.setEntry(firstRow, 1, -d1);
                jacobian.setEntry(firstRow, 2, -d2);
                jacobian.setEntry(firstRow, 3, -d3);
                jacobian.setEntry(firstRow, 4, -d4);
                jacobian.setEntry(firstRow, U_INDEX, 1.0);

                jacobian.setEntry(secondRow, 0, d0);
                jacobian.setEntry(secondRow, 1, d1);
                jacobian.setEntry(secondRow, 2, d2);
                jacobian.setEntry(secondRow, 3, d3);
                jacobian.setEntry(secondRow, 4, d4);
                jacobian.setEntry(secondRow, U_INDEX, 1.0);
            }

            return jacobian;
        }
    }

    /**
     * Evaluate the rational approximation
     *
     * <pre>
     * (x1 + y x2) /
     * (1 + y x3 + y^2 x4 + y^3 x5).
     * </pre>
     *
     * @param point coefficient vector
     * @param y sample abscissa
     * @return rational approximation
     */
    private static double rationalValue(final RealVector point,
                                        final double y) {

        final double y2 = y * y;
        final double y3 = y2 * y;

        final double numerator =
                point.getEntry(0) +
                y * point.getEntry(1);

        final double denominator =
                1.0 +
                y  * point.getEntry(2) +
                y2 * point.getEntry(3) +
                y3 * point.getEntry(4);

        return numerator / denominator;
    }

    @Test
    public void testHALDMADS() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(5000),
                        new InitialGuess(START),
                        new ObjectiveFunction(new Objective()),
                        new MinimaxConstraints(),
                        new SimpleBounds(LOWER, UPPER),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}