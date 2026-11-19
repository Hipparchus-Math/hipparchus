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
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/**
 * Hock-Schittkowski test problem TP364.
 *
 * <p>The problem contains six variables, two linear inequalities and
 * two nonlinear inequalities. All inequalities follow the Hipparchus
 * convention:</p>
 *
 * <pre>
 *     g(x) >= 0
 * </pre>
 *
 * <p>The objective is the root-mean-square fitting error of a planar
 * mechanism evaluated at 31 angular positions.</p>
 */
public class HS364Test {

    /** Problem dimension. */
    private static final int DIM = 6;

    /** Number of sampled angular positions. */
    private static final int SAMPLE_COUNT = 31;

    /** Original Fortran value of pi. */
    private static final double PI = 3.141592654;

    /** First fixed angle used by the constraints. */
    private static final double XMU1 = 0.7853981633;

    /** Second fixed angle used by the constraints. */
    private static final double XMU2 = 2.356194491;

    /** Validated objective value. */
    private static final double EXPECTED_OBJECTIVE = 0.0606002;

    /** TP364 objective function. */
    private static final class HS364Obj
        extends TwiceDifferentiableFunction {

        /** {@inheritDoc} */
        @Override
        public int dim() {
            return DIM;
        }

        /** {@inheritDoc} */
        @Override
        public double value(final RealVector x) {
            return tp364a(x);
        }

        /** {@inheritDoc} */
        @Override
        public RealVector gradient(final RealVector x) {
            /*
             * MODE=3 in the original Fortran routine returns without
             * providing analytical objective derivatives.
             */
            throw new UnsupportedOperationException(
                    "The analytical objective gradient is not provided for TP364.");
        }

        /** {@inheritDoc} */
        @Override
        public RealMatrix hessian(final RealVector x) {
            /*
             * The original Fortran problem does not provide an exact
             * objective Hessian.
             */
            throw new UnsupportedOperationException(
                    "The analytical objective Hessian is not provided for TP364.");
        }
    }

    /**
     * Four inequality constraints:
     *
     * <pre>
     * g1 = -x1 + x2 + x3 - x4
     *
     * g2 = -x1 - x2 + x3 + x4
     *
     * g3 = -x2² - x3² + (x4 - x1)²
     *      + 2 x2 x3 cos(mu1)
     *
     * g4 =  x2² + x3² - (x4 + x1)²
     *      - 2 x2 x3 cos(mu2)
     * </pre>
     *
     * <p>The first two constraints are linear and the final two are
     * nonlinear. They are kept in one generic inequality constraint
     * object, as expected by the SQP solver.</p>
     */
    private static final class HS364Ineq
        extends InequalityConstraint {

        /** Simple constructor. */
        HS364Ineq() {
            super(new ArrayRealVector(4));
        }

        /** {@inheritDoc} */
        @Override
        public int dim() {
            return DIM;
        }

        /** {@inheritDoc} */
        @Override
        public RealVector value(final RealVector x) {

            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);

            final double difference41 = x4 - x1;
            final double sum41 = x4 + x1;

            final double g1 =
                    -x1 + x2 + x3 - x4;

            final double g2 =
                    -x1 - x2 + x3 + x4;

            final double g3 =
                    -x2 * x2 -
                    x3 * x3 +
                    difference41 * difference41 +
                    2.0 * x2 * x3 * FastMath.cos(XMU1);

            final double g4 =
                    x2 * x2 +
                    x3 * x3 -
                    sum41 * sum41 -
                    2.0 * x2 * x3 * FastMath.cos(XMU2);

            return new ArrayRealVector(
                    new double[] {
                        g1,
                        g2,
                        g3,
                        g4
                    },
                    false);
        }

        /** {@inheritDoc} */
        @Override
        public RealMatrix jacobian(final RealVector x) {

            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);

            final double cosine1 = FastMath.cos(XMU1);
            final double cosine2 = FastMath.cos(XMU2);

            final RealMatrix jacobian =
                    MatrixUtils.createRealMatrix(4, DIM);

            /*
             * g1 = -x1 + x2 + x3 - x4
             */
            jacobian.setEntry(0, 0, -1.0);
            jacobian.setEntry(0, 1,  1.0);
            jacobian.setEntry(0, 2,  1.0);
            jacobian.setEntry(0, 3, -1.0);

            /*
             * g2 = -x1 - x2 + x3 + x4
             */
            jacobian.setEntry(1, 0, -1.0);
            jacobian.setEntry(1, 1, -1.0);
            jacobian.setEntry(1, 2,  1.0);
            jacobian.setEntry(1, 3,  1.0);

            /*
             * g3 = -x2² - x3² + (x4 - x1)²
             *      + 2 x2 x3 cos(mu1)
             */
            jacobian.setEntry(2, 0,
                              2.0 * (x1 - x4));

            jacobian.setEntry(2, 1,
                              -2.0 * x2 +
                               2.0 * x3 * cosine1);

            jacobian.setEntry(2, 2,
                              -2.0 * x3 +
                               2.0 * x2 * cosine1);

            jacobian.setEntry(2, 3,
                              2.0 * (x4 - x1));

            /*
             * g4 = x2² + x3² - (x4 + x1)²
             *      - 2 x2 x3 cos(mu2)
             */
            jacobian.setEntry(3, 0,
                              -2.0 * (x4 + x1));

            jacobian.setEntry(3, 1,
                               2.0 * x2 -
                               2.0 * x3 * cosine2);

            jacobian.setEntry(3, 2,
                               2.0 * x3 -
                               2.0 * x2 * cosine2);

            jacobian.setEntry(3, 3,
                              -2.0 * (x4 + x1));

            /*
             * Neither g1, g2, g3 nor g4 depends on x5 or x6.
             * Their Jacobian entries therefore remain zero.
             */

            return jacobian;
        }
    }

    /**
     * Original TP364 starting point.
     *
     * @return initial point
     */
    private static double[] start() {
        return new double[] {
            1.0,
            4.5,
            4.0,
            5.0,
            3.0,
            3.0
        };
    }

    /**
     * Original TP364 bounds.
     *
     * <pre>
     * 0.5  <= x1 <= 3
     * 0.01 <= x2
     * 0    <= x3
     * 2    <= x4 <= 10
     *
     * x5 and x6 are unbounded.
     * </pre>
     *
     * @return variable bounds
     */
    private static SimpleBounds bounds() {

        final double[] lower = {
            0.5,
            0.01,
            0.0,
            2.0,
            Double.NEGATIVE_INFINITY,
            Double.NEGATIVE_INFINITY
        };

        final double[] upper = {
            3.0,
            Double.POSITIVE_INFINITY,
            Double.POSITIVE_INFINITY,
            10.0,
            Double.POSITIVE_INFINITY,
            Double.POSITIVE_INFINITY
        };

        return new SimpleBounds(lower, upper);
    }

    /** Test TP364 using its original starting point and bounds. */
    @Test
    public void testHS364() {

        final SQPOptimizerS2 optimizer =
                HSProblemTestUtils.newOptimizer();

        final LagrangeSolution solution = optimizer.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS364Obj()),
                new HS364Ineq(),
                bounds());

        HSProblemTestUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }

    /**
     * Objective function TP364A.
     *
     * @param x current point
     * @return objective value
     */
    private static double tp364a(final RealVector x) {

        final double[] phi = new double[SAMPLE_COUNT];
        final double[] x1 = new double[SAMPLE_COUNT];
        final double[] y1 = new double[SAMPLE_COUNT];

        final double angularIncrement =
                2.0 * PI / 30.0;

        for (int i = 0; i < SAMPLE_COUNT; ++i) {
            phi[i] = angularIncrement * i;
        }

        tp364b(phi, x1, y1);

        double sum = 0.0;

        for (int i = 0; i < SAMPLE_COUNT; ++i) {

            final double cosineS =
                    tp364c(x, phi[i]);

            /*
             * Exact translation of:
             *
             *     WP = DABS(1.D0 - COSS * COSS)
             *     SINS = DSQRT(WP)
             *
             * The absolute value is intentional and must not be replaced
             * by max(0, 1 - COSS²).
             */
            final double wp =
                    FastMath.abs(1.0 - cosineS * cosineS);

            final double sineS =
                    wp > 0.0 ? FastMath.sqrt(wp) : 0.0;

            final double cosineY =
                    (x.getEntry(3) +
                     x.getEntry(2) * cosineS -
                     x.getEntry(0) * FastMath.cos(phi[i])) /
                    x.getEntry(1);

            final double sineY =
                    (x.getEntry(2) * sineS -
                     x.getEntry(0) * FastMath.sin(phi[i])) /
                    x.getEntry(1);

            final double x1Approx =
                    x.getEntry(0) * FastMath.cos(phi[i]) +
                    x.getEntry(4) * cosineY -
                    x.getEntry(5) * sineY;

            final double y1Approx =
                    x.getEntry(0) * FastMath.sin(phi[i]) +
                    x.getEntry(4) * sineY +
                    x.getEntry(5) * cosineY;

            final double deltaX =
                    x1Approx - x1[i];

            final double deltaY =
                    y1Approx - y1[i];

            sum += deltaX * deltaX +
                   deltaY * deltaY;
        }

        final double mean =
                sum / SAMPLE_COUNT;

        return mean > 0.0 ?
               FastMath.sqrt(mean) :
               0.0;
    }

    /**
     * Generate the target curve, corresponding to TP364B.
     *
     * @param phi angular samples
     * @param x1 target x coordinates
     * @param y1 target y coordinates
     */
    private static void tp364b(final double[] phi,
                               final double[] x1,
                               final double[] y1) {

        for (int i = 0; i < SAMPLE_COUNT; ++i) {

            x1[i] =
                    0.4 +
                    FastMath.sin(
                            2.0 * PI *
                            ((PI - phi[i]) /
                             (2.0 * PI) -
                             0.16));

            y1[i] =
                    2.0 +
                    0.9 * FastMath.sin(PI - phi[i]);
        }
    }

    /**
     * Compute COSS, corresponding to TP364C.
     *
     * @param x current point
     * @param phi angular position
     * @return cosine-like mechanism quantity
     */
    private static double tp364c(final RealVector x,
                                 final double phi) {

        final double x1 = x.getEntry(0);
        final double x2 = x.getEntry(1);
        final double x3 = x.getEntry(2);
        final double x4 = x.getEntry(3);

        final double m =
                2.0 * x1 * x3 * FastMath.sin(phi);

        final double l =
                2.0 * x3 * x4 -
                2.0 * x1 * x3 * FastMath.cos(phi);

        final double k =
                x1 * x1 -
                x2 * x2 +
                x3 * x3 +
                x4 * x4 -
                2.0 * x4 * x1 * FastMath.cos(phi);

        final double a =
                l * l + m * m;

        final double b =
                2.0 * k * l;

        final double c =
                k * k - m * m;

        double term =
                FastMath.sqrt(
                        FastMath.abs(
                                b * b -
                                4.0 * a * c));

        if (PI - phi < 0.0) {
            term = -term;
        }

        return (-b + term) /
               (2.0 * a);
    }
}