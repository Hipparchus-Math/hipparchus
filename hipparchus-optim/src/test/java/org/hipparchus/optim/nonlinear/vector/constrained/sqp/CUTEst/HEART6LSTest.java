/*
 * Licensed to the Hipparchus project under one or more contributor
 * license agreements. See the NOTICE file distributed with this work
 * for additional information regarding copyright ownership.
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
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.CUTEst;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.OpenMapRealMatrix;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.MaxIter;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

/**
 * CUTEst problem {@code HEART6LS}.
 *
 * <p>Least-squares form of the 6 by 6 dipole model of the heart.
 * The variables are ordered as {@code a, c, t, u, v, w}.</p>
 */
public class HEART6LSTest {

    /** Number of variables and residuals. */
    private static final int N = 6;

    /** Dipole-moment constants. */
    private static final double SUM_MX = -0.816;
    private static final double SUM_MY = -0.017;

    /** Target values for the six nonlinear equations. */
    private static final double[] TARGET = {
        -1.826,
        -0.754,
        -4.839,
        -3.259,
        -14.023,
        15.467
    };

    /** Exact least-squares lower bound reported by the SIF. */
    private static final double EXPECTED_OBJECTIVE = 0.0;

    /** Official starting point: a = c = 0 and t = u = v = w = 1. */
    private static final double[] START = {
        0.0, 0.0, 1.0, 1.0, 1.0, 1.0
    };

    /** All variables are free. */
    private static final double[] LOWER = {
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY
    };

    /** All variables are free. */
    private static final double[] UPPER = {
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY
    };

    /** Nonlinear least-squares objective. */
    private static final class Objective
            extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {

            final double[] residuals = residuals(point);

            double objective = 0.0;
            for (int i = 0; i < N; ++i) {
                objective += residuals[i] * residuals[i];
            }

            return objective;
        }

        @Override
        public RealVector gradient(final RealVector point) {

            final double[] residuals = residuals(point);
            final RealMatrix jacobian = residualJacobian(point);
            final double[] gradient = new double[N];

            for (int row = 0; row < N; ++row) {
                final double factor = 2.0 * residuals[row];

                for (int column = 0; column < N; ++column) {
                    gradient[column] +=
                            factor *
                            jacobian.getEntry(row, column);
                }
            }

            return new ArrayRealVector(gradient, false);
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            /*
             * SQPOptimizerS2 uses its BFGS approximation in this test.
             */
            return new OpenMapRealMatrix(N, N);
        }
    }

    /**
     * Evaluate the six HEART6 residuals.
     *
     * @param point current point
     * @return residual vector
     */
    private static double[] residuals(final RealVector point) {

        final double a = point.getEntry(0);
        final double c = point.getEntry(1);
        final double t = point.getEntry(2);
        final double u = point.getEntry(3);
        final double v = point.getEntry(4);
        final double w = point.getEntry(5);

        final double t2 = t * t;
        final double u2 = u * u;
        final double v2 = v * v;
        final double w2 = w * w;

        final double mxMinusA = SUM_MX - a;
        final double myMinusC = SUM_MY - c;

        return new double[] {
            a * t +
            u * mxMinusA -
            c * v -
            w * myMinusC -
            TARGET[0],

            a * v +
            w * mxMinusA +
            c * t +
            u * myMinusC -
            TARGET[1],

            a * (t2 - v2) -
            2.0 * c * t * v +
            mxMinusA * (u2 - w2) -
            2.0 * myMinusC * u * w -
            TARGET[2],

            c * (t2 - v2) +
            2.0 * a * t * v +
            myMinusC * (u2 - w2) +
            2.0 * mxMinusA * u * w -
            TARGET[3],

            a * t * (t2 - 3.0 * v2) +
            c * v * (v2 - 3.0 * t2) +
            mxMinusA * u * (u2 - 3.0 * w2) +
            myMinusC * w * (w2 - 3.0 * u2) -
            TARGET[4],

            c * t * (t2 - 3.0 * v2) -
            a * v * (v2 - 3.0 * t2) +
            myMinusC * u * (u2 - 3.0 * w2) -
            mxMinusA * w * (w2 - 3.0 * u2) -
            TARGET[5]
        };
    }

    /**
     * Evaluate the Jacobian of the six residuals.
     *
     * @param point current point
     * @return residual Jacobian
     */
    private static RealMatrix residualJacobian(
            final RealVector point) {

        final double a = point.getEntry(0);
        final double c = point.getEntry(1);
        final double t = point.getEntry(2);
        final double u = point.getEntry(3);
        final double v = point.getEntry(4);
        final double w = point.getEntry(5);

        final double t2 = t * t;
        final double u2 = u * u;
        final double v2 = v * v;
        final double w2 = w * w;

        final double mxMinusA = SUM_MX - a;
        final double myMinusC = SUM_MY - c;

        final RealMatrix jacobian =
                new OpenMapRealMatrix(N, N);

        jacobian.setEntry(0, 0, t - u);
        jacobian.setEntry(0, 1, -v + w);
        jacobian.setEntry(0, 2, a);
        jacobian.setEntry(0, 3, mxMinusA);
        jacobian.setEntry(0, 4, -c);
        jacobian.setEntry(0, 5, c - SUM_MY);

        jacobian.setEntry(1, 0, v - w);
        jacobian.setEntry(1, 1, t - u);
        jacobian.setEntry(1, 2, c);
        jacobian.setEntry(1, 3, myMinusC);
        jacobian.setEntry(1, 4, a);
        jacobian.setEntry(1, 5, mxMinusA);

        jacobian.setEntry(
                2, 0,
                t2 - v2 - u2 + w2);
        jacobian.setEntry(
                2, 1,
                -2.0 * t * v + 2.0 * u * w);
        jacobian.setEntry(
                2, 2,
                2.0 * a * t - 2.0 * c * v);
        jacobian.setEntry(
                2, 3,
                2.0 * mxMinusA * u -
                2.0 * myMinusC * w);
        jacobian.setEntry(
                2, 4,
                -2.0 * a * v - 2.0 * c * t);
        jacobian.setEntry(
                2, 5,
                -2.0 * mxMinusA * w -
                2.0 * myMinusC * u);

        jacobian.setEntry(
                3, 0,
                2.0 * t * v - 2.0 * u * w);
        jacobian.setEntry(
                3, 1,
                t2 - v2 - u2 + w2);
        jacobian.setEntry(
                3, 2,
                2.0 * c * t + 2.0 * a * v);
        jacobian.setEntry(
                3, 3,
                2.0 * myMinusC * u +
                2.0 * mxMinusA * w);
        jacobian.setEntry(
                3, 4,
                -2.0 * c * v + 2.0 * a * t);
        jacobian.setEntry(
                3, 5,
                -2.0 * myMinusC * w +
                2.0 * mxMinusA * u);

        jacobian.setEntry(
                4, 0,
                t * (t2 - 3.0 * v2) -
                u * (u2 - 3.0 * w2));
        jacobian.setEntry(
                4, 1,
                v * (v2 - 3.0 * t2) -
                w * (w2 - 3.0 * u2));
        jacobian.setEntry(
                4, 2,
                3.0 * a * (t2 - v2) -
                6.0 * c * v * t);
        jacobian.setEntry(
                4, 3,
                3.0 * mxMinusA * (u2 - w2) -
                6.0 * myMinusC * u * w);
        jacobian.setEntry(
                4, 4,
                -6.0 * a * t * v +
                3.0 * c * (v2 - t2));
        jacobian.setEntry(
                4, 5,
                -6.0 * mxMinusA * u * w +
                3.0 * myMinusC * (w2 - u2));

        jacobian.setEntry(
                5, 0,
                -v * (v2 - 3.0 * t2) +
                w * (w2 - 3.0 * u2));
        jacobian.setEntry(
                5, 1,
                t * (t2 - 3.0 * v2) -
                u * (u2 - 3.0 * w2));
        jacobian.setEntry(
                5, 2,
                3.0 * c * (t2 - v2) +
                6.0 * a * v * t);
        jacobian.setEntry(
                5, 3,
                3.0 * myMinusC * (u2 - w2) +
                6.0 * mxMinusA * u * w);
        jacobian.setEntry(
                5, 4,
                -6.0 * c * t * v +
                3.0 * a * (t2 - v2));
        jacobian.setEntry(
                5, 5,
                -6.0 * myMinusC * u * w +
                3.0 * mxMinusA * (u2 - w2));

        return jacobian;
    }

    @Test
    public void testHEART6LS() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();
        option.setMaxIteration(3000);
        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(5000),
                        new InitialGuess(START),
                        new ObjectiveFunction(new Objective()),
                        new SimpleBounds(LOWER, UPPER),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}
