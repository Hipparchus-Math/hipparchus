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
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

/**
 * CUTEst problem {@code AIRCRFTA}.
 *
 * <p>Aircraft stability equations with fixed elevator, aileron and rudder
 * controls. The eight variables are ordered as roll rate, pitch rate,
 * yaw rate, attack angle, sideslip angle, elevator, aileron and rudder.</p>
 */
public class AIRCRFTATest {

    /** Number of variables. */
    private static final int N = 8;

    /** Number of nonlinear equality equations. */
    private static final int M = 5;

    /** Expected constant objective. */
    private static final double EXPECTED_OBJECTIVE = 0.0;

    /** Official starting point. */
    private static final double[] START = {
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.1, 0.0, 0.0
    };

    /** Bounds fixing the three control variables. */
    private static final double[] LOWER = {
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        0.1, 0.0, 0.0
    };

    /** Bounds fixing the three control variables. */
    private static final double[] UPPER = {
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        0.1, 0.0, 0.0
    };

    /** Constant zero objective. */
    private static final class Objective
            extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {
            return 0.0;
        }

        @Override
        public RealVector gradient(final RealVector point) {
            return new ArrayRealVector(N);
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            return new OpenMapRealMatrix(N, N);
        }
    }

    /** Five aircraft stability equations. */
    private static final class Equalities
            extends EqualityConstraint {

        Equalities() {
            super(new ArrayRealVector(M));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector point) {

            final double p = point.getEntry(0);
            final double q = point.getEntry(1);
            final double r = point.getEntry(2);
            final double alpha = point.getEntry(3);
            final double beta = point.getEntry(4);
            final double elevator = point.getEntry(5);
            final double aileron = point.getEntry(6);
            final double rudder = point.getEntry(7);

            final double[] value = new double[M];

            value[0] =
                    -3.933 * p +
                     0.107 * q +
                     0.126 * r -
                     9.99 * beta -
                    45.83 * aileron -
                     7.64 * rudder -
                     0.727 * q * r +
                     8.39 * r * alpha -
                   684.4 * alpha * beta +
                    63.5 * q * alpha;

            value[1] =
                    -0.987 * q -
                    22.95 * alpha -
                    28.37 * elevator +
                     0.949 * p * r +
                     0.173 * p * beta;

            value[2] =
                     0.002 * p -
                     0.235 * r +
                     5.67 * beta -
                     0.921 * aileron -
                     6.51 * rudder -
                     0.716 * p * q -
                     1.578 * p * alpha +
                     1.132 * q * alpha;

            value[3] =
                     q -
                     alpha -
                     1.168 * elevator -
                     p * beta;

            value[4] =
                    -r -
                     0.196 * beta -
                     0.0071 * aileron +
                     p * alpha;

            return new ArrayRealVector(value, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector point) {

            final double p = point.getEntry(0);
            final double q = point.getEntry(1);
            final double r = point.getEntry(2);
            final double alpha = point.getEntry(3);
            final double beta = point.getEntry(4);

            final RealMatrix jacobian =
                    new OpenMapRealMatrix(M, N);

            jacobian.setEntry(0, 0, -3.933);
            jacobian.setEntry(0, 1,
                    0.107 - 0.727 * r + 63.5 * alpha);
            jacobian.setEntry(0, 2,
                    0.126 - 0.727 * q + 8.39 * alpha);
            jacobian.setEntry(0, 3,
                    8.39 * r - 684.4 * beta + 63.5 * q);
            jacobian.setEntry(0, 4,
                    -9.99 - 684.4 * alpha);
            jacobian.setEntry(0, 6, -45.83);
            jacobian.setEntry(0, 7, -7.64);

            jacobian.setEntry(1, 0,
                    0.949 * r + 0.173 * beta);
            jacobian.setEntry(1, 1, -0.987);
            jacobian.setEntry(1, 2, 0.949 * p);
            jacobian.setEntry(1, 3, -22.95);
            jacobian.setEntry(1, 4, 0.173 * p);
            jacobian.setEntry(1, 5, -28.37);

            jacobian.setEntry(2, 0,
                    0.002 - 0.716 * q - 1.578 * alpha);
            jacobian.setEntry(2, 1,
                    -0.716 * p + 1.132 * alpha);
            jacobian.setEntry(2, 2, -0.235);
            jacobian.setEntry(2, 3,
                    -1.578 * p + 1.132 * q);
            jacobian.setEntry(2, 4, 5.67);
            jacobian.setEntry(2, 6, -0.921);
            jacobian.setEntry(2, 7, -6.51);

            jacobian.setEntry(3, 0, -beta);
            jacobian.setEntry(3, 1, 1.0);
            jacobian.setEntry(3, 3, -1.0);
            jacobian.setEntry(3, 4, -p);
            jacobian.setEntry(3, 5, -1.168);

            jacobian.setEntry(4, 0, alpha);
            jacobian.setEntry(4, 2, -1.0);
            jacobian.setEntry(4, 3, p);
            jacobian.setEntry(4, 4, -0.196);
            jacobian.setEntry(4, 6, -0.0071);

            return jacobian;
        }
    }

    @Test
    public void testAIRCRFTA() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(5000),
                        new InitialGuess(START),
                        new ObjectiveFunction(new Objective()),
                        new Equalities(),
                        new SimpleBounds(LOWER, UPPER),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}