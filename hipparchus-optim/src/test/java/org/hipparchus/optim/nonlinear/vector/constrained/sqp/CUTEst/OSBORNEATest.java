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
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/**
 * CUTEst problem {@code OSBORNEA}.
 *
 * <p>Osborne's first nonlinear least-squares problem. The model contains
 * five variables and is fitted to thirty-three observations.</p>
 */
public class OSBORNEATest {

    /** Number of variables. */
    private static final int N = 5;

    /** Number of observations. */
    private static final int M = 33;

    /** Expected objective reported by the SIF. */
    private static final double EXPECTED_OBJECTIVE = 5.46489e-5;

    /** Official starting point. */
    private static final double[] START = {
        0.5, 1.5, -1.0, 0.01, 0.02
    };

    /** All variables are free. */
    private static final double[] LOWER = {
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
        Double.POSITIVE_INFINITY
    };

    /** Observed data from the SIF file. */
    private static final double[] OBSERVATIONS = {
        0.844, 0.908, 0.932, 0.936, 0.925, 0.908, 0.881,
        0.850, 0.818, 0.784, 0.751, 0.718, 0.685, 0.658,
        0.628, 0.603, 0.580, 0.558, 0.538, 0.522, 0.506,
        0.490, 0.478, 0.467, 0.457, 0.448, 0.438, 0.431,
        0.424, 0.420, 0.414, 0.411, 0.406
    };

    /**
     * Nonlinear least-squares objective.
     *
     * <pre>
     * r_i = x1
     *       + x2 exp(-10 i x4)
     *       + x3 exp(-10 i x5)
     *       - y_i,
     *
     * f(x) = sum(i = 0,...,32) r_i^2.
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

            double objective = 0.0;

            for (int i = 0; i < M; ++i) {
                final double residual = residual(point, i);
                objective += residual * residual;
            }

            return objective;
        }

        @Override
        public RealVector gradient(final RealVector point) {

            final double x2 = point.getEntry(1);
            final double x3 = point.getEntry(2);
            final double x4 = point.getEntry(3);
            final double x5 = point.getEntry(4);

            final double[] gradient = new double[N];

            for (int i = 0; i < M; ++i) {

                final double t = -10.0 * i;
                final double exp4 = FastMath.exp(t * x4);
                final double exp5 = FastMath.exp(t * x5);
                final double residual = residual(point, i);
                final double scale = 2.0 * residual;

                gradient[0] += scale;
                gradient[1] += scale * exp4;
                gradient[2] += scale * exp5;
                gradient[3] += scale * t * x2 * exp4;
                gradient[4] += scale * t * x3 * exp5;
            }

            return new ArrayRealVector(gradient, false);
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            /*
             * SQPOptimizerS2 uses its BFGS approximation in this test.
             * The exact Hessian is therefore intentionally not supplied.
             */
            return new OpenMapRealMatrix(N, N);
        }

        /**
         * Evaluate one least-squares residual.
         *
         * @param point current point
         * @param observationIndex zero-based observation index
         * @return residual
         */
        private static double residual(final RealVector point,
                                       final int observationIndex) {

            final double t = -10.0 * observationIndex;

            return point.getEntry(0) +
                   point.getEntry(1) *
                       FastMath.exp(t * point.getEntry(3)) +
                   point.getEntry(2) *
                       FastMath.exp(t * point.getEntry(4)) -
                   OBSERVATIONS[observationIndex];
        }
    }

    @Test
    public void testOSBORNEA() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

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