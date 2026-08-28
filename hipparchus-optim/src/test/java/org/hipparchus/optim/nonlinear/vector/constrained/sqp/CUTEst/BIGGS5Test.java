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
 * CUTEst problem {@code BIGGS5}.
 *
 * <p>Nonlinear exponential least-squares problem with six variables and
 * thirteen residuals. The sixth variable is fixed at {@code 3.0}.</p>
 */
public class BIGGS5Test {

    /** Number of variables. */
    private static final int N = 6;

    /** Number of residuals. */
    private static final int M = 13;

    /** Exact minimum: all residuals vanish simultaneously. */
    private static final double EXPECTED_OBJECTIVE = 0.0;

    /** Official starting point. */
    private static final double[] START = {
        1.0, 2.0, 1.0, 1.0, 4.0, 3.0
    };

    /** Lower bounds; x6 is fixed at 3. */
    private static final double[] LOWER = {
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        3.0
    };

    /** Upper bounds; x6 is fixed at 3. */
    private static final double[] UPPER = {
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        3.0
    };

    /**
     * BIGGS5 objective.
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

            for (int sample = 1; sample <= M; ++sample) {
                final double residual = residual(point, sample);
                objective += residual * residual;
            }

            return objective;
        }

        @Override
        public RealVector gradient(final RealVector point) {

            final double[] gradient = new double[N];

            for (int sample = 1; sample <= M; ++sample) {

                final double t = 0.1 * sample;

                final double e1 =
                        FastMath.exp(-t * point.getEntry(0));
                final double e2 =
                        FastMath.exp(-t * point.getEntry(1));
                final double e5 =
                        FastMath.exp(-t * point.getEntry(4));

                final double residual = residual(point, sample);
                final double factor = 2.0 * residual;

                gradient[0] += factor *
                        (-t * point.getEntry(2) * e1);

                gradient[1] += factor *
                        (t * point.getEntry(3) * e2);

                gradient[2] += factor * e1;
                gradient[3] -= factor * e2;

                gradient[4] += factor *
                        (-t * point.getEntry(5) * e5);

                gradient[5] += factor * e5;
            }

            return new ArrayRealVector(gradient, false);
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            return null;
        }

        /**
         * Evaluate one residual.
         *
         * @param point current point
         * @param sample one-based sample index
         * @return residual
         */
        private static double residual(final RealVector point,
                                       final int sample) {

            final double t = 0.1 * sample;

            final double target =
                    FastMath.exp(-t) -
                    5.0 * FastMath.exp(-sample) +
                    3.0 * FastMath.exp(-0.4 * sample);

            final double model =
                    point.getEntry(2) *
                    FastMath.exp(-t * point.getEntry(0)) -
                    point.getEntry(3) *
                    FastMath.exp(-t * point.getEntry(1)) +
                    point.getEntry(5) *
                    FastMath.exp(-t * point.getEntry(4));

            return model - target;
        }
    }

    @Test
    public void testBIGGS5() {

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