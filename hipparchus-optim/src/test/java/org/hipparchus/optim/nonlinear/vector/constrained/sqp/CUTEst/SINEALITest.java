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
 * CUTEst problem {@code SINEALI}.
 *
 * <p>This problem is a sine-based variation of the extended Rosenbrock
 * function. Bounds are imposed to exclude the additional minima introduced
 * by the periodic sine function.</p>
 *
 * <p>The S2MPJ test instance uses {@code N = 10} variables.</p>
 */
public class SINEALITest {

    /** Number of variables selected by the SIF test instance. */
    private static final int N = 10;

    /** Value of pi used literally by the SIF file. */
    private static final double PI = 3.1415926535;

    /** Scale factor produced by the SIF group scale 0.01. */
    private static final double CHAIN_WEIGHT = 100.0;

    /** Expected objective reported by the SIF for N = 10. */
    private static final double EXPECTED_OBJECTIVE = -901.0;

    /** Official starting point: all variables are zero. */
    private static final double[] START = new double[N];

    /** Variable lower bounds. */
    private static final double[] LOWER = buildLowerBounds();

    /** Variable upper bounds. */
    private static final double[] UPPER = buildUpperBounds();

    /**
     * SINEALI objective
     *
     * <pre>
     * sin(x1 - 1)
     * + 100 sum(i=2,...,N) sin(xi - x(i-1)^2).
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

            double objective =
                    FastMath.sin(point.getEntry(0) - 1.0);

            for (int i = 1; i < N; ++i) {
                final double previous =
                        point.getEntry(i - 1);
                final double argument =
                        point.getEntry(i) -
                        previous * previous;

                objective +=
                        CHAIN_WEIGHT *
                        FastMath.sin(argument);
            }

            return objective;
        }

        @Override
        public RealVector gradient(final RealVector point) {

            final double[] gradient = new double[N];

            gradient[0] =
                    FastMath.cos(point.getEntry(0) - 1.0);

            for (int i = 1; i < N; ++i) {

                final double previous =
                        point.getEntry(i - 1);
                final double argument =
                        point.getEntry(i) -
                        previous * previous;
                final double cosine =
                        CHAIN_WEIGHT *
                        FastMath.cos(argument);

                gradient[i] += cosine;
                gradient[i - 1] -=
                        2.0 * previous * cosine;
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
    }

    /**
     * Build the recursive upper bounds exactly as specified by the SIF.
     *
     * @return upper bounds
     */
    private static double[] buildUpperBounds() {

        final double[] upper = new double[N];

        upper[0] = 0.5 * PI;

        for (int i = 1; i < N; ++i) {
            upper[i] =
                    FastMath.sqrt(upper[i - 1] +
                                  0.5 * PI);
        }

        return upper;
    }

    /**
     * Build the lower bounds from the corresponding upper bounds.
     *
     * @return lower bounds
     */
    private static double[] buildLowerBounds() {

        final double[] upper = buildUpperBounds();
        final double[] lower = new double[N];

        for (int i = 0; i < N; ++i) {
            lower[i] =
                    upper[i] - 2.0 * PI;
        }

        return lower;
    }

    @Test
    public void testSINEALI() {

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