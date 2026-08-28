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
 * CUTEst problem {@code BDQRTIC}, using the original instance {@code N = 100}.
 *
 * <p>The problem is unconstrained, smooth and globally convex:</p>
 *
 * <pre>
 * f(x) = sum(i = 1,...,N-4) [
 *            (-4 x(i) + 3)^2
 *          + (x(i)^2
 *             + 2 x(i+1)^2
 *             + 3 x(i+2)^2
 *             + 4 x(i+3)^2
 *             + 5 x(N)^2)^2
 *        ].
 * </pre>
 *
 * <p>The official starting point is {@code x(i) = 1} for every variable.
 * The objective gradient is supplied analytically, while the objective
 * Hessian is intentionally omitted so that {@link SQPOptimizerS2} exercises
 * its BFGS update.</p>
 */
public class BDQRTICTest {

    /** Original CUTEst problem dimension. */
    private static final int N = 100;

    /** Number of partially separable groups. */
    private static final int NUMBER_OF_GROUPS = N - 4;

    /**
     * Objective at the global minimizer for the N = 100 instance.
     *
     * <p>This value was independently refined until the infinity norm of the
     * analytic gradient was approximately 5.3e-15.</p>
     */
    private static final double EXPECTED_OBJECTIVE =
            378.7691918086844;

    /** Official starting point: all variables equal to one. */
    private static final double[] START =
            constantArray(1.0);

    /** All variables are free. */
    private static final double[] LOWER =
            constantArray(Double.NEGATIVE_INFINITY);

    /** All variables are free. */
    private static final double[] UPPER =
            constantArray(Double.POSITIVE_INFINITY);

    /** BDQRTIC objective with analytic gradient. */
    private static final class Objective
            extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {

            final double last =
                    point.getEntry(N - 1);
            final double lastSquared =
                    last * last;

            double objective = 0.0;

            for (int i = 0; i < NUMBER_OF_GROUPS; ++i) {

                final double xi =
                        point.getEntry(i);
                final double xi1 =
                        point.getEntry(i + 1);
                final double xi2 =
                        point.getEntry(i + 2);
                final double xi3 =
                        point.getEntry(i + 3);

                final double linearResidual =
                        -4.0 * xi + 3.0;

                final double quarticResidual =
                        xi * xi +
                        2.0 * xi1 * xi1 +
                        3.0 * xi2 * xi2 +
                        4.0 * xi3 * xi3 +
                        5.0 * lastSquared;

                objective +=
                        linearResidual * linearResidual +
                        quarticResidual * quarticResidual;
            }

            return objective;
        }

        @Override
        public RealVector gradient(final RealVector point) {

            final double[] gradient =
                    new double[N];

            final double last =
                    point.getEntry(N - 1);
            final double lastSquared =
                    last * last;

            for (int i = 0; i < NUMBER_OF_GROUPS; ++i) {

                final double xi =
                        point.getEntry(i);
                final double xi1 =
                        point.getEntry(i + 1);
                final double xi2 =
                        point.getEntry(i + 2);
                final double xi3 =
                        point.getEntry(i + 3);

                final double linearResidual =
                        -4.0 * xi + 3.0;

                final double quarticResidual =
                        xi * xi +
                        2.0 * xi1 * xi1 +
                        3.0 * xi2 * xi2 +
                        4.0 * xi3 * xi3 +
                        5.0 * lastSquared;

                /*
                 * Derivative of (-4 xi + 3)^2.
                 */
                gradient[i] +=
                        -8.0 * linearResidual;

                /*
                 * Derivatives of quarticResidual^2.
                 */
                gradient[i] +=
                        4.0 * quarticResidual * xi;

                gradient[i + 1] +=
                        8.0 * quarticResidual * xi1;

                gradient[i + 2] +=
                        12.0 * quarticResidual * xi2;

                gradient[i + 3] +=
                        16.0 * quarticResidual * xi3;

                gradient[N - 1] +=
                        20.0 * quarticResidual * last;
            }

            return new ArrayRealVector(
                    gradient,
                    false);
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
     * Build an array filled with one constant value.
     *
     * @param value constant value
     * @return filled array
     */
    private static double[] constantArray(
            final double value) {

        final double[] array =
                new double[N];

        for (int i = 0; i < N; ++i) {
            array[i] = value;
        }

        return array;
    }

    @Test
    public void testBDQRTIC() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils
                        .newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(5000),
                        new InitialGuess(START),
                        new ObjectiveFunction(
                                new Objective()),
                        new SimpleBounds(
                                LOWER,
                                UPPER),
                        option);

        CUTEstProblemUtils
                .assertExpectedObjective(
                        EXPECTED_OBJECTIVE,
                        solution);
    }
}