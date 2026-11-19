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
 * CUTEst problem {@code DECONVU}.
 *
 * <p>Unconstrained deconvolution problem. The fixed zero coefficients
 * {@code C(-9), ..., C(0)} are eliminated, leaving the 51 free variables
 * used by the CUTEst benchmark: {@code C(1), ..., C(40)} followed by
 * {@code SG(1), ..., SG(11)}.</p>
 */
public class DECONVUTest {

    /** Number of measured trace values. */
    private static final int TRACE_LENGTH = 40;

    /** Number of signal coefficients. */
    private static final int SIGNAL_LENGTH = 11;

    /** Total number of free variables. */
    private static final int N =
            TRACE_LENGTH + SIGNAL_LENGTH;

    /** First index of the SG variables. */
    private static final int SIGNAL_OFFSET =
            TRACE_LENGTH;

    /** Sum-of-squares lower bound and expected optimum. */
    private static final double EXPECTED_OBJECTIVE = 0.0;

    /** Measured trace values TR(1), ..., TR(40). */
    private static final double[] TRACE = {
        0.0,
        0.0,
        0.0016,
        0.0054,
        0.0702,
        0.1876,
        0.3320,
        0.7640,
        0.9320,
        0.8120,
        0.3464,
        0.2064,
        0.0830,
        0.0340,
        0.06179999,
        1.2,
        1.8,
        2.4,
        9.0,
        2.4,
        1.801,
        1.325,
        0.0762,
        0.2104,
        0.2680,
        0.5520,
        0.9960,
        0.3600,
        0.2400,
        0.1510,
        0.0248,
        0.2432,
        0.3602,
        0.4800,
        1.8,
        0.4800,
        0.3600,
        0.2640,
        0.0060,
        0.0060
    };

    /** Official initial SG coefficients. */
    private static final double[] INITIAL_SIGNAL = {
        0.01,
        0.02,
        0.40,
        0.60,
        0.80,
        3.00,
        0.80,
        0.60,
        0.44,
        0.01,
        0.01
    };

    /** Official starting point. */
    private static final double[] START =
            buildStartPoint();

    /** All reduced variables are free. */
    private static final double[] LOWER =
            constantArray(Double.NEGATIVE_INFINITY);

    /** All reduced variables are free. */
    private static final double[] UPPER =
            constantArray(Double.POSITIVE_INFINITY);

    /**
     * Nonlinear least-squares objective.
     *
     * <pre>
     * r_k =
     *     sum(i = 1,...,min(11,k))
     *         SG_i C_(k-i+1)
     *     - TR_k,
     *
     * f = sum(k = 1,...,40) r_k^2.
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

            for (int k = 0; k < TRACE_LENGTH; ++k) {
                final double residual =
                        residual(point, k);
                objective += residual * residual;
            }

            return objective;
        }

        @Override
        public RealVector gradient(final RealVector point) {

            final double[] gradient = new double[N];

            for (int k = 0; k < TRACE_LENGTH; ++k) {

                final double residual =
                        residual(point, k);
                final double factor =
                        2.0 * residual;

                final int maximumSignalIndex =
                        Math.min(SIGNAL_LENGTH - 1, k);

                for (int signalIndex = 0;
                     signalIndex <= maximumSignalIndex;
                     ++signalIndex) {

                    final int coefficientIndex =
                            k - signalIndex;

                    final double coefficient =
                            point.getEntry(coefficientIndex);

                    final double signal =
                            point.getEntry(
                                    SIGNAL_OFFSET +
                                    signalIndex);

                    gradient[coefficientIndex] +=
                            factor * signal;

                    gradient[SIGNAL_OFFSET + signalIndex] +=
                            factor * coefficient;
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

        /**
         * Evaluate one convolution residual.
         *
         * @param point current point
         * @param traceIndex zero-based trace index
         * @return residual
         */
        private static double residual(
                final RealVector point,
                final int traceIndex) {

            double convolution = 0.0;

            final int maximumSignalIndex =
                    Math.min(
                            SIGNAL_LENGTH - 1,
                            traceIndex);

            for (int signalIndex = 0;
                 signalIndex <= maximumSignalIndex;
                 ++signalIndex) {

                final int coefficientIndex =
                        traceIndex - signalIndex;

                convolution +=
                        point.getEntry(
                                SIGNAL_OFFSET +
                                signalIndex) *
                        point.getEntry(
                                coefficientIndex);
            }

            return convolution -
                   TRACE[traceIndex];
        }
    }

    /**
     * Build the official CUTEst starting point.
     *
     * @return starting point
     */
    private static double[] buildStartPoint() {

        final double[] start = new double[N];

        for (int i = 0; i < SIGNAL_LENGTH; ++i) {
            start[SIGNAL_OFFSET + i] =
                    INITIAL_SIGNAL[i];
        }

        return start;
    }

    /**
     * Build a constant array.
     *
     * @param value constant value
     * @return filled array
     */
    private static double[] constantArray(
            final double value) {

        final double[] array = new double[N];

        for (int i = 0; i < N; ++i) {
            array[i] = value;
        }

        return array;
    }

    @Test
    public void testDECONVU() {

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