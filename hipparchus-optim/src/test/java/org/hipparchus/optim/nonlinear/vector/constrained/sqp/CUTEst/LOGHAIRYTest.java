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
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

/**
 * CUTEst problem {@code LOGHAIRY}.
 *
 * <p>This is a difficult two-variable, unconstrained and nonconvex
 * logarithmic transformation of the HAIRY problem. The surface contains
 * many sharp hills and saddle points, while a narrow valley leads to the
 * global minimizer.</p>
 */
@Disabled
public class LOGHAIRYTest {

    /** Number of variables. */
    private static final int N = 2;

    /** Density of the oscillatory HAIRY term. */
    private static final double DENSITY = 7.0;

    /** Weight of the oscillatory term. */
    private static final double HAIR_WEIGHT = 30.0;

    /** Weight of each smoothed cup term. */
    private static final double CUP_WEIGHT = 100.0;

    /** Smoothing parameter of both cup terms. */
    private static final double SMOOTHING = 0.01;

    /** Constant appearing in the logarithmic group function. */
    private static final double LOG_SCALE = 100.0;

    /**
     * The global minimum is attained at x = (0, 0):
     *
     * <pre>
     * g(0, 0) = 100 sqrt(0.01) + 100 sqrt(0.01) = 20
     * f(0, 0) = log((100 + 20) / 100) = log(1.2).
     * </pre>
     */
    private static final double EXPECTED_OBJECTIVE =
            0.1823215567939546;

    /** Official CUTEst starting point. */
    private static final double[] START = {
        -500.0, -700.0
    };

    /** Both variables are free. */
    private static final double[] LOWER = {
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY
    };

    /** Both variables are free. */
    private static final double[] UPPER = {
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY
    };

    /** LOGHAIRY objective with analytic gradient. */
    private static final class Objective
            extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {

            final double x1 = point.getEntry(0);
            final double x2 = point.getEntry(1);

            final double sin1 =
                    FastMath.sin(DENSITY * x1);
            final double cos2 =
                    FastMath.cos(DENSITY * x2);

            final double hairy =
                    sin1 * sin1 * cos2 * cos2;

            final double difference =
                    x1 - x2;

            final double diagonalCup =
                    FastMath.sqrt(
                            SMOOTHING +
                            difference * difference);

            final double oneDimensionalCup =
                    FastMath.sqrt(
                            SMOOTHING +
                            x1 * x1);

            final double groupArgument =
                    HAIR_WEIGHT * hairy +
                    CUP_WEIGHT * diagonalCup +
                    CUP_WEIGHT * oneDimensionalCup;

            return FastMath.log(
                    (LOG_SCALE + groupArgument) /
                    LOG_SCALE);
        }

        @Override
        public RealVector gradient(final RealVector point) {

            final double x1 = point.getEntry(0);
            final double x2 = point.getEntry(1);

            final double densityX1 =
                    DENSITY * x1;
            final double densityX2 =
                    DENSITY * x2;

            final double sin1 =
                    FastMath.sin(densityX1);
            final double cos2 =
                    FastMath.cos(densityX2);

            final double sin1Squared =
                    sin1 * sin1;
            final double cos2Squared =
                    cos2 * cos2;

            final double hairy =
                    sin1Squared * cos2Squared;

            final double difference =
                    x1 - x2;

            final double diagonalCup =
                    FastMath.sqrt(
                            SMOOTHING +
                            difference * difference);

            final double oneDimensionalCup =
                    FastMath.sqrt(
                            SMOOTHING +
                            x1 * x1);

            final double groupArgument =
                    HAIR_WEIGHT * hairy +
                    CUP_WEIGHT * diagonalCup +
                    CUP_WEIGHT * oneDimensionalCup;

            /*
             * Derivatives of
             *
             * sin(7 x1)^2 cos(7 x2)^2.
             */
            final double hairyDerivativeX1 =
                    DENSITY *
                    FastMath.sin(2.0 * densityX1) *
                    cos2Squared;

            final double hairyDerivativeX2 =
                    -DENSITY *
                    sin1Squared *
                    FastMath.sin(2.0 * densityX2);

            final double groupDerivativeX1 =
                    HAIR_WEIGHT * hairyDerivativeX1 +
                    CUP_WEIGHT *
                    difference / diagonalCup +
                    CUP_WEIGHT *
                    x1 / oneDimensionalCup;

            final double groupDerivativeX2 =
                    HAIR_WEIGHT * hairyDerivativeX2 -
                    CUP_WEIGHT *
                    difference / diagonalCup;

            final double denominator =
                    LOG_SCALE + groupArgument;

            return new ArrayRealVector(
                    new double[] {
                        groupDerivativeX1 / denominator,
                        groupDerivativeX2 / denominator
                    },
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

    @Test
    public void testLOGHAIRY() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils
                        .newExternalOption();

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