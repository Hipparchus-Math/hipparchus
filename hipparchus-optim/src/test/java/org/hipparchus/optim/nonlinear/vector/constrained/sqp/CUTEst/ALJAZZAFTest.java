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
 * CUTEst problem {@code ALJAZZAF}.
 *
 * <p>This test uses the original three-variable instance from the SIF,
 * corresponding to {@code N = 3} and {@code N1 = 2}.</p>
 */
public class ALJAZZAFTest {

    /** Number of variables in the original instance. */
    private static final int N = 3;

    /** Expected objective reported by the SIF. */
    private static final double EXPECTED_OBJECTIVE = 75.004996;

    /** Official starting point. */
    private static final double[] START = {
        0.0, 0.0, 0.0
    };

    /** All variables are free. */
    private static final double[] LOWER = {
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY,
        Double.NEGATIVE_INFINITY
    };

    /** All variables are free. */
    private static final double[] UPPER = {
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY
    };

    /*
     * For N = 3:
     *
     * F  = (100^2 - 1) / (N - 1) = 4999.5
     * F2 = F / 100                 = 49.995
     *
     * A1 = 100
     * A2 = 50.005
     * A3 = 0.01
     *
     * B1 = 1
     * B2 = 5000.5
     * B3 = 10000
     */
    private static final double A1 = 100.0;
    private static final double A2 = 50.005;
    private static final double A3 = 0.01;

    private static final double B1 = 1.0;
    private static final double B2 = 5000.5;
    private static final double B3 = 10000.0;

    /** Quadratic objective. */
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
            final double x3 = point.getEntry(2);

            final double d1 = x1 - 0.5;
            final double d2 = x2 + 1.0;
            final double d3 = x3 - 1.0;

            return A1 * d1 * d1 +
                   A2 * d2 * d2 +
                   A3 * d3 * d3;
        }

        @Override
        public RealVector gradient(final RealVector point) {

            return new ArrayRealVector(
                    new double[] {
                        2.0 * A1 * (point.getEntry(0) - 0.5),
                        2.0 * A2 * (point.getEntry(1) + 1.0),
                        2.0 * A3 * (point.getEntry(2) - 1.0)
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

    /** Single nonlinear equality constraint. */
    private static final class Equality
            extends EqualityConstraint {

        Equality() {
            super(new ArrayRealVector(1));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector point) {

            final double x1 = point.getEntry(0);
            final double x2 = point.getEntry(1);
            final double x3MinusOne =
                    point.getEntry(2) - 1.0;

            final double value =
                    B1 -
                    B1 * x1 +
                    B2 * x2 * x2 +
                    B3 * x3MinusOne * x3MinusOne;

            return new ArrayRealVector(
                    new double[] { value },
                    false);
        }

        @Override
        public RealMatrix jacobian(final RealVector point) {

            final RealMatrix jacobian =
                    new OpenMapRealMatrix(1, N);

            jacobian.setEntry(0, 0, -B1);
            jacobian.setEntry(
                    0,
                    1,
                    2.0 * B2 * point.getEntry(1));
            jacobian.setEntry(
                    0,
                    2,
                    2.0 * B3 *
                    (point.getEntry(2) - 1.0));

            return jacobian;
        }
    }

    @Test
    public void testALJAZZAF() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(5000),
                        new InitialGuess(START),
                        new ObjectiveFunction(new Objective()),
                        new Equality(),
                        new SimpleBounds(LOWER, UPPER),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}