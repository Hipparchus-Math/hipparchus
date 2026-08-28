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

import static org.junit.jupiter.api.Assertions.assertTrue;

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
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

/**
 * CUTEst problem {@code HIMMELBD}.
 *
 * <p>A two-variable nonlinear equation problem by Himmelblau. CUTEst
 * represents it as a feasibility problem with a constant zero objective
 * and two quadratic equality constraints.</p>
 */
@Disabled
public class HIMMELBDTest {

    /** Number of variables and equality constraints. */
    private static final int N = 2;

    /** Constant objective value. */
    private static final double EXPECTED_OBJECTIVE = 0.0;

    /** Official CUTEst starting point. */
    private static final double[] START = {
        1.0, 1.0
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

    /** Two quadratic equality constraints. */
    private static final class Equalities
            extends EqualityConstraint {

        Equalities() {
            super(new ArrayRealVector(N));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector point) {

            final double x1 = point.getEntry(0);
            final double x2 = point.getEntry(1);

            final double[] value = {
                x1 * x1 +
                12.0 * x2 -
                1.0,

                49.0 * x1 * x1 +
                49.0 * x2 * x2 +
                84.0 * x1 +
                2324.0 * x2 -
                681.0
            };

            return new ArrayRealVector(value, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector point) {

            final double x1 = point.getEntry(0);
            final double x2 = point.getEntry(1);

            final RealMatrix jacobian =
                    new OpenMapRealMatrix(N, N);

            jacobian.setEntry(0, 0, 2.0 * x1);
            jacobian.setEntry(0, 1, 12.0);

            jacobian.setEntry(
                    1, 0,
                    98.0 * x1 + 84.0);

            jacobian.setEntry(
                    1, 1,
                    98.0 * x2 + 2324.0);

            return jacobian;
        }
    }

    @Test
    public void testHIMMELBD() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final Equalities equalities = new Equalities();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(5000),
                        new InitialGuess(START),
                        new ObjectiveFunction(new Objective()),
                        equalities,
                        new SimpleBounds(LOWER, UPPER),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);

        /*
         * The objective is identically zero, so feasibility must be checked
         * explicitly; an objective-only assertion would be meaningless.
         */
        final RealVector residual =
                equalities.value(solution.getX());

        double maximumViolation = 0.0;

        for (int i = 0; i < residual.getDimension(); ++i) {
            maximumViolation =
                    FastMath.max(
                            maximumViolation,
                            FastMath.abs(residual.getEntry(i)));
        }

        assertTrue(
                maximumViolation < 1.0e-2,
                "maximum equality violation: " +
                maximumViolation);
    }
}