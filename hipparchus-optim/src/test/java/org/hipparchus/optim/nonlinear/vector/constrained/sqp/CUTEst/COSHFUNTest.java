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
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/**
 * CUTEst problem {@code COSHFUN}.
 *
 * <p>This test uses the original instance with 20 minimax functions.
 * There are 60 independent variables, one minimax variable {@code F},
 * and 20 nonlinear inequalities.</p>
 */
public class COSHFUNTest {

    /** Number of minimax component functions. */
    private static final int M = 20;

    /** Number of independent variables. */
    private static final int INDEPENDENT_VARIABLES = 3 * M;

    /** Total number of variables, including the minimax variable F. */
    private static final int N = INDEPENDENT_VARIABLES + 1;

    /** Index of the minimax variable F. */
    private static final int F_INDEX = INDEPENDENT_VARIABLES;

    /** Objective obtained by SNOPT on the original CUTEst instance. */
    private static final double EXPECTED_OBJECTIVE =
            -0.7732665517400054;

    /** Official all-zero starting point. */
    private static final double[] START =
            constantArray(0.0);

    /** All variables are free. */
    private static final double[] LOWER =
            constantArray(Double.NEGATIVE_INFINITY);

    /** All variables are free. */
    private static final double[] UPPER =
            constantArray(Double.POSITIVE_INFINITY);

    /**
     * Linear minimax objective f(x) = F.
     */
    private static final class Objective
            extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {
            return point.getEntry(F_INDEX);
        }

        @Override
        public RealVector gradient(final RealVector point) {

            final double[] gradient = new double[N];
            gradient[F_INDEX] = 1.0;

            return new ArrayRealVector(gradient, false);
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            return new OpenMapRealMatrix(N, N);
        }
    }

    /**
     * Minimax inequalities converted to the Hipparchus convention
     * {@code g(x) >= 0}.
     *
     * <p>For block {@code i}, let</p>
     *
     * <pre>
     * a_i = x(3 i)
     * b_i = x(3 i + 1)
     * c_i = x(3 i + 2).
     * </pre>
     *
     * <p>The nonlinear part is</p>
     *
     * <pre>
     * c_i^2 + cosh(b_i) + 2 a_i^2 c_i.
     * </pre>
     */
    private static final class MinimaxConstraints
            extends InequalityConstraint {

        MinimaxConstraints() {
            super(new ArrayRealVector(M));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector point) {

            final double[] constraints = new double[M];
            final double f = point.getEntry(F_INDEX);

            for (int i = 0; i < M; ++i) {

                final int first = 3 * i;
                final double a = point.getEntry(first);
                final double b = point.getEntry(first + 1);
                final double c = point.getEntry(first + 2);

                double component =
                        c * c +
                        FastMath.cosh(b) +
                        2.0 * a * a * c -
                        2.0 * c;

                if (i == 0) {
                    component -= point.getEntry(5);
                } else {
                    component += point.getEntry(first - 3);

                    if (i < M - 1) {
                        component -= point.getEntry(first + 5);
                    }
                }

                /*
                 * CUTEst defines component - F <= 0.
                 * Hipparchus uses g(x) >= 0.
                 */
                constraints[i] = f - component;
            }

            return new ArrayRealVector(constraints, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector point) {

            final RealMatrix jacobian =
                    new OpenMapRealMatrix(M, N);

            for (int i = 0; i < M; ++i) {

                final int first = 3 * i;
                final double a = point.getEntry(first);
                final double b = point.getEntry(first + 1);
                final double c = point.getEntry(first + 2);

                jacobian.setEntry(
                        i,
                        first,
                        -4.0 * a * c);

                jacobian.setEntry(
                        i,
                        first + 1,
                        -FastMath.sinh(b));

                jacobian.setEntry(
                        i,
                        first + 2,
                        2.0 - 2.0 * c - 2.0 * a * a);

                if (i == 0) {
                    jacobian.setEntry(i, 5, 1.0);
                } else {
                    jacobian.addToEntry(i, first - 3, -1.0);

                    if (i < M - 1) {
                        jacobian.setEntry(i, first + 5, 1.0);
                    }
                }

                jacobian.setEntry(i, F_INDEX, 1.0);
            }

            return jacobian;
        }
    }

    /**
     * Build an array filled with a constant.
     *
     * @param value constant value
     * @return filled array
     */
    private static double[] constantArray(final double value) {

        final double[] array = new double[N];

        for (int i = 0; i < N; ++i) {
            array[i] = value;
        }

        return array;
    }

    @Test
    public void testCOSHFUN() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(5000),
                        new InitialGuess(START),
                        new ObjectiveFunction(new Objective()),
                        new MinimaxConstraints(),
                        new SimpleBounds(LOWER, UPPER),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}