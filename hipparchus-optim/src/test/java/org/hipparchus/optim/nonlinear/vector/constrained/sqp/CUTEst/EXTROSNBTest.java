/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements. See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.CUTEst;

import java.util.Arrays;

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
 * CUTEst problem {@code EXTROSNB}.
 *
 * <p>Extended Rosenbrock function, nonseparable version:</p>
 *
 * <pre>
 * f(x) = (1 - x[0])^2
 *        + 100 * sum((x[i] - x[i - 1]^2)^2, i = 1,...,N-1).
 * </pre>
 *
 * <p>The active SIF instance has 1000 free variables, no constraints and
 * starting point {@code x[i] = -1}. The global minimum is zero at the vector
 * of all ones.</p>
 */
public class EXTROSNBTest {

    /** Active SIF dimension. */
    private static final int N = 1000;

    /** Exact global minimum. */
    private static final double EXPECTED_OBJECTIVE = 0.0;

    /** Extended Rosenbrock nonseparable objective. */
    private static final class Objective
            extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {

            final double firstResidual =
                    1.0 - point.getEntry(0);

            double objective =
                    firstResidual * firstResidual;

            for (int i = 1; i < N; ++i) {

                final double previous =
                        point.getEntry(i - 1);

                final double residual =
                        point.getEntry(i) - previous * previous;

                objective +=
                        100.0 * residual * residual;
            }

            return objective;
        }

        @Override
        public RealVector gradient(final RealVector point) {
            throw new UnsupportedOperationException(
                    "Gradient is evaluated by finite differences.");
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            throw new UnsupportedOperationException(
                    "Hessian is not required by the finite-difference option.");
        }
    }

    /** Official SIF starting point: every variable is -1. */
    private static double[] initialPoint() {
        final double[] initial = new double[N];
        Arrays.fill(initial, -1.0);
        return initial;
    }

    /** EXTROSNB is completely unbounded. */
    private static SimpleBounds bounds() {

        final double[] lower = new double[N];
        final double[] upper = new double[N];

        Arrays.fill(lower, Double.NEGATIVE_INFINITY);
        Arrays.fill(upper, Double.POSITIVE_INFINITY);

        return new SimpleBounds(lower, upper);
    }

    @Test
    public void testEXTROSNB() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(10000),
                        new InitialGuess(initialPoint()),
                        new ObjectiveFunction(new Objective()),
                        bounds(),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}