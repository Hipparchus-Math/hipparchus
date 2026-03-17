/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

/** HS201 (TP201): unconstrained quadratic. */
public class HS201Test {

    /** f(x) = 4 (x1-5)^2 + (x2-6)^2. */
    static final class HS201Objective extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }

        @Override
        public double value(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            return 4.0 * (x1 - 5.0) * (x1 - 5.0) + (x2 - 6.0) * (x2 - 6.0);
        }

        @Override
        public RealVector gradient(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            // ∂f/∂x1 = 8 (x1-5), ∂f/∂x2 = 2 (x2-6)
            return new ArrayRealVector(new double[] {
                8.0 * (x1 - 5.0),
                2.0 * (x2 - 6.0)
            }, false);
        }

        @Override
        public org.hipparchus.linear.RealMatrix hessian(final RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    /** Solve utility (no constraints, only bounds ±∞). */
    static LagrangeSolution solve(final double[] start) {
        final double SUP = Double.POSITIVE_INFINITY;
        final double INF = Double.NEGATIVE_INFINITY;

        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        return optimizer.optimize(
            new InitialGuess(start),
            new ObjectiveFunction(new HS201Objective()),
            new SimpleBounds(new double[] { INF, INF }, new double[] { SUP, SUP })
        );
    }

    // ---------- Minimal JUnit test: check only objective value ----------
    @Test
    public void testHS201() {
        final double[] x0 = { 8.0, 9.0 };     // Fortran start
        final LagrangeSolution sol = solve(x0);

        final double fEx = 0.0;               // reference optimum
        final double f   = sol.getValue();
        assertEquals(fEx, f, 1.0e-6 * (fEx + 1.0), "objective mismatch");
    }
}
