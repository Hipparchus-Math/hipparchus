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

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

/**
 * HS209 / TP209 — Unconstrained Rosenbrock with C = 10000.
 * Fortran reference:
 *   f(x) = C*(x2 - x1^2)^2 + (1 - x1)^2,  C = 10000
 *   x0 = (-1.2, 1)
 *   x* = (1, 1),  f* = 0
 */
public class HS209Test {

    /** Objective f(x) = C*(x2 - x1^2)^2 + (1 - x1)^2. */
    private static class HS209Obj extends TwiceDifferentiableFunction {
        private static final double C = 10000.0;
        @Override public int dim() { return 2; }
        @Override public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            return C * FastMath.pow(x2 - x1 * x1, 2) + FastMath.pow(1.0 - x1, 2);
        }
        // Grad/Hess not required for the test — left unimplemented.
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testHS209_Unconstrained() {
        // Initial guess as in Fortran
        InitialGuess guess = new InitialGuess(new double[]{ -1.2, 1.0 });

        // Optimizer instance
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        // Unconstrained optimize: only objective
        LagrangeSolution sol = optimizer.optimize(
                guess,
                new ObjectiveFunction(new HS209Obj())
        );

        // Expect optimum at (1,1) and f* = 0
        HSProblemTestUtils.assertExpectedObjective(0.0, sol);
        
    }
}
