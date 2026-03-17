/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
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

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * HS210 (TP210) – Scaled Rosenbrock variant.
 *
 * N = 2
 * No constraints.
 *
 * Objective:
 *   f(x) = (C*(x2 - x1^2)^2 + (1 - x1)^2) / C
 * with C = 1e6.
 *
 * Gradient:
 *   df/dx1 = (-4*C*(x2 - x1^2)*x1 - 2*(1 - x1)) / C
 *   df/dx2 =  2*C*(x2 - x1^2) / C
 *
 * Minimum:
 *   x* = (1, 1)
 *   f* = 0
 */
public class HS210Test {

    private static final int DIM = 2;
    private static final double C = 1_000_000.0;

    // -----------------------------------------
    // Objective function
    // -----------------------------------------
    private static class HS210Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() { return DIM; }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            double t = x2 - x1 * x1;
            return (C * t * t + (1.0 - x1) * (1.0 - x1)) / C;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            double t = x2 - x1 * x1;

            double g1 = (-4.0 * C * t * x1 - 2.0 * (1.0 - x1)) / C;
            double g2 = ( 2.0 * C * t ) / C;

            return new ArrayRealVector(new double[]{ g1, g2 }, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Start with zero Hessian; BFGS updates it.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -----------------------------------------
    // Test
    // -----------------------------------------
    @Test
    public void testHS210_optimization() {

        // Initial guess (from Fortran mode 1)
        double[] x0 = new double[DIM];
        x0[0] = -1.2;
        x0[1] =  1.0;

        // Bounds: Fortran has no bounds → use infinite
        double[] lower = { Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY };
        double[] upper = { Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY };
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS210Obj()),
                null,        // no eq constraints
                null,        // no inequality constraints
                bounds
        );

       
        double f = sol.getValue();

        // Expected minimum
        double[] expected = { 1.0, 1.0 };
        double fExpected = 0.0;
        final double tol = 1.0e-5 * (FastMath.abs(fExpected) + 1.0);
       
        assertEquals(fExpected, f, tol);
    }
}
