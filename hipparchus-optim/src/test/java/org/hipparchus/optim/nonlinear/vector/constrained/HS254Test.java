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

/*
 * HS254 (TP254)
 *
 * N    = 3
 * NENL = 2  (two nonlinear equalities)
 *
 * Objective:
 *   f(x) = log10(x3) - x2
 *
 * Equalities (Fortran G(i) used as equality constraints):
 *   h1(x) = x2^2 + x3^2 - 4         = 0
 *   h2(x) = x3 - 1 - x1^2           = 0
 *
 * Bounds:
 *   x1 free
 *   x2 free
 *   x3 >= 1
 *
 * Initial guess (MODE=1):
 *   x = (1, 1, 1)
 *
 * Reference solution:
 *   x* = (0, sqrt(3), 1)
 *   f* = -sqrt(3)
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

public class HS254Test {

    private static final int DIM     = 3;
    private static final int NUM_EQ  = 2;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS254Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {

            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            // f(x) = log10(x3) - x2
            return FastMath.log10(x3) - x2;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x3 = x.getEntry(2);

            double df1 = 0.0;                              // ∂f/∂x1
            double df2 = -1.0;                             // ∂f/∂x2
            double df3 = 1.0 / (x3 * FastMath.log(10.0));  // ∂f/∂x3

            return new ArrayRealVector(new double[]{df1, df2, df3}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {

            double x3 = x.getEntry(2);

            RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);

            // Second derivatives:
            // ∂²f/∂x3² = -1 / (x3^2 * ln(10))
            H.setEntry(2, 2, -1.0 / (x3 * x3 * FastMath.log(10.0)));

            // others are zero
            return H;
        }
    }

    // -------------------------------------------------------------------------
    // Equality constraints h(x) = 0
    //   h1(x) = x2^2 + x3^2 - 4
    //   h2(x) = x3 - 1 - x1^2
    // -------------------------------------------------------------------------
    private static class HS254Eq extends EqualityConstraint {

        HS254Eq() {
            super(new ArrayRealVector(new double[NUM_EQ])); // RHS = 0
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double h1 = x2 * x2 + x3 * x3 - 4.0;
            double h2 = x3 - 1.0 - x1 * x1;

            return new ArrayRealVector(new double[]{h1, h2}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            RealMatrix J = new Array2DRowRealMatrix(NUM_EQ, DIM);

            // h1 = x2^2 + x3^2 - 4
            // ∂h1/∂x1 = 0
            // ∂h1/∂x2 = 2*x2
            // ∂h1/∂x3 = 2*x3
            J.setEntry(0, 0, 0.0);
            J.setEntry(0, 1, 2.0 * x2);
            J.setEntry(0, 2, 2.0 * x3);

            // h2 = x3 - 1 - x1^2
            // ∂h2/∂x1 = -2*x1
            // ∂h2/∂x2 = 0
            // ∂h2/∂x3 = 1
            J.setEntry(1, 0, -2.0 * x1);
            J.setEntry(1, 1, 0.0);
            J.setEntry(1, 2, 1.0);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS254() {

        // Initial guess (MODE=1): x = (1,1,1)
        double[] x0 = new double[]{1.0, 1.0, 1.0};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // Bounds:
        // x1 free, x2 free, x3 ≥ 1
        SimpleBounds bounds = new SimpleBounds(
                new double[]{
                        Double.NEGATIVE_INFINITY, // x1 lower
                        Double.NEGATIVE_INFINITY, // x2 lower
                        1.0                       // x3 lower
                },
                new double[]{
                        Double.POSITIVE_INFINITY, // x1 upper
                        Double.POSITIVE_INFINITY, // x2 upper
                        Double.POSITIVE_INFINITY  // x3 upper
                }
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS254Obj()),
                new HS254Eq(),   // 2 equalities
                null,            // no inequalities
                bounds
        );

        double f = sol.getValue();

        double fExpected = -FastMath.sqrt(3.0);
        double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
