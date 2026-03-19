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
 * HS56 (TP56) – 7D non-linear equality constrained problem.
 *
 * Fortran TP56 summary:
 *
 *   N    = 7
 *   NILI = 0
 *   NINL = 0
 *   NELI = 0
 *   NENL = 4   (4 nonlinear equality constraints)
 *
 * Variables:
 *   x = (x1..x7)
 *
 * Objective:
 *
 *   f(x) = -x1 * x2 * x3
 *
 * Constraints (G(i) = 0):
 *
 *   G1(x) = x1 - 4.2 * sin(x4)^2
 *   G2(x) = x2 - 4.2 * sin(x5)^2
 *   G3(x) = x3 - 4.2 * sin(x6)^2
 *   G4(x) = x1 + 2*x2 + 2*x3 - 7.2 * sin(x7)^2
 *
 * Jacobian (dG/dx) from Fortran GG:
 *
 *   ∂G1/∂x1 = 1
 *   ∂G1/∂x4 = -8.4 * sin(x4) * cos(x4)
 *
 *   ∂G2/∂x2 = 1
 *   ∂G2/∂x5 = -8.4 * sin(x5) * cos(x5)
 *
 *   ∂G3/∂x3 = 1
 *   ∂G3/∂x6 = -8.4 * sin(x6) * cos(x6)
 *
 *   ∂G4/∂x1 = 1
 *   ∂G4/∂x2 = 2
 *   ∂G4/∂x3 = 2
 *   ∂G4/∂x7 = -14.4 * sin(x7) * cos(x7)
 *
 * No bounds: LXL/LXU are .FALSE. for all components.
 *
 * Initial guess (mode 1):
 *
 *   X(1) = 1
 *   X(2) = 1
 *   X(3) = 1
 *   X(4) = asin(sqrt(1 / 4.2))
 *   X(5) = X(4)
 *   X(6) = X(4)
 *   X(7) = asin(sqrt(5 / 7.2))
 *
 * Exact solution:
 *
 *   XEX(1) = 2.4
 *   XEX(2) = 1.2
 *   XEX(3) = 1.2
 *   XEX(4) = asin(sqrt(4/7))
 *   XEX(5) = asin(sqrt(2/7))
 *   XEX(6) = XEX(5)
 *   XEX(7) = 2*atan(1)   ( = π/2 )
 *   FEX    = -3.456
 *
 * LEX = .TRUE. (we check that f(x*) == FEX).
 */
public class HS056Test {

    private static final int DIM    = 7;
    private static final int NUM_EQ = 4;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------

    private static class HS56Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {

            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);

            // FX = -X(1)*X(2)*X(3)
            return -x1 * x2 * x3;
        }

        @Override
        public RealVector gradient(RealVector x) {

            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);

            // GF(1) = -X(2)*X(3)
            // GF(2) = -X(1)*X(3)
            // GF(3) = -X(1)*X(2)
            // GF(4..7) = 0
            final double[] g = new double[DIM];
            g[0] = -x2 * x3;
            g[1] = -x1 * x3;
            g[2] = -x1 * x2;
            g[3] = 0.0;
            g[4] = 0.0;
            g[5] = 0.0;
            g[6] = 0.0;

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Start SQP with zero Hessian; BFGS will update it.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Equality constraints G(x) = 0
    // -------------------------------------------------------------------------

    private static class HS56Eq extends EqualityConstraint {

        HS56Eq() {
            // RHS = 0 for all 4 constraints:
            super(new ArrayRealVector(new double[NUM_EQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);
            final double x5 = x.getEntry(4);
            final double x6 = x.getEntry(5);
            final double x7 = x.getEntry(6);

            final double s4 = FastMath.sin(x4);
            final double s5 = FastMath.sin(x5);
            final double s6 = FastMath.sin(x6);
            final double s7 = FastMath.sin(x7);

            // G1 = X(1) - 4.2 * sin(X(4))^2
            // G2 = X(2) - 4.2 * sin(X(5))^2
            // G3 = X(3) - 4.2 * sin(X(6))^2
            // G4 = X(1) + 2*X(2) + 2*X(3) - 7.2 * sin(X(7))^2
            final double g1 = x1 - 4.2 * s4 * s4;
            final double g2 = x2 - 4.2 * s5 * s5;
            final double g3 = x3 - 4.2 * s6 * s6;
            final double g4 = x1 + 2.0 * x2 + 2.0 * x3 - 7.2 * s7 * s7;

            return new ArrayRealVector(new double[] { g1, g2, g3, g4 }, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            final double x4 = x.getEntry(3);
            final double x5 = x.getEntry(4);
            final double x6 = x.getEntry(5);
            final double x7 = x.getEntry(6);

            final double s4 = FastMath.sin(x4);
            final double c4 = FastMath.cos(x4);
            final double s5 = FastMath.sin(x5);
            final double c5 = FastMath.cos(x5);
            final double s6 = FastMath.sin(x6);
            final double c6 = FastMath.cos(x6);
            final double s7 = FastMath.sin(x7);
            final double c7 = FastMath.cos(x7);

            // From Fortran mode 1 + mode 5:
            //
            // GG(1,1) = 1
            // GG(2,2) = 1
            // GG(3,3) = 1
            // GG(4,1) = 1
            // GG(4,2) = 2
            // GG(4,3) = 2
            //
            // GG(1,4) = -8.4 * sin(X(4)) * cos(X(4))
            // GG(2,5) = -8.4 * sin(X(5)) * cos(X(5))
            // GG(3,6) = -8.4 * sin(X(6)) * cos(X(6))
            // GG(4,7) = -14.4 * sin(X(7)) * cos(X(7))

            RealMatrix J = new Array2DRowRealMatrix(NUM_EQ, DIM);

            // Row 0: G1
            J.setEntry(0, 0, 1.0);                            // dG1/dx1
            J.setEntry(0, 3, -8.4 * s4 * c4);                 // dG1/dx4

            // Row 1: G2
            J.setEntry(1, 1, 1.0);                            // dG2/dx2
            J.setEntry(1, 4, -8.4 * s5 * c5);                 // dG2/dx5

            // Row 2: G3
            J.setEntry(2, 2, 1.0);                            // dG3/dx3
            J.setEntry(2, 5, -8.4 * s6 * c6);                 // dG3/dx6

            // Row 3: G4
            J.setEntry(3, 0, 1.0);                            // dG4/dx1
            J.setEntry(3, 1, 2.0);                            // dG4/dx2
            J.setEntry(3, 2, 2.0);                            // dG4/dx3
            J.setEntry(3, 6, -14.4 * s7 * c7);                // dG4/dx7

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // JUnit test
    // -------------------------------------------------------------------------

    @Test
    public void testHS56_optimization() {

        // Initial guess from Fortran mode 1:
        double[] x0 = new double[DIM];

        x0[0] = 1.0;  // X(1)
        x0[1] = 1.0;  // X(2)
        x0[2] = 1.0;  // X(3)

        // X(4) = asin(sqrt(1 / 4.2))
        x0[3] = FastMath.asin(FastMath.sqrt(1.0 / 4.2));

        // X(5) = X(4), X(6) = X(4)
        x0[4] = x0[3];
        x0[5] = x0[3];

        // X(7) = asin(sqrt(5 / 7.2))
        x0[6] = FastMath.asin(FastMath.sqrt(5.0 / 7.2));
        
        // Optimizer
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        SQPOption option=new SQPOption();
        option.setGradientMode(GradientMode.FORWARD);
        LagrangeSolution sol = opt.optimize(
                option,
            new InitialGuess(x0),
            new ObjectiveFunction(new HS56Obj()),
            new HS56Eq() // 4 equality constraints G(x) = 0
            
        );

        final double f = sol.getValue();

        // Fortran: LEX = .TRUE., FEX = -3.456D0
        final double fExpected = -3.456;
        

        HSProblemTestUtils.assertExpectedObjective(fExpected, sol);
    }
}
