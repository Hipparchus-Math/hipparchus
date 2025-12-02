/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 * https://www.apache.org/licenses/LICENSE-2.0
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
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Test case for the HS360 optimization problem (Hock & Schittkowski).
 * Problem: Nonlinear Programming (NLP) with 5 variables and 2 nonlinear inequality constraints,
 * plus box constraints.
 * The problem definition follows the logic of the Fortran subroutine TP360.
 */
public class HS360Test {

    private static final int DIM = 5;
    private static final int NUM_INEQUALITIES = 2;

    // Coefficients C(1) to C(10) from Fortran DATA block
    private static final double[] C = {
        -0.8720288849e7,  // C(1)
         0.1505125233e6,  // C(2)
        -0.1566950325e3,  // C(3)
         0.4764703222e6,  // C(4)
         0.7294828271e6,  // C(5)
        -0.3266695104e6,  // C(6)
         0.739068412e4,   // C(7)
        -0.278986976e2,   // C(8)
         0.16643076e5,    // C(9)
         0.30988146e5     // C(10)
    };

    /**
     * Objective Function F(X) (Fortran MODE 2).
     * F(X) = (C(1) + C(2)*X(2) + C(3)*X(3) + C(4)*X(4) + C(5)*X(5)) * X(1) + 24345.0
     * Note: The Fortran code has a leading minus sign in the first term, which is included here.
     */
    private static class HS360Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }

        @Override public double value(RealVector X) {
            double term2_5 = C[1] * X.getEntry(1) + C[2] * X.getEntry(2) + 
                             C[3] * X.getEntry(3) + C[4] * X.getEntry(4);
            
            // F(X) = (-C(1) - term2_5) * X(1) + 24345.0
            return (-C[0] - term2_5) * X.getEntry(0) + 24345.0;
        }

        @Override public RealVector gradient(RealVector X) {
            final double x1 = X.getEntry(0);
            final double[] g = new double[DIM];

            // GF(1) = -C(1) - C(2)*X(2) - C(3)*X(3) - C(4)*X(4) - C(5)*X(5)
            g[0] = -C[0] - C[1] * X.getEntry(1) - C[2] * X.getEntry(2) - 
                   C[3] * X.getEntry(3) - C[4] * X.getEntry(4);

            // GF(I) = -C(I)*X(1) for I=2 to 5 (indices 1 to 4 in Java)
            for (int i = 1; i < DIM; i++) {
                g[i] = -C[i] * x1;
            }
            return new ArrayRealVector(g);
        }

        @Override public RealMatrix hessian(RealVector X) {
            // Hessian H_ij = d^2F / (dX_i dX_j)
            final RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);

            // H_1i = d/dX_i (GF(1))
            // d/dX_1 (GF(1)) = 0
            // d/dX_i (GF(1)) = -C(i) for i=2 to 5 (indices 1 to 4 in Java)
            for (int i = 1; i < DIM; i++) {
                H.setEntry(0, i, -C[i]);
            }

            // H_i1 = d/dX_1 (GF(i)) for i=2 to 5
            // d/dX_1 (GF(i)) = -C(i)
            for (int i = 1; i < DIM; i++) {
                H.setEntry(i, 0, -C[i]);
            }
            
            // Other terms are zero (GF(i) for i > 1 depends only on X(1))
            return H;
        }
    }

    /**
     * Nonlinear Inequality Constraints G(X) >= 0 (Fortran MODE 4).
     * Two constraints: 0 <= H(X) <= 277200.0, where H(X) is defined below.
     */
    private static class HS360Ineq extends InequalityConstraint {
        // RHS = 0 for both constraints
        HS360Ineq() { super(new ArrayRealVector(new double[NUM_INEQUALITIES])); } 

        // Helper function to calculate H(X)
        private double calculateH(RealVector X) {
            final double x1 = X.getEntry(0);
            // D_i = C(i+5)
            double innerSum = C[5] + C[6] * X.getEntry(1) + C[7] * X.getEntry(2) + 
                              C[8] * X.getEntry(3) + C[9] * X.getEntry(4);
            // H = X(1) * (C(6) + C(7)*X(2) + C(8)*X(3) + C(9)*X(4) + C(10)*X(5))
            return x1 * innerSum;
        }
        
        // Helper function to calculate the gradient of H(X)
        private RealVector calculateGradH(RealVector X) {
            final double x1 = X.getEntry(0);
            final double[] gH = new double[DIM];

            // HH(1) = C(6)+C(7)*X(2)+C(8)*X(3)+C(9)*X(4)+C(10)*X(5) (Derivative w.r.t X1)
            gH[0] = C[5] + C[6] * X.getEntry(1) + C[7] * X.getEntry(2) + 
                    C[8] * X.getEntry(3) + C[9] * X.getEntry(4);

            // HH(I) = C(I+5)*X(1) for I=2 to 5 (indices 1 to 4 in Java)
            for (int i = 1; i < DIM; i++) {
                gH[i] = C[i+5] * x1;
            }
            return new ArrayRealVector(gH);
        }

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector X) {
            final double H = calculateH(X);
            double[] g = new double[NUM_INEQUALITIES];
            
            // G(1): H >= 0
            g[0] = H;                     
            // G(2): 277200.0 - H >= 0 (i.e., H <= 277200.0)
            g[1] = 277200.0 - H;          

            return new ArrayRealVector(g);
        }

        @Override public RealMatrix jacobian(RealVector X) {
            final RealVector gradH = calculateGradH(X);
            final RealMatrix GG = new Array2DRowRealMatrix(NUM_INEQUALITIES, DIM);

            // Row 1 (G(1) = H): Jacobian is grad(H)
            GG.setRowVector(0, gradH);
            
            // Row 2 (G(2) = 277200.0 - H): Jacobian is -grad(H)
            GG.setRowVector(1, gradH.mapMultiply(-1.0));

            return GG;
        }
    }
    
    // Initial guess (Fortran MODE 1)
    private static final double[] X_START = {
        2.52, 2.0, 37.5, 9.25, 6.8
    };

    // Expected solution (FEX and XEX from Fortran MODE 1)
    private static final double F_EXPECTED = -5280335.1;
   


    @Test
    public void testHS360Optimization() {
        // Initialize the SQP optimizer
        final SQPOptimizerS2 opt = new SQPOptimizerS2();
        
        // Enable debug output if the system property is set.
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
             opt.setDebugPrinter(System.out::println);
        }
        
        // Define Box Constraints (Fortran MODE 1)
        final double[] lowerBounds = { 0.0, 1.2, 20.0, 9.0, 6.5 };
        final double[] upperBounds = { 
            Double.POSITIVE_INFINITY, // X(1) is unbounded above
            2.4, 
            60.0, 
            9.3, 
            7.0 
        };
        SimpleBounds bounds = new SimpleBounds(lowerBounds, upperBounds);
        

        final LagrangeSolution sol = opt.optimize(
            new InitialGuess(X_START),
            new ObjectiveFunction(new HS360Obj()),
            new HS360Ineq(),
            bounds
        );

        double f = sol.getValue();
        
        // Assertions
        // Use relative tolerance for objective value
        final double toleranceF = 1.0e-4 * (FastMath.abs(F_EXPECTED) + 1.0);
        // Use a reasonable distance tolerance for the X vector
        final double toleranceX = 1.0e-2;
        
        // Verify the objective function value (comparing FEX)
        assertEquals(F_EXPECTED, f, toleranceF, "Discrepancy in the final objective value.");
        
      
    }
}