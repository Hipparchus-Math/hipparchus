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
 * Test case for the HS359 optimization problem (Hock & Schittkowski).
 * Problem: Linear Programming (LP) with 5 variables and 14 linear inequality constraints.
 * The problem definition follows the logic of the Fortran subroutine TP359.
 */
public class HS359Test {

    private static final int DIM = 5;
    private static final int NUM_INEQUALITIES = 14;

    // Coefficients as defined in the Fortran DATA block (MODE 1)
    private static final double[] A = {
        -0.8720288849e7, 0.1505125253e6, -0.1566950325e3, 0.4764703222e6, 0.7294828271e6
    };
    private static final double[] B = {
        -0.145421402e6, 0.29311506e4, -0.40427932e2, 0.5106192e4, 0.1571136e5
    };
    private static final double[] C = {
        -0.1550111084e6, 0.436053352e4, 0.129492344e2, 0.10236884e5, 0.13176786e5
    };
    private static final double[] D = {
        -0.3266695104e6, 0.739068412e4, -0.278986976e2, 0.16643076e5, 0.30988146e5
    };

    /**
     * Objective Function F(X).
     * F(X) = 24345.0 - SUM(A(I) * X(I)) (Fortran MODE 2)
     */
    private static class HS359Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }

        @Override public double value(RealVector X) {
            double s = 0.0;
            for (int i = 0; i < DIM; i++) s += A[i] * X.getEntry(i);
            
            // Simplified logic from Fortran MODE 2: FX = 24345.0 - SUM(A*X)
            return 24345.0 - s;
        }

        @Override public RealVector gradient(RealVector X) {
            // Gradient dF/dX_i = -A_i (Fortran GF)
            double[] g = new double[DIM];
            for (int i = 0; i < DIM; i++) g[i] = -A[i];
            return new ArrayRealVector(g);
        }

        // The Hessian of a linear function is the zero matrix
        @Override public RealMatrix hessian(RealVector x) { return new Array2DRowRealMatrix(DIM, DIM); }
    }

    /**
     * Linear Inequality Constraints G(X) >= 0. (Fortran MODE 4)
     */
    private static class HS359Ineq extends InequalityConstraint {
        // RHS = 0 for all 14 inequalities, as we move constants to the left side
        HS359Ineq() { super(new ArrayRealVector(new double[NUM_INEQUALITIES])); } 

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector X) {
            final double x1 = X.getEntry(0), x2 = X.getEntry(1), x3 = X.getEntry(2),
                         x4 = X.getEntry(3), x5 = X.getEntry(4);

            double[] g = new double[NUM_INEQUALITIES];

            // G(1) to G(8): Inter-variable constraints
            g[0] =  2.4*x1 - x2;      // G(1): 2.4*x1 - x2 >= 0
            g[1] = -1.2*x1 + x2;      // G(2): x2 - 1.2*x1 >= 0
            g[2] = 60.0*x1 - x3;      // G(3): 60.0*x1 - x3 >= 0
            g[3] = -20.0*x1 + x3;     // G(4): x3 - 20.0*x1 >= 0
            g[4] = 9.3*x1 - x4;       // G(5): 9.3*x1 - x4 >= 0
            g[5] = -9.0*x1 + x4;      // G(6): x4 - 9.0*x1 >= 0
            g[6] = 7.0*x1 - x5;       // G(7): 7.0*x1 - x5 >= 0
            g[7] = -6.5*x1 + x5;      // G(8): x5 - 6.5*x1 >= 0

            // Calculate sums S_B, S_C, S_D
            double sB = 0.0, sC = 0.0, sD = 0.0;
            for (int i = 0; i < DIM; i++) {
                final double xi = X.getEntry(i);
                sB += B[i] * xi;
                sC += C[i] * xi;
                sD += D[i] * xi;
            }

            // G(9) to G(11): Lower bound constraints (SUM(...) >= 0)
            g[8]  = sB;                     // G(9): sum(B*x) >= 0
            g[9]  = sC;                     // G(10): sum(C*x) >= 0
            g[10] = sD;                     // G(11): sum(D*x) >= 0

            // G(12) to G(14): Upper bound constraints (SUM(...) <= C) -> C - SUM(...) >= 0
            g[11] = 294000.0 - sB;          // G(12): sum(B*x) <= 294000
            g[12] = 294000.0 - sC;          // G(13): sum(C*x) <= 294000
            g[13] = 277200.0 - sD;          // G(14): sum(D*x) <= 277200

            return new ArrayRealVector(g);
        }

        // Jacobian is constant for linear constraints.
        // We throw UnsupportedOperationException to rely on the optimizer's finite differences 
        // as per the clean example provided.
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }
    
    // Initial guess (Fortran MODE 1)
    private static final double[] X_START = {
        2.52, 5.04, 94.5, 23.31, 17.136
    };

    // Expected solution (FEX and XEX from Fortran MODE 1)
    private static final double F_EXPECTED = -5280416.8;
   

    @Test
    public void testHS359() {
        // Initialize the SQP optimizer
        final SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        
        // Enable debug output if the system property is set.
        
        // Define Box Constraints.
        SimpleBounds bounds = new SimpleBounds(
            new double[] { 0.0, 0.0, 0.0, 0.0, 0.0 },
            new double[] { 
                Double.POSITIVE_INFINITY, 
                Double.POSITIVE_INFINITY, 
                Double.POSITIVE_INFINITY, 
                Double.POSITIVE_INFINITY, 
                Double.POSITIVE_INFINITY 
            }
        );
        

        final LagrangeSolution sol = opt.optimize(
            new InitialGuess(X_START),
            new ObjectiveFunction(new HS359Obj()),
            new HS359Ineq(),
            bounds
        );

        double f = sol.getValue();
        
        // Assertions
        final double toleranceF = 1.0e-4 * (FastMath.abs(F_EXPECTED) + 1.0);
        final double toleranceX = 1.0e-2;
        
        // Verify the objective function value (comparing FEX)
        assertEquals(F_EXPECTED, f, toleranceF, "Discrepancy in the final objective value.");
        
       
    }
}
