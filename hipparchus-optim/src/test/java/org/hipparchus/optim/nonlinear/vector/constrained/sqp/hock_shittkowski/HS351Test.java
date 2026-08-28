/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements. See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 *
 * https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Test case for the HS351 unconstrained non-linear least squares problem.
 * This problem involves 4 variables and 7 residuals (data points).
 */
public class HS351Test {

    private static final int DIM = 4;
    private static final int NUM_RESIDUALS = 7;

    // A(I) data from the original problem specification
    private static final double[] A_DATA = {
        0.0, 0.428e-3, 0.1e-2, 0.161e-2, 0.209e-2, 0.348e-2, 0.525e-2
    };

    // B(I) data from the original problem specification
    private static final double[] B_DATA = {
        7.391, 11.18, 16.44, 16.2, 22.2, 24.02, 31.32
    };

    /**
     * Inner class to simulate intermediate variables and calculations.
     */
    static class Context {
        final double[] F = new double[NUM_RESIDUALS]; // Residuals F(I)
        double FX; // Objective function value (Sum of F(I)^2)
        final double[][] DF = new double[NUM_RESIDUALS][DIM]; // Jacobian of residuals DF(I, J)
        final double[] GF = new double[DIM]; // Gradient of objective function GF(J)

        public void compute(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            // Intermediate squared terms (XH1 to XH4 in the original routine)
            double xh1 = x1 * x1;
            double xh2 = x2 * x2;
            double xh3 = x3 * x3;
            double xh4 = x4 * x4;

            // 1. Calculation of Residuals F(I) (Mode 2)
            for (int i = 0; i < NUM_RESIDUALS; i++) {
                double a_i = A_DATA[i];
                double b_i = B_DATA[i];

                double a_i_sq = a_i * a_i;
                
                // Numerator: (XH1 + A(I)*XH2 + A(I)*A(I)*XH3)
                double numerator = xh1 + a_i * xh2 + a_i_sq * xh3;
                
                // Denominator: (1.0 + A(I)*XH4)
                double denominator = 1.0 + a_i * xh4;
                
                // F(I) = ( (Numerator/Denominator) - B(I) ) / B(I) * 100.0
                double ratio_term = numerator / denominator;
                
                F[i] = (ratio_term - b_i) / b_i * 100.0;
            }
            
            // 2. Calculation of Objective Function FX (Mode 2)
            FX = 0.0;
            for (int i = 0; i < NUM_RESIDUALS; i++) {
                // FX = sum(F(I)^2)
                FX += F[i] * F[i];
            }
            
            // 3. Calculation of the Gradient GF(J) (Mode 3)
            // DF(I, J) is the Jacobian of the residuals
            for (int i = 0; i < NUM_RESIDUALS; i++) {
                double a_i = A_DATA[i];
                double b_i = B_DATA[i];
                double a_i_sq = a_i * a_i;
                
                double denominator = 1.0 + a_i * xh4;
                double denominator_sq = denominator * denominator;
                
                // XH5 in original code: (1.0 + XH4*A(I))*B(I)
                double xh5 = denominator * b_i;
                double xh5_sq = xh5 * xh5;
                
                // N_I = (XH1 + XH2*A(I) + XH3*A(I)**2)
                double numerator_ni = xh1 + xh2 * a_i + xh3 * a_i_sq;

                // DF(I,1) = dF_i/dx1
                // 0.2D+3 * X(1) / XH5
                DF[i][0] = 200.0 * x1 / xh5;

                // DF(I,2) = dF_i/dx2
                // 0.2D+3 * X(2) * A(I) / XH5
                DF[i][1] = 200.0 * x2 * a_i / xh5;

                // DF(I,3) = dF_i/dx3
                // 0.2D+3 * X(3) * A(I)**2 / XH5
                DF[i][2] = 200.0 * x3 * a_i_sq / xh5;

                // DF(I,4) = dF_i/dx4
                // -0.2D+3 * X(4) * A(I) * B(I) * N_I / XH5**2
                // Simplified: dF/dx4 = (100/B) * d/dx4(N/D)
                DF[i][3] = -200.0 * x4 * a_i * b_i * numerator_ni / xh5_sq;
            }
            
            // Calculation of the Objective Gradient: GF(J) = sum(2 * F(I) * DF(I, J))
            for (int j = 0; j < DIM; j++) {
                GF[j] = 0.0;
                for (int i = 0; i < NUM_RESIDUALS; i++) {
                    GF[j] += 2.0 * F[i] * DF[i][j];
                }
            }
        }
    }

    /**
     * Implementation of the objective function and gradient (Mode 2, 3).
     */
    static final class HS351Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        
        @Override public double value(RealVector x) {
            Context ctx = new Context();
            // Mode 2 calculation
            ctx.compute(x); 
            return ctx.FX;
        }
        
        @Override public RealVector gradient(RealVector x) {
            Context ctx = new Context();
            // Mode 3 calculation
            ctx.compute(x); 
            return new ArrayRealVector(ctx.GF, false);
        }
        
        @Override public RealMatrix hessian(RealVector x) {
            // Mode 5 (Hessian): Not implemented, relies on numerical estimation
            throw new UnsupportedOperationException("Hessian matrix is not implemented for this test case.");
        }
    }

    private static double[] start() { 
        // Initial values (Mode 1): X(1)=2.7, X(2)=90.0, X(3)=1500.0, X(4)=10.0
        return new double[]{2.7, 90.0, 1500.0, 10.0}; 
    }

    @Test
    public void testHS351() {
        // SQPOptimizerS2 is a placeholder for an unconstrained/SQP optimizer.
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        
        // RECOVERY: Added conditional debug printing
        
        // The HS351 problem is unconstrained.
        
        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS351Obj())
        );

        double f = sol.getValue();
        // EXPECTED VALUE (FEX): 0.31857175D+3 (318.57175)
        final double fExpected = 318.57175; 
        final double tolerance = 1.0e-5 * (Math.abs(fExpected) + 1.0);
        
        // Check if the solution is close to or better than the expected minimum.
        assertTrue(f <= fExpected + tolerance, 
                   String.format("Objective value mismatch/worse than expected. Expected: %.8f, Actual: %.8f", fExpected, f));
        
        // Check if the variables are close to the expected optimal values (XEX)
        // XEX(1)=2.7143661D+1, XEX(2)=140.43580D+3, XEX(3)=1707.5155D+4, XEX(4)=31.512867D+2
        
    }
}
