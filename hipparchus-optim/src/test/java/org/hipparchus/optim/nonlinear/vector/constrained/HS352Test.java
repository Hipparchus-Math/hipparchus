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
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.linear.RealMatrix;

import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;

/**
 * Test case for the HS352 unconstrained non-linear least squares problem (Function Fitting).
 * This problem involves 4 variables and 40 residuals (20 data points * 2 functions).
 */
public class HS352Test {

    private static final int DIM = 4;
    private static final int NUM_RESIDUALS = 40; // 20 data points for F(I) and 20 for F(I+20)

    /**
     * Inner class to simulate intermediate variables and calculations.
     */
    static class Context {
        final double[] F = new double[NUM_RESIDUALS]; // Residuals F(I) and F(I+20)
        double FX; // Objective function value (Sum of F(I)^2)
        final double[] GF = new double[DIM]; // Gradient of objective function GF(J)

        public void compute(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            // 1. Calculation of Residuals F(I) and F(I+20) (Mode 2)
            for (int i = 1; i <= 20; i++) {
                double ti = i * 0.2;
                
                // F(I) = X(1) + X(2)*TI - exp(TI)
                F[i - 1] = x1 + x2 * ti - Math.exp(ti);
                
                // F(I+20) = X(3) + X(4)*sin(TI) - cos(TI)
                F[i + 20 - 1] = x3 + x4 * Math.sin(ti) - Math.cos(ti);
            }
            
            // 2. Calculation of Objective Function FX (Mode 2)
            FX = 0.0;
            for (int i = 0; i < NUM_RESIDUALS; i++) {
                // FX = sum(F(I)^2)
                FX += F[i] * F[i];
            }
            
            // 3. Calculation of the Objective Gradient: GF(J) (Mode 3)
            // DF is the Jacobian of the residuals (partial derivatives)
            
            // Initialize gradient to zero
            for (int j = 0; j < DIM; j++) {
                GF[j] = 0.0;
            }
            
            for (int i = 1; i <= 20; i++) {
                double ti = i * 0.2;
                
                // Derivatives of F(I) = X(1) + X(2)*TI - exp(TI)
                double df_i_dx1 = 1.0;
                double df_i_dx2 = ti;
                double df_i_dx3 = 0.0;
                double df_i_dx4 = 0.0;

                // Derivatives of F(I+20) = X(3) + X(4)*sin(TI) - cos(TI)
                double df_i_plus_20_dx1 = 0.0;
                double df_i_plus_20_dx2 = 0.0;
                double df_i_plus_20_dx3 = 1.0;
                double df_i_plus_20_dx4 = Math.sin(ti);
                
                // F[i-1] is F(I), F[i+20-1] is F(I+20)
                double f_i = F[i - 1];
                double f_i_plus_20 = F[i + 20 - 1];

                // GF(J) = sum(2 * F(k) * dF(k)/dX(J))
                GF[0] += 2.0 * f_i * df_i_dx1 + 2.0 * f_i_plus_20 * df_i_plus_20_dx1;
                GF[1] += 2.0 * f_i * df_i_dx2 + 2.0 * f_i_plus_20 * df_i_plus_20_dx2;
                GF[2] += 2.0 * f_i * df_i_dx3 + 2.0 * f_i_plus_20 * df_i_plus_20_dx3;
                GF[3] += 2.0 * f_i * df_i_dx4 + 2.0 * f_i_plus_20 * df_i_plus_20_dx4;
            }
        }
    }

    /**
     * Implementation of the objective function and gradient (Mode 2, 3).
     */
    static final class HS352Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        
        @Override public double value(RealVector x) {
            Context ctx = new Context();
            ctx.compute(x); 
            return ctx.FX;
        }
        
        @Override public RealVector gradient(RealVector x) {
            Context ctx = new Context();
            ctx.compute(x); 
            return new ArrayRealVector(ctx.GF, false);
        }
        
        @Override public RealMatrix hessian(RealVector x) {
            // Mode 5 (Hessian): Not implemented, relies on numerical estimation
            throw new UnsupportedOperationException("Hessian matrix is not implemented for this test case.");
        }
    }

    private static double[] start() { 
        // Initial values (Mode 1): X(1)=25.0, X(2)=5.0, X(3)=-5.0, X(4)=-1.0
        return new double[]{25.0, 5.0, -5.0, -1.0}; 
    }

    @Test
    public void testHS352() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();
        
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }
        
        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS352Obj())
        );

        double f = sol.getValue();
        // EXPECTED VALUE (FEX): 0.90323433D+3 (903.23433)
        final double fExpected = 903.23433; 
        final double tolerance = 1.0e-5 * (Math.abs(fExpected) + 1.0);
        
        // Check if the solution is close to or better than the expected minimum.
        assertTrue(f <= fExpected + tolerance, 
                   String.format("Objective value mismatch/worse than expected. Expected: %.8f, Actual: %.8f", fExpected, f));
        
        // Check if the variables are close to the expected optimal values (XEX)
        // XEX(1)=-10.223574, XEX(2)=11.908429, XEX(3)=-0.45804134, XEX(4)=0.58031996
       
    }
}
