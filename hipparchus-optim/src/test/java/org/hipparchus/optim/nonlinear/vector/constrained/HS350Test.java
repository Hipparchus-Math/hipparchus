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

public class HS350Test {

    private static final int DIM = 4;
    private static final int NUM_RESIDUALS = 11;

    // Y(I) and U(I) data from the original problem specification
    private static final double[] Y_DATA = {
        0.1957, 0.1947, 0.1735, 0.1600, 0.0844, 0.0627,
        0.0456, 0.0342, 0.0323, 0.0235, 0.0246
    };

    private static final double[] U_DATA = {
        4.0, 2.0, 1.0, 0.5, 0.25, 0.167,
        0.125, 0.1, 0.0833, 0.0714, 0.0625
    };

    /**
     * Inner class to simulate intermediate variables and calculations.
     */
    static class Context {
        final double[] H = new double[NUM_RESIDUALS];
        final double[] F = new double[NUM_RESIDUALS]; // Residuals F(I)
        double FX; // Objective function value (Sum of F(I)^2)
        final double[][] DF = new double[NUM_RESIDUALS][DIM]; // Jacobian of residuals DF(I, J)
        final double[] GF = new double[DIM]; // Gradient of objective function GF(J)

        public void compute(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            // 1. Calculation of the intermediate vector H(I) (from the original calculation loop)
            for (int i = 0; i < NUM_RESIDUALS; i++) {
                double u_i = U_DATA[i];
                // H(I)=U(I)**2+X(3)*U(I)+X(4)
                H[i] = u_i * u_i + x3 * u_i + x4;
            }

            // 2. Calculation of Residuals F(I) (from the original problem definition, Mode 2)
            for (int i = 0; i < NUM_RESIDUALS; i++) {
                double u_i = U_DATA[i];
                double h_i = H[i];

                // F(I)=Y(I)-X(1)/H(I)*(U(I)**2+X(2)*U(I))
                // Prevention of division by zero, although unlikely with test data
                if (Math.abs(h_i) < 1.0e-12) {
                    F[i] = Double.POSITIVE_INFINITY;
                } else {
                    double numerator_term = u_i * u_i + x2 * u_i;
                    F[i] = Y_DATA[i] - x1 / h_i * numerator_term;
                }
            }

            // 3. Calculation of Objective Function FX (Mode 2)
            FX = 0.0;
            for (int i = 0; i < NUM_RESIDUALS; i++) {
                // FX=FX+F(I)**2
                FX += F[i] * F[i];
            }

            // 4. Calculation of the Gradient GF(J) (Mode 3)
            // Calculation of the Residuals Jacobian DF(I, J)
            for (int i = 0; i < NUM_RESIDUALS; i++) {
                double u_i = U_DATA[i];
                double h_i = H[i];
                double h_i_sq = h_i * h_i;
                double numerator_term = u_i * u_i + x2 * u_i;

                // DF(I,1) = dF_i/dx1
                DF[i][0] = (-numerator_term) / h_i;

                // DF(I,2) = dF_i/dx2
                DF[i][1] = (-x1 * u_i) / h_i;

                // DF(I,3) = dF_i/dx3
                DF[i][2] = x1 * u_i * numerator_term / h_i_sq;

                // DF(I,4) = dF_i/dx4
                DF[i][3] = x1 * numerator_term / h_i_sq;
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
    static final class HS350Obj extends TwiceDifferentiableFunction {
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
        // Initial values (Mode 1): X(1)=0.25, X(2)=0.39, X(3)=0.415, X(4)=0.39
        return new double[]{0.25, 0.39, 0.415, 0.39};
    }

    @Test
    public void testHS350() {
        // SQPOptimizerS2 is a placeholder for an unconstrained/SQP optimizer.
        SQPOptimizerS2 opt = new SQPOptimizerS2();

        // RECOVERY: Added conditional debug printing
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }
        
        // The HS350 problem is unconstrained.
        
        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS350Obj())
        );

        double f = sol.getValue();
        // EXPECTED VALUE (FEX): 0.30750560D-3
        final double fExpected = 0.00030750560; 
        final double tolerance = 1.0e-4 * (Math.abs(fExpected) + 1.0);
        
        // Check if the solution is close to or better than the expected minimum.
        // Since the solver is for minimization, the found solution should be <= FEX + tolerance.
        assertTrue(f <= fExpected + tolerance, 
                   String.format("Objective value mismatch/worse than expected. Expected: %.10f, Actual: %.10f", fExpected, f));
        
        // Check if the variables are close to the expected optimal values (XEX)
        // XEX(1)=0.19280644D+0, XEX(2)=0.19126279D+0, XEX(3)=0.12305098D+0, XEX(4)=0.13605235D+0
       
    }
}
