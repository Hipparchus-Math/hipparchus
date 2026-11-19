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
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

public class HS328Test {

    private static final int DIM = 2;

    static final class HS328Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            
            // Reusable powers
            double x1_sq = x1 * x1;
            double x2_sq = x2 * x2;
            double x1x2_sq = x1_sq * x2_sq;
            double x1x2_pow4 = x1x2_sq * x1x2_sq;
            
            // A = (1 + X2^2) / X1^2
            double a = (1.0 + x2_sq) / x1_sq;
            
            // B = ((X1*X2)^2 + 100) / (X1*X2)^4
            double b = (x1x2_sq + 100.0) / x1x2_pow4;
            
            // F(X) = (12 + X1^2 + A + B) / 10
            return (12.0 + x1_sq + a + b) / 10.0;
        }
        
        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            
            // Reusable powers
            double x1_sq = x1 * x1;
            double x2_sq = x2 * x2;
            double x1_cub = x1_sq * x1;
            double x1_pow4 = x1_cub * x1;
            double x1_pow5 = x1_pow4 * x1;
            double x2_cub = x2_sq * x2;
            double x2_pow4 = x2_cub * x2;
            double x2_pow5 = x2_pow4 * x2;
            
            // --- Partial derivative w.r.t X1 (dF/dX1) ---
            
            // d(A)/dX1 = -2*(1 + X2^2) / X1^3
            double dA_dx1 = -2.0 * (1.0 + x2_sq) / x1_cub;
            
            // d(B)/dX1 = d/dX1 [ X1^-2 * X2^-2 + 100 * X1^-4 * X2^-4 ] 
            // d(B)/dX1 = -2*X1^-3*X2^-2 - 400*X1^-5*X2^-4
            double dB_dx1 = -2.0 / (x1_cub * x2_sq) - 400.0 / (x1_pow5 * x2_pow4);
            
            // d(F*10)/dX1 = 2*X1 + dA/dX1 + dB/dX1
            double g1 = (2.0 * x1 + dA_dx1 + dB_dx1) / 10.0;
            
            // --- Partial derivative w.r.t X2 (dF/dX2) ---

            // d(A)/dX2 = 2*X2 / X1^2
            double dA_dx2 = 2.0 * x2 / x1_sq;
            
            // d(B)/dX2 = d/dX2 [ X1^-2 * X2^-2 + 100 * X1^-4 * X2^-4 ]
            // d(B)/dX2 = -2*X1^-2*X2^-3 - 400*X1^-4*X2^-5
            double dB_dx2 = -2.0 / (x1_sq * x2_cub) - 400.0 / (x1_pow4 * x2_pow5);
            
            // d(F*10)/dX2 = 0 + dA/dX2 + dB/dX2
            double g2 = (dA_dx2 + dB_dx2) / 10.0;
            
            return new ArrayRealVector(new double[]{g1, g2}, false);
        }
        
        @Override public RealMatrix hessian(RealVector x) {
            // Hessian matrix is not implemented for this test case.
            throw new UnsupportedOperationException("Hessian matrix is not implemented for this test case.");
        }
    }

    private static double[] start() { 
        return new double[]{0.5, 0.5}; 
    }
//    @Disabled
    @Test
    public void testHS328() {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }
        SQPOption option=new SQPOption();
        option.setGradientMode(GradientMode.CENTRAL);
        // Box constraints: 0.1 <= X1, X2 <= 3.0
        SimpleBounds bounds = new SimpleBounds(
            new double[]{0.1, 0.1}, 
            new double[]{3.0, 3.0}
        );

        LagrangeSolution sol = opt.optimize(
                option,
                new InitialGuess(start()),
                new ObjectiveFunction(new HS328Obj()),
                bounds
        );

        double f = sol.getValue();
        final double fExpected = 1.744152;
        
        assertEquals(fExpected, f, 1.0e-5 * (Math.abs(fExpected) + 1.0), "objective mismatch");
       
    }
}