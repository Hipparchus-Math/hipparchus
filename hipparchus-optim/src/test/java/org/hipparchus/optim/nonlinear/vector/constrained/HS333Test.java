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

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS333Test {

    private static final int DIM = 3;
    private static final int DATA_SIZE = 8;
    
    // Independent variable data A
    private static final double[] A_DATA = {
        4.0, 5.75, 7.5, 24.0, 32.0, 48.0, 72.0, 96.0
    };
    
    // Dependent variable data Y
    private static final double[] Y_DATA = {
        72.1, 65.6, 55.9, 17.1, 9.8, 4.5, 1.3, 0.6
    };

    static final class HS333Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        
        // F(X) = sum( ( (Y(I) - X1*exp(-X2*A(I)) - X3) / Y(I) )^2 )
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double fx = 0.0;
            
            for (int i = 0; i < DATA_SIZE; i++) {
                double exp_term = Math.exp(-x2 * A_DATA[i]);
                // Residual F(I)
                double fi = (Y_DATA[i] - x1 * exp_term - x3) / Y_DATA[i];
                fx += fi * fi;
            }
            return fx;
        }
        
        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            
            double g1 = 0.0;
            double g2 = 0.0;
            double g3 = 0.0;
            
            for (int i = 0; i < DATA_SIZE; i++) {
                double exp_term = Math.exp(-x2 * A_DATA[i]);
                
                // Residual F(I)
                double fi = (Y_DATA[i] - x1 * exp_term - x3) / Y_DATA[i];

                // DF(I,1) = (-exp(-X2*A(I)))/Y(I)
                double dfi_dx1 = -exp_term / Y_DATA[i];
                
                // DF(I,2) = (X1*A(I)*exp(-X2*A(I)))/Y(I)
                double dfi_dx2 = (x1 * A_DATA[i] * exp_term) / Y_DATA[i];
                
                // DF(I,3) = -1.0/Y(I)
                double dfi_dx3 = -1.0 / Y_DATA[i];
                
                // Gradient GF(j) = 2 * sum(F(I) * DF(I, j))
                g1 += dfi_dx1 * fi * 2.0;
                g2 += dfi_dx2 * fi * 2.0;
                g3 += dfi_dx3 * fi * 2.0;
            }
            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }
        
        @Override public RealMatrix hessian(RealVector x) {
            // Hessian matrix is not implemented for this test case.
            throw new UnsupportedOperationException("Hessian matrix is not implemented for this test case.");
        }
    }

    private static double[] start() { 
        return new double[]{30.0, 0.04, 3.0}; 
    }

    @Test
    public void testHS333() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();

        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }
        
       
        SimpleBounds bounds = new SimpleBounds(
            new double[]{Double.NEGATIVE_INFINITY, 0.0, Double.NEGATIVE_INFINITY}, 
            new double[]{Double.POSITIVE_INFINITY, 0.07, Double.POSITIVE_INFINITY}
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS333Obj()),
                bounds
        );

        double f = sol.getValue();
        final double fExpected = 0.0432;
        
        assertEquals(fExpected, f, 1.0e-4 * (Math.abs(fExpected) + 1.0), "objective mismatch");
        
        
       
    }
}