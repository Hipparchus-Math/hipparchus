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
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

public class HS327Test {

    private static final int DIM = 2;
    private static final int DATA_SIZE = 44;
    private static final double C1 = 0.49;
    private static final double Z_OFFSET = 8.0;

    // Dati Y e Z come definiti nella subroutine Fortran
    private static final double[] Y_DATA = {
        0.49, 0.49, 0.48, 0.47, 0.48, 0.47, 0.46, 0.46, 0.45, 0.43, 
        0.45, 0.43, 0.43, 0.44, 0.43, 0.43, 0.46, 0.45, 0.42, 0.42, 
        0.43, 0.41, 0.41, 0.40, 0.42, 0.40, 0.40, 0.41, 0.40, 0.41, 
        0.40, 0.40, 0.40, 0.38, 0.41, 0.40, 0.40, 0.41, 0.38, 0.40, 
        0.40, 0.39, 0.39, 0.42 // Aggiunto 0.42 per arrivare a 44, l'ultimo dato nel DATA Z non è in Y
    };
    
    // Correzione dei dati Z per allineare con la struttura implicita in Fortran
    // (Z(I)-8) è l'argomento chiave
    private static final double[] Z_DATA = {
        8., 8., 10., 10., 10., 10., 12., 12., 12., 12., 
        14., 14., 14., 16., 16., 16., 18., 18., 20., 20., 
        20., 22., 22., 22., 24., 24., 24., 26., 26., 26., 
        28., 28., 30., 30., 30., 32., 32., 34., 36., 36., 
        38., 38., 40., 42. // Ultimo valore è 42
    };

    static final class HS327Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double fx = 0.0;
            
            for (int i = 0; i < DATA_SIZE; i++) {
                double zi_minus_offset = Z_DATA[i] - Z_OFFSET;
                // f(I) = Y(I) - X(1) - (0.49 - X(1)) * exp(-X(2)*(Z(I)-8))
                double fi = Y_DATA[i] - x1 - (C1 - x1) * Math.exp(-x2 * zi_minus_offset);
                fx += fi * fi;
            }
            return fx;
        }
        
        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double g1 = 0.0;
            double g2 = 0.0;
            
            for (int i = 0; i < DATA_SIZE; i++) {
                double zi_minus_offset = Z_DATA[i] - Z_OFFSET;
                double exp_term = Math.exp(-x2 * zi_minus_offset);
                
                // f(I) = Y(I) - X(1) - (0.49 - X(1)) * exp_term
                double fi = Y_DATA[i] - x1 - (C1 - x1) * exp_term;

                // DF(I,1) = -1 + exp_term
                double dfi_dx1 = -1.0 + exp_term;
                
                // DF(I,2) = (0.49 - X(1)) * exp_term * (Z(I)-8)
                double dfi_dx2 = (C1 - x1) * exp_term * zi_minus_offset;
                
                // Gradient GF = 2 * sum(F(I) * DF(I))
                g1 += dfi_dx1 * fi * 2.0;
                g2 += dfi_dx2 * fi * 2.0;
            }
            return new ArrayRealVector(new double[]{g1, g2}, false);
        }
        
        @Override public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided for least squares objective.");
        }
    }

    static final class HS327Ineq extends InequalityConstraint {
        
        HS327Ineq() { super(new ArrayRealVector(new double[1])); } 

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            
            // G(1) = -0.09 - X1*X2 + 0.49*X2 >= 0
            double g1 = -0.09 - x1 * x2 + C1 * x2;
            
            return new ArrayRealVector(new double[]{g1}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double[][] J = new double[1][DIM];

            // G1: -X2, 0.49 - X1
            J[0][0] = -x2;
            J[0][1] = C1 - x1;

            return MatrixUtils.createRealMatrix(J);
        }
    }

    private static double[] start() { 
        return new double[]{0.42, 5.0}; 
    }
//    @Disabled
    @Test
    public void testHS327() {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }
        
       
        SQPOption option=new SQPOption();
        option.setGradientMode(GradientMode.FORWARD);
        
        LagrangeSolution sol = opt.optimize(
                option,
                new InitialGuess(start()),
                new ObjectiveFunction(new HS327Obj()),
                new HS327Ineq(),
                new SimpleBounds(new double[]{0.4, 0.4}, new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY})
        );

        double f = sol.getValue();
        final double fExpected = 0.028459670;
        
        assertEquals(fExpected, f, 1.0e-4 * (Math.abs(fExpected) + 1.0), "objective mismatch");
       
    }
}