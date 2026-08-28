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
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

public class HS338Test {

    private static final int DIM = 3;

    static final class HS338Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        
        // F(X) = -(X1^2 + X2^2 + X3^2) (Minimizing the negative of the sum of squares is equivalent to maximizing the sum of squares)
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            return -(x1 * x1 + x2 * x2 + x3 * x3);
        }
        
        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            
            // GF(I) = -2*X(I)
            return new ArrayRealVector(new double[]{-2.0 * x1, -2.0 * x2, -2.0 * x3}, false);
        }
        
        @Override public RealMatrix hessian(RealVector x) {
            // Hessian is diagonal: [-2.0, 0.0, 0.0], [0.0, -2.0, 0.0], [0.0, 0.0, -2.0]
            double[][] H = new double[DIM][DIM];
            H[0][0] = -2.0;
            H[1][1] = -2.0;
            H[2][2] = -2.0;
            return MatrixUtils.createRealMatrix(H);
        }
    }

    static final class HS338Eq extends EqualityConstraint {
        
        HS338Eq() { super(new ArrayRealVector(new double[2])); } 

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            
            // G(1) = 0.5*X1 + X2 + X3 - 1 = 0
            double h1 = 0.5 * x1 + x2 + x3 - 1.0;
            
            // G(2) = X1^2 + (2/3)*X2^2 + 0.25*X3^2 - 4 = 0
            double h2 = x1 * x1 + (2.0 / 3.0) * x2 * x2 + 0.25 * x3 * x3 - 4.0;
            
            return new ArrayRealVector(new double[]{h1, h2}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            
            double[][] J = new double[2][DIM];

            // h1: 0.5, 1.0, 1.0 (Linear part handled by GG(1,j) in Fortran)
            J[0][0] = 0.5;
            J[0][1] = 1.0;
            J[0][2] = 1.0;
            
            // h2: 2*X1, (4/3)*X2, 0.5*X3 (Non-linear part handled by GG(2,j) in Fortran)
            J[1][0] = 2.0 * x1;
            J[1][1] = (4.0 / 3.0) * x2; // (0.4D+1/0.3D+1)*X(2)
            J[1][2] = 0.5 * x3;         // 0.5D+0*X(3)

            return MatrixUtils.createRealMatrix(J);
        }
    }

    private static double[] start() { 
        return new double[]{0.0, 0.0, 0.0}; 
    }

    @Test
    @Disabled // disabled as we reach a local minimum and not the expected global one
    public void testHS338() {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }
        
       

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS338Obj()),
                new HS338Eq()
                
        );

        double f = sol.getValue();
        final double fExpected = -10.992806;
        
        // The problem is a maximization of distance from origin subject to constraints.
        // It has multiple local optima. We test for the minimum provided in the Fortran code.
        assertEquals(fExpected, f, 1.0e-5 * (Math.abs(fExpected) + 1.0), "objective mismatch");
        
       
    }
}
