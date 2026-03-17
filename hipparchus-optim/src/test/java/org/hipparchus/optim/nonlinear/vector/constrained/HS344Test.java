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
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS344Test {

    private static final int DIM = 3;
    private static final double CONST_TERM = 4.0 + 3.0 * Math.sqrt(2.0);

    static final class HS344Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        
        // F(X) = (X1-1)^2 + (X1-X2)^2 + (X2-X3)^4
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            
            double d1 = x1 - 1.0;
            double d2 = x1 - x2;
            double d3 = x2 - x3;

            return d1 * d1 + d2 * d2 + Math.pow(d3, 4);
        }
        
        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            
            double d3_3 = Math.pow(x2 - x3, 3);

            // GF(1) = 2*(X1-1) + 2*(X1-X2)
            double g1 = 2.0 * (x1 - 1.0) + 2.0 * (x1 - x2);
            
            // GF(2) = -2*(X1-X2) + 4*(X2-X3)^3
            double g2 = -2.0 * (x1 - x2) + 4.0 * d3_3;
            
            // GF(3) = -4*(X2-X3)^3
            double g3 = -4.0 * d3_3;

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }
        
        @Override public RealMatrix hessian(RealVector x) {
            // H11 = 2 + 2 = 4
            // H12 = -2
            // H13 = 0
            // H22 = 2 + 12*(X2-X3)^2
            // H23 = -12*(X2-X3)^2
            // H33 = 12*(X2-X3)^2
            double d3_2 = Math.pow(x.getEntry(1) - x.getEntry(2), 2);

            double[][] H = new double[DIM][DIM];

            H[0][0] = 4.0;
            H[0][1] = -2.0;
            H[1][0] = -2.0;
            H[1][1] = 2.0 + 12.0 * d3_2;
            H[1][2] = -12.0 * d3_2;
            H[2][1] = -12.0 * d3_2;
            H[2][2] = 12.0 * d3_2;

            return MatrixUtils.createRealMatrix(H);
        }
    }

    static final class HS344Eq extends EqualityConstraint {
        
        HS344Eq() { super(new ArrayRealVector(new double[1])); } 

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            
            // G(1) = X1*(1 + X2^2) + X3^4 - (4 + 3*sqrt(2)) = 0
            double h1 = x1 * (1.0 + x2 * x2) + Math.pow(x3, 4) - CONST_TERM;
            
            return new ArrayRealVector(new double[]{h1}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            
            double x2_2 = x2 * x2;
            double x3_3 = Math.pow(x3, 3);

            double[][] J = new double[1][DIM];

            // dH1/dX1 = 1 + X2^2
            J[0][0] = 1.0 + x2_2;
            
            // dH1/dX2 = 2*X1*X2
            J[0][1] = 2.0 * x1 * x2;
            
            // dH1/dX3 = 4*X3^3
            J[0][2] = 4.0 * x3_3;

            return MatrixUtils.createRealMatrix(J);
        }
    }

    private static double[] start() { 
        return new double[]{2.0, 2.0, 2.0}; 
    }

    @Test
    public void testHS344() {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        
       

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS344Obj()),
                new HS344Eq()
                
        );

        double f = sol.getValue();
        final double fExpected = 0.032568200;
        
        assertEquals(fExpected, f, 1.0e-5 * (Math.abs(fExpected) + 1.0), "objective mismatch");
        
       
    }
}
