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
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS329Test {

    private static final int DIM = 2;

    static final class HS329Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            // F(X) = (X1 - 10)^3 + (X2 - 20)^3
            return Math.pow(x1 - 10.0, 3) + Math.pow(x2 - 20.0, 3);
        }
        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            // Gradient: 3*(X1 - 10)^2, 3*(X2 - 20)^2
            double g1 = 3.0 * Math.pow(x1 - 10.0, 2);
            double g2 = 3.0 * Math.pow(x2 - 20.0, 2);
            return new ArrayRealVector(new double[]{g1, g2}, false);
        }
        @Override public RealMatrix hessian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double[][] H = new double[DIM][DIM];
            // Hessian: H11 = 6*(X1 - 10), H22 = 6*(X2 - 20)
            H[0][0] = 6.0 * (x1 - 10.0);
            H[1][1] = 6.0 * (x2 - 20.0);
            return MatrixUtils.createRealMatrix(H);
        }
    }

    static final class HS329Ineq extends InequalityConstraint {
        
        HS329Ineq() { super(new ArrayRealVector(new double[3])); } 

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            
            // G(1) = (X1 - 5)^2 + (X2 - 5)^2 - 100 >= 0
            double term1_sq = Math.pow(x1 - 5.0, 2);
            double term2_sq = Math.pow(x2 - 5.0, 2);
            double g1 = term1_sq + term2_sq - 100.0;
            
            // G(2) = (X1 - 6)^2 + (X2 - 5)^2 >= 0 (Always satisfied, but included for completeness)
            double term3_sq = Math.pow(x1 - 6.0, 2);
            double g2 = term3_sq + term2_sq; 
            
            // G(3) = 82.81 - (X1 - 6)^2 - (X2 - 5)^2 >= 0
            double g3 = 82.81 - term3_sq - term2_sq;
            
            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double[][] J = new double[3][DIM];

            // G1: 2*(X1 - 5), 2*(X2 - 5)
            J[0][0] = 2.0 * (x1 - 5.0);
            J[0][1] = 2.0 * (x2 - 5.0);
            
            // G2: 2*(X1 - 6), 2*(X2 - 5)
            J[1][0] = 2.0 * (x1 - 6.0);
            J[1][1] = 2.0 * (x2 - 5.0);
            
            // G3: -2*(X1 - 6), -2*(X2 - 5)
            J[2][0] = -2.0 * (x1 - 6.0);
            J[2][1] = -2.0 * (x2 - 5.0);

            return MatrixUtils.createRealMatrix(J);
        }
    }

    private static double[] start() { 
        return new double[]{14.35, 8.6}; 
    }

    @Test
    public void testHS329() {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        
        // Box constraints: 13.0 <= X1 <= 16.0, 0.0 <= X2 <= 15.0
        SimpleBounds bounds = new SimpleBounds(
            new double[]{13.0, 0.0}, 
            new double[]{16.0, 15.0}
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS329Obj()),
                new HS329Ineq(),
                bounds
        );

        double f = sol.getValue();
        final double fExpected = -6961.8139;
        
        assertEquals(fExpected, f, 1.0e-3 * (Math.abs(fExpected) + 1.0), "objective mismatch");
        
       
    }
}
