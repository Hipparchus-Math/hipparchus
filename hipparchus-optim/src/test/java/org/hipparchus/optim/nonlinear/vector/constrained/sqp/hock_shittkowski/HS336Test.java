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
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS336Test {

    private static final int DIM = 3;

    static final class HS336Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        
        // F(X) = 7*X1 - 6*X2 + 4*X3
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            return 7.0 * x1 - 6.0 * x2 + 4.0 * x3;
        }
        
        @Override public RealVector gradient(RealVector x) {
            // Gradient GF: [7.0, -6.0, 4.0]
            return new ArrayRealVector(new double[]{7.0, -6.0, 4.0}, false);
        }
        
        @Override public RealMatrix hessian(RealVector x) {
            // Hessian is all zeros (linear function)
            return MatrixUtils.createRealMatrix(DIM, DIM);
        }
    }

    static final class HS336Eq extends EqualityConstraint {
        
        HS336Eq() { super(new ArrayRealVector(new double[2])); } 

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            
            // G(1) = 5*X1 + 5*X2 - 3*X3 - 6 = 0
            double h1 = 5.0 * x1 + 5.0 * x2 - 3.0 * x3 - 6.0;
            
            // G(2) = X1^2 + 2*X2^2 + 3*X3^2 - 1 = 0
            double h2 = x1 * x1 + 2.0 * x2 * x2 + 3.0 * x3 * x3 - 1.0;
            
            return new ArrayRealVector(new double[]{h1, h2}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            
            double[][] J = new double[2][DIM];

            // h1: 5.0, 5.0, -3.0 (Linear part handled by GG(1,j) )
            J[0][0] = 5.0;
            J[0][1] = 5.0;
            J[0][2] = -3.0;
            
            // h2: 2*X1, 4*X2, 6*X3 (Non-linear part handled by GG(2,j) )
            J[1][0] = 2.0 * x1;
            J[1][1] = 4.0 * x2;
            J[1][2] = 6.0 * x3;

            return MatrixUtils.createRealMatrix(J);
        }
    }

    private static double[] start() { 
        return new double[]{0.0, 0.0, 0.0}; 
    }

    @Test
    public void testHS336() {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        
        

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS336Obj()),
                new HS336Eq()
                
        );

        double f = sol.getValue();
        final double fExpected = -0.33789573;
        
        assertEquals(fExpected, f, 1.0e-5 * (Math.abs(fExpected) + 1.0), "objective mismatch");
        
       
    }
}
