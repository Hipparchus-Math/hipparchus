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

public class HS335Test {

    private static final int DIM = 3;

    static final class HS335Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        
        // F(X) = -(0.001*X1 + X2)
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            return -(0.001 * x1 + x2);
        }
        
        @Override public RealVector gradient(RealVector x) {
            // Gradient GF: [-0.001, -1.0, 0.0]
            return new ArrayRealVector(new double[]{-0.001, -1.0, 0.0}, false);
        }
        
        @Override public RealMatrix hessian(RealVector x) {
            // Hessian is all zeros (linear function)
            return MatrixUtils.createRealMatrix(DIM, DIM);
        }
    }

    static final class HS335Ineq extends InequalityConstraint {
        
        HS335Ineq() { super(new ArrayRealVector(new double[2])); } 

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            
            double x1_sq = x1 * x1;
            double x2_sq = x2 * x2;
            
            // G(1) = 1000*X1^2 + 100*X2^2 - X3 >= 0
            double g1 = 1000.0 * x1_sq + 100.0 * x2_sq - x3;
            
            // G(2) = 100*X1^2 + 400*X2^2 + X3 - 0.01 >= 0
            double g2 = 100.0 * x1_sq + 400.0 * x2_sq + x3 - 0.01;
            
            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            
            double[][] J = new double[2][DIM];

            // G1: 2000*X1, 200*X2, -1
            J[0][0] = 2000.0 * x1;
            J[0][1] = 200.0 * x2;
            J[0][2] = -1.0;
            
            // G2: 200*X1, 800*X2, 1
            J[1][0] = 200.0 * x1;
            J[1][1] = 800.0 * x2;
            J[1][2] = 1.0;

            return MatrixUtils.createRealMatrix(J);
        }
    }

    private static double[] start() { 
        return new double[]{1.0, 1.0, 1.0}; 
    }

//    @Test
//    public void testHS335() {
//        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer()
//        SQPOption sqpOpt=new SQPOption();
//        sqpOpt.setGradientMode(GradientMode.EXTERNAL);
//        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
//            opt.setDebugPrinter(System.out::println);
//        }
//        
//        
//       
//
//        LagrangeSolution sol = opt.optimize(
//                new InitialGuess(start()),
//                new ObjectiveFunction(new HS335Obj()),
//                new HS335Ineq()
//                
//        );
//
//        double f = sol.getValue();
//        final double fExpected = -0.0044721370;
//        
//        assertEquals(fExpected, f, 1.0e-5 * (Math.abs(fExpected) + 1.0), "objective mismatch");
//       
//    }
}
