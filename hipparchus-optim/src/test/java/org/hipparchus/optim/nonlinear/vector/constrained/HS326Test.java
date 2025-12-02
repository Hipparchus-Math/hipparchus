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

public class HS326Test {

    private static final int DIM = 2;

    static final class HS326Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            // FX = X1^2 + X2^2 - 16*X1 - 10*X2
            return x1 * x1 + x2 * x2 - 16.0 * x1 - 10.0 * x2;
        }
        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            // GF(1) = 2*X1 - 16, GF(2) = 2*X2 - 10
            double g1 = 2.0 * x1 - 16.0;
            double g2 = 2.0 * x2 - 10.0;
            return new ArrayRealVector(new double[]{g1, g2}, false);
        }
        @Override public RealMatrix hessian(RealVector x) {
            // Constant Hessian H = [[2, 0], [0, 2]]
            double[][] H = new double[DIM][DIM];
            H[0][0] = 2.0;
            H[1][1] = 2.0;
            return MatrixUtils.createRealMatrix(H);
        }
    }

    static final class HS326Ineq extends InequalityConstraint {
        
        HS326Ineq() { super(new ArrayRealVector(new double[2])); } 

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            
            // G(1) = 11 - X1^2 + 6*X1 - 4*X2 >= 0
            double g1 = 11.0 - x1 * x1 + 6.0 * x1 - 4.0 * x2;
            
            // G(2) = X1*X2 - 3*X2 - exp(X1 - 3) + 1 >= 0
            double g2 = x1 * x2 - 3.0 * x2 - Math.exp(x1 - 3.0) + 1.0; 
            
            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double[][] J = new double[2][DIM];
            
            // G1: -2*X1 + 6, -4
            J[0][0] = -2.0 * x1 + 6.0;
            J[0][1] = -4.0;
            
            // G2: X2 - exp(X1 - 3), X1 - 3
            J[1][0] = x2 - Math.exp(x1 - 3.0);
            J[1][1] = x1 - 3.0;

            return MatrixUtils.createRealMatrix(J);
        }
    }

    private static double[] start() { 
        return new double[]{4.0, 3.0}; 
    }

    @Test
    public void testHS326() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();

        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }
        
       

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS326Obj()),
                new HS326Ineq(),
                new SimpleBounds(new double[]{0.0, 0.0}, new double[]{10.0, 10.0})
        );

        double f = sol.getValue();
        final double fExpected = -79.807821;
        
        assertEquals(fExpected, f, 1.0e-6 * (Math.abs(fExpected) + 1.0), "objective mismatch");
        
        
    }
}