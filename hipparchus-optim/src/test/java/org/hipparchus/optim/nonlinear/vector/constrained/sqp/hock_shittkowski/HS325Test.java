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
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS325Test {

    private static final int DIM = 2;

    static final class HS325Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            // FX = X1^2 + X2
            return x1 * x1 + x2;
        }
        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            // GF(1) = 2*X1, GF(2) = 1
            double g1 = 2.0 * x1;
            double g2 = 1.0;
            return new ArrayRealVector(new double[]{g1, g2}, false);
        }
        @Override public RealMatrix hessian(RealVector x) {
            // Constant Hessian H = [[2, 0], [0, 0]]
            double[][] H = new double[DIM][DIM];
            H[0][0] = 2.0;
            return MatrixUtils.createRealMatrix(H);
        }
    }

    static final class HS325Ineq extends InequalityConstraint {
        
        HS325Ineq() { super(new ArrayRealVector(new double[2])); } 

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            
            // G(1) = -(X1 + X2) + 1 >= 0
            double g1 = -(x1 + x2) + 1.0;
            // G(2) = -(X1 + X2^2) + 1 >= 0
            double g2 = -(x1 + x2 * x2) + 1.0; 
            
            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x2 = x.getEntry(1);
            double[][] J = new double[2][DIM];

            // G1: -1, -1
            J[0][0] = -1.0;
            J[0][1] = -1.0;
            
            // G2: -1, -2*X2
            J[1][0] = -1.0;
            J[1][1] = -2.0 * x2;

            return MatrixUtils.createRealMatrix(J);
        }
    }

    static final class HS325Eq extends EqualityConstraint {
        
        HS325Eq() { super(new ArrayRealVector(new double[1])); } 

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            // G(3) = X1^2 + X2^2 - 9 = 0
            double g3 = x1 * x1 + x2 * x2 - 9.0;
            return new ArrayRealVector(new double[]{g3}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double[][] J = new double[1][DIM];

            // G3: 2*X1, 2*X2
            J[0][0] = 2.0 * x1;
            J[0][1] = 2.0 * x2;

            return MatrixUtils.createRealMatrix(J);
        }
    }
    
    private static double[] start() { 
        return new double[]{-3.0, 0.0}; 
    }

    @Test
    public void testHS325() {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS325Obj()),
                new HS325Ineq(),
                new HS325Eq()
        );

        double f = sol.getValue();
        final double fExpected = 3.7913414;
        
        assertEquals(fExpected, f, 1.0e-6 * (Math.abs(fExpected) + 1.0), "objective mismatch");
       
    }
}
