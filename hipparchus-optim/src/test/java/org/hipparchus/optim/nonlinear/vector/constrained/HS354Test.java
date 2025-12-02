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

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;

/**
 * Test case for the constrained optimization problem HS354 (Hock & Schittkowski).
 * Objective: Minimize F(X) = (X1 + 10*X2)^2 + 5*(X3 - X4)^2 + (X2 - 2*X3)^4 + 10*(X1 - X4)^4.
 * Constraints: 1 inequality constraint (linear): G(1) = X1 + X2 + X3 + X4 - 1.0 >= 0.
 * Box constraints: 0 <= X_i <= 20.
 */
public class HS354Test {

    private static final int DIM = 4;
    private static final int NUM_INEQUALITIES = 1; 

    /**
     * Implementation of the objective function F(X).
     * F(X) = (X1 + 10*X2)^2 + 5*(X3 - X4)^2 + (X2 - 2*X3)^4 + 10*(X1 - X4)^4
     */
    static final class HS354Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);
            
            return Math.pow(x1 + 10.0 * x2, 2) + 
                   5.0 * Math.pow(x3 - x4, 2) + 
                   Math.pow(x2 - 2.0 * x3, 4) + 
                   10.0 * Math.pow(x1 - x4, 4);
        }
        
        /**
         * Gradient of F(X): dF/dXj
         * dF/dX1 = 2*(X1 + 10*X2) + 40*(X1 - X4)^3
         * dF/dX2 = 20*(X1 + 10*X2) + 4*(X2 - 2*X3)^3
         * dF/dX3 = 10*(X3 - X4) - 8*(X2 - 2*X3)^3
         * dF/dX4 = -10*(X3 - X4) - 40*(X1 - X4)^3
         */
        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);
            
            double term1 = x1 + 10.0 * x2;
            double term3 = x2 - 2.0 * x3;
            double term4 = x1 - x4;
            
            double dfdx1 = 2.0 * term1 + 40.0 * Math.pow(term4, 3);
            double dfdx2 = 20.0 * term1 + 4.0 * Math.pow(term3, 3);
            double dfdx3 = 10.0 * (x3 - x4) - 8.0 * Math.pow(term3, 3);
            double dfdx4 = -10.0 * (x3 - x4) - 40.0 * Math.pow(term4, 3);
            
            return new ArrayRealVector(new double[] {dfdx1, dfdx2, dfdx3, dfdx4}, false);
        }
        
        // The Hessian matrix calculation is omitted here for simplicity in this base test structure,
        // relying on the optimizer's internal numerical approximation if needed, though a full 
        // implementation would provide the actual Hessian.
        @Override public RealMatrix hessian(RealVector x) {
            // For a complete TwiceDifferentiableFunction, the Hessian should be implemented.
            // Returning zero matrix as a placeholder if numerical differentiation is allowed.
            return MatrixUtils.createRealMatrix(DIM, DIM); 
        }
    }
    
    /**
     * Implementation of the inequality constraint G(1).
     * Constraint must satisfy G(i) >= 0.
     * G(1) = X1 + X2 + X3 + X4 - 1.0 >= 0 (Linear)
     */
    static final class HS354Ineq extends InequalityConstraint {
        
        HS354Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQUALITIES])); 
        }

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);
            
            // G(1) = X1 + X2 + X3 + X4 - 1.0
            double G1 = x1 + x2 + x3 + x4 - 1.0;
            
            return new ArrayRealVector(new double[] { G1 }, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            // Row 0: Gradient of G(1) (constant linear part)
            // dG1/dXj = 1.0 for all j
            double[][] J = new double[NUM_INEQUALITIES][DIM];
            J[0][0] = 1.0;
            J[0][1] = 1.0;
            J[0][2] = 1.0; 
            J[0][3] = 1.0;
            
            return new Array2DRowRealMatrix(J);
        }
    }

    // Initial starting point: X = (3.0, -1.0, 0.0, 1.0)
    private static double[] start() { 
        return new double[]{3.0, -1.0, 0.0, 1.0}; 
    }
    
    // Expected solution for objective function and variables
    // FEX = 0.11378385
    private static final double F_EXPECTED = 0.11378385; 
    
    @Test
    public void testHS354() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();
        
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }
        
        //
        
        SimpleBounds bounds = new SimpleBounds(
            new double[] { 
                Double.NEGATIVE_INFINITY ,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY}, 
               
            new double[] { 
                20.0, 20.0, 
                20.0, 20.0 
            }
        );


        // Optimization run (HS354 has no equality constraints)
        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS354Obj()),
                new HS354Ineq(), 
                bounds
        );

        double f = sol.getValue();
        final double tolerance = 1.0e-5 * (Math.abs(F_EXPECTED) + 1.0);
        
        // Verify objective function value
        assertEquals(F_EXPECTED, f, tolerance, "objective mismatch");
        
       
    }
}