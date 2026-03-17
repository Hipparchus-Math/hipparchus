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
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;

/**
 * Test case for the constrained optimization problem HS353 (Hock & Schittkowski).
 * Objective: Minimize F(X) = -(24.55*X1 + 26.75*X2 + 39.0*X3 + 40.5*X4).
 * Constraints: 2 inequality constraints (1 linear, 1 non-linear), 1 equality constraint (linear), and X_i >= 0.
 * The inequality constraints are of the form G(X) >= 0.
 */
public class HS353Test {

    private static final int DIM = 4;
    private static final int NUM_INEQUALITIES = 2; 
    private static final int NUM_EQUALITIES = 1;   

    /**
     * Calculates the non-linear term Q, used in the second inequality constraint.
     * Q = (0.53*X1)^2 + (0.44*X2)^2 + (4.5*X3)^2 + (0.79*X4)^2
     * @param x the current optimization variables vector.
     * @return the value of Q.
     */
    private static double calculateQ(RealVector x) {
        double x1 = x.getEntry(0);
        double x2 = x.getEntry(1);
        double x3 = x.getEntry(2);
        double x4 = x.getEntry(3);
        
        return Math.pow(0.53 * x1, 2) + Math.pow(0.44 * x2, 2) +
               Math.pow(4.5 * x3, 2) + Math.pow(0.79 * x4, 2);
    }

    /**
     * Implementation of the objective function F(X).
     */
    static final class HS353Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        
        // F(X) = -(24.55*X1 + 26.75*X2 + 39.0*X3 + 40.5*X4)
        @Override public double value(RealVector x) {
            double fx = 24.55 * x.getEntry(0) + 26.75 * x.getEntry(1) + 39.0 * x.getEntry(2) + 40.5 * x.getEntry(3);
            return -fx; // Minimization
        }
        
        // Gradient of F(X) (constant vector)
        @Override public RealVector gradient(RealVector x) {
            return new ArrayRealVector(new double[] {-24.55, -26.75, -39.0, -40.5}, false);
        }
        
        // Hessian is zero (linear function)
        @Override public RealMatrix hessian(RealVector x) {
            return MatrixUtils.createRealMatrix(DIM, DIM);
        }
    }
    
    /**
     * Implementation of the inequality constraints G(1) and G(2).
     * These constraints must satisfy G(i) >= 0.
     */
    static final class HS353Ineq extends InequalityConstraint {
        
        HS353Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQUALITIES])); 
        }

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);
            
            double Q = calculateQ(x);

            // G(1) >= 0: 2.3*X1 + 5.6*X2 + 11.1*X3 + 1.3*X4 - 5.0 >= 0
            double G1 = 2.3 * x1 + 5.6 * x2 + 11.1 * x3 + 1.3 * x4 - 5.0;

            // G(2) >= 0: 12.0*X1 + 11.9*X2 + 41.8*X3 + 52.1*X4 - 1.645*sqrt(Q) - 12.0 >= 0
            double G2 = 12.0 * x1 + 11.9 * x2 + 41.8 * x3 + 52.1 * x4 - 
                        1.645 * Math.sqrt(Q) - 12.0;
            
            // Return G(i)
            return new ArrayRealVector(new double[] { G1, G2 }, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1 = x.getEntry(0);
            
            double[][] J = new double[NUM_INEQUALITIES][DIM];
            
            // Row 0: Gradient of G(1) (linear part)
            J[0][0] = 2.3;
            J[0][1] = 5.6;
            J[0][2] = 11.1; 
            J[0][3] = 1.3;

            // Row 1: Gradient of G(2)
            double Q = calculateQ(x);
            double sqrtQ = Math.sqrt(Q);
            
            if (sqrtQ > 1.0e-10) { 
                double termFactor = 1.645 / sqrtQ;

                // Calculate dG2/dXj. The derivative of sqrt(Q) w.r.t Xj is (1/2*sqrt(Q)) * dQ/dXj.
                // dQ/dXj = 2 * (Cj*Xj) * Cj, where Cj is the coefficient of Xj in Q's term.
                // Final term derivative: -1.645 * (1/sqrt(Q)) * (Cj^2 * Xj)
                
                // dG2/dX1
                J[1][0] = 12.0 - termFactor * Math.pow(0.53, 2) * x1;
                // dG2/dX2
                J[1][1] = 11.9 - termFactor * Math.pow(0.44, 2) * x.getEntry(1);
                // dG2/dX3
                J[1][2] = 41.8 - termFactor * Math.pow(4.5, 2) * x.getEntry(2); 
                // dG2/dX4
                J[1][3] = 52.1 - termFactor * Math.pow(0.79, 2) * x.getEntry(3);

            } else {
                // Corner case Q=0. Only the linear part remains.
                J[1][0] = 12.0;
                J[1][1] = 11.9;
                J[1][2] = 41.8;
                J[1][3] = 52.1;
            }
            
            return new Array2DRowRealMatrix(J);
        }
    }

    /**
     * Implementation of the equality constraint H(1).
     * H(1) = X1 + X2 + X3 + X4 - 1.0 = 0.
     */
    static final class HS353Eq extends EqualityConstraint {
        
        HS353Eq() {
            super(new ArrayRealVector(new double[NUM_EQUALITIES])); 
        }

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            // H(1) = X1 + X2 + X3 + X4 - 1.0
            double H1 = x1 + x2 + x3 + x4 - 1.0;
            
            return new ArrayRealVector(new double[] { H1 }, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            // Gradient of H(1) (constant): [1.0, 1.0, 1.0, 1.0]
            double[][] J = { { 1.0, 1.0, 1.0, 1.0 } };
            return MatrixUtils.createRealMatrix(J);
        }
    }

    // Initial starting point
    private static double[] start() { 
        return new double[]{0.0, 0.0, 0.4, 0.6}; 
    }
    
    // Expected solution for objective function and variables
    private static final double F_EXPECTED = -39.933673; 
    
    @Test
    public void testHS353() {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        
        // Box constraints: X_i >= 0
        SimpleBounds bounds = new SimpleBounds(
            new double[] { 0.0, 0.0, 0.0, 0.0 }, 
            new double[] { 
                Double.POSITIVE_INFINITY, 
                Double.POSITIVE_INFINITY, 
                Double.POSITIVE_INFINITY, 
                Double.POSITIVE_INFINITY 
            }
        );

        // Optimization run
        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS353Obj()),
                new HS353Ineq(), 
                new HS353Eq(), // Added equality constraint
                bounds
        );

        double f = sol.getValue();
        final double tolerance = 1.0e-5 * (Math.abs(F_EXPECTED) + 1.0);
        
        // Verify objective function value
        assertEquals(F_EXPECTED, f, tolerance, "objective mismatch");
        
        
    }
}
