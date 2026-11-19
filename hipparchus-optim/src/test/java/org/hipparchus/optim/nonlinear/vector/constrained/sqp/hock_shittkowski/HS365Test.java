/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;



public class HS365Test {

    private static final int DIM = 7;
    private static final int NUM_INEQ = 5;

   
    private static class HS365Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x3 = x.getEntry(2);
            return x1 * x3;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] g = new double[DIM];
            double x1 = x.getEntry(0);
            double x3 = x.getEntry(2);

            // df/dx1 = x3
            g[0] = x3;
            // df/dx3 = x1
            g[2] = x1;
            
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            RealMatrix h = new Array2DRowRealMatrix(DIM, DIM);
            // d²f/(dx1 dx3) = 1, 
            h.setEntry(0, 2, 1.0);
            h.setEntry(2, 0, 1.0);
            return h;
        }
    }

   
    private static class HS365Ineq extends InequalityConstraint {

        HS365Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQ])); // RHS = 0
        }

        @Override
        public int dim() {
            return DIM;
        }

        
        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);
            double x5 = x.getEntry(4);
            double x6 = x.getEntry(5);
            double x7 = x.getEntry(6);

            double P = Math.sqrt(x2 * x2 + x3 * x3);
            double Q = Math.sqrt(x3 * x3 + (x2 - x1) * (x2 - x1));

          
            if (P == 0.0) {
                P = 1e-16;
            }
            if (Q == 0.0) {
                Q = 1e-16;
            }

            double[] g = new double[NUM_INEQ];

            // G1 = (x4 - x6)^2 + (x5 - x7)^2 - 4
            g[0] = (x4 - x6) * (x4 - x6) + (x5 - x7) * (x5 - x7) - 4.0;

            // G2 = (x3*x4 - x2*x5) / P - 1
            g[1] = (x3 * x4 - x2 * x5) / P - 1.0;

            // G3 = (x3*x6 - x2*x7) / P - 1
            g[2] = (x3 * x6 - x2 * x7) / P - 1.0;

            // G4 = (x1*x3 + (x2 - x1)*x5 - x3*x4) / Q - 1
            g[3] = (x1 * x3 + (x2 - x1) * x5 - x3 * x4) / Q - 1.0;

            // G5 = (x1*x3 + (x2 - x1)*x7 - x3*x6) / Q - 1
            g[4] = (x1 * x3 + (x2 - x1) * x7 - x3 * x6) / Q - 1.0;

            return new ArrayRealVector(g, false);
        }

       
        @Override
        public RealMatrix jacobian(RealVector x) {
            final double eps = 1.0e-6;
            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            RealVector g0 = value(x);

            for (int j = 0; j < DIM; j++) {
                double xj = x.getEntry(j);

                // x + eps e_j
                x.setEntry(j, xj + eps);
                RealVector gp = value(x);

                // x - eps e_j
                x.setEntry(j, xj - eps);
                RealVector gm = value(x);

                // ripristina x(j)
                x.setEntry(j, xj);

                // derivata centrale
                for (int i = 0; i < NUM_INEQ; i++) {
                    double dij = (gp.getEntry(i) - gm.getEntry(i)) / (2.0 * eps);
                    J.setEntry(i, j, dij);
                }
            }

            return J;
        }
    }

    @Test
    public void testHS365_optimization() {

        
        double[] x0 = {
            3.0,    // X(1)
            0.0,    // X(2)
            2.0,    // X(3)
            -1.5,   // X(4)
            1.5,    // X(5)
            5.0,    // X(6)
            1.0     // X(7)
        };

        
        double[] lower = new double[] {
            0.0,                     // X(1) >= 0
            Double.NEGATIVE_INFINITY, // X(2) no LB
            0.0,                     // X(3) >= 0
            Double.NEGATIVE_INFINITY, // X(4)
            1.0,                     // X(5) >= 1
            Double.NEGATIVE_INFINITY, // X(6)
            1.0                      // X(7) >= 1
        };

        double[] upper = new double[] {
            Double.POSITIVE_INFINITY, // no UB
            Double.POSITIVE_INFINITY,
            Double.POSITIVE_INFINITY,
            Double.POSITIVE_INFINITY,
            Double.POSITIVE_INFINITY,
            Double.POSITIVE_INFINITY,
            Double.POSITIVE_INFINITY
        };

        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
            new InitialGuess(x0),
            new ObjectiveFunction(new HS365Obj()),
            new HS365Ineq(),
            bounds
        );

        double f = sol.getValue();

        // FEX = 0.23313708D+2, LEX = .FALSE.  → FEX >= f
        final double fExpected = 0.23313708e2; // 23.313708

        final double tolF = 1.0e-4 * (FastMath.abs(fExpected) + 1.0);
        assertTrue(fExpected+tolF >= f,
                   "HS362: expected F <= " + fExpected + " but got F = " + f);
        //assertEquals(fExpected, f, tolF, "HS361: objective mismatch");
    }
}
