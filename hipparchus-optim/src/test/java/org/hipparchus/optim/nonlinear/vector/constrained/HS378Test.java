/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;


public class HS378Test {

    
    static final class HS378Obj extends TwiceDifferentiableFunction {

        
        private static final double[] A = {
            -6.089,   // -0.6089D+1
            -17.164,  // -0.17164D+2
            -34.054,  // -0.34054D+2
            -5.914,   // -0.5914D+1
            -24.721,  // -0.24721D+2
            -14.986,  // -0.14986D+2
            -24.100,  // -0.24100D+2
            -10.708,  // -0.10708D+2
            -26.662,  // -0.26662D+2
            -22.179   // -0.22179D+2
        };

        @Override
        public int dim() {
            return 10;
        }

        @Override
        public double value(RealVector x) {
            // CON = log(sum_j exp(x_j))
            double sumExp = 0.0;
            for (int j = 0; j < 10; j++) {
                sumExp += Math.exp(x.getEntry(j));
            }
            double con = Math.log(sumExp);

            double fx = 0.0;
            for (int i = 0; i < 10; i++) {
                double xi = x.getEntry(i);
                fx += Math.exp(xi) * (A[i] + xi - con);
            }
            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            // GF(i) = exp(x_i) * (A_i + x_i - CON)
            double sumExp = 0.0;
            for (int j = 0; j < 10; j++) {
                sumExp += Math.exp(x.getEntry(j));
            }
            double con = Math.log(sumExp);

            double[] g = new double[10];
            for (int i = 0; i < 10; i++) {
                double xi = x.getEntry(i);
                g[i] = Math.exp(xi) * (A[i] + xi - con);
            }
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    
    static final class HS378Eq extends EqualityConstraint {

        HS378Eq() {
            // 3 vincoli, RHS = 0
            super(new ArrayRealVector(new double[3]));
        }

        @Override
        public int dim() {
            return 10;
        }

        @Override
        public RealVector value(RealVector x) {
            double[] g = new double[3];

            double x1  = x.getEntry(0);
            double x2  = x.getEntry(1);
            double x3  = x.getEntry(2);
            double x4  = x.getEntry(3);
            double x5  = x.getEntry(4);
            double x6  = x.getEntry(5);
            double x7  = x.getEntry(6);
            double x8  = x.getEntry(7);
            double x9  = x.getEntry(8);
            double x10 = x.getEntry(9);

            // G(1)=DEXP(X(1))+0.2D+1*DEXP(X(2))+0.2D+1*DEXP(X(3))+DEXP(X(6))+DEXP(X(10))-0.2D+1
            // 0.2D+1 = 2.0
            g[0] = Math.exp(x1)
                 + 2.0 * Math.exp(x2)
                 + 2.0 * Math.exp(x3)
                 + Math.exp(x6)
                 + Math.exp(x10)
                 - 2.0;

            // G(2)=DEXP(X(4))+0.2D+1*DEXP(X(5))+DEXP(X(6))+DEXP(X(7))-0.1D+1
            // 0.1D+1 = 1.0
            g[1] = Math.exp(x4)
                 + 2.0 * Math.exp(x5)
                 + Math.exp(x6)
                 + Math.exp(x7)
                 - 1.0;

            // G(3)=DEXP(X(3))+DEXP(X(7))+DEXP(X(8))+0.2D+1*DEXP(X(9))+DEXP(X(10))-0.1D+1
            g[2] = Math.exp(x3)
                 + Math.exp(x7)
                 + Math.exp(x8)
                 + 2.0 * Math.exp(x9)
                 + Math.exp(x10)
                 - 1.0;

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            double[][] J = new double[3][10];

            double x1  = x.getEntry(0);
            double x2  = x.getEntry(1);
            double x3  = x.getEntry(2);
            double x4  = x.getEntry(3);
            double x5  = x.getEntry(4);
            double x6  = x.getEntry(5);
            double x7  = x.getEntry(6);
            double x8  = x.getEntry(7);
            double x9  = x.getEntry(8);
            double x10 = x.getEntry(9);

            // g1
            J[0][0] = Math.exp(x1);
            J[0][1] = 2.0 * Math.exp(x2);
            J[0][2] = 2.0 * Math.exp(x3);
            J[0][5] = Math.exp(x6);
            J[0][9] = Math.exp(x10);

            // g2
            J[1][3] = Math.exp(x4);
            J[1][4] = 2.0 * Math.exp(x5);
            J[1][5] = Math.exp(x6);
            J[1][6] = Math.exp(x7);

            // g3
            J[2][2] = Math.exp(x3);
            J[2][6] = Math.exp(x7);
            J[2][7] = Math.exp(x8);
            J[2][8] = 2.0 * Math.exp(x9);
            J[2][9] = Math.exp(x10);

            return MatrixUtils.createRealMatrix(J);
        }
    }

    private LagrangeSolution solve() {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        // Start: X(i) = -0.23D+1 = -2.3
        double[] x0 = new double[10];
        java.util.Arrays.fill(x0, -2.3);

        // Bounds: XL(i)=-16.0, XU(i)=-0.1
        double[] lower = new double[10];
        double[] upper = new double[10];
        for (int i = 0; i < 10; i++) {
            lower[i] = -16.0;
            upper[i] = -0.1;
        }

        return opt.optimize(
            new InitialGuess(x0),
            new ObjectiveFunction(new HS378Obj()),
            new HS378Eq() ,
            new SimpleBounds(lower, upper) // sostituisci se hai un tuo tipo bounds
        );
    }

    @Test
    public void testHS378() {
    
        final double fExpected = -47.761091;
        LagrangeSolution sol = solve();
        double f = sol.getValue();
        assertEquals(fExpected, f, 1.0e-2 * (Math.abs(fExpected) + 1.0), "objective mismatch at optimum");
    }
}
