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
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class HS375Test {

    
    static final class HS375Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return 10;
        }

        @Override
        public double value(RealVector x) {
            double fx = 0.0;
            for (int i = 0; i < 10; i++) {
                double xi = x.getEntry(i);
                fx -= xi * xi;
            }
            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] g = new double[10];
            for (int i = 0; i < 10; i++) {
                g[i] = -2.0 * x.getEntry(i);
            }
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    
    static final class HS375Eq extends EqualityConstraint {

        HS375Eq() {
           
            super(new ArrayRealVector(new double[9]));
        }

        @Override
        public int dim() {
            return 10;
        }

        @Override
        public RealVector value(RealVector x) {
            double[] g = new double[9];

            //  1..8 (index 0..7)
            for (int j = 0; j < 8; j++) {
                double sum = 0.0;
                for (int i = 0; i < 10; i++) {
                    // TP375A(i+1, j+1) = 2.0 se i==j, 1.0
                    double aij = (i == j) ? 2.0 : 1.0;
                    sum += x.getEntry(i) / aij;
                }
                g[j] = sum - 1.0;
            }

            //  9 (index 8): sum x_i^2 / (1 + i/3) - 4 = 0
            double sumNL = 0.0;
            for (int i = 0; i < 10; i++) {
                double xi = x.getEntry(i);
                double denom = 1.0 + (double) i / 3.0; // 1 + (i)/3
                sumNL += (xi * xi) / denom;
            }
            g[8] = sumNL - 4.0;

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            double[][] J = new double[9][10];

            //  1..8: g_j = sum_i x_i / a_ij - 10
            // ∂g_j/∂x_i = 1 / a_ij
            for (int j = 0; j < 8; j++) {
                for (int i = 0; i < 10; i++) {
                    double aij = (i == j) ? 2.0 : 1.0;
                    J[j][i] = 1.0 / aij;
                }
            }

            //  9: g9 = sum_i x_i^2 / denom_i - 4
            // ∂g9/∂x_i = 2*x_i / denom_i
            for (int i = 0; i < 10; i++) {
                double xi = x.getEntry(i);
                double denom = 1.0 + (double) i / 3.0;
                J[8][i] = 2.0 * xi / denom;
            }

            return MatrixUtils.createRealMatrix(J);
        }
    }

    private LagrangeSolution solve() {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        
        double[] x0 = new double[10];
        java.util.Arrays.fill(x0, 1.0);

        return opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS375Obj()),
                 new HS375Eq() );
        
    }

    @Test
    public void testHS375() {
        
        final double fExpected = -15.161;
        double f = solve().getValue();
        assertTrue(fExpected >= f, "objective should be <= reference");
    }
}
