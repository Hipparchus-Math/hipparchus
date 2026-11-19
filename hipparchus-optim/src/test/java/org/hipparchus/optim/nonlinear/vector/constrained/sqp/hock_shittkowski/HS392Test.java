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

package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS392Test{

    private static int ix(int oneBased) { return oneBased - 1; }

    
    private static final double[][] R1 = new double[][]{
        {1000.0, 1000.0, 1000.0, 1100.0, 1100.0}, // j=1
        { 520.0,  520.0,  520.0,  600.0,  600.0}, // j=2
        { 910.0,  910.0, 1000.0, 1000.0, 1000.0}  // j=3
    };

    private static final double[][] R2 = new double[][]{
        {0.3, 0.3, 0.3, 0.3, 0.3},
        {0.1, 0.1, 0.1, 0.1, 0.1},
        {0.2, 0.2, 0.2, 0.2, 0.2}
    };

    private static final double[][] KA = new double[][]{
        {120.0, 150.0, 150.0, 170.0, 170.0},
        { 65.0,  65.0,  80.0,  80.0,  80.0},
        {105.0, 105.0, 120.0, 120.0, 120.0}
    };

    private static final double[][] K1C = new double[][]{
        {150.0, 150.0, 150.0, 170.0, 170.0},
        { 75.0,  75.0,  75.0,  90.0,  90.0},
        {140.0, 140.0, 140.0, 150.0, 150.0}
    };

    private static final double[][] KP = new double[][]{
        {160.0, 160.0, 160.0, 180.0, 180.0},
        { 75.0,  75.0,  75.0,  90.0,  90.0},
        {140.0, 140.0, 140.0, 150.0, 150.0}
    };

    private static final double[][] K3 = new double[][]{
        {0.02, 0.20, 0.25, 0.25, 0.25},
        {0.01, 0.10, 0.10, 0.15, 0.15},
        {0.015,0.15, 0.15, 0.15, 0.15}
    };

    private static final double[][] KL1 = new double[][]{
        {0.005, 0.05, 0.6, 0.6, 0.6},
        {0.005, 0.05, 0.6, 0.6, 0.6},
        {0.005, 0.05, 0.6, 0.6, 0.6}
    };

    private static final double[][] KL2 = new double[][]{
        { 80.0,  80.0, 100.0, 100.0, 100.0},
        { 45.0,  45.0,  45.0,  50.0,  50.0},
        { 75.0,  75.0,  90.0,  90.0,  90.0}
    };

    private static final double[][] H = new double[][]{
        {100.0, 180.0, 220.0, 150.0, 100.0},
        {280.0, 400.0, 450.0, 450.0, 400.0},
        {520.0, 400.0, 500.0, 630.0, 600.0}
    };

    
    private static final double[][] T = new double[][]{
        {0.6, 0.4, 0.1},
        {0.3, 0.1, 0.12},
        {0.36,0.08,0.06}
    };

    
    private static double B(int j, int i) { return (j == 3) ? 180.0 : 170.0; }

    
    private static class TP392Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 30; }

        @Override public double value(RealVector Xv) {
            final double[] x = Xv.toArray();
            double FX = 0.0;

            for (int i = 1; i <= 5; i++) {
                double SUM = 0.0;
                for (int j = 1; j <= 3; j++) {
                    final int Lin  = 3*(i-1) + j;        // 1..15
                    final int Lout = 12 + 3*i + j;      // 16..30

                    // SUM1 = Σ_{K=1..i} (x_out(K,j) − x_in(K,j))
                    double SUM1 = 0.0;
                    for (int K = 1; K <= i; K++) {
                        final int LinK  = 3*(K-1) + j;
                        final int LoutK = 12 + 3*K + j;
                        SUM1 += x[ix(LoutK)] - x[ix(LinK)];
                    }

                    final double xin  = x[ix(Lin)];
                    final double xout = x[ix(Lout)];

                    final double term =
                            xin * (R1[j-1][i-1] - KA[j-1][i-1])
                          - xin * xin * R2[j-1][i-1]
                          - xout * (K1C[j-1][i-1] + KP[j-1][i-1])
                          - (xout - xin)*(xout - xin) * (K3[j-1][i-1] + KL1[j-1][i-1])
                          - KL2[j-1][i-1] * SUM1;

                    SUM += term;
                }
                FX -= SUM;
            }
            return FX;
        }

        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

   
    private static class TP392IneqAll extends InequalityConstraint {
        TP392IneqAll() { super(new ArrayRealVector(new double[45])); }
        @Override public int dim() { return 30; }

        @Override public RealVector value(RealVector Xv) {
            final double[] x = Xv.toArray();
            final double[] g = new double[45];
            int p = 0;

            // (1) H(j,i) − x_in(i,j) ≥ 0
            for (int i = 1; i <= 5; i++) {
                for (int j = 1; j <= 3; j++) {
                    final int Lin = 3*(i-1) + j;
                    g[p++] = H[j-1][i-1] - x[ix(Lin)];
                }
            }

            // (2) B(j,i) − T(j,:)*x_out(i,:) ≥ 0
            for (int i = 1; i <= 5; i++) {
                for (int j = 1; j <= 3; j++) {
                    final int Lout1 = 12 + 3*i + 1;
                    final int Lout2 = 12 + 3*i + 2;
                    final int Lout3 = 12 + 3*i + 3;
                    final double dot =
                          T[j-1][0]*x[ix(Lout1)]
                        + T[j-1][1]*x[ix(Lout2)]
                        + T[j-1][2]*x[ix(Lout3)];
                    g[p++] = B(j,i) - dot;
                }
            }

            // (3) Σ_{K=1..i} (x_out(K,j) − x_in(K,j)) ≥ 0
            for (int i = 1; i <= 5; i++) {
                for (int j = 1; j <= 3; j++) {
                    double sum = 0.0;
                    for (int K = 1; K <= i; K++) {
                        final int LinK  = 3*(K-1) + j;
                        final int LoutK = 12 + 3*K + j;
                        sum += x[ix(LoutK)] - x[ix(LinK)];
                    }
                    g[p++] = sum;
                }
            }

            return new ArrayRealVector(g);
        }

        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testTP392_corretto() {
        
        final double[] x0 = {
            80,100,400, 100,200,200, 100,250,400, 50,200,500,
            50,200,500, 100,120,410, 120,250,250, 150,300,410,
            600,250,510, 100,250,510
        };

        final double[] lb = new double[30];
        final double[] ub = new double[30];
        for (int i = 0; i < 30; i++) {
            lb[i] = 0.0;
            ub[i] = Double.POSITIVE_INFINITY;
        }
        
        final SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        SQPOption option=new SQPOption();
        option.setGradientMode(GradientMode.FORWARD);
        final LagrangeSolution sol = opt.optimize(
            new InitialGuess(x0),
            new ObjectiveFunction(new TP392Obj()),
            new TP392IneqAll(),
            new SimpleBounds(lb, ub)
            
        );

        
        HSProblemTestUtils.assertExpectedObjective(-1.6960671e6, sol);
    }
}
