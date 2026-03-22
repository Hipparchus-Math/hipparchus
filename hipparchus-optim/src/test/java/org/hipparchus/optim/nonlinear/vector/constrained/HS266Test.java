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

/*
 * HS266 (TP266)
 *
 * N    = 5
 * NILI = 0
 * NINL = 0
 * NELI = 0
 * NENL = 0
 *
 * F_i(x) (i = 1..10) defined via TP266A:
 *   F_i(x) = A_i + sum_{k=1..5} x_k * ( C_{i,k} + 0.5 * D_i * sum_{l=1..5} B_{k,l} x_l )
 *
 * Objective:
 *   f(x) = sum_{i=1..10} F_i(x)^2
 *
 * Bounds (MODE=1):
 *   x_i >= 0, no upper bounds
 *
 * Initial guess:
 *   x0 = (0.1, 0.1, 0.1, 0.1, 0.1)
 *
 * Reference solution (Fortran):
 *   x* ≈ (0, 0, 0.029297857, 0, 0)
 *   f* ≈ 0.99597447
 */

package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS266Test {

    private static final int DIM = 5;

    /** A(10) from Fortran DATA A. */
    private static final double[] A = {
        0.0426149,
        0.0352053,
        0.0878058,
        0.0330812,
        0.0580924,
        0.649704,
        0.344144,
        -0.627443,
        0.001828,
        -0.224783
    };

    /** D(10) from Fortran DATA D. */
    private static final double[] D = {
        2.34659,
        2.84048,
        1.13888,
        3.02286,
        1.72139,
        0.153917,
        0.290577,
        -0.159378,
        54.6910,
        -0.444873
    };

    /*
     * C(10,5) from Fortran DATA C / ... /.
     *
     * Fortran fills column-major:
     *   C(1,1),...,C(10,1), C(1,2),...,C(10,2), ..., C(1,5)..C(10,5)
     *
     * Qui: C[i][k] = C(i+1, k+1).
     */
    private static final double[][] C = {
        // i = 1: C(1,1..5)
        { -0.564255,   0.0392417, -0.404979,   0.927589,  -0.0735083 },
        // i = 2: C(2,1..5)
        {  0.535493,   0.658799,  -0.0636666, -0.681091,  -0.869487  },
        // i = 3: C(3,1..5)
        {  0.586387,   0.289826,   0.854402,   0.789312,   0.949721  },
        // i = 4: C(4,1..5)
        {  0.608734,   0.984915,   0.375699,   0.239547,   0.463136  },
        // i = 5: C(5,1..5)
        {  0.774227,   0.325421,  -0.151719,   0.448051,   0.149926  },
        // i = 6: C(6,1..5)
        { -0.435033,  -0.688583,   0.0222278, -0.524653,   0.413248  },
        // i = 7: C(7,1..5)
        {  0.759468,  -0.627795,   0.0403142,  0.724666,  -0.0182537 },
        // i = 8: C(8,1..5)
        { -0.152448,  -0.546437,   0.484134,   0.353951,   0.887866  },
        // i = 9: C(9,1..5)
        { -0.821772,  -0.53412,   -0.798498,  -0.658572,   0.662362  },
        // i =10: C(10,1..5)
        {  0.819831,  -0.910632,  -0.480344,  -0.871758,  -0.978666  }
    };

    /*
     * B(5,5) from Fortran DATA B / ... /.
     *
     * Fortran column-major; qui B[k][l] = B(k+1, l+1).
     */
    private static final double[][] B = {
        {  0.354033,   -0.0230349, -0.211938,  -0.0554288,  0.220429  },
        { -0.0230349,   0.29135,   -0.00180333,-0.111141,   0.0485461 },
        { -0.211938,   -0.00180333,-0.815808,  -0.133538,  -0.38067   },
        { -0.0554288,  -0.111141,  -0.133538,   0.389198,  -0.131586  },
        {  0.220429,    0.0485461, -0.38067,   -0.131586,   0.534706  }
    };

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------

    private static class HS266Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double[] X = x.toArray();
            double fx = 0.0;

            // F(i) = A_i + sum_k X_k * ( C_{i,k} + 0.5 * D_i * sum_l B_{k,l} X_l )
            for (int i = 0; i < 10; i++) {
                double fi = A[i];

                for (int k = 0; k < 5; k++) {
                    double hf = 0.0;
                    for (int l = 0; l < 5; l++) {
                        hf += B[k][l] * X[l];
                    }
                    fi += X[k] * (C[i][k] + 0.5 * D[i] * hf);
                }

                fx += fi * fi;
            }

            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] X = x.toArray();
            double[] grad = new double[DIM];

            // Prima ricalcoliamo F(i)
            double[] F = new double[10];
            for (int i = 0; i < 10; i++) {
                double fi = A[i];
                for (int k = 0; k < 5; k++) {
                    double hf = 0.0;
                    for (int l = 0; l < 5; l++) {
                        hf += B[k][l] * X[l];
                    }
                    fi += X[k] * (C[i][k] + 0.5 * D[i] * hf);
                }
                F[i] = fi;
            }

            // Fortran:
            // HF = sum_l (B(K,L)+B(L,K)) * X(L)
            // DF(I,K) = C(I,K) + 0.5*D(I)*HF
            // GF(K)   = 2 * sum_I F(I) * DF(I,K)
            for (int k = 0; k < 5; k++) {
                double gk = 0.0;

                for (int i = 0; i < 10; i++) {
                    double hf = 0.0;
                    for (int l = 0; l < 5; l++) {
                        hf += (B[k][l] + B[l][k]) * X[l];
                    }
                    double dfik = C[i][k] + 0.5 * D[i] * hf;
                    gk += 2.0 * F[i] * dfik;
                }

                grad[k] = gk;
            }

            return new ArrayRealVector(grad, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Hessiano reale complicato; lasciamo che il BFGS lo costruisca.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS266() {

        // Initial guess: X(I) = 0.1
        double[] x0 = new double[]{0.1, 0.1, 0.1, 0.1, 0.1};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // Bounds: x_i >= 0 (XL(I)=0, LXL(I)=.TRUE., LXU(I)=.FALSE.)
        SimpleBounds bounds = new SimpleBounds(
            new double[]{0.0, 0.0, 0.0, 0.0, 0.0},
            new double[]{
                Double.POSITIVE_INFINITY,
                Double.POSITIVE_INFINITY,
                Double.POSITIVE_INFINITY,
                Double.POSITIVE_INFINITY,
                Double.POSITIVE_INFINITY
            }
        );

        LagrangeSolution sol = opt.optimize(
            new InitialGuess(x0),
            new ObjectiveFunction(new HS266Obj()),
            null,   // no equalities
            null,   // no inequalities
            bounds
        );

        double f = sol.getValue();

        double fExpected = 0.99597447;
        double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
