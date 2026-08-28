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
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * HS379 (TP379) – Nonlinear least-squares fit with 11 parameters and 65 data points.
 *
 * From TP379:
 *
 *   N     = 11
 *   NILI  = 0
 *   NINL  = 0
 *   NELI  = 0
 *   NENL  = 0
 *   LSUM  = 65
 *
 * Parameters: x1..x11.
 *
 * Data Y(i), i = 1..65, at t_i = 0.1 * (i-1).
 *
 * Model:
 *
 *   M(t) = x1 * exp(-x5 * t)
 *        + x2 * exp(-x6 * (t - x9 )^2)
 *        + x3 * exp(-x7 * (t - x10)^2)
 *        + x4 * exp(-x8 * (t - x11)^2)
 *
 * Residuals:
 *   F(i) = Y(i) - M(t_i)
 *
 * Objective:
 *   FX = sum_{i=1..65} F(i)^2
 *
 * Bounds:
 *   In MODE=1:
 *     LXL(i) = .TRUE., XL(i) = 0.0
 *     LXU(i) = .FALSE.
 *   → 0 <= x_i, no explicit upper bound.
 *
 * Fortran in MODE=2 re-clamps x: x(i) = min(max(x(i), XL(i)), XU(i)), ma XU non è impostato qui;
 * nel test Java imponiamo solo il lower bound 0 tramite SimpleBounds.
 *
 * Gradient (MODE=3 in TP379):
 *
 *   DF(i,1)  = -exp(-x5 * t)
 *   DF(i,2)  = -exp(-x6 * (t - x9 )^2)
 *   DF(i,3)  = -exp(-x7 * (t - x10)^2)
 *   DF(i,4)  = -exp(-x8 * (t - x11)^2)
 *   DF(i,5)  =  x1 * t              * exp(-x5 * t)
 *   DF(i,6)  =  x2 * (t - x9 )^2    * exp(-x6 * (t - x9 )^2)
 *   DF(i,7)  =  x3 * (t - x10)^2    * exp(-x7 * (t - x10)^2)
 *   DF(i,8)  =  x4 * (t - x11)^2    * exp(-x8 * (t - x11)^2)
 *   DF(i,9)  = -2*x2*x6*(t - x9 )   * exp(-x6 * (t - x9 )^2)
 *   DF(i,10) = -2*x3*x7*(t - x10)   * exp(-x7 * (t - x10)^2)
 *   DF(i,11) = -2*x4*x8*(t - x11)   * exp(-x8 * (t - x11)^2)
 *
 *   GF(j) = 2 * sum_{i=1..65} F(i) * DF(i,j)
 *
 * Reference:
 *   LEX = .FALSE.
 *   FEX = 0.401377D-1 = 0.0401377
 *   → usiamo FEX come upper bound su f.
 */
public class HS379Test {

    private static final int DIM   = 11;
    private static final int NDATA = 65;

    // Y(1..65) from DATA statement
    private static final double[] Y = {
        1.366, 1.191, 1.112, 1.013, 0.991,
        0.885, 0.831, 0.847, 0.786, 0.725, 0.746,
        0.679, 0.608, 0.655, 0.616, 0.606, 0.602,
        0.626, 0.651, 0.724, 0.649, 0.649, 0.694,
        0.644, 0.624, 0.661, 0.612, 0.558, 0.533,
        0.495, 0.500, 0.423, 0.395, 0.375, 0.372,
        0.391, 0.396, 0.405, 0.428, 0.429, 0.523,
        0.562, 0.607, 0.653, 0.672, 0.708, 0.633,
        0.668, 0.645, 0.632, 0.591, 0.559, 0.597,
        0.625, 0.739, 0.710, 0.729, 0.720, 0.636,
        0.581, 0.428, 0.292, 0.162, 0.098, 0.054
    };

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS379Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        /**
         * Compute residuals F(i) = Y(i) - M(t_i).
         * t_i = 0.1 * (i-1), i=0..64.
         */
        private double[] computeResiduals(double[] x) {
            double[] F = new double[NDATA];

            final double x1  = x[0];
            final double x2  = x[1];
            final double x3  = x[2];
            final double x4  = x[3];
            final double x5  = x[4];
            final double x6  = x[5];
            final double x7  = x[6];
            final double x8  = x[7];
            final double x9  = x[8];
            final double x10 = x[9];
            final double x11 = x[10];

            for (int i = 0; i < NDATA; i++) {
                double t = 0.1 * i;

                double term1 = x1 * FastMath.exp(-x5 * t);
                double dt9   = t - x9;
                double dt10  = t - x10;
                double dt11  = t - x11;

                double term2 = x2 * FastMath.exp(-x6  * dt9  * dt9);
                double term3 = x3 * FastMath.exp(-x7  * dt10 * dt10);
                double term4 = x4 * FastMath.exp(-x8  * dt11 * dt11);

                double model = term1 + term2 + term3 + term4;

                F[i] = Y[i] - model;
            }

            return F;
        }

        @Override
        public double value(RealVector xVec) {
            double[] x = xVec.toArray();
            // Bound handling in Fortran (MODE=2) clippa X a [XL,XU];
            // qui imponiamo solo >=0 via SimpleBounds nel test.

            double[] F = computeResiduals(x);

            double fx = 0.0;
            for (int i = 0; i < NDATA; i++) {
                fx += F[i] * F[i];
            }
            return fx;
        }

        @Override
        public RealVector gradient(RealVector xVec) {
            double[] x = xVec.toArray();

            final double x1  = x[0];
            final double x2  = x[1];
            final double x3  = x[2];
            final double x4  = x[3];
            final double x5  = x[4];
            final double x6  = x[5];
            final double x7  = x[6];
            final double x8  = x[7];
            final double x9  = x[8];
            final double x10 = x[9];
            final double x11 = x[10];

            double[] F = new double[NDATA];
            double[][] DF = new double[NDATA][DIM];

            // Reproduce exactly MODE=2+3 of Fortran:
            for (int i = 0; i < NDATA; i++) {
                double t = 0.1 * i;

                double dt9   = t - x9;
                double dt10  = t - x10;
                double dt11  = t - x11;

                double e1 = FastMath.exp(-x5 * t);
                double e2 = FastMath.exp(-x6 * dt9  * dt9);
                double e3 = FastMath.exp(-x7 * dt10 * dt10);
                double e4 = FastMath.exp(-x8 * dt11 * dt11);

                double model = x1 * e1 +
                               x2 * e2 +
                               x3 * e3 +
                               x4 * e4;

                F[i] = Y[i] - model;

                // DF(i,j) = ∂F/∂x_j
                DF[i][0] = -e1;                             // dF/dx1
                DF[i][1] = -e2;                             // dF/dx2
                DF[i][2] = -e3;                             // dF/dx3
                DF[i][3] = -e4;                             // dF/dx4

                DF[i][4] =  x1 * t      * e1;               // dF/dx5
                DF[i][5] =  x2 * dt9*dt9  * e2;             // dF/dx6
                DF[i][6] =  x3 * dt10*dt10 * e3;            // dF/dx7
                DF[i][7] =  x4 * dt11*dt11 * e4;            // dF/dx8

                DF[i][8]  = -x2 * x6 * 2.0 * dt9  * e2;     // dF/dx9
                DF[i][9]  = -x3 * x7 * 2.0 * dt10 * e3;     // dF/dx10
                DF[i][10] = -x4 * x8 * 2.0 * dt11 * e4;     // dF/dx11
            }

            double[] g = new double[DIM];
            // GF(j) = 2 * sum_i F(i) * DF(i,j)
            for (int j = 0; j < DIM; j++) {
                double gj = 0.0;
                for (int i = 0; i < NDATA; i++) {
                    gj += 2.0 * F[i] * DF[i][j];
                }
                g[j] = gj;
            }

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Non fornita nel Fortran; usiamo Hessian nullo.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS379_optimization() {

        // Initial point from MODE=1:
        double[] x0 = new double[DIM];
        x0[0]  = 1.3;
        x0[1]  = 0.65;
        x0[2]  = 0.65;
        x0[3]  = 0.7;
        x0[4]  = 0.6;
        x0[5]  = 3.0;
        x0[6]  = 5.0;
        x0[7]  = 7.0;
        x0[8]  = 2.0;
        x0[9]  = 4.5;
        x0[10] = 5.5;

        // Bounds: XL(i)=0, LXL(i)=TRUE, LXU(i)=FALSE → x_i >= 0, no upper bound.
        double[] lower = new double[DIM];
        double[] upper = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            lower[i] = 0.0;
            upper[i] = Double.POSITIVE_INFINITY;
        }
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS379Obj()),
                bounds
        );

        double f = sol.getValue();

        // LEX = .FALSE., FEX = 0.401377D-1 → FEX as upper bound
        final double fExpected = 0.401377e-1;
        final double tolF      = 1.0e-4 * (FastMath.abs(fExpected) + 1.0);

        assertTrue(fExpected + tolF >= f,
                   "HS379: expected F <= " + (fExpected + tolF) + " but got F = " + f);
    }
}
