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
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

/**
 * HS383 (TP383) – Linear-fractional type problem (14 variables) with
 * one linear equality constraint and simple bounds.
 *
 * From TP383:
 *
 * N    = 14
 * NILI = 0
 * NINL = 0
 * NELI = 1
 * NENL = 0
 *
 * Objective:
 *   f(x) = sum_{i=1..14} A(i) / x(i)
 *
 * where:
 *
 *   A = {
 *     12842.275, 634.25, 634.25, 634.125,
 *     1268.0,   633.875, 633.75, 1267.0,
 *     760.05,   633.25,  1266.25, 632.875,
 *     394.46,   940.838
 *   }
 *
 * Bounds (MODE=1):
 *   X(i)  initial = 0.1D-1 = 0.01
 *   LXL(i) = TRUE, LXU(i) = TRUE
 *   XL(i)  = 0.1D-3 = 1.0e-4
 *   XU(i)  = 0.1D+1 / B(i) = 1.0 / B(i)
 *
 *   with
 *   B = {
 *     25, 26, 26, 27, 28, 29, 30, 32,
 *     33, 34, 35, 37, 38, 36
 *   }
 *
 * Equality constraint (MODE=4):
 *
 *   G1(x) = sum_{i=1..14} C(i) * X(i) - 0.1D+1 = 0
 *         = sum C(i) * X(i) - 1.0 = 0
 *
 *   with
 *   C = {
 *     5.47934, 0.83234, 0.94749, 1.11082,
 *     2.64824, 1.55868, 1.73215, 3.90896,
 *     2.74284, 2.60541, 5.96184, 3.29522,
 *     1.83517, 2.81372
 *   }
 *
 * Reference:
 *   LEX = .FALSE.
 *   FEX = 0.728566D+6 = 728566.0
 *   → we use FEX as an upper bound on f.
 */
public class HS383Test {

    private static final int DIM    = 14;
    private static final int NUM_EQ = 1;

    // XL(i) = 0.1D-3 = 1e-4
    private static final double XL = 1.0e-4;

    private static final double[] A = {
        0.12842275e5, 0.63425e3, 0.63425e3, 0.634125e3,
        0.1268e4,     0.633875e3, 0.63375e3, 0.1267e4,
        0.76005e3,    0.63325e3,  0.126625e4, 0.632875e3,
        0.39446e3,    0.940838e3
    };

    private static final double[] B = {
        0.25e2, 0.26e2, 0.26e2, 0.27e2,
        0.28e2, 0.29e2, 0.30e2, 0.32e2,
        0.33e2, 0.34e2, 0.35e2, 0.37e2,
        0.38e2, 0.36e2
    };

    private static final double[] C = {
        0.547934e1, 0.83234e0, 0.94749e0, 0.111082e1,
        0.264824e1, 0.155868e1, 0.173215e1, 0.390896e1,
        0.274284e1, 0.260541e1, 0.596184e1, 0.329522e1,
        0.183517e1, 0.281372e1
    };

    private static class HS383Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double fx = 0.0;
            for (int i = 0; i < DIM; i++) {
                final double xi = FastMath.max(x.getEntry(i), XL);
                fx += A[i] / xi;
            }
            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] g = new double[DIM];
            for (int i = 0; i < DIM; i++) {
                final double xi = x.getEntry(i);
                g[i] = xi <= XL ? 0.0 : -A[i] / (xi * xi);
            }
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException();
        }
    }

    private static class HS383Eq extends EqualityConstraint {

        HS383Eq() {
            super(new ArrayRealVector(new double[NUM_EQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {
            double sum = 0.0;
            for (int i = 0; i < DIM; i++) {
                sum += C[i] * x.getEntry(i);
            }
            double g1 = sum - 1.0;
            return new ArrayRealVector(new double[] { g1 }, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            RealMatrix J = new Array2DRowRealMatrix(NUM_EQ, DIM);
            for (int j = 0; j < DIM; j++) {
                J.setEntry(0, j, C[j]);
            }
            return J;
        }
    }

    @Test
    @Disabled // disabled as we reach a local minimum and not the expected global one
    public void testHS383_optimization() {

        double[] x0 = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            x0[i] = 0.01;
        }

        double[] lower = new double[DIM];
        double[] upper = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            lower[i] = XL;
            upper[i] = 1.0 / B[i];
        }
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        
        SQPOption option = new SQPOption();
        option.setGradientMode(GradientMode.FORWARD);
        LagrangeSolution sol = opt.optimize(
                option,
                new InitialGuess(x0),
                new ObjectiveFunction(new HS383Obj()),
                new HS383Eq(),
                bounds
        );

        HSProblemTestUtils.assertBetterObjective(0.728566e6, sol);
    }
}
