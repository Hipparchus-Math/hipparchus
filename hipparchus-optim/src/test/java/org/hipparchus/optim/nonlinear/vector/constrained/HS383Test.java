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

package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;


import org.junit.jupiter.api.Test;

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
                final double xi = x.getEntry(i);
                fx += A[i] / xi;
            }
            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] g = new double[DIM];
            for (int i = 0; i < DIM; i++) {
                final double xi = x.getEntry(i);
                g[i] = -A[i] / (xi * xi);
            }
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            return new Array2DRowRealMatrix(DIM, DIM);
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
        
        SQPOption option =new SQPOption();
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