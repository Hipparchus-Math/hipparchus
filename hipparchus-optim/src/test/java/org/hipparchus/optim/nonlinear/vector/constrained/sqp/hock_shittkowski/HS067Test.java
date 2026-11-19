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
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/** TP67 / HS67 transcription from PROB.FOR (3 vars, 14 nonlinear inequalities). */
public class HS067Test {

    private static final int DIM = 3;

    /**
     * Compute Y(2)..Y(8) exactly as TP67 does (fixed-point inner loops).
     */
    private static double[] computeY(final RealVector x) {
        final double x1 = x.getEntry(0);
        final double x2 = x.getEntry(1);
        final double x3 = x.getEntry(2);

        final double[] y = new double[9]; // use Fortran-like 1-based indexing

        final double rx = 1.0 / x1;

        // First fixed-point loop for Y2 (labels 100..102)
        y[2] = 1.6 * x1;
        for (int irep = 0; irep <= 100; ++irep) {
            y[3] = 1.22 * y[2] - x1;
            y[6] = (x2 + y[3]) * rx;

            final double v2 = (112.0 + (13.167 - 0.6667 * y[6]) * y[6]) * 0.01;
            final double y2c = x1 * v2;

            if (FastMath.abs(y2c - y[2]) <= 1.0e-3) {
                break;
            }
            y[2] = y2c;
        }

        // Second fixed-point loop for Y4 (labels 105..109)
        y[4] = 93.0;
        for (int irep = 0; irep <= 100; ++irep) {
            y[5] = 86.35 + 1.098 * y[6] - 0.038 * y[6] * y[6] + 0.325 * (y[4] - 89.0);
            y[8] = -133.0 + 3.0 * y[5];
            y[7] = 35.82 - 0.222 * y[8];

            final double denom = y[2] * y[7] + 1.0e3 * x3;
            final double y4c = 9.8e4 * x3 / denom;

            if (FastMath.abs(y4c - y[4]) <= 1.0e-4) {
                break;
            }
            y[4] = y4c;
        }

        return y;
    }

    private static class HS067Obj extends TwiceDifferentiableFunction {
        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(final RealVector x) {
            final double[] y = computeY(x);
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);

            // TP67 MODE=2
            return -(0.063 * y[2] * y[5] - 5.04 * x1 - 3.36 * y[3] - 0.035 * x2 - 10.0 * x3);
        }

        @Override
        public RealVector gradient(final RealVector x) {
            throw new UnsupportedOperationException();
        }

        @Override
        public RealMatrix hessian(final RealVector x) {
            throw new UnsupportedOperationException();
        }
    }

    private static class HS067Ineq extends InequalityConstraint {
        HS067Ineq() {
            super(new ArrayRealVector(new double[14]));
        }

        @Override
        public RealVector value(final RealVector x) {
            final double[] y = computeY(x);
            return new ArrayRealVector(new double[] {
                y[2],
                y[3],
                y[4] - 85.0,
                y[5] - 90.0,
                y[6] - 3.0,
                y[7] - 0.01,
                y[8] - 145.0,
                5000.0 - y[2],
                2000.0 - y[3],
                93.0 - y[4],
                95.0 - y[5],
                12.0 - y[6],
                4.0 - y[7],
                162.0 - y[8]
            }, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {
            throw new UnsupportedOperationException();
        }

        @Override
        public int dim() {
            return DIM;
        }
    }

    @Test
    public void testHS067() {
        final InitialGuess guess = new InitialGuess(new double[] {1745.0, 12000.0, 110.0});
        final SimpleBounds bounds = new SimpleBounds(
                new double[] {1.0e-5, 1.0e-5, 1.0e-5},
                new double[] {2000.0, 16000.0, 120.0}
        );

        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        final LagrangeSolution sol = optimizer.optimize(
                guess,
                new ObjectiveFunction(new HS067Obj()),
                new HS067Ineq(),
                bounds
        );

        HSProblemTestUtils.assertExpectedObjective(-0.116203650728e4, sol);
    }
}