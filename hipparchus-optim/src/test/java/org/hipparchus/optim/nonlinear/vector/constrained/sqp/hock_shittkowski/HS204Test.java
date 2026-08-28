/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 */
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS204Test {

    /** Objective: sum_i Fi(x)^2, with Fi = Ai + p + 0.5 * Di * p^2 and p = Hi·x. */
    static final class HS204Obj extends TwiceDifferentiableFunction {
        private static final double[] A = { 0.13294, -0.244378, 0.325895 };
        private static final double[] D = { 2.5074, -1.36401, 1.02282 };
        // Correct mapping from Fortran DATA (column-major) to Java rows:
        private static final double[][] H = {
                { -0.564255,  0.392417 },
                { -0.404979,  0.927589 },
                { -0.0735084, 0.535493 }
        };

        @Override public int dim() { return 2; }

        @Override public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            double f = 0.0;
            for (int i = 0; i < 3; i++) {
                double p  = H[i][0] * x1 + H[i][1] * x2;
                double Fi = A[i] + p + 0.5 * D[i] * p * p;
                f += Fi * Fi;
            }
            return f;
        }

        @Override public RealVector gradient(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            double g1 = 0.0, g2 = 0.0;
            for (int i = 0; i < 3; i++) {
                double p  = H[i][0] * x1 + H[i][1] * x2;
                double Fi = A[i] + p + 0.5 * D[i] * p * p;
                double c  = 1.0 + D[i] * p;             // dFi/dp
                g1 += 2.0 * Fi * c * H[i][0];
                g2 += 2.0 * Fi * c * H[i][1];
            }
            return new ArrayRealVector(new double[]{ g1, g2 }, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            return new Array2DRowRealMatrix(2, 2); // not used
        }
    }

    private static LagrangeSolution solve() {
        final SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        final double[] x0 = { 0.1, 0.1 };
        return opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS204Obj())
        );
    }

    @Test
    public void testHS204() {
        final LagrangeSolution sol = solve();
        final double expected = 0.183601;
        final double f = sol.getValue();
        assertEquals(expected, f, 1.0e-6 * (Math.abs(expected) + 1.0), "objective mismatch at optimum");
    }
}
