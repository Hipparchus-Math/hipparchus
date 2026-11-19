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

public class HS203Test {

    /** Objective: sum_{i=1..3} (c_i - x1*(1 - x2^i))^2 with exact gradient. */
    static final class HS203Obj extends TwiceDifferentiableFunction {
        private static final double[] C = {1.5, 2.25, 2.625};
        @Override public int dim() { return 2; }

        private static double f(int i, double x1, double x2) {
            // F_i = c_i - x1 * (1 - x2^i)
            return C[i] - x1 * (1.0 - Math.pow(x2, i + 1)); // i: 0..2 -> power 1..3
        }

        @Override public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            double F1 = f(0, x1, x2);
            double F2 = f(1, x1, x2);
            double F3 = f(2, x1, x2);
            return F1*F1 + F2*F2 + F3*F3;
        }

        @Override public RealVector gradient(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            // F1, F2, F3
            final double x2_1 = x2;                 // x2^1
            final double x2_2 = x2 * x2;           // x2^2
            final double x2_3 = x2_2 * x2;         // x2^3
            final double F1 = C[0] - x1 * (1.0 - x2_1);
            final double F2 = C[1] - x1 * (1.0 - x2_2);
            final double F3 = C[2] - x1 * (1.0 - x2_3);

            // Jacobian of F wrt (x1, x2):
            // dF1/dx1 = -(1 - x2),     dF1/dx2 = x1
            // dF2/dx1 = -(1 - x2^2),   dF2/dx2 = 2*x1*x2
            // dF3/dx1 = -(1 - x2^3),   dF3/dx2 = 3*x1*x2^2
            final double dF1dx1 = -(1.0 - x2_1);
            final double dF2dx1 = -(1.0 - x2_2);
            final double dF3dx1 = -(1.0 - x2_3);
            final double dF1dx2 = x1;
            final double dF2dx2 = 2.0 * x1 * x2_1;
            final double dF3dx2 = 3.0 * x1 * x2_2;

            final double g1 = 2.0 * (F1*dF1dx1 + F2*dF2dx1 + F3*dF3dx1);
            final double g2 = 2.0 * (F1*dF1dx2 + F2*dF2dx2 + F3*dF3dx2);
            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    private static LagrangeSolution solve() {
        final SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        // No bounds, no constraints for HS203.
        final double[] start = {2.0, 0.2}; // Fortran start

        return opt.optimize(
                new InitialGuess(start),
                new ObjectiveFunction(new HS203Obj())
        );
    }

    @Test
    public void testHS203() {
        final LagrangeSolution sol = solve();
        final double f = sol.getValue();
        // Best-known minimum f* = 0 at x* = (3, 0.5)
        assertEquals(0.0, f, 1.0e-6, "objective mismatch at optimum");
    }
}
