/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 *
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0.
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

/**
 * HS242 (TP242)
 *
 * N    = 3
 * NILI = 0  (no linear inequalities)
 * NINL = 0  (no nonlinear inequalities)
 * NELI = 0  (no linear equalities)
 * NENL = 0  (no nonlinear equalities)
 *
 * Fortran defines residuals F(i), i=1..10, for
 *   t_i = 0.1 * i
 *
 *   F_i(x) = exp(-x1 * t_i) - exp(-x2 * t_i)
 *            - x3 * (exp(-t_i) - exp(-10 * t_i))
 *
 * Objective (MODE = 2):
 *   f(x) = sum_{i=1}^{10} F_i(x)^2
 *
 * Gradient (MODE = 3):
 *   g = 2 * Σ_i F_i * ∇F_i
 * where
 *   ∂F_i/∂x1 = -t_i * exp(-x1 * t_i)
 *   ∂F_i/∂x2 =  t_i * exp(-x2 * t_i)
 *   ∂F_i/∂x3 =  exp(-10 * t_i) - exp(-t_i)
 *
 * Bounds (MODE = 1):
 *   0 <= x_j <= 10, j = 1..3
 *
 * Reference solution (MODE = 1):
 *   x*  = (1, 10, 1)
 *   f*  = 0
 */
public class HS242Test {

    private static final int DIM      = 3;
    private static final int NUM_INEQ = 0;
    private static final int NUM_EQ   = 0;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS242Obj extends TwiceDifferentiableFunction {

        private static final int N_POINTS = 10;

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double fx = 0.0;

            for (int i = 1; i <= N_POINTS; i++) {
                double ti = 0.1 * i;

                double e1 = FastMath.exp(-x1 * ti);
                double e2 = FastMath.exp(-x2 * ti);
                double eT = FastMath.exp(-ti);
                double e10T = FastMath.exp(-10.0 * ti);

                double fi = e1 - e2 - x3 * (eT - e10T);

                fx += fi * fi;
            }

            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double g1 = 0.0;
            double g2 = 0.0;
            double g3 = 0.0;

            for (int i = 1; i <= 10; i++) {
                double ti = 0.1 * i;

                double e1 = FastMath.exp(-x1 * ti);
                double e2 = FastMath.exp(-x2 * ti);
                double eT = FastMath.exp(-ti);
                double e10T = FastMath.exp(-10.0 * ti);

                double fi = e1 - e2 - x3 * (eT - e10T);

                double dfi_dx1 = -ti * e1;
                double dfi_dx2 =  ti * e2;
                double dfi_dx3 =  e10T - eT;

                g1 += 2.0 * fi * dfi_dx1;
                g2 += 2.0 * fi * dfi_dx2;
                g3 += 2.0 * fi * dfi_dx3;
            }

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Let SQP/BFGS build an approximation; start from zero matrix.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test (bound-constrained, no inequalities/equalities)
    // -------------------------------------------------------------------------
    @Test
    public void testHS242_optimization() {

        // Initial guess (MODE=1): X = (2.5, 10, 10)
        double[] x0 = new double[]{2.5, 10.0, 10.0};

        // Bounds: 0 <= x_j <= 10
        double[] lower = new double[]{0.0, 0.0, 0.0};
        double[] upper = new double[]{10.0, 10.0, 10.0};
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS242Obj()),
                null,        // no equalities
                null,        // no inequalities
                bounds       // bounds 0 <= x_j <= 10
        );

        double f = sol.getValue();

        // Reference optimum from Fortran: FEX = 0.0 at x* = (1, 10, 1)
        final double fExpected = 0.0;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
