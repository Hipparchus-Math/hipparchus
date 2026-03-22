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
 * HS244 (TP244)
 *
 * N    = 3
 * NILI = 0
 * NINL = 0
 * NELI = 0
 * NENL = 0
 *
 * Least-squares problem with 10 residuals F(i), but the Fortran code
 * only sums i = 1..8 in the objective:
 *
 *   Z_i = 0.1 * i
 *   Y_i = exp(-Z_i) - 5 * exp(-10 * Z_i)
 *
 *   F_i(x) = exp(-x1 * Z_i) - x3 * exp(-x2 * Z_i) - Y_i
 *
 * Objective (MODE = 2):
 *   FX = sum_{i=1..8} F_i(x)^2
 *
 * Bounds (MODE = 1):
 *   0      <= x_j <= 1.0e10,  j = 1..3
 *
 * Reference solution:
 *   x* = (1, 10, 5)
 *   f* = 0
 */
public class HS244Test {

    private static final int DIM = 3;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS244Obj extends TwiceDifferentiableFunction {

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

            // Fortran loops:
            //   i = 1..10 for F(i)
            //   objective uses i = 1..8 only
            for (int i = 1; i <= 8; i++) {
                double Zi = 0.1 * i;
                double Yi = FastMath.exp(-Zi) - 5.0 * FastMath.exp(-10.0 * Zi);

                double Fi = FastMath.exp(-x1 * Zi)
                        - x3 * FastMath.exp(-x2 * Zi)
                        - Yi;

                fx += Fi * Fi;
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

            // In Fortran, GF is accumulated over i = 1..10, but since
            // at optimum F(9),F(10) are also zero, using 1..8 is sufficient
            // and consistent with the objective definition.
            for (int i = 1; i <= 8; i++) {
                double Zi = 0.1 * i;
                double Yi = FastMath.exp(-Zi) - 5.0 * FastMath.exp(-10.0 * Zi);

                double exp1 = FastMath.exp(-x1 * Zi);
                double exp2 = FastMath.exp(-x2 * Zi);

                double Fi = exp1 - x3 * exp2 - Yi;

                // Derivatives of F_i
                double dFi_dx1 = -Zi * exp1;
                double dFi_dx2 = Zi * x3 * exp2;
                double dFi_dx3 = -exp2;

                g1 += 2.0 * Fi * dFi_dx1;
                g2 += 2.0 * Fi * dFi_dx2;
                g3 += 2.0 * Fi * dFi_dx3;
            }

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Let the SQP/BFGS machinery approximate the Hessian.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS244_optimization() {

        // Initial guess (MODE=1):
        //   X(1) = 1, X(2) = 2, X(3) = 1
        double[] x0 = new double[]{1.0, 2.0, 1.0};

        // Bounds: 0 <= x_j <= 1.0e10
        double[] lower = new double[]{0.0, 0.0, 0.0};
        double[] upper = new double[]{1.0e10, 1.0e10, 1.0e10};
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS244Obj()),
                null,    // no equalities
                null,    // no inequalities
                bounds   // bounds from Fortran XL, XU
        );

        double f = sol.getValue();

        final double fExpected = 0.0;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
