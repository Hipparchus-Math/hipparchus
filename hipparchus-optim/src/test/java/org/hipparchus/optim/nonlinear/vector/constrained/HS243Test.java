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
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * HS243 (TP243)
 *
 * N = 3
 * No constraints (unconstrained least-squares).
 *
 * Residuals:
 *   For i = 1..4:
 *   F_i(x) = A[i] + E[i,1] x1 + E[i,2] x2 + E[i,3] x3 + 0.5 * XBX * D[i]
 *
 * where
 *   XBX = xᵀ B x
 *
 * Objective:
 *   f(x) = sum_{i=1..4} F_i(x)^2
 *
 * Reference optimum:
 *   x* = (0,0,0)
 *   f* = 0.79657853
 */
public class HS243Test {

    private static final int DIM = 3;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS243Obj extends TwiceDifferentiableFunction {

        // A(1..4)
        private static final double[] A = {
                0.14272,
                -0.184981,
                -0.521869,
                -0.685306
        };

        // D(1..4)
        private static final double[] D = {
                1.75168,
                -1.35195,
                -0.479048,
                -0.3648
        };

        // B 3×3
        private static final double[][] B = {
                {  2.95137,  4.87407, -2.0506 },
                {  4.87407,  9.39321, -3.93189 },
                { -2.0506,  -3.93189,  2.64745 }
        };

        // E 4×3
        private static final double[][] E = {
                {-0.564255,   0.392417,  -0.404979},
                { 0.927589,  -0.0735083,  0.535493},
                { 0.658799,  -0.636666,  -0.681091},
                {-0.869487,   0.586387,   0.289826}
        };

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            // Compute XBX = xᵀ B x
            double XBX =
                    (x1 * B[0][0] + x2 * B[1][0] + x3 * B[2][0]) * x1 +
                    (x1 * B[0][1] + x2 * B[1][1] + x3 * B[2][1]) * x2 +
                    (x1 * B[0][2] + x2 * B[1][2] + x3 * B[2][2]) * x3;

            double fx = 0.0;

            for (int i = 0; i < 4; i++) {
                double Fi = A[i]
                        + E[i][0] * x1
                        + E[i][1] * x2
                        + E[i][2] * x3
                        + 0.5 * XBX * D[i];

                fx += Fi * Fi;
            }

            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            // XBX derivatives
            double dXBXdx1 = 2.0 * (x1 * B[0][0] + x2 * B[1][0] + x3 * B[2][0]);
            double dXBXdx2 = 2.0 * (x1 * B[0][1] + x2 * B[1][1] + x3 * B[2][1]);
            double dXBXdx3 = 2.0 * (x1 * B[0][2] + x2 * B[1][2] + x3 * B[2][2]);

            double g1 = 0.0;
            double g2 = 0.0;
            double g3 = 0.0;

            for (int i = 0; i < 4; i++) {

                double Fi = A[i]
                        + E[i][0] * x1
                        + E[i][1] * x2
                        + E[i][2] * x3
                        + 0.5 * (x.dotProduct(new ArrayRealVector(B[0]))
                                  * x1 +
                                  x.dotProduct(new ArrayRealVector(B[1]))
                                  * x2 +
                                  x.dotProduct(new ArrayRealVector(B[2]))
                                  * x3) * D[i]; // replaced by XBX later if needed

                double dFi_dx1 = E[i][0] + 0.5 * dXBXdx1 * D[i];
                double dFi_dx2 = E[i][1] + 0.5 * dXBXdx2 * D[i];
                double dFi_dx3 = E[i][2] + 0.5 * dXBXdx3 * D[i];

                g1 += 2.0 * Fi * dFi_dx1;
                g2 += 2.0 * Fi * dFi_dx2;
                g3 += 2.0 * Fi * dFi_dx3;
            }

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Let SQP/BFGS approximate it.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS243_optimization() {

        // Initial guess (MODE=1): X = (0.1, 0.1, 0.1)
        double[] x0 = new double[]{0.1, 0.1, 0.1};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS243Obj()),
                null,   // no equalities
                null,   // no inequalities
                null    // no bounds
        );

        double f = sol.getValue();

        final double fExpected = 0.79657853;
        final double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
