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
 * HS252 (TP252)
 *
 * N    = 3
 * NILI = 0
 * NINL = 0
 * NELI = 0
 * NENL = 1  (one nonlinear inequality)
 *
 * Objective (MODE=2):
 *   f(x) = 0.01 (x1 - 1)^2 + (x2 - x1^2)^2
 *
 * Fortran inequality (G(i) ≥ 0):
 *   G1(x) = x1 + x3^2 + 1 ≥ 0
 *
 * We keep the same ≥ 0 convention used by the original Fortran
 * and by the SQP optimizer: the constraint function returned by
 * {@link InequalityConstraint#value(RealVector)} is G1 itself.
 *
 * Bounds (from MODE=1):
 *   LXU(1) = .TRUE., XU(1) = -1
 *   all other LXL/LXU = .FALSE.
 *
 * So in our wrapper:
 *   -∞ < x1 ≤ -1
 *   x2, x3 free
 *
 * Reference solution (MODE=1):
 *   x*  = (-1, 1, 0)
 *   f*  = 0.04
 */
public class HS252Test {

    private static final int DIM      = 3;
    private static final int NUM_INEQ = 1;
    private static final int NUM_EQ   = 0;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS252Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            // x3 does not appear in the objective
            return 0.01 * FastMath.pow(x1 - 1.0, 2)
                 + FastMath.pow(x2 - x1 * x1, 2);
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // From Fortran:
            // GF(1) = 0.02*(x1 - 1) - 4*(x2 - x1^2)*x1
            // GF(2) = 2*(x2 - x1^2)
            // GF(3) = 0
            double g1 = 0.02 * (x1 - 1.0) - 4.0 * (x2 - x1 * x1) * x1;
            double g2 = 2.0 * (x2 - x1 * x1);
            double g3 = 0.0;

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // f(x) = 0.01 (x1 - 1)^2 + (x2 - x1^2)^2
            //
            // ∂f/∂x1 = 0.02 (x1 - 1) - 4 x1 (x2 - x1^2)
            // ∂f/∂x2 = 2 (x2 - x1^2)
            //
            // Second derivatives:
            // d²f/dx1² = 0.02 + d/dx1[-4 x1 (x2 - x1^2)]
            //          = 0.02 + (-4 x2 + 12 x1^2)
            //          = 0.02 - 4 x2 + 12 x1²
            //
            // d²f/dx1dx2 = d/dx2[∂f/∂x1] = -4 x1
            // d²f/dx2²   = d/dx2[2 (x2 - x1^2)] = 2
            //
            // All derivatives w.r.t x3 are zero.
            RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);

            double h11 = 0.02 - 4.0 * x2 + 12.0 * x1 * x1;
            double h12 = -4.0 * x1;
            double h22 = 2.0;

            H.setEntry(0, 0, h11);
            H.setEntry(0, 1, h12);
            H.setEntry(1, 0, h12);
            H.setEntry(1, 1, h22);
            // x3-related second derivatives are zero
            H.setEntry(2, 2, 0.0);

            return H;
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints (G >= 0, Fortran convention)
    // -------------------------------------------------------------------------
    private static class HS252Ineq extends InequalityConstraint {

        HS252Ineq() {
            // RHS = 0 for all inequalities (G(x) >= 0)
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x3 = x.getEntry(2);

            // Fortran:
            // G(1) = x1 + x3^2 + 1  (>= 0)
            double g1 = x1 + x3 * x3 + 1.0;

            // SQP optimizer convention: constraints are returned as G(x) >= 0
            return new ArrayRealVector(new double[]{g1}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            double x3 = x.getEntry(2);

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // g1 = x1 + x3^2 + 1 → grad = [1, 0, 2*x3]
            J.setEntry(0, 0, 1.0);
            J.setEntry(0, 1, 0.0);
            J.setEntry(0, 2, 2.0 * x3);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS252_optimization() {

        // Initial guess (MODE=1): X(1) = -1, X(2) = 2, X(3) = 2
        double[] x0 = new double[]{-1.0, 2.0, 2.0};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // Bounds:
        // From TP252:
        //   LXU(1) = .TRUE., XU(1) = -1.0
        //   all other LXL/LXU are .FALSE.
        //
        // So:
        //   -∞ < x1 ≤ -1
        //   x2, x3 free
        SimpleBounds bounds = new SimpleBounds(
                new double[]{
                        Double.NEGATIVE_INFINITY, // x1 lower
                        Double.NEGATIVE_INFINITY, // x2 lower
                        Double.NEGATIVE_INFINITY  // x3 lower
                },
                new double[]{
                        -1.0,                     // x1 upper
                        Double.POSITIVE_INFINITY, // x2 upper
                        Double.POSITIVE_INFINITY  // x3 upper
                }
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS252Obj()),
                null,            // no equalities
                new HS252Ineq(), // 1 inequality
                bounds
        );

        double f = sol.getValue();

        // Reference optimum from Fortran: FEX = 0.4D-1 = 0.04
        final double fExpected = 0.04;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
