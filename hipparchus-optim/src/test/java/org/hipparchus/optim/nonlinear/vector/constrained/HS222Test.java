/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
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
 * HS222 (TP222)
 *
 * N    = 2
 * NILI = 0
 * NINL = 1
 * NELI = 0
 * NENL = 0
 *
 * Inequality:
 *   G1(x) = 0.125 - (x1 - 1)^3 - x2 <= 0
 *
 * Objective:
 *   f(x) = -x1
 *
 * Solution (LEX = .TRUE.):
 *   x* = (1.5, 0)
 *   f* = -1.5
 */
public class HS222Test {

    private static final int DIM      = 2;
    private static final int NUM_INEQ = 1;

    // -------------------------------------------------------------------------
    // Objective f = -x1
    // -------------------------------------------------------------------------
    private static class HS222Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            return -x.getEntry(0);
        }

        @Override
        public RealVector gradient(RealVector x) {
            return new ArrayRealVector(new double[]{-1.0, 0.0}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Zero Hessian – BFGS will build curvature
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality: G1 = 0.125 - (x1 - 1)^3 - x2 <= 0
    // -------------------------------------------------------------------------
    private static class HS222Ineq extends InequalityConstraint {

        HS222Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQ])); // RHS = 0
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            double g = 0.125 - FastMath.pow(x1 - 1.0, 3.0) - x2;
            return new ArrayRealVector(new double[]{g}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            final double x1 = x.getEntry(0);

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // dG/dx1 = -3 * (x1 - 1)^2
            // dG/dx2 = -1
            J.setEntry(0, 0, -3.0 * FastMath.pow(x1 - 1.0, 2.0));
            J.setEntry(0, 1, -1.0);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS222_optimization() {

        // Initial guess from Fortran:
        // X(1) = 1.3, X(2) = 0.2
        double[] x0 = new double[]{1.3, 0.2};

        // Bounds: x1 >= 0, x2 >= 0, no upper bounds
        double[] lower = new double[]{0.0, 0.0};
        double[] upper = new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS222Obj()),
                null,              // no equalities
                new HS222Ineq(),   // 1 inequality
                bounds
        );

        double f = sol.getValue();

        final double fExpected = -1.5;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
