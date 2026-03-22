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


class HS249Test {

    private static final int DIM      = 3;
    private static final int NUM_INEQ = 1;
    private static final int NUM_EQ   = 0;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS249Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            return x1 * x1 + x2 * x2 + x3 * x3;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            return new ArrayRealVector(new double[]{
                    2.0 * x1,
                    2.0 * x2,
                    2.0 * x3
            }, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Hessian of x1^2 + x2^2 + x3^2 is 2 * I
            RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);
            for (int i = 0; i < DIM; ++i) {
                H.setEntry(i, i, 2.0);
            }
            return H;
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraint (c(x) ≥ 0)
    // -------------------------------------------------------------------------
    private static class HS249Ineq extends InequalityConstraint {

        HS249Ineq() {
            // Single inequality, RHS = 0
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // Fortran G1 = x1^2 + x2^2 - 1  (≥ 0)
            double c1 = x1 * x1 + x2 * x2 - 1.0;

            // Optimizer convention: c(x) ≥ 0
            return new ArrayRealVector(new double[]{c1}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // c1 = x1^2 + x2^2 - 1 → grad = [2*x1, 2*x2, 0]
            J.setEntry(0, 0, 2.0 * x1);
            J.setEntry(0, 1, 2.0 * x2);
            J.setEntry(0, 2, 0.0);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS249_optimization() {

        // Initial guess from Fortran (MODE=1): X(i)=1 for all i
        double[] x0 = new double[]{1.0, 1.0, 1.0};

        // Bounds: x1 ≥ 1, x2 and x3 free
        double[] lower = new double[]{
                1.0,                     // x1 lower
                Double.NEGATIVE_INFINITY,
                Double.NEGATIVE_INFINITY
        };
        double[] upper = new double[]{
                Double.POSITIVE_INFINITY,
                Double.POSITIVE_INFINITY,
                Double.POSITIVE_INFINITY
        };

        SimpleBounds bounds =
                new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS249Obj()),
                null,            // no equalities
                new HS249Ineq(), // one inequality
                bounds           // bound on x1
        );

        double f = sol.getValue();

        final double fExpected = 1.0; // FEX
        final double tol = 1.0e-4 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}











