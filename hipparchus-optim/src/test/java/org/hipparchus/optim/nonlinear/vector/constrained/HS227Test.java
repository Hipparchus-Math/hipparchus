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
 * HS227 (TP227)
 *
 * N    = 2
 * NINL = 2 (two nonlinear inequality constraints)
 *
 * Objective:
 *   f(x) = (x1 - 2)^2 + (x2 - 1)^2
 *
 * Constraints:
 *   G1(x) = -x1^2 + x2
 *   G2(x) = x1 - x2^2
 *
 * Reference optimum:
 *   x* = (1, 1)
 *   f* = 1.0
 */
public class HS227Test {

    private static final int DIM = 2;
    private static final int NUM_INEQ = 2;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS227Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            return (x1 - 2.0) * (x1 - 2.0) + (x2 - 1.0) * (x2 - 1.0);
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            return new ArrayRealVector(new double[]{
                2.0 * (x1 - 2.0),
                2.0 * (x2 - 1.0)
            }, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Hessian is constant:
            // d²f/dx1² = 2, d²f/dx2² = 2, cross terms = 0
            RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);
            H.setEntry(0, 0, 2.0);
            H.setEntry(1, 1, 2.0);
            return H;
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints
    // -------------------------------------------------------------------------
    private static class HS227Ineq extends InequalityConstraint {

        HS227Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQ])); // RHS=0
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            double g1 = -x1 * x1 + x2;      // From MODE=4
            double g2 =  x1 - x2 * x2;

            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // From MODE=5:
            // GG(1,1) = -2*x1
            // GG(2,2) = -2*x2

            // G1 gradient: [-2*x1, 1]
            J.setEntry(0, 0, -2.0 * x1);
            J.setEntry(0, 1, 1.0);

            // G2 gradient: [1, -2*x2]
            J.setEntry(1, 0, 1.0);
            J.setEntry(1, 1, -2.0 * x2);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS227_optimization() {

        double[] x0 = new double[]{0.5, 0.5}; // From MODE=1

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS227Obj()),
                null,
                new HS227Ineq(),
                null
        );

        double f = sol.getValue();

        // Reference optimum
        final double fExpected = 1.0;
        final double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
