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
 * HS228 (TP228)
 *
 * N    = 2
 * NILI = 1  (one linear inequality)
 * NINL = 1  (one nonlinear inequality)
 *
 * Objective (MODE=2):
 *   f(x) = x1^2 + x2
 *
 * Fortran constraints (G(i) ≥ 0):
 *   G1(x) = -x1 - x2 + 1
 *   G2(x) = -(x1^2 + x2^2) + 9
 *
 * In this wrapper we use g(x) <= 0, so:
 *   g1(x) = x1 + x2 - 1       <= 0
 *   g2(x) = x1^2 + x2^2 - 9   <= 0
 *
 * Reference solution (MODE=1):
 *   x*  = (0, -3)
 *   f*  = -3
 */
public class HS228Test {

    private static final int DIM      = 2;
    private static final int NUM_INEQ = 2;
    private static final int NUM_EQ   = 0;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS228Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            return x1 * x1 + x2;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            return new ArrayRealVector(new double[]{
                    2.0 * x1, // df/dx1
                    1.0       // df/dx2
            }, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);
            H.setEntry(0, 0, 2.0); // d²f/dx1²
            H.setEntry(1, 1, 0.0); // d²f/dx2²
            return H;
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints (g <= 0)
    // -------------------------------------------------------------------------
    private static class HS228Ineq extends InequalityConstraint {

        HS228Ineq() {
            // RHS = 0 for all inequalities
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

            // Fortran G1 = -x1 - x2 + 1  (>= 0)
            // wrapper g1 = x1 + x2 - 1   (<= 0)
            double g1 = x1 + x2 - 1.0;

            // Fortran G2 = -(x1^2 + x2^2) + 9  (>= 0)
            // wrapper g2 = x1^2 + x2^2 - 9      (<= 0)
            double g2 = x1 * x1 + x2 * x2 - 9.0;

            return new ArrayRealVector(new double[]{-g1, -g2}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // g1 = x1 + x2 - 1 → grad = [1, 1]
            J.setEntry(0, 0, -1.0);
            J.setEntry(0, 1, -1.0);

            // g2 = x1^2 + x2^2 - 9 → grad = [2*x1, 2*x2]
            J.setEntry(1, 0, -2.0 * x1);
            J.setEntry(1, 1, -2.0 * x2);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS228_optimization() {

        // Initial guess (MODE=1): X(1)=0, X(2)=0
        double[] x0 = new double[]{0.0, 0.0};

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS228Obj()),
                null,              // no equalities
                new HS228Ineq(),   // 2 inequalities
                null               // no bounds
        );

        double f = sol.getValue();

        // Reference optimum from Fortran: FEX = -3.0
        final double fExpected = -3.0;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
