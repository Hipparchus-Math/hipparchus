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
 * HS234 (TP234)
 *
 * N    = 2
 * NILI = 0
 * NINL = 1 (one nonlinear inequality)
 * NELI = 0
 * NENL = 0
 *
 * Objective (MODE=2):
 *   f(x) = (x2 - x1)^4 - (1 - x1)
 *
 * Fortran constraint (G(x) ≥ 0):
 *   G1(x) = -x1^2 - x2^2 + 1
 *
 * Bounds (MODE=1):
 *   0.2 ≤ x1 ≤ 2
 *   0.2 ≤ x2 ≤ 2
 *
 * Reference solution (MODE=1):
 *   x*  = (0.2, 0.2)
 *   f*  = -0.8
 */
public class HS234Test {

    private static final int DIM      = 2;
    private static final int NUM_INEQ = 1;
    private static final int NUM_EQ   = 0;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS234Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double t  = x2 - x1;
            return FastMath.pow(t, 4) - (1.0 - x1);
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double t  = x2 - x1;
            double t3 = t * t * t;

            // From Fortran:
            // GF(1) = -4*(x2 - x1)^3 + 1
            // GF(2) =  4*(x2 - x1)^3
            double df1 = -4.0 * t3 + 1.0;
            double df2 =  4.0 * t3;

            return new ArrayRealVector(new double[]{df1, df2}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Let SQP/BFGS update the Hessian
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraint G(x) ≥ 0
    // -------------------------------------------------------------------------
    private static class HS234Ineq extends InequalityConstraint {

        HS234Ineq() {
            // RHS = 0 for the inequality
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

            // Fortran: G1 = -x1^2 - x2^2 + 1  (≥ 0)
            double G1 = -x1 * x1 - x2 * x2 + 1.0;

            return new ArrayRealVector(new double[]{G1}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // dG1/dx1 = -2*x1
            // dG1/dx2 = -2*x2
            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);
            J.setEntry(0, 0, -2.0 * x1);
            J.setEntry(0, 1, -2.0 * x2);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS234_optimization() {

        // Initial guess (MODE=1): X(1)=1, X(2)=1
        double[] x0 = new double[]{1.0, 1.0};

        // Bounds: 0.2 ≤ xi ≤ 2.0
        double[] lower = new double[]{0.2, 0.2};
        double[] upper = new double[]{2.0, 2.0};
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS234Obj()),
                null,              // no equalities
                new HS234Ineq(),   // 1 inequality
                bounds
        );

        double f = sol.getValue();

        final double fExpected = -0.8;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
