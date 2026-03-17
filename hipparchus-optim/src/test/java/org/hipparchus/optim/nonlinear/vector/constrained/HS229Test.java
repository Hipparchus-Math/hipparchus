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
 * HS229 (TP229)
 *
 * N    = 2
 * NILI = 0
 * NINL = 0
 * NELI = 0
 * NENL = 0
 *
 * Objective:
 *   f(x) = 100 (x2 - x1^2)^2 + (1 - x1)^2
 *
 * Bounds:
 *   -2 <= x1 <= 2
 *   -2 <= x2 <= 2
 *
 * Reference solution (MODE=1):
 *   x*  = (1, 1)
 *   f*  = 0
 */
public class HS229Test {

    private static final int DIM      = 2;
    private static final int NUM_INEQ = 0;
    private static final int NUM_EQ   = 0;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS229Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double t  = x2 - x1 * x1;
            return 100.0 * t * t + (1.0 - x1) * (1.0 - x1);
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double t  = x2 - x1 * x1;

            double df1 = -400.0 * x1 * t - 2.0 * (1.0 - x1);
            double df2 =  200.0 * t;

            return new ArrayRealVector(new double[]{df1, df2});
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Simple zero Hessian; SQP/BFGS will build curvature.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // No inequality constraints for TP229 → no InequalityConstraint wrapper

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS229_optimization() {

        // Initial guess (MODE=1): X(1) = -1.2, X(2) = 1
        double[] x0 = new double[]{-1.2, 1.0};

        // Bounds: -2 <= x_i <= 2
        double[] lower = new double[]{-2.0, -2.0};
        double[] upper = new double[]{ 2.0,  2.0};
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS229Obj()),
                null,   // no equalities
                null,   // no inequalities
                bounds  // box constraints
        );

        double f = sol.getValue();

        final double fExpected = 0.0;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
