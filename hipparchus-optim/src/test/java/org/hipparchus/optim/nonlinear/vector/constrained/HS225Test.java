/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements…
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
 * HS225 (TP225)
 *
 * N    = 2
 * NILI = 0
 * NINL = 5  (nonlinear inequalities)
 * NELI = 0
 * NENL = 0
 *
 * No bounds: both variables are free.
 *
 * Objective:
 *   f(x) = x1^2 + x2^2
 *
 * Constraints:
 *   G1 = x1 + x2 - 1
 *   G2 = x1^2 + x2^2 - 1
 *   G3 = 9*x1^2 + x2^2 - 9
 *   G4 = x1^2 - x2
 *   G5 = x2^2 - x1
 *
 * Optimum (LEX = TRUE):
 *   x* = (1,1)
 *   f* = 2
 */
public class HS225Test {

    private static final int DIM      = 2;
    private static final int NUM_INEQ = 5;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS225Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() { return DIM; }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            return x1*x1 + x2*x2;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            return new ArrayRealVector(new double[]{
                2.0*x1,
                2.0*x2
            }, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);
            H.setEntry(0,0,2.0);
            H.setEntry(1,1,2.0);
            return H;
        }
    }

    // -------------------------------------------------------------------------
    // Inequalities G1..G5
    // -------------------------------------------------------------------------
    private static class HS225Ineq extends InequalityConstraint {

        HS225Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQ])); // RHS = 0
        }

        @Override
        public int dim() { return DIM; }

        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            double g1 = x1 + x2 - 1.0;
            double g2 = x1*x1 + x2*x2 - 1.0;
            double g3 = 9.0*x1*x1 + x2*x2 - 9.0;
            double g4 = x1*x1 - x2;
            double g5 = x2*x2 - x1;

            return new ArrayRealVector(new double[]{g1,g2,g3,g4,g5}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // G1 = x1 + x2 - 1
            J.setEntry(0,0,1.0);
            J.setEntry(0,1,1.0);

            // G2 = x1^2 + x2^2 -1
            J.setEntry(1,0,2.0*x1);
            J.setEntry(1,1,2.0*x2);

            // G3 = 9*x1^2 + x2^2 - 9
            J.setEntry(2,0,18.0*x1);
            J.setEntry(2,1, 2.0*x2);

            // G4 = x1^2 - x2
            J.setEntry(3,0,2.0*x1);
            J.setEntry(3,1,-1.0);

            // G5 = x2^2 - x1
            J.setEntry(4,0,-1.0);
            J.setEntry(4,1, 2.0*x2);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS225_optimization() {

        double[] x0 = new double[]{3.0, 1.0};

        // No bounds → free variables
        double[] lower = new double[]{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
        double[] upper = new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS225Obj()),
                null,
                new HS225Ineq(),
                bounds
        );

        double f = sol.getValue();

        final double fExpected = 2.0;
        final double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
