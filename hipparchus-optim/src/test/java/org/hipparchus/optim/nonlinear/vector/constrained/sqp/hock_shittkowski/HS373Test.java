/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  
 */

package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * HS373 (TP373) – 9 variables, 6 nonlinear *equality* constraints.
 *
 * From TP373:
 *
 *   f(x) = sum_{i=4..9} x(i)^2
 *
 * Constraints (EQUALITY):
 *   For i = 1..6, let k_i = 2*i - 7 = [-5,-3,-1,1,3,5]
 *                     c_i = [127,151,379,421,460,426]
 *
 *   Gi(x) = x1 + x2 * exp(k_i x3) + x_{i+3} - c_i = 0
 *
 * Bounds:
 *   x1,x2,x4,x5,x6,x7,x8,x9 : free
 *   x3: -1 <= x3 <= 0
 *
 * Reference solution:
 *   FEX = 13390.093
 *   XEX(1..9) = provided below
 *   LEX = TRUE → strict matching
 */
public class HS373Test {

    private static final int DIM = 9;
    private static final int NUM_EQ = 6;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS373Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        /** f(x) = sum_{i=4..9} x_i^2 */
        @Override
        public double value(RealVector x) {
            double fx = 0.0;
            for (int i = 3; i < DIM; i++) {
                double v = x.getEntry(i);
                fx += v * v;
            }
            return fx;
        }

        /** Gradient: df/dx_i = 0 for i<=3, = 2*x_i for i>=4 */
        @Override
        public RealVector gradient(RealVector x) {
            double[] g = new double[DIM];

            // i = 1..3 → 0
            g[0] = 0.0;
            g[1] = 0.0;
            g[2] = 0.0;

            // i = 4..9 → 2*x_i
            for (int i = 3; i < DIM; i++) {
                g[i] = 2.0 * x.getEntry(i);
            }
            return new ArrayRealVector(g, false);
        }

        /** Hessian: diagonal 2's for x4..x9 */
        @Override
        public RealMatrix hessian(RealVector x) {
            RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);
            for (int i = 3; i < DIM; i++) {
                H.setEntry(i, i, 2.0);
            }
            return H;
        }
    }

    // -------------------------------------------------------------------------
    // Equality constraints (G1..G6 = 0)
    // -------------------------------------------------------------------------
    private static class HS373Eq extends EqualityConstraint {

        HS373Eq() {
            super(new ArrayRealVector(new double[NUM_EQ])); // RHS=0
        }

        @Override
        public int dim() { return DIM; }

        /**
         * Gi(x) = x1 + x2 * exp(k_i * x3) + x_{i+3} - c_i
         *
         * k_i = [-5,-3,-1,1,3,5]
         * c_i = [127,151,379,421,460,426]
         */
        @Override
        public RealVector value(RealVector x) {

            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);

            final int[] k = {-5, -3, -1, 1, 3, 5};
            final double[] c = {127, 151, 379, 421, 460, 426};

            double[] g = new double[NUM_EQ];

            for (int i = 0; i < 6; i++) {
                double expTerm = FastMath.exp(k[i] * x3);
                double xip3 = x.getEntry(3 + i); // x4..x9
                g[i] = x1 + x2 * expTerm + xip3 - c[i];
            }
            return new ArrayRealVector(g, false);
        }

        /** Jacobian of equality constraints */
        @Override
        public RealMatrix jacobian(RealVector x) {

            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);

            final int[] k = {-5, -3, -1, 1, 3, 5};

            RealMatrix J = new Array2DRowRealMatrix(NUM_EQ, DIM);

            for (int i = 0; i < 6; i++) {
                double expTerm = FastMath.exp(k[i] * x3);
                double d_dx3 = x2 * k[i] * expTerm;
                int colXiP3 = 3 + i;

                // dG_i/dx1
                J.setEntry(i, 0, 1.0);

                // dG_i/dx2
                J.setEntry(i, 1, expTerm);

                // dG_i/dx3
                J.setEntry(i, 2, d_dx3);

                // dG_i/dx_{i+3}
                J.setEntry(i, colXiP3, 1.0);
            }

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // TEST
    // -------------------------------------------------------------------------
    @Test
    public void testHS373_optimization() {

        // Initial point from MODE=1:
        double[] x0 = new double[DIM];
        x0[0] = 300.0;
        x0[1] = -100.0;
        x0[2] = -0.1997;
        x0[3] = -127.0;
        x0[4] = -151.0;
        x0[5] = 379.0;
        x0[6] = 421.0;
        x0[7] = 460.0;
        x0[8] = 426.0;

        // Bounds:
        // x3 ∈ [-1, 0], others free
        double[] lower = new double[DIM];
        double[] upper = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            lower[i] = Double.NEGATIVE_INFINITY;
            upper[i] = Double.POSITIVE_INFINITY;
        }
        lower[2] = -1.0;
        upper[2] = 0.0;

        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS373Obj()),
                new HS373Eq(),   // 6 equality constraints
                bounds
        );

        final double f = sol.getValue();

        // Reference (LEX = TRUE → exact match):
        final double fExpected = 0.13390093e5;
        final double tol = 1e-4 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol,
                "HS373: objective mismatch with reference FEX");
    }
}
