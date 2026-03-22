/*
 * HS270 (TP270)
 *
 * N    = 5
 * NINL = 1  (1 nonlinear inequality)
 *
 * Objective (polynomial in x1..x5):
 *   as in Fortran TP270
 *
 * Inequality:
 *   G1 = 34 - x1^2 - x2^2 - x3^2 - x4^2 - x5^2 ≥ 0
 *
 * Bounds:
 *   x1 ≥ 1, x2 ≥ 2, x3 ≥ 3, x4 ≥ 4, x5 free
 *
 * Initial guess:
 *   x0 = (1.1, 2.1, 3.1, 4.1, -1.0)
 *
 * Reference solution (Fortran):
 *   x* = (1, 2, 3, 4, 2)
 *   f* = -1
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

public class HS270Test {

    private static final int DIM      = 5;
    private static final int NUM_INEQ = 1;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS270Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);
            double x5 = x.getEntry(4);

            double fx =
                    x1 * x2 * x3 * x4
                    - 3.0 * x1 * x2 * x4
                    - 4.0 * x1 * x2 * x3
                    + 12.0 * x1 * x2
                    - x2 * x3 * x4
                    + 3.0 * x2 * x4
                    + 4.0 * x2 * x3
                    - 12.0 * x2
                    - 2.0 * x1 * x3 * x4
                    + 6.0 * x1 * x4
                    + 8.0 * x1 * x3
                    - 24.0 * x1
                    + 2.0 * x3 * x4
                    - 6.0 * x4
                    - 8.0 * x3
                    + 24.0
                    + 1.5 * FastMath.pow(x5, 4)
                    - 5.75 * FastMath.pow(x5, 3)
                    + 5.25 * x5 * x5;

            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);
            double x5 = x.getEntry(4);

            // From Fortran GF(1..5)
            double g1 =
                    x2 * x3 * x4
                    - 3.0 * x2 * x4
                    - 4.0 * x2 * x3
                    + 12.0 * x2
                    - 2.0 * x3 * x4
                    + 6.0 * x4
                    + 8.0 * x3
                    - 24.0;

            double g2 =
                    x1 * x3 * x4
                    - 3.0 * x1 * x4
                    - 4.0 * x1 * x3
                    + 12.0 * x1
                    - x3 * x4
                    + 3.0 * x4
                    + 4.0 * x3
                    - 12.0;

            double g3 =
                    x1 * x2 * x4
                    - 4.0 * x1 * x2
                    - x2 * x4
                    + 4.0 * x2
                    - 2.0 * x1 * x4
                    + 8.0 * x1
                    + 2.0 * x4
                    - 8.0;

            double g4 =
                    x1 * x2 * x3
                    - 3.0 * x1 * x2
                    - x2 * x3
                    + 3.0 * x2
                    - 2.0 * x1 * x3
                    + 6.0 * x1
                    + 2.0 * x3
                    - 6.0;

            double g5 =
                    10.5 * x5
                    - 17.25 * x5 * x5
                    + 6.0 * FastMath.pow(x5, 3);

            return new ArrayRealVector(new double[]{g1, g2, g3, g4, g5}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Nonlinear, but we let BFGS approximate Hessian.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Nonlinear inequality G(x) ≥ 0:
    //   G1 = 34 - sum_i x_i^2
    // -------------------------------------------------------------------------
    private static class HS270Ineq extends InequalityConstraint {

        HS270Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {
            double sumSq = 0.0;
            for (int i = 0; i < DIM; i++) {
                double xi = x.getEntry(i);
                sumSq += xi * xi;
            }
            double g1 = 34.0 - sumSq;
            return new ArrayRealVector(new double[]{g1}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);
            // ∂G1/∂xi = -2*xi
            for (int i = 0; i < DIM; i++) {
                J.setEntry(0, i, -2.0 * x.getEntry(i));
            }
            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS270() {

        double[] x0 = new double[]{1.1, 2.1, 3.1, 4.1, -1.0};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // Bounds: x1 ≥1, x2 ≥2, x3 ≥3, x4 ≥4, x5 free
        SimpleBounds bounds = new SimpleBounds(
                new double[]{1.0, 2.0, 3.0, 4.0, -Double.POSITIVE_INFINITY},
                new double[]{
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY
                }
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS270Obj()),
                null,              // no equalities
                new HS270Ineq(),   // 1 inequality
                bounds
        );

        double f = sol.getValue();

        double fExpected = -1.0;
        double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
