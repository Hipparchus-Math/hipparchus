/*
 * HS264 (TP264)
 *
 * N    = 4
 * NILI = 0
 * NINL = 3  (3 nonlinear inequalities)
 * NELI = 0
 * NENL = 0
 *
 * Objective:
 *   f(x) = x1^2 + x2^2 + 2 x3^2 + x4^2
 *          - 5 x1 - 5 x2 - 21 x3 + 7 x4
 *
 * Inequalities (G >= 0):
 *
 *   G1 = 8  - x1^2 - x2^2 - x3^2 - x4^2 - x1 + x2 - x3 + x4
 *   G2 = 9  - x1^2 - 2 x2^2 - x3^2 - 2 x4^2 + x1 + x4
 *   G3 = 5  - 2 x1^2 - x2^2 - x3^2 - 2 x1 + x2 + x4
 *
 * Initial guess: x = (0, 0, 0, 0)
 *
 * Reference solution (Fortran):
 *   x*  = (0, 10, 20, -10)
 *   f*  = -44
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

public class HS264Test {

    private static final int DIM       = 4;
    private static final int NUM_INEQ  = 3;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS264Obj extends TwiceDifferentiableFunction {

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

            return (x1 * x1 + x2 * x2 + 2.0 * x3 * x3 + x4 * x4
                    - 5.0 * x1 - 5.0 * x2 - 21.0 * x3 + 7.0 * x4);
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            double g1 = 2.0 * x1 - 5.0;
            double g2 = 2.0 * x2 - 5.0;
            double g3 = 4.0 * x3 - 21.0;
            double g4 = 2.0 * x4 + 7.0;

            return new ArrayRealVector(new double[]{g1, g2, g3, g4}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);

            H.setEntry(0, 0, 2.0); // d²/dx1²
            H.setEntry(1, 1, 2.0); // d²/dx2²
            H.setEntry(2, 2, 4.0); // d²/dx3² (2 * 2)
            H.setEntry(3, 3, 2.0); // d²/dx4²

            return H;
        }
    }

    // -------------------------------------------------------------------------
    // Nonlinear inequalities (G >= 0)
    // -------------------------------------------------------------------------
    private static class HS264Ineq extends InequalityConstraint {

        HS264Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQ])); // RHS = 0
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            double g1 = 8.0 - x1 * x1 - x2 * x2 - x3 * x3 - x4 * x4
                        - x1 + x2 - x3 + x4;

            double g2 = 9.0 - x1 * x1 - 2.0 * x2 * x2 - x3 * x3 - 2.0 * x4 * x4
                        + x1 + x4;

            double g3 = 5.0 - 2.0 * x1 * x1 - x2 * x2 - x3 * x3
                        - 2.0 * x1 + x2 + x4;

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // g1 = 8 - x1^2 - x2^2 - x3^2 - x4^2 - x1 + x2 - x3 + x4
            J.setEntry(0, 0, -2.0 * x1 - 1.0); // ∂g1/∂x1
            J.setEntry(0, 1, -2.0 * x2 + 1.0); // ∂g1/∂x2
            J.setEntry(0, 2, -2.0 * x3 - 1.0); // ∂g1/∂x3
            J.setEntry(0, 3, -2.0 * x4 + 1.0); // ∂g1/∂x4

            // g2 = 9 - x1^2 - 2 x2^2 - x3^2 - 2 x4^2 + x1 + x4
            J.setEntry(1, 0, -2.0 * x1 + 1.0);    // ∂g2/∂x1
            J.setEntry(1, 1, -4.0 * x2);          // ∂g2/∂x2
            J.setEntry(1, 2, -2.0 * x3);          // ∂g2/∂x3
            J.setEntry(1, 3, -4.0 * x4 + 1.0);    // ∂g2/∂x4

            // g3 = 5 - 2 x1^2 - x2^2 - x3^2 - 2 x1 + x2 + x4
            J.setEntry(2, 0, -4.0 * x1 - 2.0);    // ∂g3/∂x1
            J.setEntry(2, 1, -2.0 * x2 + 1.0);    // ∂g3/∂x2
            J.setEntry(2, 2, -2.0 * x3);          // ∂g3/∂x3
            J.setEntry(2, 3,  1.0);               // ∂g3/∂x4

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS264() {

        // Initial guess: x = (0, 0, 0, 0)
        double[] x0 = new double[]{0.0, 0.0, 0.0, 0.0};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS264Obj()),
                null,              // no equalities
                new HS264Ineq(),   // 3 inequalities
                null               // no bounds
        );

        double f = sol.getValue();

        double fExpected = -44.0;
        double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
