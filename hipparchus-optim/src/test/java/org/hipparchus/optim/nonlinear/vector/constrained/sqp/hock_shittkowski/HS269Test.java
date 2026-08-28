/*
 * HS269 (TP269)
 *
 * N    = 5
 * NELI = 3  (3 linear equalities)
 *
 * F1 = x1 - x2
 * F2 = x2 + x3 - 2
 * F3 = x4 - 1
 * F4 = x5 - 1
 *
 * Objective:
 *   f(x) = sum_{i=1..4} F_i(x)^2
 *
 * Equalities:
 *   G1 = x1 + 3 x2              = 0
 *   G2 = x3 + x4 - 2 x5         = 0
 *   G3 = x2 - x5                = 0
 *
 * Initial guess:
 *   x0 = (2, 2, 2, 2, 2)
 *
 * Reference solution (Fortran):
 *   x* = ( -33/43,  11/43, 27/43, -5/43, 11/43 )
 *   f* = 176/43 ≈ 4.0930232558
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

public class HS269Test {

    private static final int DIM      = 5;
    private static final int NUM_EQ   = 3;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS269Obj extends TwiceDifferentiableFunction {

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

            double f1 = x1 - x2;
            double f2 = x2 + x3 - 2.0;
            double f3 = x4 - 1.0;
            double f4 = x5 - 1.0;

            return f1 * f1 + f2 * f2 + f3 * f3 + f4 * f4;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);
            double x5 = x.getEntry(4);

            // From Fortran:
            // GF(1)=2*(x1-x2)
            // GF(2)=2*(2*x2 + x3 - x1 - 2)
            // GF(3)=2*(x2 + x3 - 2)
            // GF(4)=2*(x4 - 1)
            // GF(5)=2*(x5 - 1)
            double g1 = 2.0 * (x1 - x2);
            double g2 = 2.0 * (2.0 * x2 + x3 - x1 - 2.0);
            double g3 = 2.0 * (x2 + x3 - 2.0);
            double g4 = 2.0 * (x4 - 1.0);
            double g5 = 2.0 * (x5 - 1.0);

            return new ArrayRealVector(new double[]{g1, g2, g3, g4, g5}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Quadratic objective; Hessian is constant, but BFGS can estimate it.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Equality constraints G(x) = 0 (NELI = 3)
    //
    // G1 = x1 + 3 x2
    // G2 = x3 + x4 - 2 x5
    // G3 = x2 - x5
    // -------------------------------------------------------------------------
    private static class HS269Eq extends EqualityConstraint {

        HS269Eq() {
            super(new ArrayRealVector(new double[NUM_EQ]));
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
            double x5 = x.getEntry(4);

            double g1 = x1 + 3.0 * x2;
            double g2 = x3 + x4 - 2.0 * x5;
            double g3 = x2 - x5;

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            RealMatrix J = new Array2DRowRealMatrix(NUM_EQ, DIM);

            // G1 = x1 + 3 x2
            J.setEntry(0, 0, 1.0);
            J.setEntry(0, 1, 3.0);

            // G2 = x3 + x4 - 2 x5
            J.setEntry(1, 2, 1.0);
            J.setEntry(1, 3, 1.0);
            J.setEntry(1, 4, -2.0);

            // G3 = x2 - x5
            J.setEntry(2, 1, 1.0);
            J.setEntry(2, 4, -1.0);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS269() {

        double[] x0 = new double[]{2.0, 2.0, 2.0, 2.0, 2.0};

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS269Obj()),
                new HS269Eq(),   // 3 equalities
                null,           // no inequalities
                null            // no bounds
        );

        double f = sol.getValue();

        double fExpected = 176.0 / 43.0;  // 0.176D+3 / 0.43D+2
        double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
