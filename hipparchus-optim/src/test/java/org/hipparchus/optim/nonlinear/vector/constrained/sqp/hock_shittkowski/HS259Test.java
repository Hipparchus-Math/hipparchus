/*
 * HS259 (TP259)
 *
 * N     = 4
 * NILI  = 0
 * NINL  = 0
 * NELI  = 0
 * NENL  = 0
 *
 * Fortran TP259:
 *
 * MODE=2:
 *   f(x) = 100 (x2 - x1^2)^2 + (1 - x1)^2
 *        + 90 (x4 - x3^2)^2 + (1 - x3)^3
 *        + 10.1 (x2 - 1)^2 + (x4 - 1)^2
 *        + 19.8 (x2 - 1)(x4 - 1)
 *
 * MODE=3:
 *   GF(1) = 4e2 (x1^3 - x1 x2) + 2 x1 - 2
 *   GF(2) = -2e2 x1^2 + 220.2 x2 + 19.8 x4 - 40
 *   GF(3) = 3.6e2 (x3^3 - x3 x4) - 3 (1 - x3)^2
 *   GF(4) = -1.8e2 x3^2 + 182 x4 + 19.8 x2 - 21.8
 *
 * Bounds (MODE=1):
 *   x4 ≤ 1, others free.
 *
 * Initial guess (MODE=1):
 *   x0 = (0, 0, 0, 0)
 *
 * Reference solution (one of them, NEX=2 in Fortran):
 *   x* ≈ (1.4358451, 2.0631635, 0.069002268, -0.099963939)
 *   f* ≈ -8.5446210
 */

package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS259Test {

    private static final int DIM = 4;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS259Obj extends TwiceDifferentiableFunction {

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

            double t1   = x2 - x1 * x1;
            double t3   = x4 - x3 * x3;
            double t2m1 = x2 - 1.0;
            double t4m1 = x4 - 1.0;

            double term1 = 100.0 * t1 * t1;
            double term2 = FastMath.pow(1.0 - x1, 2);
            double term3 = 90.0 * t3 * t3;
            double term4 = FastMath.pow(1.0 - x3, 3); // (1 - x3)^3
            double term5 = 10.1 * t2m1 * t2m1;
            double term6 = FastMath.pow(t4m1, 2);
            double term7 = 19.8 * t2m1 * t4m1;

            return term1 + term2 + term3 + term4 + term5 + term6 + term7;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            double g1 = 4.0e2 * (FastMath.pow(x1, 3) - x1 * x2)
                      + 2.0 * x1
                      - 2.0;

            double g2 = -2.0e2 * x1 * x1
                      + 220.2 * x2
                      + 19.8 * x4
                      - 40.0;

            double g3 = 3.6e2 * (FastMath.pow(x3, 3) - x3 * x4)
                      - 3.0 * FastMath.pow(1.0 - x3, 2);

            double g4 = -1.8e2 * x3 * x3
                      + 182.0 * x4
                      + 19.8 * x2
                      - 21.8;

            return new ArrayRealVector(new double[]{g1, g2, g3, g4}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Let BFGS approximate the Hessian
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS259() {

        // Initial guess from Fortran:
        double[] x0 = new double[]{0.0, 0.0, 0.0, 0.0};

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        // Only x4 has an upper bound: x4 ≤ 1
        SimpleBounds bounds = new SimpleBounds(
                new double[]{
                        Double.NEGATIVE_INFINITY,
                        Double.NEGATIVE_INFINITY,
                        Double.NEGATIVE_INFINITY,
                        Double.NEGATIVE_INFINITY
                },
                new double[]{
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY,
                        1.0
                }
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS259Obj()),
                null,   // no equalities
                null,   // no inequalities
                bounds  // bound on x4
        );

        double f = sol.getValue();

        double fExpected = -8.5446210;
        double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
