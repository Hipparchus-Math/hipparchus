/*
 * HS257 (TP257)
 *
 * N     = 4
 * NILI  = 0
 * NINL  = 0
 * NELI  = 0
 * NENL  = 0
 *
 * Fortran TP257:
 *
 *   MODE=2:
 *     FX = 100*(x1^2 - x2)^2 + (x1 - 1)^2
 *        + 90*(x3^2 - x4)^2 + (x3 - 1)^2
 *        + 10.1*((x2 - 1)^2 + (x4 - 1)^2)
 *        + 19.8*(x1 - 1)*(x4 - 1)
 *
 *   MODE=3:
 *     GF(1)=  4.0E2*(x1^3 - x1*x2) + 2*x1 + 19.8*x4 - 21.8
 *     GF(2)= -2.0E2*x1^2 + 220.2*x2 - 20.2
 *     GF(3)=  3.6E2*(x3^3 - x3*x4) + 2*x3 - 2
 *     GF(4)= -1.8E2*x3^2 + 200.2*x4 + 19.8*x1 - 40
 *
 * Bounds (MODE=1):
 *   x_i ≥ 0, no explicit upper bound in Fortran.
 *
 * Initial guess (MODE=1):
 *   x0 = (0, 0, 0, 0)
 *
 * Reference solution:
 *   x*  = (1, 1, 1, 1)
 *   f*  = 0
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

public class HS257Test {

    private static final int DIM = 4;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS257Obj extends TwiceDifferentiableFunction {

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

            double term1 = 100.0 * FastMath.pow(x1 * x1 - x2, 2);     // 100*(x1^2 - x2)^2
            double term2 = FastMath.pow(x1 - 1.0, 2);                  // (x1 - 1)^2
            double term3 = 90.0 * FastMath.pow(x3 * x3 - x4, 2);      // 90*(x3^2 - x4)^2
            double term4 = FastMath.pow(x3 - 1.0, 2);                  // (x3 - 1)^2

            double t2m1 = x2 - 1.0;
            double t4m1 = x4 - 1.0;

            double term5 = 10.1 * (t2m1 * t2m1 + t4m1 * t4m1);         // 10.1*((x2-1)^2 + (x4-1)^2)
            double term6 = 19.8 * (x1 - 1.0) * (x4 - 1.0);             // 19.8*(x1-1)*(x4-1)

            return term1 + term2 + term3 + term4 + term5 + term6;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            // From Fortran GF(...)
            double g1 = 4.0e2 * (FastMath.pow(x1, 3) - x1 * x2)
                      + 2.0 * x1
                      + 19.8 * x4
                      - 21.8;

            double g2 = -2.0e2 * x1 * x1
                      + 220.2 * x2
                      - 20.2;

            double g3 = 3.6e2 * (FastMath.pow(x3, 3) - x3 * x4)
                      + 2.0 * x3
                      - 2.0;

            double g4 = -1.8e2 * x3 * x3
                      + 200.2 * x4
                      + 19.8 * x1
                      - 40.0;

            return new ArrayRealVector(new double[]{g1, g2, g3, g4}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Let BFGS build an approximation
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS257() {

        // Initial guess: X(i) = 0
        double[] x0 = new double[]{0.0, 0.0, 0.0, 0.0};

        // Bounds: x_i ≥ 0, no upper bounds in the original problem
        SimpleBounds bounds = new SimpleBounds(
                new double[]{0.0, 0.0, 0.0, 0.0},
                new double[]{
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY
                }
        );

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS257Obj()),
                null,   // no equalities
                null,   // no inequalities
                bounds
        );

        double f = sol.getValue();

        double fExpected = 0.0;
        double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
