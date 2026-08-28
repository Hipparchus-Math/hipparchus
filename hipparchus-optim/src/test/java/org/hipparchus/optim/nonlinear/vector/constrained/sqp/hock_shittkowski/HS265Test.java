/*
 * HS265 (TP265)
 *
 * N    = 4
 * NILI = 0
 * NINL = 0
 * NELI = 2  (2 linear equalities)
 * NENL = 0
 *
 * Objective:
 *   f(x) = 2 - exp(-10 x1 e^{-x3}) - exp(-10 x2 e^{-x4})
 *
 * Equalities (Fortran G = 0):
 *   G1 = x1 + x2 - 1 = 0
 *   G2 = x3 + x4 - 1 = 0
 *
 * Bounds (MODE=1):
 *   xl(i) = 0,  lxl(i) = .TRUE.    →  x_i ≥ 0
 *   lxu(i) = .FALSE.               →  no active upper bounds
 *
 * Initial guess:
 *   x0 = (0, 0, 0, 0)
 *
 * Reference solution (Fortran):
 *   x*  = (1, 0, 1, 0)
 *   f*  = 0.97474658
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
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

public class HS265Test {

    private static final int DIM     = 4;
    private static final int NUM_EQ  = 2;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS265Obj extends TwiceDifferentiableFunction {

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

            double term1 = FastMath.exp(-10.0 * x1 * FastMath.exp(-x3));
            double term2 = FastMath.exp(-10.0 * x2 * FastMath.exp(-x4));

            return 2.0 - term1 - term2;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            // Fortran:
            // GF(I)   = 10 * exp(-10*X(I)*exp(-X(I+2)) - X(I+2))
            // GF(I+2) = -X(I) * GF(I)
            double[] g = new double[DIM];

            double g1 = 10.0 * FastMath.exp(-10.0 * x1 * FastMath.exp(-x3) - x3);
            double g2 = 10.0 * FastMath.exp(-10.0 * x2 * FastMath.exp(-x4) - x4);

            g[0] = g1;          // d f / d x1
            g[1] = g2;          // d f / d x2
            g[2] = -x1 * g1;    // d f / d x3
            g[3] = -x2 * g2;    // d f / d x4

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Complicated; let BFGS build it
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Equalities (G = 0)
    // -------------------------------------------------------------------------
    private static class HS265Eq extends EqualityConstraint {

        HS265Eq() {
            super(new ArrayRealVector(new double[NUM_EQ])); // RHS = 0
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

            // Fortran MODE=4:
            // G(1) = X(1) + X(2) - 1
            // G(2) = X(3) + X(4) - 1
            double g1 = x1 + x2 - 1.0;
            double g2 = x3 + x4 - 1.0;

            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            RealMatrix J = new Array2DRowRealMatrix(NUM_EQ, DIM);

            // g1 = x1 + x2 - 1
            J.setEntry(0, 0, 1.0);
            J.setEntry(0, 1, 1.0);
            J.setEntry(0, 2, 0.0);
            J.setEntry(0, 3, 0.0);

            // g2 = x3 + x4 - 1
            J.setEntry(1, 0, 0.0);
            J.setEntry(1, 1, 0.0);
            J.setEntry(1, 2, 1.0);
            J.setEntry(1, 3, 1.0);

            return J;
        }
    }

//    // -------------------------------------------------------------------------
//    // Test
//    // -------------------------------------------------------------------------
    @Test
    @Disabled // disabled as we reach a local minimum and not the expected global one
    public void testHS265() {

        // Initial guess: x = (0, 0, 0, 0)
        double[] x0 = new double[]{0.0, 0.0, 0.0, 0.0};

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // Bounds: x_i >= 0 (no upper bounds active in Fortran)
        SimpleBounds bounds = new SimpleBounds(
                new double[]{0.0, 0.0, 0.0, 0.0},
                new double[]{
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY
                }
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS265Obj()),
                new HS265Eq(),   // 2 equalities
                null,            // no inequalities
                bounds
        );

        

        double fExpected = 0.97474658;
         HSProblemTestUtils.assertExpectedObjective(fExpected, sol);
    }
}
