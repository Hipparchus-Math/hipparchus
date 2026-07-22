/*
 * HS272 (TP272)
 *
 * N    = 6
 * No constraints
 *
 * For i = 1..13, with H = 0.1 * i:
 *
 *   F(i) = x4*exp(-x1*H) - x5*exp(-x2*H) + x6*exp(-x3*H)
 *          - exp(-H) + 5*exp(-10*H) - 3*exp(-4*H)
 *
 * Objective:
 *   f(x) = sum_{i=1..13} F(i)^2
 *
 * Bounds:
 *   XL(i)=0, LXL(i)=TRUE → x_i >= 0
 *   LXU(i)=FALSE        → no upper bounds
 *
 * Initial guess:
 *   x0 = (1, 2, 1, 1, 1, 1)
 *
 * Reference solution (Fortran):
 *   x* = (1, 10, 4, 1, 5, 3)
 *   f* = 0
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

public class HS272Test {

    private static final int DIM = 6;
    private static final int M   = 13; // number of residuals

    private static class HS272Obj extends TwiceDifferentiableFunction {

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
            double x6 = x.getEntry(5);

            double fx = 0.0;

            for (int i = 1; i <= M; i++) {
                double h  = 0.1 * i;
                double eh = FastMath.exp(-h);

                double term =
                        x4 * FastMath.exp(-x1 * h)
                      - x5 * FastMath.exp(-x2 * h)
                      + x6 * FastMath.exp(-x3 * h)
                      - eh
                      + 5.0 * FastMath.exp(-10.0 * h)
                      - 3.0 * FastMath.exp(-4.0 * h);

                fx += term * term;
            }

            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);
            double x5 = x.getEntry(4);
            double x6 = x.getEntry(5);

            double g1 = 0.0;
            double g2 = 0.0;
            double g3 = 0.0;
            double g4 = 0.0;
            double g5 = 0.0;
            double g6 = 0.0;

            // DF(J,4) = exp(-x1*h)
            // DF(J,5) = -exp(-x2*h)
            // DF(J,6) = exp(-x3*h)
            // DF(J,1) = -h*x4*DF(J,4)
            // DF(J,2) = -h*x5*DF(J,5)
            // DF(J,3) = -h*x6*DF(J,6)
            for (int i = 1; i <= M; i++) {
                double h = 0.1 * i;

                double e1 = FastMath.exp(-x1 * h);
                double e2 = FastMath.exp(-x2 * h);
                double e3 = FastMath.exp(-x3 * h);
                double eh = FastMath.exp(-h);

                double r =
                        x4 * e1
                      - x5 * e2
                      + x6 * e3
                      - eh
                      + 5.0 * FastMath.exp(-10.0 * h)
                      - 3.0 * FastMath.exp(-4.0 * h);

                double dFdx1 = -h * x4 * e1;
                double dFdx2 =  h * x5 * e2;
                double dFdx3 = -h * x6 * e3;
                double dFdx4 = e1;
                double dFdx5 = -e2;
                double dFdx6 = e3;

                g1 += 2.0 * r * dFdx1;
                g2 += 2.0 * r * dFdx2;
                g3 += 2.0 * r * dFdx3;
                g4 += 2.0 * r * dFdx4;
                g5 += 2.0 * r * dFdx5;
                g6 += 2.0 * r * dFdx6;
            }

            return new ArrayRealVector(new double[]{g1, g2, g3, g4, g5, g6}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // True Hessian is complicated; let BFGS build it
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    @Test
    public void testHS272() {

        // Initial guess: X(1..6) = 1, then X(2) = 2
        double[] x0 = new double[]{1.0, 2.0, 1.0, 1.0, 1.0, 1.0};

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        // Bounds: x_i >= 0
        SimpleBounds bounds = new SimpleBounds(
                new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                new double[]{
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY
                }
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS272Obj()),
                null,   // no equalities
                null,   // no inequalities
                bounds
        );

        double f = sol.getValue();
        double fExpected = 0.0;
        double tol = 1e-2 * (FastMath.abs(fExpected) + 1.0);
        assertEquals(fExpected, f, tol);
    }
}
