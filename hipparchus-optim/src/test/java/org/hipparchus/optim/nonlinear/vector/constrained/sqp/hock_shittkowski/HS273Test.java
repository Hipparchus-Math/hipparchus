/*
 * HS273 (TP273)
 *
 * N    = 6
 * No constraints
 *
 * H(X) = sum_{i=1..6} (16 - i) * (x_i - 1)^2
 *
 * Fortran:
 *   HX = TP273A(X)
 *   FX = 10 * HX * (1 + HX)
 *
 * So:
 *   f(x) = 10 * H * (1 + H)
 *
 * Gradient (Fortran):
 *   GF(i) = 20 * (16 - i) * (x_i - 1) * (1 + 2*H)
 *
 * Initial guess:
 *   x0 = (0,0,0,0,0,0)
 *
 * Reference solution:
 *   x* = (1,1,1,1,1,1)
 *   f* = 0
 */

package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS273Test {

    private static final int DIM = 6;

    private static class HS273Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        private double computeH(double[] X) {
            double h = 0.0;
            for (int i = 0; i < DIM; i++) {
                int idx = i + 1;
                double w = 16.0 - idx;
                double diff = X[i] - 1.0;
                h += w * diff * diff;
            }
            return h;
        }

        @Override
        public double value(RealVector x) {
            double[] X = x.toArray();
            double h = computeH(X);
            return 10.0 * h * (1.0 + h);
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] X = x.toArray();
            double h = computeH(X);

            double[] g = new double[DIM];
            for (int i = 0; i < DIM; i++) {
                int idx = i + 1;
                double w = 16.0 - idx;
                g[i] = 20.0 * w * (X[i] - 1.0) * (1.0 + 2.0 * h);
            }
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Hessian is messy; let BFGS approximate it
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    @Test
    public void testHS273() {

        double[] x0 = new double[DIM]; // all zeros

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS273Obj()),
                null,   // no equalities
                null,   // no inequalities
                null    // no bounds
        );

        double f = sol.getValue();
        double fExpected = 0.0;
        double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);
        assertEquals(fExpected, f, tol);
    }
}
