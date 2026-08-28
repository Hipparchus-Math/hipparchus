/*
 * HS271 (TP271)
 *
 * N    = 6
 * No constraints
 *
 * In Fortran:
 *   F(i) = sqrt(10*(16 - i)) * (x_i - 1),    i = 1..6
 *   FX   = sum_i F(i)^2
 *
 * So:
 *   f(x) = sum_{i=1..6} 10*(16 - i) * (x_i - 1)^2
 *   grad_i = 20*(16 - i)*(x_i - 1)
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

public class HS271Test {

    private static final int DIM = 6;

    private static class HS271Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double fx = 0.0;
            for (int i = 0; i < DIM; i++) {
                int idx = i + 1; // Fortran index
                double wi = 10.0 * (16.0 - idx);
                double diff = x.getEntry(i) - 1.0;
                fx += wi * diff * diff;
            }
            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] g = new double[DIM];
            for (int i = 0; i < DIM; i++) {
                int idx = i + 1;
                double coeff = 20.0 * (16.0 - idx);
                g[i] = coeff * (x.getEntry(i) - 1.0);
            }
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Diagonal Hessian with 2*wi, but we let BFGS build it
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    @Test
    public void testHS271() {

        double[] x0 = new double[DIM]; // all zeros

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS271Obj()),
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
