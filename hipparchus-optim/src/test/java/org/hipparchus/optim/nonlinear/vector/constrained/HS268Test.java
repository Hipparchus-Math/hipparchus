/*
 * HS268 (TP268)
 *
 * N    = 5
 * NILI = 5  (5 linear inequalities)
 * NINL = 0
 * NELI = 0
 * NENL = 0
 *
 * Objective (from Fortran):
 *   FX = DVDV + sum_{i=1..5} x_i * ( (DD x)_i - 2 * DDVEKT_i )
 *
 * Bounds: none
 *
 * Initial guess:
 *   x0 = (1, 1, 1, 1, 1)
 *
 * Reference solution (Fortran):
 *   x*  = (1, 2, -1, 3, -4)
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

public class HS268Test {

    private static final int DIM      = 5;
    private static final int NUM_INEQ = 5;

    // DVDV = DVEKT' * DVEKT
    private static final double DVDV = 14463.0;

    // DD(5,5) from Fortran DATA DD / ... / (symmetric)
    // DD(i,j) = DD(i+1,j+1)
    private static final double[][] DD = {
            { 10197.0, -12454.0, -1013.0,  1948.0,   329.0 },
            { -12454.0, 20909.0, -1733.0, -4914.0,  -186.0 },
            { -1013.0,  -1733.0,  1755.0,  1089.0,  -174.0 },
            {  1948.0,  -4914.0,  1089.0,  1515.0,   -22.0 },
            {   329.0,   -186.0,  -174.0,   -22.0,    27.0 }
    };

    // DDVEKT(5) from Fortran DATA DDVEKT
    private static final double[] DDVEKT = {
            -9170.0, 17099.0, -2271.0, -4336.0, -43.0
    };

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS268Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double[] X = x.toArray();
            double fx = DVDV;

            // FX = DVDV + sum_i X_i * ( (DD X)_i - 2*DDVEKT_i )
            for (int i = 0; i < DIM; i++) {
                double hf = 0.0;
                for (int j = 0; j < DIM; j++) {
                    hf += DD[i][j] * X[j];
                }
                fx += X[i] * (hf - 2.0 * DDVEKT[i]);
            }

            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] X = x.toArray();
            double[] g = new double[DIM];

            // Fortran:
            // GF(I) = - 2*DDVEKT(I) + sum_j (DD(I,J)+DD(J,I))*X(J)
            for (int i = 0; i < DIM; i++) {
                double gi = -2.0 * DDVEKT[i];
                for (int j = 0; j < DIM; j++) {
                    gi += (DD[i][j] + DD[j][i]) * X[j];
                }
                g[i] = gi;
            }

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // True Hessian = (DD+DD^T), but we let BFGS build it.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints G(x) ≥ 0 (NILI = 5)
    //
    // Fortran:
    // G(1) = -x1 -x2 -x3 -x4 -x5 + 5
    // G(2) = 10*x1 +10*x2 -3*x3 +5*x4 +4*x5 -20
    // G(3) = -8*x1 + x2 -2*x3 -5*x4 +3*x5 +40
    // G(4) =  8*x1 - x2 +2*x3 +5*x4 -3*x5 -11     <-- 0.11D+2 = 11.0
    // G(5) = -4*x1 -2*x2 +3*x3 -5*x4 + x5 +30
    // -------------------------------------------------------------------------
    private static class HS268Ineq extends InequalityConstraint {

        HS268Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQ]));
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

            double g1 = -x1 - x2 - x3 - x4 - x5 + 5.0;
            double g2 = 10.0 * x1 + 10.0 * x2 - 3.0 * x3 + 5.0 * x4 + 4.0 * x5 - 20.0;
            double g3 = -8.0 * x1 + x2 - 2.0 * x3 - 5.0 * x4 + 3.0 * x5 + 40.0;
            double g4 =  8.0 * x1 - x2 + 2.0 * x3 + 5.0 * x4 - 3.0 * x5 - 11.0; // FIXED
            double g5 = -4.0 * x1 - 2.0 * x2 + 3.0 * x3 - 5.0 * x4 + x5 + 30.0;

            return new ArrayRealVector(new double[]{g1, g2, g3, g4, g5}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // G1
            J.setEntry(0, 0, -1.0);
            J.setEntry(0, 1, -1.0);
            J.setEntry(0, 2, -1.0);
            J.setEntry(0, 3, -1.0);
            J.setEntry(0, 4, -1.0);

            // G2
            J.setEntry(1, 0, 10.0);
            J.setEntry(1, 1, 10.0);
            J.setEntry(1, 2, -3.0);
            J.setEntry(1, 3, 5.0);
            J.setEntry(1, 4, 4.0);

            // G3
            J.setEntry(2, 0, -8.0);
            J.setEntry(2, 1, 1.0);
            J.setEntry(2, 2, -2.0);
            J.setEntry(2, 3, -5.0);
            J.setEntry(2, 4, 3.0);

            // G4  (same derivatives, only constant term changed)
            J.setEntry(3, 0, 8.0);
            J.setEntry(3, 1, -1.0);
            J.setEntry(3, 2, 2.0);
            J.setEntry(3, 3, 5.0);
            J.setEntry(3, 4, -3.0);

            // G5
            J.setEntry(4, 0, -4.0);
            J.setEntry(4, 1, -2.0);
            J.setEntry(4, 2, 3.0);
            J.setEntry(4, 3, -5.0);
            J.setEntry(4, 4, 1.0);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS268() {

        double[] x0 = new double[]{1.0, 1.0, 1.0, 1.0, 1.0};

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        // No box bounds in this problem
        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS268Obj()),
                null,              // no equalities
                new HS268Ineq(),   // 5 inequalities
                null               // no bounds
        );

        double f = sol.getValue();

        double fExpected = 0.0;
        double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
