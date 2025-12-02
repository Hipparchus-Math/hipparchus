/*
 * HS274–HS276 (TP274/TP275/TP276)
 *
 * Same model with different dimensions N = 2, 4, 6:
 *
 * COMMON /D274/ A(6,6)
 *
 * Initialization (MODE=1):
 *   N = 2  (TP274), 4 (TP275), 6 (TP276)
 *   NILI = NINL = NELI = NENL = 0
 *   X(i)   = -4 / i,  i = 1..N
 *   A(i,j) = 1 / (i + j - 1)   (Hilbert matrix)
 *   XEX(i) = 0
 *   FEX    = 0
 *
 * Objective (MODE=2):
 *   FX = sum_{i=1..N} sum_{j=1..N} A(i,j) * X(i) * X(j)
 *
 * Gradient (MODE=3):
 *   GF(i) = sum_{j=1..N} X(j) * (A(i,j) + A(j,i))
 *
 * --> In Java:
 *   A(i,j) = 1 / (i + j - 1)  (we recompute on the fly)
 *   f(x)   = x^T A x
 *   grad   = (A + A^T) x
 * For Hilbert, A is symmetric, so grad = 2 A x, Hessian = 2 A.
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

public class HS274ToHS276Test {

    /**
     * Quadratic objective with Hilbert-like matrix:
     * A(i,j) = 1 / (i + j - 1), i,j = 1..N
     *
     * f(x)   = x^T A x
     * grad_i = sum_j x_j (A(i,j) + A(j,i))
     */
    private static class HilbertQuadratic extends TwiceDifferentiableFunction {

        private final int n;

        HilbertQuadratic(int n) {
            this.n = n;
        }

        @Override
        public int dim() {
            return n;
        }

        /** Fortran A(i,j) = 1.0 / (i + j - 1), with i,j = 1..N. */
        private double A(int i0, int j0) {
            // i0,j0 are 0-based; Fortran uses i=j=1..N
            int ip1 = i0 + 1;
            int jp1 = j0 + 1;
            return 1.0 / (ip1 + jp1 - 1.0);
        }

        @Override
        public double value(RealVector x) {
            double[] X = x.toArray();
            double fx = 0.0;

            for (int i = 0; i < n; i++) {
                double xi = X[i];
                for (int j = 0; j < n; j++) {
                    fx += A(i, j) * xi * X[j];
                }
            }

            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] X = x.toArray();
            double[] g = new double[n];

            // GF(i) = sum_j X(j) * (A(i,j) + A(j,i))
            for (int i = 0; i < n; i++) {
                double gi = 0.0;
                for (int j = 0; j < n; j++) {
                    gi += X[j] * (A(i, j) + A(j, i));
                }
                g[i] = gi;
            }

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Constant Hessian: H(i,j) = A(i,j) + A(j,i) = 2*A(i,j) (symmetric)
            double[][] H = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    H[i][j] = A(i, j) + A(j, i);
                }
            }
            return new Array2DRowRealMatrix(H, false);
        }
    }

    // -------------------------------------------------------------------------
    // TP274: N = 2
    // -------------------------------------------------------------------------
    @Test
    public void testHS274_N2() {

        int n = 2;

        // Initial guess X(i) = -4 / i (Fortran i = 1..N)
        double[] x0 = new double[n];
        for (int i = 0; i < n; i++) {
            int ip1 = i + 1;
            x0[i] = -4.0 / ip1;
        }

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HilbertQuadratic(n)),
                null,   // no equalities
                null,   // no inequalities
                null    // no bounds
        );

        double f = sol.getValue();
        double fExpected = 0.0;
        double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }

    // -------------------------------------------------------------------------
    // TP275: N = 4
    // -------------------------------------------------------------------------
    @Test
    public void testHS275_N4() {

        int n = 4;

        double[] x0 = new double[n];
        for (int i = 0; i < n; i++) {
            int ip1 = i + 1;
            x0[i] = -4.0 / ip1;
        }

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HilbertQuadratic(n)),
                null,
                null,
                null
        );

        double f = sol.getValue();
        double fExpected = 0.0;
        double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }

    // -------------------------------------------------------------------------
    // TP276: N = 6
    // -------------------------------------------------------------------------
    @Test
    public void testHS276_N6() {

        int n = 6;

        double[] x0 = new double[n];
        for (int i = 0; i < n; i++) {
            int ip1 = i + 1;
            x0[i] = -4.0 / ip1;
        }

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HilbertQuadratic(n)),
                null,
                null,
                null
        );

        double f = sol.getValue();
        double fExpected = 0.0;
        double tol = 1e-4 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
