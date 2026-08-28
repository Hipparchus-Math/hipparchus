package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

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

/**
 * HS119 — Traduzione fedele del problema TP119 (Schittkowski).
 * N = 16 variabili, 8 vincoli di uguaglianza, 0 <= x[i] <= 5.
 * Obiettivo: FX = Σ_i Σ_j A(i,j)*(x_i^2 + x_i + 1)*(x_j^2 + x_j + 1)
 * con A definita come nel codice Fortran originale.
 */
public class HS119Test {

    private static final double[] C = { 2.5, 1.1, -3.1, -3.5, 1.3, 2.1, 2.3, -1.5 };

    /** q(x) = x^2 + x + 1 */
    private static double q(double x) { return x * x + x + 1.0; }

    /** Obiettivo: FX = sum_{i,j} A(i,j)*q(x_i)*q(x_j). */
    private static class HS119Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 16; }

        @Override
        public double value(RealVector X) {
            final double[] x = new double[17];
            for (int i = 1; i <= 16; i++) x[i] = X.getEntry(i - 1);

            // Costruisci A come nel Fortran (solo dove A(i,j)=1)
            boolean[][] A = new boolean[17][17];
            for (int i = 1; i <= 16; i++) A[i][i] = true;
            // termini extra di A(i,j)=1
            int[][] extra = {
                    {1,4},{1,7},{1,8},{1,16},
                    {2,3},{2,7},{2,10},
                    {3,7},{3,9},{3,10},{3,14},
                    {4,7},{4,11},{4,15},
                    {5,6},{5,10},{5,12},{5,16},
                    {6,8},{6,15},
                    {7,11},{7,13},
                    {8,10},{8,15},
                    {9,12},{9,16},
                    {10,14},
                    {11,13},
                    {12,14},
                    {13,14}
            };
            for (int[] e : extra) A[e[0]][e[1]] = true;

            double T = 0.0;
            for (int i = 1; i <= 16; i++) {
                for (int j = 1; j <= 16; j++) {
                    if (A[i][j])
                        T += q(x[i]) * q(x[j]);
                }
            }
            return T;
        }

        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    /** Vincoli s[i] = c[i] come in Fortran (G(i)=Σ B(i,j)x_j - C(i)=0). */
    private static class HS119Eq extends EqualityConstraint {
        HS119Eq() { super(new ArrayRealVector(new double[8])); }
        @Override public int dim() { return 16; }

        @Override
        public RealVector value(RealVector X) {
            final double[] x = new double[17];
            for (int i = 1; i <= 16; i++) x[i] = X.getEntry(i - 1);

            double[] g = new double[8];

            g[0] =  0.22*x[1] + 0.20*x[2] + 0.19*x[3] + 0.25*x[4] + 0.15*x[5]
                   +0.11*x[6] + 0.12*x[7] + 0.13*x[8] + 1.0*x[9]  - C[0];

            g[1] = -1.46*x[1] -1.30*x[3] + 1.82*x[4] -1.15*x[5] + 0.80*x[7]
                   +1.0*x[10] - C[1];

            g[2] =  1.29*x[1] -0.89*x[2] -1.16*x[5] -0.96*x[6] -0.49*x[8]
                   +1.0*x[11] - C[2];

            g[3] = -1.10*x[1] -1.06*x[2] +0.95*x[3] -0.54*x[4] -1.78*x[6]
                   -0.41*x[7] +1.0*x[12] - C[3];

            g[4] = -1.43*x[4] +1.51*x[5] +0.59*x[6] -0.33*x[7] -0.43*x[8]
                   +1.0*x[13] - C[4];

            g[5] = -1.72*x[2] -0.33*x[3] +1.62*x[5] +1.24*x[6] +0.21*x[7]
                   -0.26*x[8] +1.0*x[14] - C[5];

            g[6] =  1.12*x[1] +0.31*x[4] +1.12*x[7] -0.36*x[9] +1.0*x[15]
                   - C[6];

            g[7] =  0.45*x[2] +0.26*x[3] -1.10*x[4] +0.58*x[5] -1.03*x[7]
                   +0.10*x[8] +1.0*x[16] - C[7];

            return new ArrayRealVector(g);
        }

        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testHS119_reference() {
        // soluzione nota dal Fortran
        double[] xEx = {
            0.0398473514099, 0.791983155694, 0.202870330224, 0.844357916347,
            1.26990645286, 0.934738707827, 1.68196196924, 0.155300877490,
            1.56787033356, 0.0, 0.0, 0.0, 0.660204066000, 0.0, 0.674255926901, 0.0
        };
        double fEx = 244.899697515;

        HS119Obj f = new HS119Obj();
        double val = f.value(new ArrayRealVector(xEx));

        assertEquals(fEx, val, 1e-6, "Objective mismatch at known solution");
    }

    @Test
    public void solveHS119() {
        double[] x0 = new double[16];
        for (int i = 0; i < 16; i++) x0[i] = 5.0;
        double[] lb = new double[16];
        double[] ub = new double[16];
        for (int i = 0; i < 16; i++) { lb[i] = 0.0; ub[i] = 5.0; }

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS119Obj()),
                new HS119Eq(),
                new SimpleBounds(lb, ub)
        );

        // controllo puntuale
        HSProblemTestUtils.assertExpectedObjective(244.8996975, sol);
    }
}
