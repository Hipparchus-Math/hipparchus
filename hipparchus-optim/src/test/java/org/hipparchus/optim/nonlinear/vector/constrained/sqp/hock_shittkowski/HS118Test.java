package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * HS118 — Traduzione fedele del problema TP118 (Schittkowski).
 * N = 15 variabili, 29 vincoli di disuguaglianza lineari (forma g(x) >= 0), bounds espliciti.
 * f(x) = somma, per m=0..4, dei blocchi:
 *   2.3*x[3m] + 1e-4*x[3m]^2 + 1.7*x[3m+1] + 1e-4*x[3m+1]^2 + 2.2*x[3m+2] + 1.5e-4*x[3m+2]^2
 */
public class HS118Test {

    /** Obiettivo con gradiente e Hessiana (diagonale) come nel Fortran TP118. */
    private static class HS118Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 15; }

        @Override
        public double value(RealVector X) {
            final double[] x = X.toArray();
            double T = 0.0;
            for (int m = 0; m < 5; m++) {
                int i = 3 * m;
                T += 2.3 * x[i]     + 1.0e-4 * x[i]     * x[i];
                T += 1.7 * x[i + 1] + 1.0e-4 * x[i + 1] * x[i + 1];
                T += 2.2 * x[i + 2] + 1.5e-4 * x[i + 2] * x[i + 2];
            }
            return T;
        }

        @Override
        public RealVector gradient(RealVector X) {
            final double[] x = X.toArray();
            final double[] g = new double[15];
            for (int m = 0; m < 5; m++) {
                int i = 3 * m;
                g[i]     = 2.3 + 2.0e-4 * x[i];          // df/dx_{3m}
                g[i + 1] = 1.7 + 2.0e-4 * x[i + 1];      // df/dx_{3m+1}
                g[i + 2] = 2.2 + 3.0e-4 * x[i + 2];      // df/dx_{3m+2}
            }
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector X) {
            // Hessiana diagonale, costante a blocchi:
            // d2f/dx_{3m}^2   = 2e-4
            // d2f/dx_{3m+1}^2 = 2e-4
            // d2f/dx_{3m+2}^2 = 3e-4
            double[][] H = new double[15][15];
            for (int m = 0; m < 5; m++) {
                int i = 3 * m;
                H[i][i]         = 2.0e-4;
                H[i + 1][i + 1] = 2.0e-4;
                H[i + 2][i + 2] = 3.0e-4;
            }
            return new Array2DRowRealMatrix(H, false);
        }
    }

    /**
     * Vincoli G(x) >= 0 come in TP118:
     * Per i = 1..4 (qui indicizzato 0..3):
     *  g1:  x(3i+1) - x(3i-2) + 7  >= 0
     *  g2:  x(3i+2) - x(3i-1) + 7  >= 0
     *  g3:  x(3i+3) - x(3i)   + 7  >= 0
     *  g4:  x(3i-2) - x(3i+1) + 6  >= 0
     *  g5:  x(3i-1) - x(3i+2) + 7  >= 0
     *  g6:  x(3i)   - x(3i+3) + 6  >= 0
     *
     * e poi le 5 somme:
     *  x1+x2+x3   - 60 >= 0
     *  x4+x5+x6   - 50 >= 0
     *  x7+x8+x9   - 70 >= 0
     *  x10+x11+x12- 85 >= 0
     *  x13+x14+x15-100 >= 0
     */
    /** Vincoli “0 ≤ expr ≤ cap” scritti come coppie di disuguaglianze g(x) ≥ 0. */
private static class HS118Ineq extends InequalityConstraint {
    HS118Ineq() { super(new ArrayRealVector(new double[29])); }
    @Override public int dim() { return 15; }

    @Override
    public RealVector value(RealVector X) {
        final double[] x = X.toArray();
        final double[] g = new double[29];
        int k = 0;

        // 1) Catena 1: 0 ≤ x4−x1+7 ≤ 13  (upper ⇒ x1−x4+6 ≥ 0)
        g[k++] = x[3]  - x[0] + 7.0;
        g[k++] = x[0]  - x[3] + 6.0;

        // 2) Catena 2: 0 ≤ x7−x4+7 ≤ 13
        g[k++] = x[6]  - x[3] + 7.0;
        g[k++] = x[3]  - x[6] + 6.0;

        // 3) Catena 3: 0 ≤ x10−x7+7 ≤ 13
        g[k++] = x[9]  - x[6] + 7.0;
        g[k++] = x[6]  - x[9] + 6.0;

        // 4) Catena 4: 0 ≤ x13−x10+7 ≤ 13
        g[k++] = x[12] - x[9] + 7.0;
        g[k++] = x[9]  - x[12] + 6.0;

        // 5) Catena 5: 0 ≤ x5−x2+7 ≤ 14  (upper ⇒ x2−x5+7 ≥ 0)
        g[k++] = x[4]  - x[1] + 7.0;
        g[k++] = x[1]  - x[4] + 7.0;

        // 6) Catena 6: 0 ≤ x8−x5+7 ≤ 14
        g[k++] = x[7]  - x[4] + 7.0;
        g[k++] = x[4]  - x[7] + 7.0;

        // 7) Catena 7: 0 ≤ x11−x8+7 ≤ 14
        g[k++] = x[10] - x[7] + 7.0;
        g[k++] = x[7]  - x[10] + 7.0;

        // 8) Catena 8: 0 ≤ x14−x11+7 ≤ 14
        g[k++] = x[13] - x[10] + 7.0;
        g[k++] = x[10] - x[13] + 7.0;

        // 9) Catena 9: 0 ≤ x6−x3+7 ≤ 13  (upper ⇒ x3−x6+6 ≥ 0)
        g[k++] = x[5]  - x[2] + 7.0;
        g[k++] = x[2]  - x[5] + 6.0;

        // 10) Catena 10: 0 ≤ x9−x6+7 ≤ 13
        g[k++] = x[8]  - x[5] + 7.0;
        g[k++] = x[5]  - x[8] + 6.0;

        // 11) Catena 11: 0 ≤ x12−x9+7 ≤ 13
        g[k++] = x[11] - x[8] + 7.0;
        g[k++] = x[8]  - x[11] + 6.0;

        // 12) Catena 12: 0 ≤ x15−x12+7 ≤ 13
        g[k++] = x[14] - x[11] + 7.0;
        g[k++] = x[11] - x[14] + 6.0;

        // Somme (identiche al modello)
        g[k++] = x[0] + x[1] + x[2]   - 60.0;
        g[k++] = x[3] + x[4] + x[5]   - 50.0;
        g[k++] = x[6] + x[7] + x[8]   - 70.0;
        g[k++] = x[9] + x[10] + x[11] - 85.0;
        g[k++] = x[12] + x[13] + x[14] - 100.0;

        return new ArrayRealVector(g, false);
    }

    @Override
    public RealMatrix jacobian(RealVector X) {
        // Jacobiano (29x15) costante, derivato dalle righe qui sopra.
        final double[][] J = new double[29][15];
        int r = 0;

        // Per ciascuna riga sopra, riempi le colonne con ±1 come in value().
        // 1)  x4 - x1 + 7
        J[r][3] = +1; J[r][0] = -1; r++;
        //     x1 - x4 + 6
        J[r][0] = +1; J[r][3] = -1; r++;

        // 2)  x7 - x4 + 7
        J[r][6] = +1; J[r][3] = -1; r++;
        //     x4 - x7 + 6
        J[r][3] = +1; J[r][6] = -1; r++;

        // 3)  x10 - x7 + 7
        J[r][9] = +1; J[r][6] = -1; r++;
        //     x7 - x10 + 6
        J[r][6] = +1; J[r][9] = -1; r++;

        // 4)  x13 - x10 + 7
        J[r][12] = +1; J[r][9]  = -1; r++;
        //     x10 - x13 + 6
        J[r][9]  = +1; J[r][12] = -1; r++;

        // 5)  x5 - x2 + 7
        J[r][4] = +1; J[r][1] = -1; r++;
        //     x2 - x5 + 7
        J[r][1] = +1; J[r][4] = -1; r++;

        // 6)  x8 - x5 + 7
        J[r][7] = +1; J[r][4] = -1; r++;
        //     x5 - x8 + 7
        J[r][4] = +1; J[r][7] = -1; r++;

        // 7)  x11 - x8 + 7
        J[r][10] = +1; J[r][7]  = -1; r++;
        //     x8 - x11 + 7
        J[r][7]  = +1; J[r][10] = -1; r++;

        // 8)  x14 - x11 + 7
        J[r][13] = +1; J[r][10] = -1; r++;
        //     x11 - x14 + 7
        J[r][10] = +1; J[r][13] = -1; r++;

        // 9)  x6 - x3 + 7
        J[r][5] = +1; J[r][2] = -1; r++;
        //     x3 - x6 + 6
        J[r][2] = +1; J[r][5] = -1; r++;

        // 10) x9 - x6 + 7
        J[r][8] = +1; J[r][5] = -1; r++;
        //     x6 - x9 + 6
        J[r][5] = +1; J[r][8] = -1; r++;

        // 11) x12 - x9 + 7
        J[r][11] = +1; J[r][8]  = -1; r++;
        //     x9 - x12 + 6
        J[r][8]  = +1; J[r][11] = -1; r++;

        // 12) x15 - x12 + 7
        J[r][14] = +1; J[r][11] = -1; r++;
        //     x12 - x15 + 6
        J[r][11] = +1; J[r][14] = -1; r++;

        // somme (righe finali)
        J[r][0]=1; J[r][1]=1; J[r][2]=1; r++;
        J[r][3]=1; J[r][4]=1; J[r][5]=1; r++;
        J[r][6]=1; J[r][7]=1; J[r][8]=1; r++;
        J[r][9]=1; J[r][10]=1; J[r][11]=1; r++;
        J[r][12]=1; J[r][13]=1; J[r][14]=1; r++;

        return new Array2DRowRealMatrix(J, false);
    }
}


    @Test
    public void testHS118() {
        final double[] x0 = new double[15];
        for (int i = 0; i < 15; i++) x0[i] = 20.0;
        x0[1]=44.0;
        x0[2]=4.0;
        
        final double[] lb = {
            8, 43, 3,
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,
            0, 0, 0
        };
        final double[] ub = {
            21, 57, 16,
            90, 120, 60,
            90, 120, 60,
            90, 120, 60,
            90, 120, 60
        };

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
            new InitialGuess(x0),
            new ObjectiveFunction(new HS118Obj()),
            new HS118Ineq(),
            new SimpleBounds(lb, ub)
        );
//FEX=0.664820449993D+03 
        HSProblemTestUtils.assertExpectedObjective(664.820449993, sol);
    }
}
