package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.*;
import org.hipparchus.util.FastMath;

/**
 * Implementazione rigorosa dell'algoritmo di Gill-Murray-Wright 
 * per la fattorizzazione di Cholesky modificata (Pure Cholesky L).
 */
public class GillMurrayWrightCholesky {



    public static RealMatrix getRegulatedMatrix(RealMatrix A, double delta) {
        int n = A.getRowDimension();
        double[][] a = A.getData(); 
        int[] piv = new int[n];
        for (int i = 0; i < n; i++) piv[i] = i;

        // 1. Calcolo rigoroso di betaSqr (Pag. 88 modchol.c)
        double gamma = 0, xi = 0;
        for (int i = 0; i < n; i++) {
            gamma = FastMath.max(gamma, FastMath.abs(a[i][i]));
            for (int j = i + 1; j < n; j++) {
                xi = FastMath.max(xi, FastMath.abs(a[i][j]));
            }
        }
        double eps = 2.22e-16; // DBL_EPSILON
        double betaSq = FastMath.max(gamma, FastMath.max(xi / FastMath.sqrt(n * n - 1.0), eps));

        double[] d = new double[n];
        double[][] l = new double[n][n];

        // 2. Fattorizzazione con Diagonal Pivoting (GMW)
        for (int j = 0; j < n; j++) {
            int p = j;
            for (int i = j + 1; i < n; i++) {
                if (FastMath.abs(a[i][i]) > FastMath.abs(a[p][p])) p = i;
            }

            if (p != j) {
                swapSym(a, j, p);
                int t = piv[j]; piv[j] = piv[p]; piv[p] = t;
            }

            l[j][j] = 1.0;
            double cjj = a[j][j];
            for (int s = 0; s < j; s++) cjj -= l[j][s] * l[j][s] * d[s];

            double thetaJ = 0;
            for (int i = j + 1; i < n; i++) {
                double cij = a[i][j];
                for (int s = 0; s < j; s++) cij -= l[i][s] * l[j][s] * d[s];
                a[i][j] = cij;
                thetaJ = FastMath.max(thetaJ, FastMath.abs(cij));
            }

            // Regola di Gill-Murray-Wright
            d[j] = FastMath.max(FastMath.abs(cjj), FastMath.max(delta, (thetaJ * thetaJ) / betaSq));

            for (int i = j + 1; i < n; i++) l[i][j] = a[i][j] / d[j];
        }

        // 3. Ricostruzione Matrice Regolarizzata (Permutata)
        RealMatrix L = MatrixUtils.createRealMatrix(l);
        RealMatrix D = MatrixUtils.createRealDiagonalMatrix(d);
        RealMatrix M = L.multiply(D).multiply(L.transpose());
        
        // 4. UN-PIVOTING (Fondamentale: riporta all'ordine originale)
        double[][] result = new double[n][n];
        double[][] mData = M.getData();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                // Riordiniamo gli indici usando il vettore piv
                result[piv[i]][piv[j]] = mData[i][j];
            }
        }
        return MatrixUtils.createRealMatrix(result);
    }

    private static void swapSym(double[][] a, int i, int j) {
        double[] temp = a[i]; a[i] = a[j]; a[j] = temp;
        for (int k = 0; k < a.length; k++) {
            double t = a[k][i]; a[k][i] = a[k][j]; a[k][j] = t;
        }
    }
}