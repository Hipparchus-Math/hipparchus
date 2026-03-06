package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.CholeskyDecomposition;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

public final class CholeskyDiagThresholdReg {

    private CholeskyDiagThresholdReg() { }

    public static final class Result {
        public final RealMatrix L;
        public final boolean regularized;
        public final double mu;
        public final double s2;

        private Result(final RealMatrix L, final boolean regularized, double mu, double s2) {
            this.L = L;
            this.regularized = regularized;
            this.mu = mu;
            this.s2 = s2;
        }
    }

    public static Result regularizeByDiagonalThreshold(final RealMatrix Lin) {
    final int n = Lin.getRowDimension();
    final double eps = Precision.EPSILON;

    // 1. Calcolo dei parametri di scala (Standard in SNOPT/MINOS)
    // Usiamo la norma infinito della diagonale per definire la "tolleranza di vita"
    double maxHii = 0.0;
    for (int i = 0; i < n; i++) {
        maxHii = FastMath.max(maxHii, FastMath.pow(Lin.getEntry(i, i), 2));
    }

    // Soglia di "morte numerica": se un pivot è sotto questa, il solver duale esplode.
    // In produzione si usa sqrt(eps) per garantire che L^-1 sia stabile.
    final double delta = FastMath.max(1e-12, FastMath.sqrt(eps) * maxHii);

    // 2. Ricostruzione H = L*L^T (Necessaria perché L è "viva" ma potenzialmente sporca)
    RealMatrix H = Lin.multiply(Lin.transpose());
    
    // 3. ALGORITMO DI GILL-MURRAY (Modified Cholesky)
    // Non aggiungiamo mu a tutti. Modifichiamo solo i pivot che falliscono.
    RealMatrix L = MatrixUtils.createRealMatrix(n, n);
    double muMax = 0.0;
    boolean regularized = false;

    for (int j = 0; j < n; j++) {
        double d = H.getEntry(j, j);
        for (int k = 0; k < j; k++) {
            double ljk = L.getEntry(j, k);
            d -= ljk * ljk;
        }

        // --- IL CUORE DEL SOLVER COMMERCIALE ---
        // Invece di far fallire la fattorizzazione, se il pivot d è troppo piccolo
        // lo "alziamo" forzatamente alla soglia delta.
        if (d < delta) {
            double mu = delta - d;
            d = delta;
            muMax = FastMath.max(muMax, mu);
            regularized = true;
        }

        double ljj = FastMath.sqrt(d);
        L.setEntry(j, j, ljj);

        for (int i = j + 1; i < n; i++) {
            double s = H.getEntry(i, j);
            for (int k = 0; k < j; k++) {
                s -= L.getEntry(i, k) * L.getEntry(j, k);
            }
            L.setEntry(i, j, s / ljj);
        }
    }

    return new Result(L, regularized, muMax, delta);
}
}