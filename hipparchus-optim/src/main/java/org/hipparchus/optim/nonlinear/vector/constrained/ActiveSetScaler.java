package org.hipparchus.optim.nonlinear.vector.constrained;

import java.util.Arrays;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * Active-set style scaling (variable diagonal scaling + row scaling on constraints).
 * Allocates new matrices/vectors and fills them, without using copy().
 *
 * Intended usage:
 *   - Keep originals as G1/g01/CE1/... (untouched)
 *   - Build scaled working data G/g0/CE/... with scaleFill(...)
 *   - Solve on scaled data
 *   - Unscale x and lambda with unscalePrimal/unscaleDual
 */
public final class ActiveSetScaler {

    private static final double MIN_DIAG      = 1e-12;
    private static final double VAR_SCALE_MIN = 1e-6;
    private static final double VAR_SCALE_MAX = 1e6;
    private static final double ROW_NORM_MIN  = 1e-12;

    /** Variable scaling (from diag(G)), size n. */
    private double[] varScale;

    /** Equality constraint row scaling factors (column norms), size p. */
    private double[] eqRowScale;

    /** Inequality constraint row scaling factors (column norms), size m. */
    private double[] ineqRowScale;

    /** Number of equality constraints. */
    private int pEq;

    /** Number of inequality constraints. */
    private int mIneq;

    public static final class ScaledData {
        public final RealMatrix G;
        public final RealVector g0;
        public final RealMatrix CE;
        public final RealVector ce0;
        public final RealMatrix CI;
        public final RealVector ci0;

        private ScaledData(RealMatrix G, RealVector g0,
                           RealMatrix CE, RealVector ce0,
                           RealMatrix CI, RealVector ci0) {
            this.G = G;
            this.g0 = g0;
            this.CE = CE;
            this.ce0 = ce0;
            this.CI = CI;
            this.ci0 = ci0;
        }
    }

    /**
     * Build a scaled problem into freshly allocated matrices/vectors (no copy()).
     * Originals are not modified.
     *
     * Scaling is applied in the mathematically consistent order:
     *  1) Variable scaling derived from diag(G): x = D * x_hat
     *     -> G_hat = D^T G D, g_hat = D^T g, A_hat = D^T A
     *  2) Row scaling of BOTH equality and inequality constraints:
     *     each constraint column is normalized by its 2-norm (after step 1).
     */
    public ScaledData scaleFill(final RealMatrix G1,
                                final RealVector g01,
                                final RealMatrix CE1,
                                final RealVector ce01,
                                final RealMatrix CI1,
                                final RealVector ci01) {

        final int n = G1.getRowDimension();
        this.pEq   = (CE1 != null) ? CE1.getColumnDimension() : 0;
        this.mIneq = (CI1 != null) ? CI1.getColumnDimension() : 0;

        // Allocate outputs (empty)
        final RealMatrix G  = MatrixUtils.createRealMatrix(n, n);
        final RealVector g0 = new ArrayRealVector(n);

        final RealMatrix CE = (CE1 != null) ? MatrixUtils.createRealMatrix(n, pEq) : null;
        final RealVector ce0 = (ce01 != null) ? new ArrayRealVector(pEq) : null;

        final RealMatrix CI = (CI1 != null) ? MatrixUtils.createRealMatrix(n, mIneq) : null;
        final RealVector ci0 = (ci01 != null) ? new ArrayRealVector(mIneq) : null;

        // ---------------------------------------------------------------------
        // 0) Copy rhs vectors (will be further scaled later)
        // ---------------------------------------------------------------------
        if (ce0 != null) {
            for (int j = 0; j < pEq; j++) {
                ce0.setEntry(j, ce01.getEntry(j));
            }
        }
        if (ci0 != null) {
            for (int j = 0; j < mIneq; j++) {
                ci0.setEntry(j, ci01.getEntry(j));
            }
        }

        // ---------------------------------------------------------------------
        // 1) Variable scaling from diag(G1)
        // ---------------------------------------------------------------------
        this.varScale = new double[n];
        for (int i = 0; i < n; i++) {
            double gii = Math.abs(G1.getEntry(i, i));
            gii = Math.max(gii, MIN_DIAG);
            double d = 1.0 / Math.sqrt(gii);
            if (d < VAR_SCALE_MIN) d = VAR_SCALE_MIN;
            if (d > VAR_SCALE_MAX) d = VAR_SCALE_MAX;
            varScale[i] = d;
        }

        // Fill G = D^{-1} * G1 * D^{-1}   (your convention; keep unchanged)
        for (int i = 0; i < n; i++) {
            final double di = varScale[i];
            for (int j = 0; j < n; j++) {
                final double dj = varScale[j];
                G.setEntry(i, j, G1.getEntry(i, j) / (di * dj));
            }
        }

        // Fill g0 = D^{-1} * g01
        for (int i = 0; i < n; i++) {
            g0.setEntry(i, g01.getEntry(i) / varScale[i]);
        }

        // Fill CE = D^{-1} * CE1
        if (CE != null) {
            for (int i = 0; i < n; i++) {
                final double di = varScale[i];
                for (int j = 0; j < pEq; j++) {
                    CE.setEntry(i, j, CE1.getEntry(i, j) / di);
                }
            }
        }

        // Fill CI = D^{-1} * CI1
        if (CI != null) {
            for (int i = 0; i < n; i++) {
                final double di = varScale[i];
                for (int j = 0; j < mIneq; j++) {
                    CI.setEntry(i, j, CI1.getEntry(i, j) / di);
                }
            }
        }

        // ---------------------------------------------------------------------
        // 2) Row scaling (column normalization) for BOTH equalities and inequalities
        //    IMPORTANT: this must be done AFTER variable scaling for coherence.
        // ---------------------------------------------------------------------

        // 2a) Equalities: normalize each CE(:,j) and scale ce0(j) accordingly
        if (CE != null) {
            this.eqRowScale = new double[pEq];
            Arrays.fill(this.eqRowScale, 1.0);

            for (int j = 0; j < pEq; j++) {
                double sumSq = 0.0;
                for (int i = 0; i < n; i++) {
                    final double aij = CE.getEntry(i, j);
                    sumSq += aij * aij;
                }
                final double norm = Math.sqrt(sumSq);

                if (norm < ROW_NORM_MIN) {
                    eqRowScale[j] = 1.0;
                    continue;
                }

                final double inv = 1.0 / norm;
                for (int i = 0; i < n; i++) {
                    CE.setEntry(i, j, CE.getEntry(i, j) * inv);
                }
                if (ce0 != null) {
                    ce0.setEntry(j, ce0.getEntry(j) * inv);
                }
                eqRowScale[j] = norm;
            }
        } else {
            this.eqRowScale = null;
        }

        // 2b) Inequalities: normalize each CI(:,j) and scale ci0(j) accordingly
        if (CI != null) {
            this.ineqRowScale = new double[mIneq];
            Arrays.fill(this.ineqRowScale, 1.0);

            for (int j = 0; j < mIneq; j++) {
                double sumSq = 0.0;
                for (int i = 0; i < n; i++) {
                    final double aij = CI.getEntry(i, j);
                    sumSq += aij * aij;
                }
                final double norm = Math.sqrt(sumSq);

                if (norm < ROW_NORM_MIN) {
                    ineqRowScale[j] = 1.0;
                    continue;
                }

                final double inv = 1.0 / norm;
                for (int i = 0; i < n; i++) {
                    CI.setEntry(i, j, CI.getEntry(i, j) * inv);
                }
                if (ci0 != null) {
                    ci0.setEntry(j, ci0.getEntry(j) * inv);
                }
                ineqRowScale[j] = norm;
            }
        } else {
            this.ineqRowScale = null;
        }

        return new ScaledData(G, g0, CE, ce0, CI, ci0);
    }

    public RealVector unscalePrimal(final RealVector xScaled) {
        final double[] x = xScaled.toArray();
        for (int i = 0; i < x.length; i++) {
            x[i] /= varScale[i];
        }
        return new ArrayRealVector(x, false);
    }

    /** Full lambda layout: eq in [0..p-1], ineq in [p..p+m-1]. */
    public RealVector unscaleDual(final RealVector lambdaScaled) {
        final double[] lam = lambdaScaled.toArray();

        // Undo equality row scaling (if applied)
        if (eqRowScale != null) {
            final int endEq = Math.min(lam.length, pEq);
            for (int idx = 0; idx < endEq; idx++) {
                lam[idx] /= eqRowScale[idx];
            }
        }

        // Undo inequality row scaling (if applied)
        if (ineqRowScale != null) {
            final int start = pEq;
            final int end = Math.min(lam.length, pEq + mIneq);
            for (int idx = start; idx < end; idx++) {
                lam[idx] /= ineqRowScale[idx - pEq];
            }
        }

        return new ArrayRealVector(lam, false);
    }
}
