package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.Precision;

public class LDL {

    private final int n;
    private final double[][] lData;
    private final double[] d;
    private final double[] w;

    public LDL(final int n) {
        this.n = n;
        this.lData = new double[n][n];
        this.d = new double[n];
        this.w = new double[n];
        reset(1.0);
    }

    public final void reset(final double gamma) {
        for (int i = 0; i < n; i++) {
            d[i] = gamma;
            for (int j = 0; j < n; j++) {
                lData[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
    }

    public int getDimension() {
        return n;
    }

    public double[][] getLData() {
        return lData;
    }

    public double[] getD() {
        return d;
    }

    /**
     * Compute H * v = (L * D * L^T) * v in O(n^2).
     */
    public RealVector operate(final RealVector vVec) {
        final double[] v = vVec.toArray();

        // work = L^T * v
        final double[] work = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = v[i];
            for (int j = i + 1; j < n; j++) {
                sum += lData[j][i] * v[j];
            }
            work[i] = sum;
        }

        // work = D * work
        for (int i = 0; i < n; i++) {
            work[i] *= d[i];
        }

        // res = L * work
        final double[] res = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = work[i];
            for (int j = 0; j < i; j++) {
                sum += lData[i][j] * work[j];
            }
            res[i] = sum;
        }

        return new ArrayRealVector(res, false);
    }

    /**
     * Factorize H = L D L^T in-place, no pivoting.
     * Assumes H is SPD-ish.
     */
    public void factorizeNoPivot(final RealMatrix H) {

        // reset L to identity first
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                lData[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }

        for (int k = 0; k < n; k++) {

            double dk = H.getEntry(k, k);
            for (int j = 0; j < k; j++) {
                final double lkj = lData[k][j];
                dk -= lkj * lkj * d[j];
            }

            // do not clamp artificially; only avoid exact zero division
            d[k] = (dk == 0.0) ? Precision.EPSILON : dk;

            for (int i = k + 1; i < n; i++) {
                double num = H.getEntry(i, k);
                for (int j = 0; j < k; j++) {
                    num -= lData[i][j] * lData[k][j] * d[j];
                }
                lData[i][k] = num / d[k];
            }

            lData[k][k] = 1.0;
        }
    }

    /**
     * Fletcher-Powell / Kraft LDL^T rank-1 update:
     *
     *     H <- H + sigma * z z^T
     *
     * z is modified in-place, exactly like the original routine.
     */
    public boolean update(final RealVector zInput, final double sigma) {

        if (Math.abs(sigma) < Precision.SAFE_MIN) {
            return true;
        }

        final double[] z = zInput.toArray();
        double t = 1.0 / sigma;
        double[] tpList = null;

        // --------------------------------------------------
        // Negative update preparation
        // --------------------------------------------------
        if (sigma < 0.0) {
            tpList = new double[n];

            // w <- z
            System.arraycopy(z, 0, w, 0, n);

            for (int i = 0; i < n; i++) {
                final double di = d[i];
                if (Math.abs(di) < Precision.SAFE_MIN) {
                    return false;
                }

                final double v = w[i];
                t += (v * v) / di;

                // IMPORTANT: update w, NOT z
                for (int j = i + 1; j < n; j++) {
                    w[j] -= v * lData[j][i];
                }

                tpList[i] = v;
            }

            if (t >= 0.0) {
                t = Precision.EPSILON / sigma;
            }

            for (int i = n - 1; i >= 0; i--) {
                final double di = d[i];
                if (Math.abs(di) < Precision.SAFE_MIN) {
                    return false;
                }

                final double v = tpList[i];
                tpList[i] = t;
                t -= (v * v) / di;
            }

            // DO NOT reset t = 1/sigma here
        }

        // --------------------------------------------------
        // Main update loop
        // --------------------------------------------------
        for (int i = 0; i < n; i++) {
            final double di = d[i];
            if (Math.abs(di) < Precision.SAFE_MIN) {
                return false;
            }

            final double vi = z[i];
            final double delta = vi / di;
            final double tp = (sigma < 0.0) ? tpList[i] : (t + delta * vi);

            if (Math.abs(tp) < Precision.SAFE_MIN || Math.abs(t) < Precision.SAFE_MIN) {
                return false;
            }

            final double alpha = tp / t;
            d[i] = alpha * di;   // NO clamp here

            if (i < n - 1) {
                final double beta = delta / tp;

                if (alpha > 4.0) {
                    final double gamma = t / tp;
                    for (int j = i + 1; j < n; j++) {
                        final double lij = lData[j][i];
                        lData[j][i] = gamma * lij + beta * z[j];
                        z[j] -= vi * lij;
                    }
                } else {
                    for (int j = i + 1; j < n; j++) {
                        final double lij = lData[j][i];
                        final double zj = z[j] - vi * lij;
                        z[j] = zj;
                        lData[j][i] = lij + beta * zj;
                    }
                }
            }

            t = tp;
        }

        return true;
    }

    /**
     * Reconstruct full Hessian H = L D L^T.
     */
    public RealMatrix getHessian() {
        final double[][] h = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int k = 0; k <= i; k++) {
                double sum = 0.0;
                for (int j = 0; j <= Math.min(i, k); j++) {
                    final double lij = (i == j) ? 1.0 : lData[i][j];
                    final double lkj = (k == j) ? 1.0 : lData[k][j];
                    sum += lij * d[j] * lkj;
                }
                h[i][k] = sum;
                h[k][i] = sum;
            }
        }

        return new Array2DRowRealMatrix(h, false);
    }
    public void copyTo(final double[][] lCopy, final double[] dCopy) {
    if (lCopy.length != n || dCopy.length != n) {
        throw new IllegalArgumentException("Dimension mismatch");
    }
    for (int i = 0; i < n; i++) {
        System.arraycopy(lData[i], 0, lCopy[i], 0, n);
    }
    System.arraycopy(d, 0, dCopy, 0, n);
}

public void restoreFrom(final double[][] lCopy, final double[] dCopy) {
    if (lCopy.length != n || dCopy.length != n) {
        throw new IllegalArgumentException("Dimension mismatch");
    }
    for (int i = 0; i < n; i++) {
        System.arraycopy(lCopy[i], 0, lData[i], 0, n);
    }
    System.arraycopy(dCopy, 0, d, 0, n);
}
}