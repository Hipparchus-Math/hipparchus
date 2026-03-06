package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/**
 * Modified Cholesky factorization for QP solves (Powell / Overton style).
 *
 * <p>
 * Computes L such that:
 * <pre>
 *     H + E = L Lᵀ
 * </pre>
 * where E is a diagonal correction introduced only if needed.
 *
 * <p>
 * Assumptions:
 * <ul>
 *   <li>H is symmetric by construction (no symmetry checks performed)</li>
 *   <li>No pivoting</li>
 *   <li>Used ONLY inside the QP solver</li>
 * </ul>
 */
public final class ModifiedCholeskyQP {

    private ModifiedCholeskyQP() { }

    /**
     * Factorize H for QP usage: returns L such that H + E = L Lᵀ.
     *
     * @param H symmetric Hessian (not modified)
     * @return lower triangular matrix L
     */
    public static RealMatrix factorize(final RealMatrix H) {

    final int n = H.getRowDimension();
    if (n == 0) {
        return MatrixUtils.createRealMatrix(0, 0);
    }

    final double eps = Precision.EPSILON;

    // gamma = max |diag(H)|, xi = max |offdiag(H)|
    double gamma = 0.0;
    double xi    = 0.0;
    for (int i = 0; i < n; ++i) {
        gamma = FastMath.max(gamma, FastMath.abs(H.getEntry(i, i)));
        for (int j = 0; j < i; ++j) {
            xi = FastMath.max(xi, FastMath.abs(H.getEntry(i, j)));
        }
    }

    final double delta = eps * FastMath.max(gamma + xi, 1.0);
    final double beta  = FastMath.sqrt(FastMath.max(FastMath.max(gamma, xi / (double) n), eps));

    final RealMatrix L = MatrixUtils.createRealMatrix(n, n);
    final double[] d   = new double[n];

    // L = I
    for (int i = 0; i < n; ++i) {
        L.setEntry(i, i, 1.0);
    }

    for (int j = 0; j < n; ++j) {

        // djtemp = H(j,j) - sum_{k<j} L(j,k)^2 * d(k)
        double djtemp = H.getEntry(j, j);
        for (int k = 0; k < j; ++k) {
            final double ljk = L.getEntry(j, k);
            djtemp -= ljk * ljk * d[k];
        }

        double theta = 0.0;

        if (j < n - 1) {

            // theta = max_i |C(i,j)|, i=j+1..n-1
            for (int i = j + 1; i < n; ++i) {
                double cij = H.getEntry(i, j);
                for (int k = 0; k < j; ++k) {
                    cij -= L.getEntry(i, k) * d[k] * L.getEntry(j, k);
                }
                theta = FastMath.max(theta, FastMath.abs(cij));
            }

            final double t = (theta / beta);
            final double bound = t * t;

            // IMPORTANT: use djtemp (not abs(djtemp))
            d[j] = FastMath.max(FastMath.max(djtemp, bound), delta);

            // L(i,j) = C(i,j) / d(j)
            for (int i = j + 1; i < n; ++i) {
                double cij = H.getEntry(i, j);
                for (int k = 0; k < j; ++k) {
                    cij -= L.getEntry(i, k) * d[k] * L.getEntry(j, k);
                }
                L.setEntry(i, j, cij / d[j]);
            }

        } else {
            // last diagonal: d(n) = max(djtemp, delta)
            d[j] = FastMath.max(djtemp, delta);
        }
    }

    // L := L * sqrt(D)
    for (int j = 0; j < n; ++j) {
        final double s = FastMath.sqrt(d[j]);
        for (int i = j; i < n; ++i) {
            L.multiplyEntry(i, j, s);
        }
    }

    return L;
}

}
