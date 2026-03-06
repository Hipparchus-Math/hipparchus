package org.hipparchus.optim.nonlinear.vector.constrained;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

public final class HessianRegularizer {
    private static double D=1.0e-4;
    private HessianRegularizer() { }

    /**
     * Paper-style Gershgorin stabilization (single call).
     *
     * Computes l = min_i (h_ii - sum_{j!=i} |h_ij|).
     * If l < 0, applies diagonal shift m = -2(1 + D) l and performs H <- H + m I in-place.
     *
     * @param H symmetric Hessian (modified in-place)
     * @param D small margin parameter (paper); use 0.0 unless you have a reason
     * @return true if H was modified (shift applied), false otherwise
     */
    public static boolean regularizeInPlacePaper(final RealMatrix H) {
        final int n = H.getRowDimension();
        if (n == 0) {
            return false;
        }

        double l = Double.POSITIVE_INFINITY;

        for (int i = 0; i < n; ++i) {
            final double hii = H.getEntry(i, i);

            double ri = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j == i) {
                    continue;
                }
                ri += FastMath.abs(H.getEntry(i, j));
            }

            final double li = hii - ri;
            if (li < l) {
                l = li;
            }
        }

        if (l >= 0.0) {
            return false; // already PD by Gershgorin bound
        }

        // Paper formula: m = -2(1 + D) l  (l < 0 => m > 0)
        double m =-(1.0 + D) * l;

        // Apply H <- H + m I
        for (int i = 0; i < n; ++i) {
            H.addToEntry(i, i, m);
        }

        return true;
    }
}
