package org.hipparchus.optim.nonlinear.vector.constrained;
import org.hipparchus.util.FastMath;
import org.hipparchus.linear.RealMatrix;

public class SchnabelEskowSE99 {
    private static final double EPS = 2.220446049250313E-16;

    public static boolean getSafeHessian(RealMatrix matrixH) {
        int n = matrixH.getRowDimension();
        // Pack: lower triangular, column-wise
        double[] sa = new double[n * (n + 1) / 2];
        int k = 0;
        for (int j = 0; j < n; j++) {
            for (int i = j; i < n; i++) sa[k++] = matrixH.getEntry(i, j);
        }

        double[] sd = new double[n];
        int[] perm = new int[n];
        for (int i = 0; i < n; i++) perm[i] = i;
        double[] modified = new double[n];

        // GMW-I: pivot_method=1, is_type1=true, nondecreasing=false, is_2phase=true, relax=0.75
        int info = mchol_gmw(n, sa, sd, perm, modified, EPS, 1, true, false, true, 0.75);

        if (info == 0) {
            boolean changed = false;
            for (int i = 0; i < n; i++) {
                if (modified[i] > 1e-16) {
                    matrixH.setEntry(i, i, matrixH.getEntry(i, i) + modified[i]);
                    changed = true;
                }
            }
            return changed;
        }
        return false;
    }

    private static int mchol_gmw(int nrc, double[] sa, double[] sd, int[] perm, double[] modified,
                                 double delta, int pivot_method, boolean is_type1, 
                                 boolean nondec, boolean is_2phase, double relax) {
        
        int steps_p1 = 0;
        if (is_2phase) {
            steps_p1 = phase_one_factorization(nrc, sa, sd, perm, delta, relax);
            for (int i = 0; i < steps_p1; i++) modified[i] = 0.0;
        }

        // Calcolo offset per iniziare la Phase II dopo steps_p1
        int sa_ptr = 0;
        for (int i = 0; i < steps_p1; i++) sa_ptr += (nrc - i);

        // Calcolo Beta2 (Test O(n^2))
        double gamma = 0.0;
        int temp_ptr = sa_ptr;
        for (int j = steps_p1; j < nrc; j++) {
            temp_ptr++; // skip diag
            for (int i = j + 1; i < nrc; i++) gamma = FastMath.max(gamma, FastMath.abs(sa[temp_ptr++]));
        }
        
        double beta2 = EPS;
        int n2 = nrc - steps_p1;
        if (n2 > 1) {
            double div = is_type1 ? (n2 * n2 - 1.0) : (n2 * n2 - n2);
            beta2 = FastMath.max(EPS, gamma / FastMath.sqrt(div));
        } else {
            beta2 = FastMath.max(EPS, gamma);
        }
        beta2 *= beta2;

        // Phase II
        for (int i = steps_p1; i < nrc; i++) {
            int idx_rel = (pivot_method == 1) ? get_idx_max_diag(nrc - i, sa, sa_ptr) : 0;
            if (idx_rel != 0) {
                interchange_row_column(nrc, sa, i, i + idx_rel);
                int tmp = perm[i]; perm[i] = perm[i + idx_rel]; perm[i + idx_rel] = tmp;
            }

            double d = sa[sa_ptr];
            if (is_type1 && d < 0) d = -d;
            d = FastMath.max(d, delta);

            double row_max_sq = 0.0;
            for (int j = 1; j < nrc - i; j++) {
                double val = sa[sa_ptr + j];
                row_max_sq = FastMath.max(row_max_sq, val * val);
            }
            d = FastMath.max(d, row_max_sq / beta2);

            modified[i] = d - sa[sa_ptr];
            if (nondec && i > steps_p1 && modified[i] < modified[i-1]) {
                d += (modified[i-1] - modified[i]);
                modified[i] = modified[i-1];
            }
            sd[i] = d;

            // LDL^T Update (Schur Complement)
            int next_col = sa_ptr + (nrc - i);
            for (int j = i + 1; j < nrc; j++) {
                double l_ji = sa[sa_ptr + (j - i)] / d;
                int target = next_col;
                for (int k = j; k < nrc; k++) sa[target++] -= l_ji * sa[sa_ptr + (k - i)];
                sa[sa_ptr + (j - i)] = l_ji;
                next_col += (nrc - j);
            }
            sa_ptr += (nrc - i);
        }

        // Ripristina ordine originale dei modificatori
        double[] final_mod = new double[nrc];
        for (int i = 0; i < nrc; i++) final_mod[perm[i]] = modified[i];
        System.arraycopy(final_mod, 0, modified, 0, nrc);
        return 0;
    }

    private static int phase_one_factorization(int nrc, double[] sa, double[] sd, int[] perm, double delta, double relax) {
        double max_diag = 0;
        for (int i = 0, p = 0; i < nrc; i++) {
            max_diag = FastMath.max(max_diag, FastMath.abs(sa[p]));
            p += (nrc - i);
        }
        double threshold = (relax > 0) ? -relax * max_diag : delta;
        
        int sa0 = 0;
        for (int i = 0; i < nrc; i++) {
            int idx = get_idx_max_diag(nrc - i, sa, sa0);
            if (idx != 0) {
                interchange_row_column(nrc, sa, i, i + idx);
                int tmp = perm[i]; perm[i] = perm[i + idx]; perm[i + idx] = tmp;
            }
            
            if (sa[sa0] < delta) return i;

            int next_diag = sa0 + (nrc - i);
            for (int j = i + 1; j < nrc; j++) {
                double off = sa[sa0 + (j - i)];
                if (sa[next_diag] - (off * off / sa[sa0]) < threshold) return i;
                next_diag += (nrc - j);
            }

            double d = sa[sa0];
            sd[i] = d;
            int update_ptr = sa0 + (nrc - i);
            for (int j = i + 1; j < nrc; j++) {
                double l_ji = sa[sa0 + (j - i)] / d;
                int target = update_ptr;
                for (int k = j; k < nrc; k++) sa[target++] -= l_ji * sa[sa0 + (k - i)];
                sa[sa0 + (j - i)] = l_ji;
                update_ptr += (nrc - j);
            }
            sa[sa0] = 1.0;
            sa0 += (nrc - i);
        }
        return nrc;
    }

    private static void interchange_row_column(int n, double[] sa, int i, int j) {
        int r1 = Math.min(i, j), r2 = Math.max(i, j);
        if (r1 == r2) return;
        int p1 = r1, p2 = r2, off = n - 1;
        for (int k = 0; k < r1; k++) {
            double t = sa[p1]; sa[p1] = sa[p2]; sa[p2] = t;
            p1 += off; p2 += off; off--;
        }
        int d1 = p1; p1++; p2 += off; off--;
        for (int k = r1 + 1; k < r2; k++) {
            double t = sa[p1]; sa[p1] = sa[p2]; sa[p2] = t;
            p1++; p2 += off; off--;
        }
        p1++; int d2 = p2; p2++;
        for (int k = r2 + 1; k < n; k++) {
            double t = sa[p1]; sa[p1] = sa[p2]; sa[p2] = t;
            p1++; p2++;
        }
        double t = sa[d1]; sa[d1] = sa[d2]; sa[d2] = t;
    }

    private static int get_idx_max_diag(int n, double[] sa, int off) {
        int best = 0; double maxv = sa[off]; int cur = off;
        for (int i = 0; i < n; i++) {
            if (sa[cur] > maxv) { maxv = sa[cur]; best = i; }
            cur += (n - i);
        }
        return best;
    }
}