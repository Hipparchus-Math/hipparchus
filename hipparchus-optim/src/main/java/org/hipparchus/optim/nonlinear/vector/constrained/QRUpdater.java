package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/**
 * Updates a QR factorization when adding or removing constraints in active-set methods.
 * <p>
 * This implementation leverages Hipparchus {@link RealMatrix} abstraction. For large
 * dimensions, {@link MatrixUtils#createRealMatrix(int, int)} automatically generates 
 * block-partitioned matrices to ensure high cache locality and JIT optimization.
 * </p>
 */
public class QRUpdater {

    /** The full working matrix J (representing the basis of the subspace). */
    private final RealMatrix J;

    /** The upper triangular factor R of the active constraints. */
    private final RealMatrix R;

    /** Number of currently active constraints. */
    private int iq;

    /** Dimension of the optimization problem. */
    private final int n;

    /** Workspace buffer for the constraint vector being added. */
    private final double[] tempDBuffer;

    /** Workspace for rotation cosines. */
    private final double[] cList;

    /** Workspace for rotation sines. */
    private final double[] sList;

    /**
     * Constructs a new QRUpdater given the lower triangular matrix L.
     * <p>
     * J is initialized as L^T using matrix-level operations to preserve
     * block-storage benefits if n is large.
     * </p>
     * @param L lower triangular matrix from Cholesky decomposition.
     */
    public QRUpdater(final RealMatrix L) {
        this.n = L.getRowDimension();
        this.iq = 0;

        // Initialize J and R through MatrixUtils to exploit BlockRealMatrix for large n
        this.J = L.transpose(); 
        this.R = MatrixUtils.createRealMatrix(n, n);

        // Pre-allocated workspace for Givens rotations and rollback
        this.tempDBuffer = new double[n];
        this.cList = new double[n];
        this.sList = new double[n];
    }

    /**
     * Adds a constraint vector and updates the QR factorization via a backward Givens cascade.
     * @param d constraint vector to add.
     * @return {@code true} if added successfully; {@code false} if degenerate.
     */
    public boolean addConstraint(final RealVector d) {
        // Load vector into workspace
        for (int i = 0; i < n; i++) {
            tempDBuffer[i] = d.getEntry(i);
        }

        // Backward Givens cascade to eliminate elements in tempDBuffer
        for (int j = n - 1; j >= iq + 1; j--) {
            double cc = tempDBuffer[j - 1];
            double ss = tempDBuffer[j];
            double h = FastMath.hypot(cc, ss);

            if (h < Precision.EPSILON) {
                cList[j] = 1.0; 
                sList[j] = 0.0;
                continue;
            }

            tempDBuffer[j] = 0.0;
            ss /= h;
            cc /= h;

            if (cc < 0.0) {
                cc = -cc; ss = -ss;
                tempDBuffer[j - 1] = -h;
            } else {
                tempDBuffer[j - 1] = h;
            }

            cList[j] = cc;
            sList[j] = ss;

            // Apply rotation to J. Using getEntry/setEntry on BlockRealMatrix is 
            // handled efficiently by tiles during JIT execution.
            final double xny = ss / (1.0 + cc);
            for (int k = 0; k < n; k++) {
                final double t1 = J.getEntry(k, j - 1);
                final double t2 = J.getEntry(k, j);
                final double newT1 = t1 * cc + t2 * ss;
                J.setEntry(k, j - 1, newT1);
                J.setEntry(k, j, xny * (t1 + newT1) - t2);
            }
        }

        // Degeneracy check
        if (FastMath.abs(tempDBuffer[iq]) < Precision.EPSILON) {
            // Rollback rotations in reverse order to restore J
            for (int j = iq + 1; j <= n - 1; j++) {
                double cc = cList[j];
                double ss = sList[j];
                if (cc == 1.0 && ss == 0.0) continue;

                final double xny = ss / (1.0 + cc);
                for (int k = 0; k < n; k++) {
                    final double t1 = J.getEntry(k, j - 1);
                    final double t2 = J.getEntry(k, j);
                    final double newT1 = t1 * cc + t2 * ss;
                    J.setEntry(k, j - 1, newT1);
                    J.setEntry(k, j, xny * (t1 + newT1) - t2);
                }
            }
            return false;
        }

        // Store new column in R factor
        for (int i = 0; i <= iq; i++) {
            R.setEntry(i, iq, tempDBuffer[i]);
        }
        iq++;
        return true;
    }

    /**
     * Deletes the active constraint at the specified index.
     * @param index index of the constraint to remove.
     */
    public void deleteConstraint(final int index) {
        if (index < 0 || index >= iq) return;

        // Shift columns of R left to close the gap
        for (int j = index; j < iq - 1; j++) {
            for (int i = 0; i < n; i++) {
                R.setEntry(i, j, R.getEntry(i, j + 1));
            }
        }
        // Nullify last column
        for (int i = 0; i < n; i++) R.setEntry(i, iq - 1, 0.0);

        iq--;
        if (iq == 0) return;

        // Restore triangular form via forward Givens rotations
        for (int j = index; j < iq; j++) {
            double cc = R.getEntry(j, j);
            double ss = R.getEntry(j + 1, j);
            double h = FastMath.hypot(cc, ss);
            if (h < Precision.EPSILON) continue;

            R.setEntry(j, j, h);
            R.setEntry(j + 1, j, 0.0);
            cc /= h; ss /= h;
            if (cc < 0.0) { cc = -cc; ss = -ss; R.setEntry(j, j, -h); }

            final double xny = ss / (1.0 + cc);
            // Apply rotations to R (active part)
            for (int k = j + 1; k < iq; k++) {
                double t1 = R.getEntry(j, k);
                double t2 = R.getEntry(j + 1, k);
                double nextT1 = t1 * cc + t2 * ss;
                R.setEntry(j, k, nextT1);
                R.setEntry(j + 1, k, xny * (t1 + nextT1) - t2);
            }
            // Apply rotations to J (full basis)
            for (int k = 0; k < n; k++) {
                double t1 = J.getEntry(k, j);
                double t2 = J.getEntry(k, j + 1);
                double nextT1 = t1 * cc + t2 * ss;
                J.setEntry(k, j, nextT1);
                J.setEntry(k, j + 1, xny * (t1 + nextT1) - t2);
            }
        }
    }

    /**
     * Computes z = J_inactive * d_inactive.
     * @param d vector containing inactive components from index iq.
     * @return resulting vector z of dimension n.
     */
    public RealVector computeZ(final RealVector d) {
        final double[] zRaw = new double[n];
        for (int col = iq; col < n; col++) {
            final double val = d.getEntry(col); 
            if (val == 0.0) continue;
            for (int row = 0; row < n; row++) {
                zRaw[row] += J.getEntry(row, col) * val;
            }
        }
        return new ArrayRealVector(zRaw, false);
    }

    /**
     * Solves R x = rhs using backward substitution.
     * @param d right-hand side vector.
     * @return solution vector of dimension iq.
     */
    public RealVector solveR(final RealVector d) {
        if (iq == 0) return new ArrayRealVector(0);
        final double[] xResult = new double[iq];
        for (int i = 0; i < iq; i++) xResult[i] = d.getEntry(i);

        for (int i = iq - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < iq; j++) sum += R.getEntry(i, j) * xResult[j];
            double rii = R.getEntry(i, i);
            if (FastMath.abs(rii) < Precision.EPSILON) return null;
            
            xResult[i] = (xResult[i] - sum) / rii;
        }
        return new ArrayRealVector(xResult, false);
    }

    /**
     * Solves R^T * x = rhs using forward substitution.
     * @param rhs right-hand side vector.
     * @return solution vector of dimension iq.
     */
    public RealVector solveRT(final RealVector rhs) {
        if (iq == 0) return new ArrayRealVector(0);
        final double[] xResult = new double[iq];
        for (int i = 0; i < iq; i++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) sum += R.getEntry(j, i) * xResult[j];
            double rii = R.getEntry(i, i);
            if (FastMath.abs(rii) < Precision.EPSILON) return null;
            xResult[i] = (rhs.getEntry(i) - sum) / rii;
        }
        return new ArrayRealVector(xResult, false);
    }

    /**
     * Applies the active part of J (J1 * coeffs).
     * @param coeffs multipliers for active constraints.
     * @return resulting vector of dimension n.
     */
    public RealVector applyJActive(final RealVector coeffs) {
        final double[] out = new double[n];
        for (int col = 0; col < iq; col++) {
            final double c = coeffs.getEntry(col);
            if (c == 0.0) continue;
            for (int row = 0; row < n; row++) {
                out[row] += J.getEntry(row, col) * c;
            }
        }
        return new ArrayRealVector(out, false);
    }

    /** @return the current active R factor (submatrix if iq < n). */
    public RealMatrix getR() { 
        return (iq == n) ? R : (iq > 0 ? R.getSubMatrix(0, iq - 1, 0, iq - 1) : null); 
    }

    /** @return the full working matrix J. */
    public RealMatrix getJ() { return J; }

    /** @return current number of active constraints. */
    public int getIq() { return iq; }

    /** Computes d = J^T * a. */
    public RealVector computeD(final RealVector a) { return J.preMultiply(a); }
}