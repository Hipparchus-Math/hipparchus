package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/**
 * Updates a QR factorization when adding or removing constraints in active-set methods.
 */
public class QRUpdaterR {

    
    /** Internal raw data for the J matrix (stored as J = L^T). */
    private final double[][] jData;

    /** Internal raw data for the upper triangular R matrix. */
    private final double[][] rData;

    /** Number of active constraints. */
    private int iq;

    /** Norm parameter of R, used to detect degeneracy. */
    private double RNorm = 1.0;

    /** Dimension of the optimization problem. */
    private final int n;

    /** Workspace buffer for the constraint vector being added. */
    private final double[] tempDBuffer;

    /** Workspace for rotation cosines to allow precise rollback. */
    private final double[] cList;

    /** Workspace for rotation sines to allow precise rollback. */
    private final double[] sList;

    /** Wrapper for the J matrix (no-copy reference). */
    private final RealMatrix J;

    /** Wrapper for the R matrix (no-copy reference). */
    private final RealMatrix R;

    /**
     * Constructs a new QRUpdater given the lower triangular matrix L.
     * <p>
     * J is initialized as L^T and R is initialized as a zero matrix.
     * Raw data is extracted once to ensure all subsequent operations use direct array access.
     * </p>
     * @param L lower triangular matrix from Cholesky decomposition.
     */
    public QRUpdaterR(final RealMatrix L) {
        this.n = L.getRowDimension();
        this.iq = 0;

        // One-time memory allocation for the solver lifecycle
        this.jData = new double[n][n];
        this.rData = new double[n][n];
        this.tempDBuffer = new double[n];
        this.cList = new double[n];
        this.sList = new double[n];

        // Direct manual transpose for initialization
        final double[][] lRaw = (L instanceof Array2DRowRealMatrix) ? 
                                ((Array2DRowRealMatrix) L).getDataRef() : L.getData();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                this.jData[j][i] = lRaw[i][j];
            }
        }

        // Create wrappers with the false flag to prevent internal copying
        this.J = new Array2DRowRealMatrix(jData, false);
        this.R = new Array2DRowRealMatrix(rData, false);
    }

    /**
     * Adds a constraint vector and updates the QR factorization via a backward Givens cascade.
     * @param d constraint vector to add.
     * @return {@code true} if added successfully; {@code false} if the constraint 
     * is linearly dependent (degenerate) and the update was rolled back.
     */
    public boolean addConstraint(final RealVector d) {
        // Extract raw data once to fill the pre-allocated buffer
        final double[] dRaw = (d instanceof ArrayRealVector) ? 
                              ((ArrayRealVector) d).getDataRef() : d.toArray();
        System.arraycopy(dRaw, 0, tempDBuffer, 0, n);

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

            // Store rotation parameters for possible rollback
            cList[j] = cc;
            sList[j] = ss;

            final double xny = ss / (1.0 + cc);
            for (int k = 0; k < n; k++) {
                final double t1 = jData[k][j - 1];
                final double t2 = jData[k][j];
                final double newT1 = t1 * cc + t2 * ss;
                jData[k][j - 1] = newT1;
                jData[k][j] = xny * (t1 + newT1) - t2;
            }
        }

        // Degeneracy check: if pivot is near zero, roll back J to its previous state
        if (FastMath.abs(tempDBuffer[iq]) < Precision.EPSILON) {
            for (int j = iq + 1; j <= n - 1; j++) {
                double cc = cList[j]; // Invert cosine sign
                double ss = sList[j]; // Invert sine sign
                if (cc == 1.0 && ss == 0.0) continue;

                final double xny = ss / (1.0 + cc);
                for (int k = 0; k < n; k++) {
                    final double t1 = jData[k][j - 1];
                    final double t2 = jData[k][j];
                    final double newT1 = t1 * cc + t2 * ss;
                    jData[k][j - 1] = newT1;
                    jData[k][j] = xny * (t1 + newT1) - t2;
                }
            }
            return false;
        }

        // Store the new transformed column into the R matrix
        for (int i = 0; i <= iq; i++) {
            rData[i][iq] = tempDBuffer[i];
        }
        iq++;
        return true;
    }

    /**
     * Deletes the active constraint at the specified index and restores 
     * upper triangular form via forward Givens rotations.
     * @param index index of the constraint to remove.
     */
    public void deleteConstraint(final int index) {
        if (index < 0 || index >= iq) return;

        // Shift columns of R left to close the gap
        for (int j = index; j < iq - 1; j++) {
            for (int i = 0; i < n; i++) {
                rData[i][j] = rData[i][j + 1];
            }
        }
        // Wipe the last active column
        for (int i = 0; i < n; i++) rData[i][iq - 1] = 0.0;

        iq--;
        if (iq == 0) return;

        // Restore triangular form
        for (int j = index; j < iq; j++) {
            double cc = rData[j][j];
            double ss = rData[j + 1][j];
            double h = FastMath.hypot(cc, ss);
            if (h < Precision.EPSILON) continue;

            rData[j][j] = h;
            rData[j + 1][j] = 0.0;
            cc /= h; ss /= h;
            if (cc < 0.0) { cc = -cc; ss = -ss; rData[j][j] = -h; }

            final double xny = ss / (1.0 + cc);
            // Apply rotations to R
            for (int k = j + 1; k < iq; k++) {
                double t1 = rData[j][k];
                double t2 = rData[j + 1][k];
                rData[j][k] = t1 * cc + t2 * ss;
                rData[j + 1][k] = xny * (t1 + rData[j][k]) - t2;
            }
            // Apply rotations to J
            for (int k = 0; k < n; k++) {
                double t1 = jData[k][j];
                double t2 = jData[k][j + 1];
                jData[k][j] = t1 * cc + t2 * ss;
                jData[k][j + 1] = xny * (t1 + jData[k][j]) - t2;
            }
        }
       
    }

     /**
     * Computes the vector z used in the primal step without explicitly
     * forming the J2 submatrix or the tail subvector of d.
     *
     * This replaces the sequence:
     *
     *   d = J^T * ai
     *   J2 = J[:, iq:n-1]
     *   z = J2 * d[iq:n-1]
     *
     * with a direct multiplication using the inactive columns of J.
     *
     * @param d vector d = J^T * ai
     * @return z vector of length n
     */
    public RealVector computeZ(final RealVector d) {
        final double[] zRaw = new double[n];
        for (int col = iq; col < n; col++) {
            final double dj = d.getEntry(col); 
            if (dj == 0.0) continue;
            for (int row = 0; row < n; row++) {
                zRaw[row] += jData[row][col] * dj;
            }
        }
        return new ArrayRealVector(zRaw, false);
    }

    /**
     * Solves R x = rhs using backward substitution on the active upper
     * triangular factor.
     *
     * <p>Only the first iq entries of d are used; remaining entries
     * are ignored.</p>
     *
     * @param d full vector whose first iq entries define the right-hand side
     * @return solution vector of dimension iq, or null if R is singular
     */
    public RealVector solveR(final RealVector d) {
        if (iq == 0) return new ArrayRealVector(0);
        final double[] x = new double[iq];
        for (int i = 0; i < iq; i++) x[i] = d.getEntry(i);

        for (int i = iq - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < iq; j++) sum += rData[i][j] * x[j];
            double rii = rData[i][i];
            if (FastMath.abs(rii) < Precision.SAFE_MIN) return null;
            x[i] = (x[i] - sum) / rii;
        }
        return new ArrayRealVector(x, false);
    }

    /**
     * Solves R^T*x = rhs using forward substitution.
     * @param rhs right-hand side vector.
     * @return solution vector of dimension iq.
     */
    public RealVector solveRT(final RealVector rhs) {
        if (iq == 0) return new ArrayRealVector(0);
        final double[] x = new double[iq];
        for (int i = 0; i < iq; i++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) sum += rData[j][i] * x[j];
            double rii = rData[i][i];
            if (FastMath.abs(rii) < Precision.SAFE_MIN) return null;
            x[i] = (rhs.getEntry(i) - sum) / rii;
        }
        return new ArrayRealVector(x, false);
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
                out[row] += jData[row][col] * c;
            }
        }
        return new ArrayRealVector(out, false);
    }

    /** @return the current active R factor. */
    public RealMatrix getR() { 
        return (iq == n) ? R : (iq > 0 ? R.getSubMatrix(0, iq - 1, 0, iq - 1) : null); 
    }

    /** @return the full working matrix J. */
    public RealMatrix getJ() { return J; }

    /** @return current number of active constraints. */
    public int getIq() { return iq; }

    /**
     * Computes d = J^T * a using Hipparchus preMultiply.
     *
     * @param a input vector
     * @return d = J^T * a
     */
    public RealVector computeD(final RealVector a) { return J.preMultiply(a); }
}