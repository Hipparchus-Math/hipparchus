/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 * https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/**
 * Updates a QR factorization when adding or removing constraints in active-set methods.
 *
 * @see "J. W. Daniel, W. B. Gragg, L. Kaufman and G. W. Stewart,
 *      Reorthogonalization and Stable Algorithms for Updating the Gram-Schmidt QR Factorization"
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

    // Load the constraint vector into the reusable workspace.
    for (int i = 0; i < n; i++) {
        tempDBuffer[i] = d.getEntry(i);
    }

    /*
     * Apply a backward cascade of orthogonal reflectors to eliminate
     * all components beyond the current active-set position.
     */
    for (int j = n - 1; j >= iq + 1; j--) {

        double cc = tempDBuffer[j - 1];
        double ss = tempDBuffer[j];

        /*
         * The component to eliminate is already exactly zero.
         * Store (0, 0) as an unambiguous marker indicating that no
         * transformation was applied.
         */
        if (ss == 0.0) {
            cList[j] = 0.0;
            sList[j] = 0.0;
            continue;
        }

        /*
         * FastMath.hypot computes the Euclidean norm without the
         * avoidable intermediate overflow/underflow of cc*cc + ss*ss.
         *
         * Give h the sign of cc. Consequently, after normalization,
         * cc is non-negative and 1 + cc cannot suffer cancellation.
         *
         * The cc == 0 case needs no special treatment:
         * it naturally produces cc = 0 and ss = +/-1.
         */
        double h = FastMath.hypot(cc, ss);

        if (cc < 0.0) {
            h = -h;
        }

        cc /= h;
        ss /= h;

        tempDBuffer[j - 1] = h;
        tempDBuffer[j] = 0.0;

        cList[j] = cc;
        sList[j] = ss;

        final double xny = ss / (1.0 + cc);

        // Apply the same reflector to columns j - 1 and j of J.
        for (int k = 0; k < n; k++) {
            final double t1 = J.getEntry(k, j - 1);
            final double t2 = J.getEntry(k, j);

            final double newT1 = t1 * cc + t2 * ss;
            final double newT2 = xny * (t1 + newT1) - t2;

            J.setEntry(k, j - 1, newT1);
            J.setEntry(k, j, newT2);
        }
    }

    /*
     * Preserve the existing degeneracy test.
     * This test concerns acceptance of the new constraint, not whether
     * an individual orthogonal transformation must be applied.
     */
    if (FastMath.abs(tempDBuffer[iq]) < Precision.EPSILON) {

        /*
         * Roll back the applied transformations in reverse order.
         * The reflector is symmetric and orthogonal, hence self-inverse.
         */
        for (int j = iq + 1; j < n; j++) {

            final double cc = cList[j];
            final double ss = sList[j];

            // No transformation was applied at this position.
            if (cc == 0.0 && ss == 0.0) {
                continue;
            }

            final double xny = ss / (1.0 + cc);

            for (int k = 0; k < n; k++) {
                final double t1 = J.getEntry(k, j - 1);
                final double t2 = J.getEntry(k, j);

                final double newT1 = t1 * cc + t2 * ss;
                final double newT2 = xny * (t1 + newT1) - t2;

                J.setEntry(k, j - 1, newT1);
                J.setEntry(k, j, newT2);
            }
        }

        return false;
    }

    // Append the new column to the active triangular factor R.
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
//
   public void deleteConstraint(final int index) {

    if (index < 0 || index >= iq) {
        return;
    }

    // Shift the columns following the deleted constraint to the left.
    for (int j = index; j < iq - 1; j++) {
        for (int i = 0; i < n; i++) {
            R.setEntry(i, j, R.getEntry(i, j + 1));
        }
    }

    // Clear the final column, which is no longer active.
    for (int i = 0; i < n; i++) {
        R.setEntry(i, iq - 1, 0.0);
    }

    iq--;

    if (iq == 0) {
        return;
    }

    /*
     * Restore the upper-triangular form of R by eliminating each
     * subdiagonal component created by the column shift.
     */
    for (int j = index; j < iq; j++) {

        double cc = R.getEntry(j, j);
        double ss = R.getEntry(j + 1, j);

        /*
         * The subdiagonal component is already exactly zero.
         * No transformation is necessary.
         */
        if (ss == 0.0) {
            continue;
        }

        /*
         * Give h the sign of cc so that the normalized cc is
         * non-negative. The cc == 0 case is handled naturally.
         */
        double h = FastMath.hypot(cc, ss);

        if (cc < 0.0) {
            h = -h;
        }

        cc /= h;
        ss /= h;

        R.setEntry(j, j, h);
        R.setEntry(j + 1, j, 0.0);

        final double xny = ss / (1.0 + cc);

        // Apply the reflector to the remaining active part of R.
        for (int k = j + 1; k < iq; k++) {
            final double t1 = R.getEntry(j, k);
            final double t2 = R.getEntry(j + 1, k);

            final double newT1 = t1 * cc + t2 * ss;
            final double newT2 = xny * (t1 + newT1) - t2;

            R.setEntry(j, k, newT1);
            R.setEntry(j + 1, k, newT2);
        }

        // Apply the same reflector to columns j and j + 1 of J.
        for (int k = 0; k < n; k++) {
            final double t1 = J.getEntry(k, j);
            final double t2 = J.getEntry(k, j + 1);

            final double newT1 = t1 * cc + t2 * ss;
            final double newT2 = xny * (t1 + newT1) - t2;

            J.setEntry(k, j, newT1);
            J.setEntry(k, j + 1, newT2);
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
            if (val == 0.0) {
                continue;
            }
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
    if (iq == 0) {
        return new ArrayRealVector(0);
    }

    final double[] xResult = new double[iq];

    for (int i = iq - 1; i >= 0; i--) {

        double sum = 0.0;

        for (int j = i + 1; j < iq; j++) {
            sum += R.getEntry(i, j) * xResult[j];
        }

        final double rii = R.getEntry(i, i);

        if (FastMath.abs(rii) < Precision.SAFE_MIN) return null;

        xResult[i] = (d.getEntry(i) - sum) / rii;
    }

    return new ArrayRealVector(xResult, false);
}

    /**
     * Solves R^T * x = rhs using forward substitution.
     * @param rhs right-hand side vector.
     * @return solution vector of dimension iq.
     */
    public RealVector solveRT(final RealVector rhs) {
        if (iq == 0) {
            return new ArrayRealVector(0);
        }
        final double[] xResult = new double[iq];
        for (int i = 0; i < iq; i++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) sum += R.getEntry(j, i) * xResult[j];
            double rii = R.getEntry(i, i);
            if (FastMath.abs(rii) < Precision.SAFE_MIN) {
                return null;
            }
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

    /** Get the current active R factor.
     * @return the current active R factor (submatrix if iq &lt; n)
     */
    public RealMatrix getR() {
        return (iq == n) ? R : (iq > 0 ? R.getSubMatrix(0, iq - 1, 0, iq - 1) : null);
    }

    /** Get the full working matrix J.
     * @return the full working matrix J
     */
    public RealMatrix getJ() {
        return J;
    }

    /** Get current number of active constraints.
     * @return current number of active constraints.
     */
    public int getIq() {
        return iq;
    }

    /** Computes d = J^T * a.
     * @param a a vector
     * @return d = J^T * a
     */
    public RealVector computeD(final RealVector a) {
        return J.preMultiply(a);
    }

}
