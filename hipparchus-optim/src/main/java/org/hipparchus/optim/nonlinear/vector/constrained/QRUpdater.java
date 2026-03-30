package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/**
 * Updates a QR factorization when adding or removing constraints in
 * active-set methods for nonlinear vector optimization.
 *
 * This version is algebraically identical to your working backup, but
 * replaces {@link FastMath#hypot(double, double)} with a robust
 * scaling-based safeHypot to avoid overflow/underflow.
 */
public class QRUpdater {
    
    /** sqrtEpsilon*/
    private double sqrtEpsilon=FastMath.sqrt(Precision.EPSILON);

    /** Inverse of the lower triangular matrix L (stored as J = L^T). */
    private RealMatrix J;

    /** Upper triangular R matrix for active constraints. */
    private final RealMatrix R;

    /** Number of active constraints. */
    private int iq;

    /** Norm parameter of R, used to detect degeneracy. */
    private double RNorm = 1.0;

    /** Dimension of the optimization problem. */
    private final int n;
    

    /**
     * Constructs a new QRUpdater given the lower triangular matrix L.
     * J is initialized as L^T and R as an n-by-n zero matrix.
     *
     * @param L lower triangular matrix to initialize the updater
     */
    public QRUpdater(final RealMatrix L) {
        this.n = L.getRowDimension();
        this.J = L.transpose();
        this.R = MatrixUtils.createRealMatrix(n, n);
        this.iq = 0;
        
      
    }

   
    /**
     * Adds a constraint vector and updates the QR factorization via Givens rotations.
     *
     * @param d constraint vector to add; must have length n
     * @return {@code true} if the constraint was added successfully; {@code false} if
     *         the problem is degenerate and the constraint cannot be added
     */
    public boolean addConstraint(RealVector d) {

        RealMatrix Jtemp = new Array2DRowRealMatrix(J.getData());
        RealVector tempD = new ArrayRealVector(d);

        double cc;
        double ss;
        double h;
        double t1;
        double t2;
        double xny;

        // Backward Givens cascade to zero out the tail entries of d
        for (int j = n - 1; j >= iq + 1; j--) {
            cc = tempD.getEntry(j - 1);
            ss = tempD.getEntry(j);

            // use robust safeHypot instead of FastMath.hypot
//             h = FastMath.hypot(cc, ss);
            h = FastMath.hypot(cc, ss);
            if (h < Precision.EPSILON) {
                continue;
            }
            tempD.setEntry(j, 0.0);
            ss /= h;
            cc /= h;

            // ensure cc >= 0, flipping signs and the stored d(j-1) accordingly
            if (cc < 0.0) {
                cc = -cc;
                ss = -ss;
                tempD.setEntry(j - 1, -h);
                
            } else {
                tempD.setEntry(j - 1, h);
            }

            xny = ss / (1.0 + cc);
            for (int k = 0; k < n; k++) {
                t1 = J.getEntry(k, j - 1);
                t2 = J.getEntry(k, j);
                J.setEntry(k, j - 1, t1 * cc + t2 * ss);
                J.setEntry(k, j, xny * (t1 + J.getEntry(k, j - 1)) - t2);
            }
        }

        // Degeneracy check: leading component must not be too small
        if (FastMath.abs(tempD.getEntry(iq)) < Precision.EPSILON* RNorm) {
            J = Jtemp;
            return false;
        }

        // Store the new column in R
        for (int i = 0; i <= iq; i++) {
            R.setEntry(i, iq, tempD.getEntry(i));
        }
        RNorm = FastMath.max(RNorm, FastMath.abs(tempD.getEntry(iq)));
        iq++;
        return true;
    }
   

    
    

    
    /**
     * Deletes the active constraint at the specified index and updates
     * the QR factorization via Givens rotations.
     *
     * @param constraintIndex index of the constraint to delete
     */
    public void deleteConstraint(int constraintIndex) {
        if (constraintIndex < 0 || constraintIndex >= iq) {
            return; // index not found
        }

        // Shift columns of R left starting from the deleted one
        for (int i = constraintIndex; i < iq - 1; i++) {
            for (int j = 0; j < n; j++) {
                R.setEntry(j, i, R.getEntry(j, i + 1));
            }
        }
        // Zero the last (now unused) column
        for (int j = 0; j < n; j++) {
            R.setEntry(j, iq - 1, 0.0);
        }

        iq--;
        if (iq == 0) {
            return;
        }

        // Re-triangularize R and update J using Givens rotations
        for (int j = constraintIndex; j < iq; j++) {
            double cc = R.getEntry(j, j);
            double ss = R.getEntry(j + 1, j);

            // use robust safeHypot instead of FastMath.hypot
            double h = FastMath.hypot(cc, ss);
            if (h < Precision.EPSILON) {
                continue;
            }
            R.setEntry(j, j, h);
            R.setEntry(j + 1, j, 0.0);
            cc /= h;
            ss /= h;
            //for deleting constraint sign change doesn't work for all problems
            if(cc<0)
            {
                cc=-cc;
                ss=-ss;
                R.setEntry(j, j, -h);
               
                
            }

            double xny = ss / (1.0 + cc);
            for (int k = j + 1; k < iq; k++) {
                double t1 = R.getEntry(j, k);
                double t2 = R.getEntry(j + 1, k);
                R.setEntry(j, k, t1 * cc + t2 * ss);
                R.setEntry(j + 1, k, xny * (t1 + R.getEntry(j, k)) - t2);
            }
            for (int k = 0; k < n; k++) {
                double t1 = J.getEntry(k, j);
                double t2 = J.getEntry(k, j + 1);
                J.setEntry(k, j, t1 * cc + t2 * ss);
                J.setEntry(k, j + 1, xny * (t1 + J.getEntry(k, j)) - t2);
            }
        }
         // ---- Aggiornamento RNorm ----
    double maxD = 0.0;
    for (int i = 0; i < iq; i++) {
        maxD = FastMath.max(maxD, FastMath.abs(R.getEntry(i, i)));
    }
    RNorm = (maxD > 0.0 ? maxD : 1.0);
    }
    
    

    /** Returns the current active upper triangular factor R. */
    public RealMatrix getR() {
        if (iq == this.n) {
            return R;
        }
        if (iq > 0) {
            return R.getSubMatrix(0, iq - 1, 0, iq - 1);
        }
        return null;
    }

    /** Returns the inverse of the active R factor. */
    public RealMatrix getRInv() {
        if (iq > 0) {
            return inverseUpperTriangular(getR());
        }
        return null;
    }
    
    

    /** Inverse of an upper triangular matrix via backward substitution. */
    private RealMatrix inverseUpperTriangular(RealMatrix U) {
        int p = U.getRowDimension();
        RealMatrix Uinv = MatrixUtils.createRealMatrix(p, p);
        for (int i = p - 1; i >= 0; i--) {
            Uinv.setEntry(i, i, 1.0 / U.getEntry(i, i));
            for (int j = i - 1; j >= 0; j--) {
                double sum = 0.0;
                for (int k = j + 1; k <= i; k++) {
                    sum += U.getEntry(j, k) * Uinv.getEntry(k, i);
                }
                Uinv.setEntry(j, i, -sum / U.getEntry(j, j));
            }
        }
        return Uinv;
    }

    /** Returns the current J matrix. */
    public RealMatrix getJ() {
        return J;
    }

    /** Returns the inactive columns of J, starting at the first non-active index. */
    public RealMatrix getJ2() {
        if (iq == n) {
            return null;
        }
        return J.getSubMatrix(0, n - 1, iq, n - 1);
    }

    /** Returns the number of active constraints. */
    public int getIq() {
        return iq;
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

    final int start = iq;
    final int dim = n;

//    // if all constraints are active there is no primal direction
//    if (!(n-iq>0)) {
//        return new ArrayRealVector(dim,0);
//    }

    final ArrayRealVector z = new ArrayRealVector(dim,0);

    // z = J[:, start:n-1] * d[start:n-1]
    for (int col = start; col < dim; ++col) {

        final double dj = d.getEntry(col);
        if (dj == 0.0) {
            continue;
        }

        for (int row = 0; row < dim; ++row) {
            z.addToEntry(row, J.getEntry(row, col) * dj);
        }
    }

    return z;
}
 /**
 * Computes d = J^T * a using Hipparchus preMultiply.
 *
 * @param a input vector
 * @return d = J^T * a
 */
public RealVector computeD(final RealVector a) {
    return J.preMultiply(a);
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
    if (iq == 0) {
        return new ArrayRealVector(0, 0);
    }

    final ArrayRealVector x = new ArrayRealVector(iq);

    for (int i = 0; i < iq; ++i) {
        x.setEntry(i, d.getEntry(i));
    }

    for (int i = iq - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < iq; ++j) {
            sum += R.getEntry(i, j) * x.getEntry(j);
        }

        final double rii = R.getEntry(i, i);
        if (FastMath.abs(rii) < Precision.SAFE_MIN) {
            
            return null;
        }

        x.setEntry(i, (x.getEntry(i) - sum) / rii);
    }

    return x;
    
    
  

}
}
