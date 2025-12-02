/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/**
 * Updates a QR factorization when adding or removing constraints in
 * <a href = "https://en.wikipedia.org/wiki/Active_set_method">active set methods</a> for
 * nonlinear vector optimization.
 * <p>
 * Maintains an inverse of the lower triangular matrix (J) and an upper triangular factor (R),
 * applying Givens rotations for efficient rank updates.
 * </p>
 *
 * @see <a href="https://en.wikipedia.org/wiki/QR_decomposition">QR decomposition</a>
 * @since 1.2
 */
public class QRUpdaterBackUp {

    /** Inverse of the lower triangular matrix L. */
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
     * <p>
     * Computes J = L^{-1} and initializes R to an n-by-n zero matrix.
     * </p>
     * @param L lower triangular matrix to initialize the updater
     */
    public QRUpdaterBackUp(final RealMatrix L) {
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
        for (int j = n - 1; j >= iq + 1; j--) {
            cc = tempD.getEntry(j - 1);
            ss = tempD.getEntry(j);
            h = FastMath.hypot(cc, ss);
            if (FastMath.abs(h) < Precision.EPSILON) {
                continue;
            }
            tempD.setEntry(j, 0.0);
            ss /= h;
            cc /= h;
//            if(FastMath.abs(ss)<Precision.EPSILON)continue;
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
        
        

        if (FastMath.abs(tempD.getEntry(iq)) <Precision.EPSILON* (RNorm)) {
            J =Jtemp;
            return false;
        }


        for (int i = 0; i <= iq; i++) {
            R.setEntry(i, iq, tempD.getEntry(i));
        }
        RNorm = FastMath.max(RNorm, FastMath.abs(tempD.getEntry(iq)));
        iq++;
        return true;
    }
    
   

    /**
     * Deletes the active constraint at the specified index and updates the QR factorization via Givens rotations.
     *
     * @param constraintIndex index of the constraint to delete
     */
    public void deleteConstraint(int constraintIndex) {
        if (constraintIndex < 0 || constraintIndex >= iq) {
            
            return; //index not found
        }
        for (int i = constraintIndex; i < iq - 1; i++) {
            for (int j = 0; j < n; j++) {
                R.setEntry(j, i, R.getEntry(j, i + 1));
            }
        }
        for (int j = 0; j < n; j++) {
            R.setEntry(j, iq - 1, 0.0);
        }
        iq--;
        if (iq == 0) {
            return;
        }
        for (int j = constraintIndex; j < iq; j++) {
            double cc = R.getEntry(j, j);
            double ss = R.getEntry(j + 1, j);
            double h = FastMath.hypot(cc, ss);
            if (FastMath.abs(h) < Precision.EPSILON) {
                continue;
            }
            R.setEntry(j, j, h);
            R.setEntry(j + 1, j, 0.0);
            cc /= h;
            ss /= h;
           
    
               //////////////// 
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
       double rNorm=0.0;
       for (int i = 0; i < iq; i++)
       {
       double norm=FastMath.abs(R.getEntry(i, i));    
       
       if(norm>rNorm)
           rNorm=norm;
    
       }
       RNorm=rNorm; 
    }
    


    /**
     * Returns the current active upper triangular factor R.
     *
     * @return submatrix of R containing active columns or {@code null} if none
     */
    public RealMatrix getR() {
        if(iq==this.n) return R;
        if (iq > 0) {
            return R.getSubMatrix(0, iq - 1, 0, iq - 1);
        }
        return null;
    }

    /**
     * Returns the inverse of the active R factor.
     *
     * @return inverse of the current R or {@code null} if no active constraints
     */
    public RealMatrix getRInv() {
        if (iq > 0) {
            return inverseUpperTriangular(getR());
        }
        return null;
    }

    /**
     * Computes the inverse of an upper triangular matrix via backward substitution.
     *
     * @param U upper triangular matrix to invert
     * @return inverse of U
     */
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

    /**
     * Returns the inverse of L used internally.
     *
     * @return current J matrix
     */
    public RealMatrix getJ() {
        return J;
    }

    /**
     * Returns the inactive columns of J, starting at the first non-active index.
     *
     * @return submatrix of J for inactive columns or {@code null} if fully occupied
     */
    public RealMatrix getJ2() {
        if (iq == n) {
            return null;
        }
        return J.getSubMatrix(0, n - 1, iq, n - 1);
    }

    /**
     * Returns the number of active constraints.
     *
     * @return count of active constraints
     */
    public int getIq() {
        return iq;
    }
}
