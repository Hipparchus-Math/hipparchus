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
import org.hipparchus.linear.CholeskyDecomposition;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/**
 * BFGS Hessian updater with dynamic damping and robustness improvements.
 * <p>
 * Manages Hessian updates for SQP solvers by:
 * </p>
 * <ul>
 * <li>Checking curvature condition</li>
 * <li>Applying dynamic damping if necessary</li>
 * <li>Skipping update if curvature still fails after damping</li>
 * <li>Soft regularization of diagonal entries on repeated failures</li>
 * <li>Automatic Hessian reset after configurable failures</li>
 * </ul>
 *
 * @since 4.1
 */
public class BFGSUpdater {

    /**
     * AutoScaling Flag.
     */
    private final boolean SCALE;

    /**
     * EPS.
     */
    private final double EPS;

    /**
     * Damping factor.
     */
    private static final double GAMMA = 0.2;

    /**
     * trigger skip update for diagonal of Hessian.
     */
    private final double sqrtEPSmachine = FastMath.sqrt(Precision.EPSILON);

    /**
     * Tolerance for symmetric matrices decomposition.
     *
     * @since 4.1
     */
    private final double decompositionEpsilon;

    /**
     * Stored initial Hessian for resets.
     */
    private final RealMatrix initialH;

    /**
     * Current Cholesky factor L such that H = L·Lᵀ.
     */
    private RealMatrix L;

    /**
     * Creates a new updater.
     *
     * @param initialHess initial positive‐definite Hessian matrix
     * @param eps treshold to apply auto scale sty<sqrt(eps)
     * @param autoSCale true apply auto hessain rescaling
     * @param decompositionEpsilon tolerance for symmetric matrices
     * decomposition
     */
    public BFGSUpdater(final RealMatrix initialHess, final double eps, final boolean autoScale, final double decompositionEpsilon) {
        this.initialH = new Array2DRowRealMatrix(initialHess.getData());
        this.EPS = eps;
        this.SCALE = autoScale;
        this.decompositionEpsilon = decompositionEpsilon;
        resetHessian();
    }

    /**
     * Returns the current Hessian matrix H = L·Lᵀ.
     *
     * @return current Hessian
     */
    public RealMatrix getHessian() {
        return L.multiplyTransposed(L);
    }

    /**
     * Updates the Hessian approximation using the BFGS formula.
     * <p>
     * If curvature condition fails, applies damping or regularization.
     * </p>
     *
     * @param s displacement vector (x_{k+1} − x_k)
     * @param y1 gradient difference (∇f_{k+1} − ∇f_k)
     * @return type of update
     * 0: update is done
     * 1:sHs too small skip update
     * 2:sty<gamma*sHs skipped update
     * 3:sty<sqrt(eps)auto scale gamma*Hini
     * 4:singulatity detected during downdate reset to gamma*Hinit
     */
    public int update(RealVector s, RealVector y1) {
        RealVector Hs = L.operate(L.preMultiply(s)); // al posto di getHessian().operate(s)
        double sHs = s.dotProduct(Hs);
        if (sHs <= 0.0) {
           
//            double sty = s.dotProduct(y1);
//            double yy = y1.dotProduct(y1);
//            double gamma =yy /sty;
//            if (sty>0) {
//                gamma = FastMath.max(sqrtEPSmachine, FastMath.min(1/sqrtEPSmachine, gamma));
//                this.resetHessian(gamma);
//                return 1;
//            }
//            this.resetHessian();
            return 1;
        }

        RealVector y = damp(s, y1, Hs, sHs);
        if (y == null ) {
            
            return 2;//damp failed
        }
        
 
        
        
        if (SCALE ) {
            double sty = s.dotProduct(y);
            double yy = y.dotProduct(y);
            double gamma =yy /sty;
            if (gamma < FastMath.sqrt(EPS)) {
                gamma = FastMath.max(sqrtEPSmachine, FastMath.min(1/sqrtEPSmachine, gamma));
                this.resetHessian(gamma);
                return 3;
            }
        }

        // Attempt rank‐one BFGS update; regularize on failure
        if (rankOneUpdate(s, y, Hs, sHs))  return 0;
       
        return 4;
    }
    
    
    /**
     * Resets the Hessian approximation to its initial value.
     */
    public void resetHessian() {
        final CholeskyDecomposition ch = new CholeskyDecomposition(initialH, decompositionEpsilon, decompositionEpsilon);
        L = ch.getL();
    }

    /**
     * Resets the Hessian approximation with information on the curvature.
     *
     * @param gamma scale factor
     */
    public void resetHessian(double gamma) {
        final CholeskyDecomposition ch = new CholeskyDecomposition(initialH.scalarMultiply(gamma), decompositionEpsilon, decompositionEpsilon);
        L = ch.getL();
    }

    /**
     * Applies dynamic damping to y to satisfy curvature condition sᵀy ≥ γ sᵀHs.
     *
     * @param s search direction
     * @param y1 raw gradient difference
     * @param Hs vector
     * @param sHs value
     * @return damped y, or null if update should be skipped
     */
    public RealVector damp(RealVector s, RealVector y1, RealVector Hs, double sHs) {
        RealVector y = new ArrayRealVector(y1);
        double sty = s.dotProduct(y1);

        if (sty < GAMMA * sHs) {
            double phi = (1.0 - GAMMA) * sHs / (sHs - sty);
            if (phi >= 1.0) return y1;
             else if (phi <= 0.0) return Hs;
//            
            y = (y1.mapMultiply(phi)).add(Hs.mapMultiply(1.0 - phi));
            sty = s.dotProduct(y);
           
           
            if (sty < GAMMA * sHs) {

                return null;
            }

        }
        return y;
    }
    
  

    /**
     * Performs a BFGS rank‐one update on L.
     *
     * @param s displacement vector
     * @param y gradient difference vector
     * @param Hs vector
     * @param sHs value
     * @return true if update succeeded, false otherwise
     */
    private boolean rankOneUpdate(RealVector s, RealVector y, RealVector Hs, double sHs) {
        RealMatrix Lcopy = new Array2DRowRealMatrix(L.getData());
        double rho = 1.0 / FastMath.sqrt(s.dotProduct(y));
        double theta = 1.0 / FastMath.sqrt(sHs);
        RealVector v = y.mapMultiply(rho);
        RealVector w = Hs.mapMultiply(theta);
        cholupdateLower(v,+1);//upgrade
        if (!cholupdateLower(w,-1)) { //downdate    
            double gamma = 1.0;
            double sty = s.dotProduct(y);
            double yy = y.dotProduct(y);
            if (sty> Precision.SAFE_MIN) {
                gamma =yy /sty;
                gamma = FastMath.max(sqrtEPSmachine, FastMath.min(1/sqrtEPSmachine, gamma));
                this.resetHessian(gamma);
            }
            else this.resetHessian();
                
            
            return false;
        }

        return true;
    }

    /**
     * Performs a rank‐one Cholesky update/downdate on L.
     * <p>
     * Updates L such that A' = A+σu uᵀ or A' = A−u uᵀ, without refactorization.
     * </p>
     *
     * @param u update vector
     * @param sigma +1 for update, -1 for downdate
     * @return true if resulting matrix remains PD, false otherwise
     */
    private boolean cholupdateLower(RealVector u, int sigma) {
        int n = u.getDimension();
        RealVector temp = new ArrayRealVector(u);
        for (int i = 0; i < n; i++) {
            double lii = L.getEntry(i, i);
            double ui = temp.getEntry(i);
            double r2 = lii * lii + sigma * ui * ui;

            //skip or update 
            if (sigma < 0 && r2 < this.sqrtEPSmachine * lii * lii) {
//                if (sigma < 0 && r2 < Precision.EPSILON) {
                return false;
            }

            double r = Math.sqrt(r2);
            double c = r / lii;
            double s = ui / lii;
            L.setEntry(i, i, r);
            for (int j = i + 1; j < n; j++) {
                double lji = L.getEntry(j, i);
                double uj = temp.getEntry(j);
                double newLji = (lji + sigma * s * uj) / c;
                double newUj = c * uj - s * newLji;
                L.setEntry(j, i, newLji);
                temp.setEntry(j, newUj);
            }
        }
        return true;
    }
    
 



}
