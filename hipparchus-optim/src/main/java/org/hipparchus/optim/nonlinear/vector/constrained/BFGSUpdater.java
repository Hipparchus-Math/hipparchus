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
import org.hipparchus.linear.MatrixUtils;
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
    private final double sqrtEPS = FastMath.sqrt(Precision.EPSILON);

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
    private boolean DAMPED;
    private final int dim;
    private double FACTOR=1.0;

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
        this.dim = initialHess.getColumnDimension();
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
    public int update(final RealVector s, final RealVector y1) {
        final RealVector Hs = L.operate(L.preMultiply(s));
        final double sHs = s.dotProduct(Hs);
        double sty = s.dotProduct(y1);
        if (sHs <= 0.0) {
            resetHessian();
            return 1;
        }

        DAMPED = false;
        RealVector y = y1;
        if (sty < GAMMA * sHs) {
            DAMPED = true;
            final double denominator = sHs - sty;
            final double phi = (sHs - GAMMA * sHs) / denominator;
            y = y1.mapMultiply(phi).add(Hs.mapMultiply(1.0 - phi));
            sty = s.dotProduct(y);
        }

        if (!(sty > 0.0)) {
            return 2;
        }
        
         
         final double yy = y.dotProduct(y);
         FACTOR=1.0;
         if (!DAMPED && sty > Precision.SAFE_MIN) {
                FACTOR= yy / sty;
            }
         final double th = 1.0e-3;
         FACTOR = FastMath.max(th, FastMath.min(1.0 / th, FACTOR));
       
        if (!rankOneUpdate(s, y, Hs, sHs, sty)) {
            return 3;
        }

        return 0;
    }

    /**
     * Resets the Hessian approximation to its initial value.
     */
    public void resetHessian() {
        L = MatrixUtils.createRealIdentityMatrix(dim);
    }

    /**
     * Resets the Hessian approximation with information on the curvature.
     *
     * @param gamma scale factor
     */
    public void resetHessian(double gamma) {
        final double sqrtGAMMA = FastMath.sqrt(gamma);
        L = MatrixUtils.createRealIdentityMatrix(dim).scalarMultiply(sqrtGAMMA);
    }
    
    /**
     * Resets the Hessian approximation with information on the curvature for esternal usage
     * 
     */
    public void resetHessianFactor() {
        final double sqrtGAMMA = FastMath.sqrt(FACTOR);
        L = MatrixUtils.createRealIdentityMatrix(dim).scalarMultiply(sqrtGAMMA);
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
    private boolean rankOneUpdate(final RealVector s, final RealVector y, final RealVector Hs,
                                  final double sHs, final double sty) {

        final double rho = 1.0 / FastMath.sqrt(sty);
        final double theta = 1.0 / FastMath.sqrt(sHs);

        cholupdateLower(y, rho, +1);

        if (!cholupdateLower(Hs, theta, -1)) {
            
            resetHessian(FACTOR);

            return false;
        }

        return true;
    }

    /**
     * Performs a rank-one Cholesky update/downdate on L.
     * <p>
     * Updates L such that A' = A + sigma * alpha² * u * uᵀ, without refactorization.
     * Uses a numerically robust variant based on Givens/hyperbolic rotations.
     * </p>
     *
     * @param u direction vector
     * @param alpha scaling factor applied to {@code u}
     * @param sigma +1 for update, -1 for downdate
     * @return true if resulting matrix remains SPD, false otherwise
     */
    private boolean cholupdateLower(final RealVector u, final double alpha, final int sigma) {
        final int n = u.getDimension();
        final double[] work = u.toArray();
        for (int i = 0; i < n; ++i) {
            work[i] *= alpha;
        }
        // --- Single global guard for downdate ---
    if (sigma < 0) {
        final double[] z = new double[n];

        // Solve L z = work
        for (int i = 0; i < n; ++i) {
            double sum = work[i];
            for (int j = 0; j < i; ++j) {
                sum -= L.getEntry(i, j) * z[j];
            }
            final double lii = L.getEntry(i, i);
            z[i] = sum / lii;
        }

        double znorm2 = 0.0;
        for (int i = 0; i < n; ++i) {
            znorm2 += z[i] * z[i];
        }

        if (znorm2 >= 1.0 - Precision.EPSILON) {
            return false;
        }
    }

        for (int k = 0; k < n; ++k) {
            final double lkk = L.getEntry(k, k);
            final double xk = work[k];

            double r;
            if (sigma > 0) {
                r = FastMath.hypot(lkk, xk);
            } else {
                
                double radicand=(lkk - xk) * (lkk + xk);
                
                if(radicand<=Precision.EPSILON) 
                    r=sqrtEPS;
                //return false;
                else
                    r = FastMath.sqrt(radicand);
                    
            }
           
          
            final double c = r / lkk;
            final double s = xk / lkk;
            final double invC = 1.0 / c;

            L.setEntry(k, k, r);

            for (int i = k + 1; i < n; ++i) {
                final double lik = L.getEntry(i, k);
                final double xi = work[i];
                final double newLik = (lik + sigma * s * xi) * invC;
                final double newXi = c * xi - s * newLik;
                L.setEntry(i, k, newLik);
                work[i] = newXi;
            }
        }
        return true;
    }
    
  

    /**
     * Scales the Hessian H = L*L^T by a positive factor gamma,
     * by scaling the Cholesky factor L in-place.
     *
     * H_new = gamma * H
     * L_new = sqrt(gamma) * L
     *
     * @param gamma scaling factor (must be > 0)
     */
    private void scaleHessian(final double gamma) {

        final double s = FastMath.sqrt(gamma);

        final int n = L.getRowDimension();
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                L.setEntry(i, j, s * L.getEntry(i, j));
            }
        }
    }

    /**
     * Returns the current lower-triangular Cholesky factor L.
     *
     * @return current lower-triangular Cholesky factor
     */
    public RealMatrix getL() {
        return L;
    }
    
    
}