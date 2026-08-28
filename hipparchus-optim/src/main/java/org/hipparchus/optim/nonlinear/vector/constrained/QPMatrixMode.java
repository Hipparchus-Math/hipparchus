package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.optim.OptimizationData;

/**
 * Marker OptimizationData indicating that the quadratic matrix "G" passed to the QP solver
 * is not the Hessian H, but the lower Cholesky factor L such that H = L*L^T or the inverse of Cholesky factor L^-1 
 *
 * When present and set to true, QPDualActiveSolver will interpret function.getP()
 * as L or L^-1 and will NOT factorize it.
 */

public enum QPMatrixMode implements OptimizationData {

    /** QP Problem with full Hessian Matrix. */
    FULL,

    /** QP Problem with Cholesky L factor. */
    CHOLESKY,

    /** QP Problem with Inverse o Cholesky factor. */
    INVCHOLESKY

}

