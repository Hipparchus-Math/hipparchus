package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.optim.OptimizationData;

/**
 * Marker OptimizationData indicating that the quadratic matrix "G" passed to the QP solver
 * is not the Hessian H, but the lower Cholesky factor L such that H = L*L^T.
 *
 * When present and set to true, QPDualActiveSolver will interpret function.getP()
 * as L and will NOT factorize it.
 */

public enum QPMatrixMode implements OptimizationData {

    /** Compute gradients from the objective and constraints functions themselves. */
    FULL,

    /** Compute gradients using forward difference. */
    CHOLESKY,

    /** Compute gradients using central differences. */
    INVCHOLESKY

}

