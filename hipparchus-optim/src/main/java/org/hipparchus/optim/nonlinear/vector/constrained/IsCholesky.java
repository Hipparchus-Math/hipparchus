package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.optim.OptimizationData;

/**
 * Marker OptimizationData indicating that the quadratic matrix "G" passed to the QP solver
 * is not the Hessian H, but the lower Cholesky factor L such that H = L*L^T.
 *
 * When present and set to true, QPDualActiveSolver will interpret function.getP()
 * as L and will NOT factorize it.
 */
public final class IsCholesky implements OptimizationData {

    private final boolean cholesky;

    public IsCholesky(final boolean cholesky) {
        this.cholesky = cholesky;
    }

    public boolean isCholesky() {
        return cholesky;
    }
}
