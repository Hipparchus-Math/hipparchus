package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.optim.OptimizationData;

/**
 * Custom convergence criterion for an SQP optimizer.
 *
 * @author rocca
 */

public interface SQPCustomCriterion extends OptimizationData {

    /**
     * Checks whether the SQP iteration has converged.
     *
     * @param alpha accepted step length
     * @param dxNorm norm of the primal step
     * @param dxHdx quadratic step measure
     * @param complSlack complementary-slackness residual
     * @param kkt KKT residual
     * @param viol constraint violation
     * @param funDiff objective-function difference
     * @return {@code true} if convergence has been reached
     */
    boolean converged(double alpha,
                      double dxNorm,
                      double dxHdx,
                      double complSlack,
                      double kkt,
                      double viol,
                      double funDiff);
}