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

import org.hipparchus.util.FastMath;

import java.util.LinkedList;
import java.util.Queue;

/** Robust Line Search strategy.
 * <p>
 * This class manages monotone and non-monotone line search. Switching is
 * automatic after repeated failures of the monotone search.
 * </p>
 *
 * <p>
 * Typical parameter values used in Schittkowski's algorithms:
 * </p>
 * <ul>
 * <li>maxHistory = 5 to 10 (memory for non-monotone search)</li>
 * <li>sigma = 1e-4 or 1e-5 (Armijo sufficient decrease parameter)</li>
 * <li>beta = 0.5 (step reduction factor)</li>
 * <li>alphaMin (minimum allowed step length)</li>
 * <li>alphaMax = 1.0 (maximum allowed step length)</li>
 * <li>maxMonotoneFailures = 20 to 50 (failures before switching to non-monotone)</li>
 * <li>maxBadSteps = 3 to 5 (allowed consecutive bad steps before Hessian reset)</li>
 * </ul>
 * @since 4.1
 */
public class LineSearch {

    /** Maximum penalty history size. */
    private final int maxHistory;

    /** Parameter for evaluation of Armijo condition for descend direction. */
    private final double sigma;

    /** Step reduction factor. */
    private final double beta;

    /** Minimum step length. */
    private  double alphaMin;

    /** Failures before switching to non-monotone. */
    private final int maxMonotoneFailures;

    /** Allowed consecutive bad steps before Hessian reset. */
    private final int maxBadSteps;

    /** Penalty history. */
    private final Queue<Double> history;

    /** Number of iterations in line search. */
    private int searchCount;

    /** Number of monotone search failures. */
    private int monotoneFailures;

    /** Number of consecutive bad steps. */
    private int badStepCount;

    /** Indicator for non-monotone search. */
    private boolean nonMonotoneEnabled;

    /** Indicator for detected bad step. */
    private boolean badStepDetected;

    /** Indicator for too many consecutive bad steps. */
    private boolean badStepFailed;

    /** Simple constructor.
     * @param eps tolerance for convergence and active constraint evaluation
     * @param maxHistory maximum penalty history size
     * @param mu parameter for evaluation of Armijo condition for descend direction
     * @param beta step reduction factor
     * @param maxMonotoneFailures failures before switching to non-monotone
     * @param maxBadSteps allowed consecutive bad steps before Hessian reset
     */
    public LineSearch(final double eps, final int maxHistory, final double mu, final double beta,
                      final int maxMonotoneFailures, final int maxBadSteps) {
        this.maxHistory          = maxHistory;
        this.sigma               = mu;
        this.beta                = beta;
        this.alphaMin            = FastMath.min(1.0e-12, eps);
        this.maxMonotoneFailures = maxMonotoneFailures;
        this.maxBadSteps         = maxBadSteps;
        this.history             = new LinkedList<>();
        this.nonMonotoneEnabled  = false;
        resetBadStepCount();
    }

    /** Check if bad step has been detected.
     * @return true if bad step has been detected
     */
    public boolean isBadStepDetected() {
        return badStepDetected;
    }

    /** Check if too many consecutive bad step have been detected.
     * @return true if too many consecutive bad step have been detected
     */
    public boolean isBadStepFailed() {
        return badStepFailed;
    }

    /** Get number of iterations in line search.
     * @return number of iterations in line search
     */
    public int getIteration() {
        return searchCount;
    }

    /** Reset bad steps count and indicators.
     */
    public void resetBadStepCount() {
        badStepCount = 0;
        monotoneFailures = 0;
        badStepDetected = false;
        badStepFailed = false;
    }

    /** Save penalty value when step is accepted for reusing in case of non-monotone research.
     * @param fx penalty
     */
    public void updateHistory(final double fx) {
        if (nonMonotoneEnabled) {
            history.add(fx);
            if (history.size() > maxHistory) {
                history.poll();
            }
        }
    }

    /** Verify Armijo condition for accept step.
     * @param fxNew penalty at candidate point x+dx*alpha
     * @param fxCurrent penalty at the current point
     * @param alpha step length
     * @param directionalDeriv penalty gradient
     * @return true or false
     */
    public boolean acceptStep(double fxNew, double fxCurrent, double alpha, double directionalDeriv) {
        double ref = fxCurrent;
        if (nonMonotoneEnabled) {
            for (double v : history) {
                ref = FastMath.max(ref, v);
            }
        }
        // alfaPenalty - currentPenalty > getSettings().getMu() * alpha * currentPenaltyGrad
        return fxNew < ref + sigma * alpha * directionalDeriv;
    }

    /** Mark Good Step if line search worked.
     */
    public void markGoodStep() {
        nonMonotoneEnabled = false;
        monotoneFailures = 0;
        badStepCount = 0;
        badStepDetected = false;
        badStepFailed = false;
    }

    /** Mark Bad Step if line search failed.
     */
    private void markBadStep() {
        nonMonotoneEnabled = false;
        badStepCount++;
        badStepDetected = true;
        monotoneFailures = 0;
        if (badStepCount > maxBadSteps) {
            badStepFailed = true;
        }
    }

    /** Update alpha qith quadratic curvature.
     * @param alpha step length
     * @param fxCurrent penalty at the current point
     * @param fxNew penalty at candidate point x+dx*alpha
     * @param directionalDeriv penalty gradient
     * @return alpha
     */
    private double updateStepLength(final double alpha, final double fxCurrent, final double fxNew,
                                    final double directionalDeriv) {
        final double numerator   = 0.5 * alpha * alpha * directionalDeriv;
        final double denominator = alpha * directionalDeriv - fxNew + fxCurrent;
        if (Math.abs(denominator) > 1e-12) {
            double alphaStar = numerator / denominator;
            return FastMath.max(alphaStar, alpha * beta);
        } else {
            return alpha * beta;
        }
    }

    /**
     * Main line search process: tries monotone first, then non-monotone.
     *
     * @param f penalty function
     * @return alpha
     */
    public double search(final MeritFunctionL2 f) {
        searchCount = 0;
        double fxCurrent = f.getPenaltyEval();
        alphaMin = FastMath.max(1e-12, 1e-6 / FastMath.max(1.0, f.getDx().getNorm()));
        //double fxCurrent = f.value(0);
        double directionalDeriv = f.getGradient();

        // Monotone Search
        nonMonotoneEnabled = false;
        double alpha = searchLoop(f, fxCurrent, directionalDeriv);
        if (!Double.isNaN(alpha)) {
            // monotone search succeeded
            return alpha;
        } else {
            // Non-monotone search
            nonMonotoneEnabled = true;
            alpha = searchLoop(f, fxCurrent, directionalDeriv);
            if (!Double.isNaN(alpha)) {
                // non-monotone search succeeded
                return alpha;
            } else {
                // last attempt before we declare bad step
                final double fxNew = f.value(alphaMin);
                searchCount++;
                if (acceptStep(fxNew, fxCurrent, alphaMin, directionalDeriv)) {
                    markGoodStep();
                    updateHistory(fxNew);
                } else {
                    markBadStep();
                }
                return alphaMin;
            }
        }

    }

    /** Iterative search (either monotone or non-monotone).
     * @param f penalty function
     * @param fxCurrent penalty function at current point
     * @param directionalDeriv directional derivative
     * @return accepted alpha, or NaN if search failed
     */
    private double searchLoop(final MeritFunctionL2 f, final double fxCurrent,
                              final double directionalDeriv) {

        // start value for step length
        double alpha = 1.0;

        // step length reduction loop
        while (alpha >= alphaMin) {

            final double fxNew = f.value(alpha);
            if (acceptStep(fxNew, fxCurrent, alpha, directionalDeriv)) {
                // we have found an acceptable step length
                markGoodStep();
                updateHistory(fxNew);
                return alpha;
            }

            // step length was not accepted, continue iteration
            alpha = updateStepLength(alpha, fxCurrent, fxNew, directionalDeriv);
            monotoneFailures++;
            if (monotoneFailures >= maxMonotoneFailures) {
                nonMonotoneEnabled = !nonMonotoneEnabled;
                break;
            }
            searchCount++;

        }

        // we failed to find an acceptable step length
        return Double.NaN;

    }

}
