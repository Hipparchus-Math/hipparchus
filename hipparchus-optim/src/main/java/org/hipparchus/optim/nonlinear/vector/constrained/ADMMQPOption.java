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

import org.hipparchus.optim.OptimizationData;

/** Container for {@link ADMMQPOptimizer} settings.
 * @since 3.1
 */
public class ADMMQPOption implements OptimizationData {

    /** Default Absolute and Relative Tolerance for convergence. */
    public static final double DEFAULT_EPS = 1.0e-5;

    /** Default Absolute and Relative Tolerance for Infeasible Criteria. */
    public static final double DEFAULT_EPS_INFEASIBLE = 1.0e-7;

    /** Default Value of regularization term sigma for Karush–Kuhn–Tucker solver. */
    public static final double DEFAULT_SIGMA = 1.0e-12;

    /** Default Value of Alpha filter for ADMM iteration. */
    public static final double DEFAULT_ALPHA = 1.6;

    /** Default Value for Enabling Problem Scaling. */
    public static final boolean DEFAULT_SCALING = true;

    /** Default Value for the Max Iteration for the scaling. */
    public static final int DEFAULT_SCALING_MAX_ITERATION = 10;

    /** Default Value for adapting the weight during iterations. */
    public static final boolean DEFAULT_RHO_UPDATE = true;

    /** Default Max Value for the Weight for ADMM iteration. */
    public static final double DEFAULT_RHO_MAX = 1.0e6;

    /** Default Min Value for the Weight for ADMM iteration. */
    public static final double DEFAULT_RHO_MIN = 1.0e-6;

    /** Default Max number of weight changes. */
    public static final int DEFAULT_MAX_RHO_ITERATION = 10;

    /** Default Value for enabling polishing the solution. */
    public static final boolean DEFAULT_POLISHING = false;

    /** Default Value for Iteration of polishing Algorithm. */
    public static final int DEFAULT_POLISHING_ITERATION = 5;

    /** Absolute and Relative Tolerance for convergence. */
    private double eps;

    /** Absolute and Relative Tolerance for Infeasible Criteria. */
    private double epsInfeasible;

    /** Value of regularization term sigma for Karush–Kuhn–Tucker solver. */
    private double sigma;

    /** Value of alpha filter for ADMM iteration. */
    private double alpha;

    /** Scaling enabling flag. */
    private boolean scaling;

    /** Value for the Max Iteration for the scaling. */
    private int scaleMaxIteration;

    /** Value for adapt the weight during iterations. */
    private boolean updateRho;

    /** Max Value for thr Weight for ADMM iteration. */
    private double rhoMax;

    /** Min Value for the Weight for ADMM iteration. */
    private double rhoMin;

    /** Max Value of changing the weight during iterations. */
    private int maxRhoIteration;

    /** Enabling flag for polishing the solution. */
    private boolean polishing;

    /** Value for Iteration of polishing Algorithm. */
    private int polishingIteration;

    /** Simple constructor.
     */
    public ADMMQPOption() {
        eps                = DEFAULT_EPS;
        epsInfeasible      = DEFAULT_EPS_INFEASIBLE;
        sigma              = DEFAULT_SIGMA;
        alpha              = DEFAULT_ALPHA;
        scaling            = DEFAULT_SCALING;
        scaleMaxIteration  = DEFAULT_SCALING_MAX_ITERATION;
        updateRho          = DEFAULT_RHO_UPDATE;
        rhoMax             = DEFAULT_RHO_MAX;
        rhoMin             = DEFAULT_RHO_MIN;
        maxRhoIteration    = DEFAULT_MAX_RHO_ITERATION;
        polishing          = DEFAULT_POLISHING;
        polishingIteration = DEFAULT_POLISHING_ITERATION;
    }

    /** Set absolute and Relative Tolerance for convergence.
     * @param eps absolute and Relative Tolerance for convergence
     */
    public void setEps(final double eps) {
        this.eps = eps;
    }

    /** Get absolute and Relative Tolerance for convergence.
     * @return absolute and Relative Tolerance for convergence
     */
    public double getEps() {
        return eps;
    }

    /** Set absolute and Relative Tolerance for infeasible criteria.
     * @param epsInfeasible absolute and Relative Tolerance for infeasible criteria
     */
    public void setEpsInfeasible(final double epsInfeasible) {
        this.epsInfeasible = epsInfeasible;
    }

    /** Get absolute and Relative Tolerance for infeasible criteria.
     * @return absolute and Relative Tolerance for infeasible criteria
     */
    public double getEpsInfeasible() {
        return epsInfeasible;
    }

    /** Set value of regularization term sigma for Karush–Kuhn–Tucker solver.
     * @param sigma value of regularization term sigma for Karush–Kuhn–Tucker solver
     */
    public void setSigma(final double sigma) {
        this.sigma = sigma;
    }

    /** Get value of regularization term sigma for Karush–Kuhn–Tucker solver.
     * @return value of regularization term sigma for Karush–Kuhn–Tucker solver
     */
    public double getSigma() {
        return sigma;
    }

    /** Set value of alpha filter for ADMM iteration.
     * @param alpha value of alpha filter for ADMM iteration
     */
    public void setAlpha(final double alpha) {
        this.alpha = alpha;
    }

    /** Get value of alpha filter for ADMM iteration.
     * @return value of alpha filter for ADMM iteration
     */
    public double getAlpha() {
        return alpha;
    }

    /** Set scaling enabling flag.
     * @param scaling if true, scaling is enabled
     */
    public void setScaling(final boolean scaling) {
        this.scaling = scaling;
    }

    /** Check if scaling is enabled.
     * @return true if scaling is enabled
     */
    public boolean isScaling() {
        return scaling;
    }

    /** Set max iteration for the scaling.
     * @param scaleMaxIteration max iteration for the scaling
     */
    public void setScaleMaxIteration(final int scaleMaxIteration) {
        this.scaleMaxIteration = scaleMaxIteration;
    }

    /** Get max iteration for the scaling.
     * @return max iteration for the scaling
     */
    public int getScaleMaxIteration() {
        return scaleMaxIteration;
    }

    /** Set weight updating flag.
     * @param updateRho if true, weight is updated during iterations
     */
    public void setUpdateRho(final boolean updateRho) {
        this.updateRho = updateRho;
    }

    /** Check if weight updating is enabled.
     * @return true if weight is updated during iterations
     */
    public boolean updateRho() {
        return updateRho;
    }

    /** Set min Value for the Weight for ADMM iteration.
     * @param rhoMin min Value for the Weight for ADMM iteration
     */
    public void setRhoMin(final double rhoMin) {
        this.rhoMin = rhoMin;
    }

    /** Get min Value for the Weight for ADMM iteration.
     * @return min Value for the Weight for ADMM iteration
     */
    public double getRhoMin() {
        return rhoMin;
    }

    /** Set max Value for the Weight for ADMM iteration.
     * @param rhoMax max Value for the Weight for ADMM iteration
     */
    public void setRhoMax(final double rhoMax) {
        this.rhoMax = rhoMax;
    }

    /** Get max Value for the Weight for ADMM iteration.
     * @return max Value for the Weight for ADMM iteration
     */
    public double getRhoMax() {
        return rhoMax;
    }

    /** Set max number of weight changes.
     * @param maxRhoIteration max number of weight changes
     */
    public void setMaxRhoIteration(final int maxRhoIteration) {
        this.maxRhoIteration = maxRhoIteration;
    }

    /** Get max number of weight changes.
     * @return max number of weight changes
     */
    public int getMaxRhoIteration() {
        return maxRhoIteration;
    }

    /** Set polishing enabling flag.
     * @param polishing if true, polishing is enabled
     */
    public void setPolishing(final boolean polishing) {
        this.polishing = polishing;
    }

    /** Check if polishing is enabled.
     * @return true if polishing is enabled
     */
    public boolean isPolishing() {
        return polishing;
    }

    /** Set number of iterations of polishing algorithm.
     * @param polishingIteration number of iterations of polishing algorithm
     */
    public void setPolishingIteration(final int polishingIteration) {
        this.polishingIteration = polishingIteration;
    }

    /** Get number of iterations of polishing algorithm.
     * @return number of iterations of polishing algorithm
     */
    public int getPolishIteration() {
        return polishingIteration;
    }

}
