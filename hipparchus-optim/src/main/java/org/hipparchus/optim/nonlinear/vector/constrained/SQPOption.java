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

/** Parameter for SQP Algorithm.
 * @since 3.1
 */
public class SQPOption implements OptimizationData {

    /** Default convergence criteria. */
    public static final int DEFAULT_CONV_CRITERIA = 1;

    /** Default tolerance for convergence and active constraint. */
    public static final double DEFAULT_EPSILON = 1.0e-7;//>0

    /** Default weight for augmented QP subproblem. */
    public static final double DEFAULT_RHO = 100.0;//rho>1

    /** Default max value admitted for additional variable in QP subproblem. */
    public static final  double DEFAULT_SIGMA_MAX = 0.90;//0<sigma<1

    /** Default max iteration admitted for QP subproblem. */
    public static final  int DEFAULT_QP_MAX_LOOP = 4;

    /** Default parameter for evaluation of Armijo condition for descend direction. */
    public static final  double DEFAULT_MU = 1.0e-4;//[0,0.5]

    /** Default parameter for quadratic line search. */
    public static final  double DEFAULT_B = 0.5;//[0;1]

    /** Default flag for using BFGS update formula. */
    public static final  boolean DEFAULT_USE_FUNCTION_HESSIAN = false;

    /** Default max iteration before reset hessian. */
    public static final  int DEFAULT_MAX_LINE_SEARCH_ITERATION = 50;

    /** Default Gradient mode. */
    public static final GradientMode DEFAULT_GRADIENT_MODE = GradientMode.FORWARD;

    /** Convergence criteria*/
    private int convCriteria;

    /** Tolerance for convergence and active constraint evaluation. */
    private double eps;

    /** Weight for augmented QP subproblem. */
    private double rhoCons;

     /** Max value admitted for the solution of the additional variable in QP subproblem. */
    private double sigmaMax;

    /** Max iteration admitted for QP subproblem evaluation.
     * (over this threshold the descend direction will be approximated using the merit function).
     */
    private int qpMaxLoop;

    /** Parameter for evaluation of Armijo condition for descend direction.
     * (fhi(alfa)-fhi(0)<=mu * alfa * fhi'(0)) */
    private double mu;

    /** Parameter for quadratic line search. */
    private double b;

    /** Max Iteration for the line search. */
    private int maxLineSearchIteration;

    /** Enable or Disable using direct the function Hessian. */
    private boolean useFunHessian;

    /** Gradient Mode. */
    private GradientMode gradientMode;

    /** Simple constructor.
     * <p>
     * This constructor uses all defaults values.
     * </p>
     */
    public SQPOption() {
        this.convCriteria           = DEFAULT_CONV_CRITERIA;
        this.eps                    = DEFAULT_EPSILON;
        this.rhoCons                = DEFAULT_RHO;
        this.sigmaMax               = DEFAULT_SIGMA_MAX;
        this.qpMaxLoop              = DEFAULT_QP_MAX_LOOP;
        this.mu                     = DEFAULT_MU;
        this.b                      = DEFAULT_B;
        this.maxLineSearchIteration = DEFAULT_MAX_LINE_SEARCH_ITERATION;
        this.useFunHessian          = DEFAULT_USE_FUNCTION_HESSIAN;
        this.gradientMode           = DEFAULT_GRADIENT_MODE;
    }

     /** Set Gradient mode
     * @param gradientMode gradient mode
     */
    public void setGradientMode(final GradientMode gradientMode) {
        this.gradientMode = gradientMode;
    }

    /** Get Gradient Mode.
     * @return Gradient Mode
     */
    public GradientMode getGradientMode() {
        return gradientMode;
    }

    /** Set convergence criteria.
     * @param convCriteria convergence criteria
     */
    public void setConvCriteria(final int convCriteria) {
        this.convCriteria = convCriteria;
    }

    /** Get convergence criteria.
     * @return convergence criteria
     */
    public int getConvCriteria() {
        return convCriteria;
    }

    /** Set tolerance for convergence and active constraint evaluation.
     * @param eps tolerance for convergence and active constraint evaluation
     */
    public void setEps(final double eps) {
        this.eps = eps;
    }

    /** Get tolerance for convergence and active constraint evaluation.
     * @return tolerance for convergence and active constraint evaluation
     */
    public double getEps() {
        return eps;
    }

    /** Set weight for augmented QP subproblem.
     * @param rhoCons weight for augmented QP subproblem
     */
    public void setRhoCons(final double rhoCons) {
        this.rhoCons = rhoCons;
    }

    /** Get weight for augmented QP subproblem.
     * @return weight for augmented QP subproblem
     */
    public double getRhoCons() {
        return rhoCons;
    }

    /** Set max value admitted for the solution of the additional variable in QP subproblem.
     * @param sigmaMax max value admitted for the solution of the additional variable in QP subproblem
     */
    public void setSigmaMax(final double sigmaMax) {
        this.sigmaMax = sigmaMax;
    }

    /** Get max value admitted for the solution of the additional variable in QP subproblem.
     * @return max value admitted for the solution of the additional variable in QP subproblem
     */
    public double getSigmaMax() {
        return sigmaMax;
    }

    /** Set max iteration admitted for QP subproblem evaluation.
     * @param qpMaxLoop max iteration admitted for QP subproblem evaluation
     */
    public void setQpMaxLoop(final int qpMaxLoop) {
        this.qpMaxLoop = qpMaxLoop;
    }

    /** Get max iteration admitted for QP subproblem evaluation.
     * @return max iteration admitted for QP subproblem evaluation
     */
    public int getQpMaxLoop() {
        return qpMaxLoop;
    }

    /** Set parameter for evaluation of Armijo condition for descend direction.
     * @param mu parameter for evaluation of Armijo condition for descend direction
     */
    public void setMu(final double mu) {
        this.mu  = mu;
    }

    /** Get parameter for evaluation of Armijo condition for descend direction.
     * @return parameter for evaluation of Armijo condition for descend direction
     */
    public double getMu() {
        return mu;
    }

    /** Set parameter for quadratic line search.
     * @param b parameter for quadratic line search
     */
    public void setB(final double b) {
        this.b = b;
    }

    /** Get parameter for quadratic line search.
     * @return parameter for quadratic line search
     */
    public double getB() {
        return b;
    }

    /** Set max Iteration for the line search
     * @param maxLineSearchIteration max Iteration for the line search
     */
    public void setMaxLineSearchIteration(final int maxLineSearchIteration) {
        this.maxLineSearchIteration = maxLineSearchIteration;
    }

    /** Get max Iteration for the line search
     * @return max Iteration for the line search
     */
    public int getMaxLineSearchIteration() {
        return maxLineSearchIteration;
    }

    /** Enable or Disable using direct the function Hessian.
     * @param useFunHessian enable or Disable using direct the function Hessian
     */
    public void setUseFunHessian(final boolean useFunHessian) {
        this.useFunHessian = useFunHessian;
    }

    /** Check if using direct the function Hessian is enabled or disabled.
     * @return true if using direct the function Hessian is enabled
     */
    public boolean useFunHessian() {
        return useFunHessian;
    }

}
