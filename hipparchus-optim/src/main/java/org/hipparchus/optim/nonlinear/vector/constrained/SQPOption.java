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
import org.hipparchus.util.Precision;

/** Parameter for SQP Algorithm.
 * @since 3.1
 */
public class SQPOption implements OptimizationData {


    /** Default tolerance for convergence and active constraint(proportional to sqrt of function an constraints evaluation accuracy) */
    public static final double DEFAULT_EPSILON = 1.0e-7;//>0

    /** Default tolerance for QP convergence(proportional to function and constraints evaluation accuracy)  */
    public static final double DEFAULT_EPSILON_QP = Precision.EPSILON;//>0

    /** Default weight for augmented QP subproblem. */
    public static final double DEFAULT_RHO = 100.0;//rho>1

    /** Default max value admitted for additional variable in QP subproblem. */
    public static final  double DEFAULT_SIGMA_MAX = 0.90;//0<sigma<1

    /** Default max iteration admitted for QP subproblem. */
    public static final  int DEFAULT_QP_MAX_LOOP = 4;

    /** Default parameter for evaluation of Armijo condition for descend direction. */
    public static final  double DEFAULT_MU = 1.0e-4;//[0,0.5]

    /** Default parameter for quadratic line search. */
    public static final  double DEFAULT_B = 0.1;//[0;1]

    /** Default parameter for alpha min during line search. */
    public static final  double DEFAULT_ALPHA_MIN = 1.0e-14;

    /** Default parameter for non monotone line search history. */
    public static final  int DEFAULT_HISTORY = 10;

    /** Default flag for using BFGS update formula. */
    public static final  boolean DEFAULT_USE_FUNCTION_HESSIAN = false;

    /** Default max iteration before reset hessian. */
    public static final  int DEFAULT_MAX_LINE_SEARCH_ITERATION = 30;

    /** Default Gradient mode. */
    public static final GradientMode DEFAULT_GRADIENT_MODE = GradientMode.FORWARD;

    /** Default Gradient EPS (proportional to function and constraints evaluation accuracy). */
    public static final double DEFAULT_GRAD_EPS = Precision.EPSILON;

    /** Default Max Iterations. */
    public static final int DEFAULT_MAX_ITERATION = 300;

    /** Convergence criteria*/
    private int convCriteria;

    /** Tolerance for convergence and active constraint evaluation. */
    private double eps;

    /** Tolerance for qp convergence */
    private double epsQP;

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

    /** Parameter line search alpha min. */
    private double alphaMin;

    /** Parameter for non monotone line search history. */
    private int history;

    /** Max Iteration for the line search. */
    private int maxLineSearchIteration;

    /** Enable or Disable using direct the function Hessian. */
    private boolean useFunHessian;

    /** Gradient Mode. */
    private GradientMode gradientMode;

    /** Gradient eps. */
    private double gradientEps;

    /** Max Iteration */
    private int maxIteration;


    /**Convergence criteria function for specific problems
     * to avoid not significative iterations
     */
    private SQPCustomCriterion customCovergence=null;

    /** Simple constructor.
     * <p>
     * This constructor uses all defaults values.
     * </p>
     */
    public SQPOption() {

        this.eps                    = DEFAULT_EPSILON;
        this.epsQP                  = DEFAULT_EPSILON_QP;
        this.rhoCons                = DEFAULT_RHO;
        this.sigmaMax               = DEFAULT_SIGMA_MAX;
        this.qpMaxLoop              = DEFAULT_QP_MAX_LOOP;
        this.mu                     = DEFAULT_MU;
        this.b                      = DEFAULT_B;
        this.alphaMin               = DEFAULT_ALPHA_MIN;
        this.maxLineSearchIteration = DEFAULT_MAX_LINE_SEARCH_ITERATION;
        this.useFunHessian          = DEFAULT_USE_FUNCTION_HESSIAN;
        this.gradientMode = DEFAULT_GRADIENT_MODE;
        this.gradientEps  = DEFAULT_GRAD_EPS;
        this.maxIteration = DEFAULT_MAX_ITERATION;
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

     /** Set Gradient Eps
      *  Choose value proportional to function and constraints evaluation accuracy)
      * @param gradientEps gradient mode
      */
    public void setGradientEps(final double gradientEps) {
        this.gradientEps = gradientEps;
    }

    /** Get Gradient Eps.
     * @return Gradient Eps
     */
    public double getGradientEps() {
        return gradientEps;
    }


    /** Set tolerance for convergence
     *  Choose value proportional to  sqrt of function and constraints evaluation accuracy
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

     /** Set tolerance for QP convergence
      *  Choose value proportional to function and constraints evaluation accuracy
      * @param epsQP tolerance for convergence
      */
    public void setEpsQP(final double epsQP) {
        this.epsQP = epsQP;
    }

    /** Get tolerance for QP convergence
     * @return tolerance for QP convergence
     */
    public double getEpsQP() {
        return epsQP;
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

    /** Set parameter for line search alpha min.
     * @param alphaMin parameter line search
     */
    public void setAlphaMin(final double alphaMin) {
        this.alphaMin = alphaMin;
    }

    /** Get parameter line search alpha min.
     * @return parameter for quadratic line search
     */
    public double getAlphaMin() {
        return this.alphaMin;
    }

    /** Set parameter for nom monotone line search history.
     * @param history parameter line search
     */
    public void setHistory(final int history) {
        this.history = history;
    }

    /** Get parameter fot non monotone line search history.
     * @return parameter for non monotone line search history
     */
    public int getHistory() {
        return this.history;
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

    /** Set max Iteration
     * @param maxIteration max Iteration
     */
    public void setMaxIteration(final int maxIteration) {
        this.maxIteration = maxIteration;
    }

    /** Get max Iteration
     * @return max Iteration
     */
    public int getMaxIteration() {
        return maxIteration;
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

     /** Set custom convergence criteria function
     * @param fun convergence criteria function
     */
    public void setConvergenceFunction(final SQPCustomCriterion fun) {
        this.customCovergence = fun;
    }

    /** Get Convergence criteria Function
     * @return the function of the custom convergenc criteria
     */
    public SQPCustomCriterion getConvergenceFunction() {
        return this.customCovergence;
    }

}
