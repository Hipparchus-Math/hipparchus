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

    /**
     * SQP_ConvCriteria default convergenzia criteria
     * dxT * H*dx<epsilon * epsilon SQP_ConvCriteria = 0
     * ||dx||<=sqrt(epsilon)*(1+||x|| SQP_ConvCriteria = 1
     * evaluation {@value SQPOption#SQP_eps}
     *
     */
    public static final int SQP_ConvCriteria =1;//

    /**
     * SQP_eps default tolerance for convergence and active constraint
     * evaluation {@value SQPOption#SQP_eps}
     *
     */
    public static final double SQP_eps = 1.0e-9;//>0
    /**
     * SQP_rhoCons default weight for augmented QP subproblem
     * {@value SQPOption#SQP_rhoCons}
     */
    public static final double SQP_rhoCons = 100.0;//rho>1

    /**
     * SQP_sigmaMax default max value admitted for additional variable in QP subproblem
     * {@value SQPOption#SQP_sigmaMax}
     */
    public static final  double SQP_sigmaMax = 0.90;//0<sigma<1
    /**
     * SQP_qpMaxLoop default max iteration admitted for QP subprobleam
     * evaluation {@value SQPOption#SQP_qpMaxLoop}
     */

    public static final  int SQP_qpMaxLoop = 4;

    /**
     * SQP_mu default parameter for evaluation of Armijo condition for descend
     * direction {@value SQPOption#SQP_mu}
     */
    public static final  double SQP_mu = 0.1;//[0,0.5]
    /**
     * SQP_b default parameter for quadratic line search
     * {@value SQPOption#SQP_b}
     */
    public static final  double SQP_b = 0.1;//[0;1]
    /**
     * SQP_useFunHessian (false use use BFGS update formula)
     * {@value SQPOption#SQP_useFunHessian}
     */
    public static final  boolean SQP_useFunHessian = false;
    /**
     * SQP_maxLineSearchIteration max iteration before reset hessian
     * {@value SQPOption#SQP_maxLineSearchIteration}
     */
    public static final  int SQP_maxLineSearchIteration = 20;





    /**
     * Convergence Criteria
     * {@link SQPOption#SQP_ConvCriteria}
     *
     *
     */

     public final int convCriteria ;//




    /**
     * Tolerance for convergence and active constraint evaluation.
     * {@link SQPOption#SQP_eps}
     *
     */
    public final double eps;


    /**
     *  Weight for augmented QP subproblem.
     * {@link SQPOption#SQP_rhoCons}
     */
    public final double rhoCons;

     /**
     *  Max value admitted for the solution of the additional variable in QP subproblem.
     * {@link SQPOption#SQP_sigmaMax}
     */
    public final double sigmaMax;

    /**
     * Max iteration admitted for QP subproblem evaluation.
     * (over this threshold the descend direction will be approximated using the merit function).
     * {@link SQPOption#SQP_qpMaxLoop}
     */
    public final int qpMaxLoop;


    /**
     * Parameter for evaluation of Armijo condition for descend direction.
     * (fhi(alfa)-fhi(0)<=mu * alfa * fhi'(0))
     * {@link SQPOption#SQP_mu}
     */
    public final double mu;

    /**
     * Parameter for quadratic line search.
     * {@value SQPOption#SQP_b}
     */
    public final double b;

    /**
     *  Enable or Disable using direct the function Hessian.
     * (false use use BFGS update formula)
     * {@link SQPOption#SQP_useFunHessian}
     */
    public final boolean useFunHessian;

    /**
     *  Max Iteration for the line search.
     * {@link SQPOption#SQP_maxLineSearchIteration}
     */
    public final int maxLineSearchIteration;

    public SQPOption() {
        this.convCriteria = SQP_ConvCriteria;
        this.eps = SQP_eps;
        this.rhoCons = SQP_rhoCons;
        this.sigmaMax = SQP_sigmaMax;
        this.qpMaxLoop = SQP_qpMaxLoop;
        this.mu = SQP_mu;
        this.b = SQP_b;
        this.useFunHessian = SQP_useFunHessian;
        this.maxLineSearchIteration = SQP_maxLineSearchIteration;
    }

    public SQPOption(int convergenceCriteria,double SQPepsAbs, double SQP_rhoCons, double SQP_sigmaMax, int SQP_qpMaxLoop, double SQP_mu, double SQP_b, int SQP_maxLineSearchIteration, boolean SQP_useFunHessian) {
         this.convCriteria = convergenceCriteria;
        this.eps = SQPepsAbs;
        this.rhoCons = SQP_rhoCons;
        this.sigmaMax = SQP_sigmaMax;
        this.qpMaxLoop = SQP_qpMaxLoop;
        this.mu = SQP_mu;
        this.b = SQP_b;
        this.useFunHessian = SQP_useFunHessian;
        this.maxLineSearchIteration = SQP_maxLineSearchIteration;

    }

}
