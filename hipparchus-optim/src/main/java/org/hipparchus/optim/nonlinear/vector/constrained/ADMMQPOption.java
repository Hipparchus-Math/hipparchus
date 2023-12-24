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

public class ADMMQPOption implements OptimizationData {

    /**
     * Default Absolute and Relative Tolerance for convergence
     * {@value ADMMQPOption#ADMM_eps}
     */
    public static final double ADMM_eps = 1.0e-5;

    /**
     * Default Absolute and Relative Tolerance for Infeasible Criteria
     * {@value ADMMQPOption#ADMM_eps}
     */
    static double ADMM_epsInf = 1.0e-7;

    /**
     * Default Value of Sigma for KKT solver{@value ADMMQPOption#ADMM_sigma}
     */
    public static final double ADMM_sigma = 1.0e-12;
    /**
     * Default Value of Alfa filter for ADMM iteration
     * {@value ADMMQPOption#ADMM_alfa}
     */
    public static final double ADMM_alfa = 1.6;

    /**
     * Default Value for Enabling Problem Scaling
     * {@value ADMMQPOption#ADMM_scale}
     */
    public static final boolean ADMM_scale =true;

    /**
     * Default Value for the Max Iteration for the scaling
     * {@value ADMMQPOption#ADMM_scale_maxIteration}
     */
    public static final int ADMM_scale_maxIteration = 10;

    /**
     * Default Value for adapt the weight during iterations
     * {@value ADMMQPOption#ADMM_rho_update}
     */
    public static final boolean ADMM_rho_update = true;

    /**
     * Default Max Value for thr Weight for ADMM iteration
     * {@value ADMMQPOption#ADMM_rho_max}
     */
    public static final double ADMM_rho_max = 10e6;

    /**
     * Default Min Value for thr Weight for ADMM iteration
     * {@value ADMMQPOption#ADMM_rho_min}
     */
    public static final double ADMM_rho_min = 1.0e-6;

    /**
     * Default Max Value of changing the weight during iterations
     * {@value ADMMQPOption#ADMM_rho_maxIteration}
     */
    public static final int ADMM_rho_maxIteration = 10;

    /**
     * Default Value for enabling the polish solution
     * {@value ADMMQPOption#ADMM_polish}
     */
    public static final boolean ADMM_polish = false;
    /**
     * Default Value for Iteration of polish Algorithm
     * {@value ADMMQPOption#ADMM_polish_iteration}
     */
    public static final int ADMM_polish_iteration = 5;

    /**
     * Absolute and Relative Tolerance for convergence
     * {@link ADMMQPOption#ADMM_eps}
     */
    public double eps = ADMMQPOption.ADMM_eps;

    /**
     * Absolute and Relative Tolerance for Infeasible Criteria
     * {@link ADMMQPOption#ADMM_eps}Inf
     */
    public double epsInf = ADMMQPOption.ADMM_epsInf;

    /**
     * Value of Sigma for KKT solver {@link ADMMQPOption#ADMM_sigma}
     */
    public double sigma = ADMMQPOption.ADMM_sigma;
    /**
     * Value of Alfa filter for ADMM iteration {@link ADMMQPOption#ADMM_alfa}
     */
    public double alfa = ADMMQPOption.ADMM_alfa;

    /**
     * Value for Enabling Problem Scaling {@link ADMMQPOption#ADMM_scale}
     */
    public boolean scale = ADMMQPOption.ADMM_scale;

    /**
     * Value for the Max Iteration for the scaling
     * {@link ADMMQPOption#ADMM_scale_maxIteration}
     */
    public int scale_maxIteration = ADMMQPOption.ADMM_scale_maxIteration;

    /**
     * Value for adapt the weight during iterations
     * {@link ADMMQPOption#ADMM_rho_update}
     */
    public boolean rho_update = ADMMQPOption.ADMM_rho_update;

    /**
     * Max Value for thr Weight for ADMM iteration
     * {@link ADMMQPOption#ADMM_rho_max}
     */
    public double rho_max = ADMMQPOption.ADMM_rho_max;

    /**
     * Min Value for thr Weight for ADMM iteration
     * {@link ADMMQPOption#ADMM_rho_min}
     */
    public double rho_min = ADMMQPOption.ADMM_rho_min;

    /**
     * Max Value of changing the weight during iterations
     * {@link ADMMQPOption#ADMM_rho_maxIteration}
     */
    public int rho_maxIteration = ADMMQPOption.ADMM_rho_maxIteration;

    /**
     * Value for enabling the polish solution {@link ADMMQPOption#ADMM_polish}
     */
    public boolean polish = ADMMQPOption.ADMM_polish;
    /**
     * Value for Iteration of polish Algorithm
     * {@link ADMMQPOption#ADMM_polish_iteration}
     */
    public int polish_iteration = ADMMQPOption.ADMM_polish_iteration;

    public ADMMQPOption() {
    }

    public void setEps(double eps) {
        this.eps = eps;
    }

    public void setEpsInf(double eps) {
        this.epsInf = eps;
    }

    public void setSigma(double sigma) {
        this.sigma = sigma;
    }

    public void setAlfa(double alfa) {
        this.alfa = alfa;
    }

    public void setScale(boolean scale) {
        this.scale = scale;
    }

    public void setScaleMaxIteration(int iteration) {
        this.scale_maxIteration = iteration;
    }

    public void setRhoUpdate(boolean rhoUpdate) {
        this.rho_update = rhoUpdate;
    }

    public void setRhoMin(double rhoMin) {
        this.rho_min = rhoMin;
    }

    public void setRhoMax(double rhoMax) {
        this.rho_max = rhoMax;
    }

    public void setRhoMaxIteration(int iteration) {
        this.rho_maxIteration = iteration;
    }

    public void setPolish(boolean polish) {
        this.polish = polish;
    }

    public void setPolishIteration(int polish) {
        this.polish_iteration = polish;
    }

}
