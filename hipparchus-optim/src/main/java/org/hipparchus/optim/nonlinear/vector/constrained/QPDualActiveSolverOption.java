/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 * https://www.apache.org/licenses/LICENSE-2.0
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

/** Options for dual active-set method.
 */
public class QPDualActiveSolverOption implements OptimizationData  {

    /** Default tolerance for constraint satisfaction */
    public static final double DEFAULT_EPSILON = Precision.EPSILON;//>0

     /** Default Matrix Mode */
    public static final QPMatrixMode DEFAULT_MATRIX_MODE = QPMatrixMode.FULL;

    /**
    *  Default Relaxed feasibility tolerance used only when both the primal and dual
     * step are not available. In that numerically stalled case, a target
     * constraint whose violation does not exceed this value is accepted as
     * numerically feasible.
     */
    public static final double DEFAULT_EPSILON_RELAXATION = 1.0e-12;//>0


    /** Tolerance for constraint satisfaction */
    private double eps;

    /**
     * Relaxed feasibility tolerance used only when both the primal and dual
     * step are not available for ill conditioning problems. In that numerically stalled case, a target
     * constraint whose violation does not exceed this value is accepted as
     * numerically feasible. For the major part of the problems choosing 10.0*eps is enough
     */
    private double epsRelaxation;

    /** Matrix Mode */
    private QPMatrixMode matrixMode;


    public QPDualActiveSolverOption()
    {
        this.eps=DEFAULT_EPSILON;
        this.epsRelaxation=DEFAULT_EPSILON_RELAXATION;
        this.matrixMode=DEFAULT_MATRIX_MODE;
    }

    /** Set tolerance for convergence
     *  typical 1.0e-16 to 1.0e-9
     * @param eps tolerance for convergence
     */
    public void setEps(final double eps) {
        this.eps = eps;
    }

    /** get tolerance for convergence and active constraint evaluation.
     * @return tolerance for convergence and active constraint evaluation
     */
    public double getEps() {
        return eps;
    }

     /** Set tolerance convergence for constraint in case of degeneracy
     *  typical 10.0*eps
     * @param epsRelaxation  tolerance for convergence
     */
    public void setEpsRelaxation(final double epsRelaxation) {
        this.epsRelaxation = epsRelaxation;
    }

    /** get tolerance for convergence for constraint in case of degeneracy.
     * @return tolerance for convergence for constraint in case of degeneracy
     */
    public double getEpsRelaxation() {
        return epsRelaxation;
    }

     /** Set MatrixMode.
      *
      * @param matrixMode  matrix mode
      */
    public void setMatrixMode(final QPMatrixMode matrixMode)
    {
        this.matrixMode=matrixMode;
    }

    /** Get MatrixMode
     * @return matric mode
     */
     public QPMatrixMode getMatrixMode()
    {
      return  this.matrixMode;
    }
}
