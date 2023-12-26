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

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.ConvergenceChecker;
import org.hipparchus.optim.OptimizationData;
import org.hipparchus.util.FastMath;

/** Convergence Checker for ADMM QP Optimizer.
 * @since 3.1
 */
public class ADMMQPConvergenceChecker implements ConvergenceChecker<LagrangeSolution>, OptimizationData  {

    /** Quadratic term matrix. */
    private  final RealMatrix h;

    /** Constraint coefficients matrix. */
    private  final RealMatrix a;

    /** Linear term matrix. */
    private  final RealVector q;

    /** Absolute tolerance for convergence. */
    private  final double epsAbs;

    /** Relative tolerance for convergence. */
    private  final double epsRel;

    /** Convergence indicator. */
    private  boolean converged;

    /** Simple constructor.
     * @param h quadratic term matrix
     * @param a constraint coefficients matrix
     * @param q linear term matrix
     * @param epsAbs
     * @param epsRel
     */
    ADMMQPConvergenceChecker(final RealMatrix h, final RealMatrix a, final RealVector q,
                             final double epsAbs, final double epsRel) {
        this.h         = h;
        this.a         = a;
        this.q         = q;
        this.epsAbs    = epsAbs;
        this.epsRel    = epsRel;
        this.converged = false;
    }

    /** {@inheritDoc} */
    @Override
    public boolean converged(final int i, final LagrangeSolution previous, final LagrangeSolution current) {
        return converged;
    }

    /** Evaluate convergence.
     * @param rp primal residual
     * @param rd dual residual
     * @param maxPrimal primal vectors max
     * @param maxDual dual vectors max
     * @return true of convergence has been reached
     */
    public boolean converged(final double rp, final double rd, final double maxPrimal, final double maxDual) {
        boolean result = false;

        if (rp <= epsPrimalDual(maxPrimal) && rd <= epsPrimalDual(maxDual)) {
            result = true;
            converged = true;
        }
        return result;
    }

    /** Compute primal residual.
     * @param x primal problem solution
     * @param z auxiliary variable
     * @return primal residual
     */
    public double residualPrime(final RealVector x, final RealVector z) {
        return a.operate(x).subtract(z).getLInfNorm();
    }

    /** Compute dual residual.
     * @param x primal problem solution
     * @param y dual problem solution
     * @return dual residual
     */
    public double residualDual(final RealVector x, final RealVector y) {
        return q.add(a.transpose().operate(y)).add(h.operate(x)).getLInfNorm();
    }

    /** Compute primal vectors max.
     * @param x primal problem solution
     * @param z auxiliary variable
     * @return primal vectors max
     */
    public double maxPrimal(final RealVector x, final RealVector z) {
        return FastMath.max(a.operate(x).getLInfNorm(), z.getLInfNorm());
    }

    /** Compute dual vectors max.
     * @param x primal problem solution
     * @param y dual problem solution
     * @return dual vectors max
     */
    public double maxDual(final RealVector x, final RealVector y) {
        return FastMath.max(FastMath.max(h.operate(x).getLInfNorm(),
                                         a.transpose().operate(y).getLInfNorm()),
                            q.getLInfNorm());
    }

    /** Combine absolute and relative tolerances.
     * @param maxPrimalDual either {@link #maxPrimal(RealVector, RealVector)}
     * or {@link #maxDual(RealVector, RealVector)}
     * @return global tolerance
     */
    private double epsPrimalDual(final double maxPrimalDual) {
        return epsAbs + epsRel * maxPrimalDual;
    }

}
