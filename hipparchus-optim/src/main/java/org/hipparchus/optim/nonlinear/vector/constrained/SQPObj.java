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

/**
 * Adapter class that exposes an {@link SQPProblem} objective function under the
 * {@link TwiceDifferentiableFunction} abstraction used by SQP optimizers.
 *
 * <p>
 * This wrapper delegates objective evaluations and gradient computations to
 * the underlying {@link SQPProblem} implementation. The Hessian is currently
 * not provided and returns {@code null}, allowing SQP solvers to supply their
 * own Hessian approximation (e.g. via BFGS or other quasi-Newton strategies).
 * </p>
 *
 * <p>
 * Instances of this class are typically constructed internally by SQP
 * optimizers when processing {@link org.hipparchus.optim.OptimizationData}
 * supplied by users.
 * </p>
 *
 * @since 1.0
 * @author rocca
 */
public class SQPObj extends TwiceDifferentiableFunction {

    /** Optimization problem providing objective and gradient information. */
    private final SQPProblem problem;

    /**
     * Creates an adapter exposing the objective of an {@link SQPProblem}
     * instance through the {@link TwiceDifferentiableFunction} API.
     *
     * @param problem SQP formulation from which objective data are obtained
     */
    public SQPObj(final SQPProblem problem) {
        this.problem = problem;
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * Delegates to {@link SQPProblem#getdim()}.
     * </p>
     */
    @Override
    public int dim() {
        return this.problem.getdim();
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * Delegates to {@link SQPProblem#getObjectiveEvaluation(RealVector)}.
     * </p>
     */
    @Override
    public double value(final RealVector rv) {
        return this.problem.getObjectiveEvaluation(rv);
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * Delegates to {@link SQPProblem#getObjectiveGradient(RealVector)}.
     * </p>
     */
    @Override
    public RealVector gradient(final RealVector rv) {
        return this.problem.getObjectiveGradient(rv);
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * Returns {@code null} because {@link SQPProblem} does not define
     * second-order information. SQP methods typically compute Hessian
     * approximations internally (e.g. BFGS, SR1, etc.).
     * </p>
     *
     * @return always {@code null}
     */
    @Override
    public RealMatrix hessian(final RealVector rv) {
        return null;
    }
}
