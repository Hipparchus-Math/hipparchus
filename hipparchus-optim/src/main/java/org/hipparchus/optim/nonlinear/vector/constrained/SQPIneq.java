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

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * Adapter representing the inequality constraint set of an {@link SQPProblem}
 * as an {@link InequalityConstraint} suitable for use by SQP solvers.
 *
 * <p>
 * This wrapper delegates the evaluation and Jacobian computation of inequality
 * constraints back to the underlying {@link SQPProblem}. The dimension of the
 * constraint mapping and lower bounds are derived directly from the problem.
 * </p>
 *
 * <p>
 * Inequality constraints are assumed to follow the canonical form:
 * </p>
 * <pre>
 *     g(x) >= lb
 * </pre>
 * <p>
 * where {@code lb} is given by {@link SQPProblem#getIneqConstraintLB()}.
 * </p>
 *
 * <p>
 * An instance of this class is normally constructed automatically by SQP
 * optimizers when processing {@link org.hipparchus.optim.OptimizationData}. It serves as a
 * structural adapter layer between user-supplied models and solver APIs.
 * </p>
 *
 * @since 5.0
 */
public class SQPIneq extends InequalityConstraint {

    /** Reference to the originating SQP problem definition. */
    private final SQPProblem problem;

    /**
     * Creates an inequality constraint adapter for the given {@link SQPProblem}.
     *
     * @param problem the problem model supplying inequality constraint values
     *                and their Jacobian
     */
    public SQPIneq(final SQPProblem problem) {
        super(problem.getIneqConstraintLB());
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
        return problem.getdim();
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * Returns the number of inequality constraint components. If no lower
     * bounds are supplied, the mapping is assumed empty.
     * </p>
     */
    @Override
    public int dimY() {
        return (problem.getIneqConstraintLB() != null) ?
                problem.getIneqConstraintLB().getDimension() :
                0;
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * Delegates evaluation to
     * {@link SQPProblem#getIneqConstraintEvaluation(RealVector)}.
     * </p>
     */
    @Override
    public RealVector value(final RealVector rv) {
        return problem.getIneqConstraintEvaluation(rv);
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * Delegates to
     * {@link SQPProblem#getIneqConstraintJacobian(RealVector)}.
     * </p>
     */
    @Override
    public RealMatrix jacobian(final RealVector rv) {
        return problem.getIneqConstraintJacobian(rv);
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * Convenience overload that evaluates constraints given a raw array input.
     * Converts the input to {@link ArrayRealVector} and delegates to
     * {@link SQPProblem#getIneqConstraintEvaluation(RealVector)}.
     * </p>
     *
     */
    @Override
    public double[] value(final double[] doubles) {
        return problem.getIneqConstraintEvaluation(new ArrayRealVector(doubles)).toArray();
    }
}
