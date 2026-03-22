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
 * Adapter exposing the equality constraint set of an {@link SQPProblem}
 * through the {@link EqualityConstraint} abstraction used by SQP solvers.
 *
 * <p>
 * This wrapper delegates equality constraint evaluation and Jacobian
 * computation to the underlying {@link SQPProblem}. The dimensionality
 * of the constraint mapping and the associated lower bounds are deduced
 * directly from the problem instance.
 * </p>
 *
 * <p>
 * Equality constraints typically model conditions of the form:
 * <pre>
 *     h(x) = 0
 * </pre>
 * but some SQP formulations internally interpret equalities through
 * a bounded form (e.g. {@code h(x) >= lb} and {@code h(x) <= ub}).
 * Here, only the lower bound is provided via
 * {@link SQPProblem#getEqCostraintLB()}, consistent with the
 * {@link EqualityConstraint} superclass.
 * </p>
 *
 * <p>
 * Instances of this class are commonly created internally by SQP
 * optimizers when parsing user-supplied {@link org.hipparchus.optim.OptimizationData}.
 * It acts as an adapter layer between domain models and solver APIs.
 * </p>
 *
 * @since 1.0
 * @author rocca
 */
public class SQPEq extends EqualityConstraint {

    /** Reference to the originating SQP model. */
    private final SQPProblem problem;

    /**
     * Constructs an equality constraint wrapper around the supplied
     * {@link SQPProblem}.
     *
     * @param problem optimization model providing equality constraint
     *                functions and gradients
     */
    public SQPEq(final SQPProblem problem) {
        super(problem.getEqCostraintLB());
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
     * Returns the number of equality constraints. If no lower bound vector
     * is provided, the set is assumed empty.
     * </p>
     */
    @Override
    public int dimY() {
        return (problem.getEqCostraintLB() != null)
                ? problem.getEqCostraintLB().getDimension()
                : 0;
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * Delegates evaluation to
     * {@link SQPProblem#getEqCostraintEvaluation(RealVector)}.
     * </p>
     */
    @Override
    public RealVector value(final RealVector rv) {
        return problem.getEqCostraintEvaluation(rv);
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * Delegates to
     * {@link SQPProblem#getEqCostraintJacobian(RealVector)}.
     * </p>
     */
    @Override
    public RealMatrix jacobian(final RealVector rv) {
        return problem.getEqCostraintJacobian(rv);
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * Convenience overload evaluating equality constraints for a raw
     * primitive array. Converts the input to {@link ArrayRealVector}
     * before delegation.
     * </p>
     *
     * @throws IllegalArgumentException if vector size mismatch occurs
     */
    @Override
    public double[] value(final double[] doubles) throws IllegalArgumentException {
        return problem.getEqCostraintEvaluation(new ArrayRealVector(doubles)).toArray();
    }
}
