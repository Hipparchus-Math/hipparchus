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
import org.hipparchus.optim.OptimizationData;

/**
 * Defines a nonlinear constrained optimization problem suitable for use
 * with Sequential Quadratic Programming (SQP) solvers.
 *
 * <p>
 * Implementations supply objective function information, constraint
 * definitions, bounds, and problem dimensionality. SQP solvers use
 * this interface to query the optimization model.
 * </p>
 *
 * <p>
 * Each component is optional: the problem may include only an objective
 * (unconstrained case), or objective plus bounds, or full equality and
 * inequality constraint sets. Boolean guard methods allow users and
 * solvers to detect which problem structures are active.
 * </p>
 *
 * <p>
 * Implementations must be deterministic and thread-safe if used across
 * multiple optimizer instances.
 * </p>
 *
 * @since 1.0
 * 
 */
public interface SQPProblem extends OptimizationData {

    /**
     * Returns the number of decision variables of the problem.
     *
     * @return problem dimension
     */
    int getdim();

    /**
     * Indicates whether this problem supplies an initial guess.
     * The optimizer may fallback to a default guess if false.
     *
     * @return true if {@link #getInitialGuess()} is defined
     */
    boolean hasInitialGuess();

    /**
     * Returns the initial guess vector for the optimization variables.
     *
     * @return initial guess values as a primitive array
     */
    double[] getInitialGuess();

    /**
     * Indicates whether explicit variable bounds are defined.
     *
     * @return true if box constraints are available
     */
    boolean hasBounds();

    /**
     * Returns lower bounds for each decision variable, if available.
     *
     * @return array of lower bound values (length equals problem dimension)
     */
    double[] getBoxConstraintLB();

    /**
     * Returns upper bounds for each decision variable, if available.
     *
     * @return array of upper bound values (length equals problem dimension)
     */
    double[] getBoxConstraintUB();

    /**
     * Evaluates the objective function at the given point.
     *
     * @param rv point at which to evaluate the objective
     * @return objective function value
     */
    double getObjectiveEvaluation(RealVector rv);

    /**
     * Evaluates the gradient of the objective function at the given point.
     *
     * @param rv point at which to evaluate the gradient
     * @return gradient vector
     */
    RealVector getObjectiveGradient(RealVector rv);

    /**
     * Indicates whether equality constraints are present.
     *
     * @return true if the problem contains equality constraints
     */
    boolean hasEquality();

    /**
     * Computes equality constraint values at the given point.
     *
     * @param rv evaluation point
     * @return vector of equality constraint residuals
     */
    RealVector getEqCostraintEvaluation(RealVector rv);

    /**
     * Computes the Jacobian of the equality constraints at the given point.
     *
     * @param rv evaluation point
     * @return Jacobian matrix of equality constraints
     */
    RealMatrix getEqCostraintJacobian(RealVector rv);

    /**
     * Returns the lower bound for equality constraints (typically zero).
     *
     * <p>
     * Some interior-point/SQP formulations treat equalities as bounded
     * two–sided constraints {@code h(x) = 0}.
     * </p>
     *
     * @return vector of equality constraint lower bounds
     */
    RealVector getEqCostraintLB();

    /**
     * Indicates whether inequality constraints are present.
     *
     * @return true if the problem contains inequality constraints
     */
    boolean hasInequality();

    /**
     * Computes inequality constraint values at the given point.
     *
     * <p>
     * Inequalities follow the canonical form {@code g(x) >= c}.
     * </p>
     *
     * @param rv evaluation point
     * @return vector of inequality constraint residuals
     */
    RealVector getIneqConstraintEvaluation(RealVector rv);

    /**
     * Computes the Jacobian of inequality constraints at the given point.
     *
     * @param rv evaluation point
     * @return Jacobian matrix of inequality constraints
     */
    RealMatrix getIneqCostraintJacobian(RealVector rv);

    /**
     * Returns the lower bound for inequality constraints.
     *
     * @return vector of lower bound values for inequality constraints
     */
    RealVector getIneqCostraintLB();
}
