/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.optim.nonlinear.scalar;

import org.hipparchus.analysis.MultivariateVectorFunction;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.optim.ConvergenceChecker;
import org.hipparchus.optim.OptimizationData;
import org.hipparchus.optim.PointValuePair;

/**
 * Base class for implementing optimizers for multivariate scalar
 * differentiable functions.
 * It contains boiler-plate code for dealing with gradient evaluation.
 *
 */
public abstract class GradientMultivariateOptimizer
    extends MultivariateOptimizer {
    /**
     * Gradient of the objective function.
     */
    private MultivariateVectorFunction gradient;

    /** Simple constructor.
     * @param checker Convergence checker.
     */
    protected GradientMultivariateOptimizer(ConvergenceChecker<PointValuePair> checker) {
        super(checker);
    }

    /**
     * Compute the gradient vector.
     *
     * @param params Point at which the gradient must be evaluated.
     * @return the gradient at the specified point.
     */
    protected double[] computeObjectiveGradient(final double[] params) {
        return gradient.value(params);
    }

    /**
     * {@inheritDoc}
     *
     * @param optData Optimization data. In addition to those documented in
     * {@link MultivariateOptimizer#parseOptimizationData(OptimizationData[])
     * MultivariateOptimizer}, this method will register the following data:
     * <ul>
     *  <li>{@link ObjectiveFunctionGradient}</li>
     * </ul>
     * @return {@inheritDoc}
     * @throws MathIllegalStateException if the maximal number of
     * evaluations (of the objective function) is exceeded.
     */
    @Override
    public PointValuePair optimize(OptimizationData... optData)
        throws MathIllegalStateException {
        // Set up base class and perform computation.
        return super.optimize(optData);
    }

    /**
     * Scans the list of (required and optional) optimization data that
     * characterize the problem.
     *
     * @param optData Optimization data.
     * The following data will be looked for:
     * <ul>
     *  <li>{@link ObjectiveFunctionGradient}</li>
     * </ul>
     */
    @Override
    protected void parseOptimizationData(OptimizationData... optData) {
        // Allow base class to register its own data.
        super.parseOptimizationData(optData);

        // The existing values (as set by the previous call) are reused if
        // not provided in the argument list.
        for (OptimizationData data : optData) {
            if  (data instanceof ObjectiveFunctionGradient) {
                gradient = ((ObjectiveFunctionGradient) data).getObjectiveFunctionGradient();
                // If more data must be parsed, this statement _must_ be
                // changed to "continue".
                break;
            }
        }
    }
}
