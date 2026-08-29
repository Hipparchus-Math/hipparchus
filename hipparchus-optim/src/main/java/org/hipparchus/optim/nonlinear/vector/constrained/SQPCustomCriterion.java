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

/**
 * Custom convergence criterion for an SQP optimizer.
 */

public interface SQPCustomCriterion extends OptimizationData {

    /**
     * Checks whether the SQP iteration has converged.
     *
     * @param alpha accepted step length
     * @param dxNorm norm of the primal step
     * @param dxHdx quadratic step measure
     * @param complSlack complementary-slackness residual
     * @param kkt KKT residual
     * @param viol constraint violation
     * @param funDiff objective-function difference
     * @param xNorm variables norm
     * @param funNorm objective-function norm
     * @return {@code true} if convergence has been reached
     */
    boolean converged(double alpha,
                      double dxNorm,
                      double dxHdx,
                      double complSlack,
                      double kkt,
                      double viol,
                      double funDiff,
                      double xNorm,
                      double funNorm);
}
