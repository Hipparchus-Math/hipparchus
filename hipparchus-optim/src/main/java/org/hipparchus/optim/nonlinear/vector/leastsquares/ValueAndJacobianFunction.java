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
package org.hipparchus.optim.nonlinear.vector.leastsquares;

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * A interface for functions that compute a vector of values and can compute their
 * derivatives (Jacobian).
 *
 */
public interface ValueAndJacobianFunction extends MultivariateJacobianFunction {
    /**
     * Compute the value.
     *
     * @param params Point.
     * @return the value at the given point.
     */
    RealVector computeValue(double[] params);

    /**
     * Compute the Jacobian.
     *
     * @param params Point.
     * @return the Jacobian at the given point.
     */
    RealMatrix computeJacobian(double[] params);
}
