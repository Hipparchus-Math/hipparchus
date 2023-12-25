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

import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.OptimizationData;

/** Generic constraint.
 * @since 3.1
 */
public interface Constraint extends VectorDifferentiableFunction, OptimizationData {

    /** Get Lower Bound.
     * @return Lower Bound Vector
     */
    RealVector getLowerBound();

    /** Get Upper Bound.
     * @return Upper Bound Vector
     */
    RealVector getUpperBound();

}
