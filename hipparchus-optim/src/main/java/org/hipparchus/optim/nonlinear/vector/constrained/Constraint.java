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

    /** Get Lower Bound for {@link #value(RealVector) value(x)}.
     * @return Lower Bound for {@link #value(RealVector) value(x)}
     */
    RealVector getLowerBound();

    /** Get Upper Bound for {@link #value(RealVector) value(x)}.
     * @return Upper Bound for {@link #value(RealVector) value(x)}
     */
    RealVector getUpperBound();

    /** Check how much a point overshoots the constraint.
     * <p>
     * The overshoot is zero if the point fulfills the constraint, and
     * positive if the {@link #value(RealVector) value} of the constraint
     * is on the wrong side of {@link #getLowerBound() lower} or {@link
     * #getUpperBound() upper} boundaries.
     * </p>
     * @param y constraint value (y = {@link #value(RealVector) value}(x))
     * @return L¹-norm of constraint overshoot
     */
    double overshoot(RealVector y);

}
