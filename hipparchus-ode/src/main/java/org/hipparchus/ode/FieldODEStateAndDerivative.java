/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.hipparchus.ode;

import org.hipparchus.RealFieldElement;

/** Container for time, main and secondary state vectors as well as their derivatives.

 * @see FieldOrdinaryDifferentialEquation
 * @see FieldSecondaryODE
 * @see FieldODEIntegrator
 * @param <T> the type of the field elements
 */

public class FieldODEStateAndDerivative<T extends RealFieldElement<T>> extends FieldODEState<T> {

    /** Derivative of the primary state at time. */
    private final T[] primaryDerivative;

    /** Derivative of the secondary state at time. */
    private final T[][] secondaryDerivative;

    /** Simple constructor.
     * <p>Calling this constructor is equivalent to call {@link
     * #FieldODEStateAndDerivative(RealFieldElement, RealFieldElement[], RealFieldElement[],
     * RealFieldElement[][], RealFieldElement[][]) FieldODEStateAndDerivative(time, state,
     * derivative, null, null)}.</p>
     * @param time time
     * @param primaryState primary state at time
     * @param primaryDerivative derivative of the primary state at time
     */
    public FieldODEStateAndDerivative(T time, T[] primaryState, T[] primaryDerivative) {
        this(time, primaryState, primaryDerivative, null, null);
    }

    /** Simple constructor.
     * @param time time
     * @param primaryState primary state at time
     * @param primaryDerivative derivative of the primary state at time
     * @param secondaryState state at time (may be null)
     * @param secondaryDerivative derivative of the state at time (may be null)
     */
    public FieldODEStateAndDerivative(T time, T[] primaryState, T[] primaryDerivative,
                                      T[][] secondaryState, T[][] secondaryDerivative) {
        super(time, primaryState, secondaryState);
        this.primaryDerivative   = primaryDerivative.clone();
        this.secondaryDerivative = copy(time.getField(), secondaryDerivative);
    }

    /** Get derivative of the primary state at time.
     * @return derivative of the primary state at time
     */
    public T[] getPrimaryDerivative() {
        return primaryDerivative.clone();
    }

    /** Get derivative of the secondary state at time.
     * @param index index of the secondary set as returned
     * by {@link FieldExpandableODE#addSecondaryEquations(FieldSecondaryODE)}
     * (beware index 0 corresponds to primary state, secondary states start at 1)
     * @return derivative of the secondary state at time
     */
    public T[] getSecondaryDerivative(final int index) {
        return index == 0 ? primaryDerivative.clone() : secondaryDerivative[index - 1].clone();
    }

}
