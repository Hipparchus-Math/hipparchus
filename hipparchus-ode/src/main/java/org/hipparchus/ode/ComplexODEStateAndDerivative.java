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

package org.hipparchus.ode;

import org.hipparchus.complex.Complex;

/** Container for time, main and secondary state vectors as well as their derivatives.

 * @see ComplexOrdinaryDifferentialEquation
 * @see ComplexSecondaryODE
 * @see ODEIntegrator
 */

public class ComplexODEStateAndDerivative extends ComplexODEState {

    /** Serializable UID. */
    private static final long serialVersionUID = 20180902L;

    /** Derivative of the primary state at time. */
    private final Complex[] primaryDerivative;

    /** Derivative of the secondary state at time. */
    private final Complex[][] secondaryDerivative;

    /** Simple constructor.
     * <p>Calling this constructor is equivalent to call {@link
     * #ComplexODEStateAndDerivative(double, Complex[], Complex[],
     * Complex[][], Complex[][]) ComplexODEStateAndDerivative(time, state,
     * derivative, null, null)}.</p>
     * @param time time
     * @param primaryState primary state at time
     * @param primaryDerivative derivative of the primary state at time
     */
    public ComplexODEStateAndDerivative(double time, Complex[] primaryState, Complex[] primaryDerivative) {
        this(time, primaryState, primaryDerivative, null, null);
    }

    /** Simple constructor.
     * @param time time
     * @param primaryState primary state at time
     * @param primaryDerivative derivative of the primary state at time
     * @param secondaryState state at time (may be null)
     * @param secondaryDerivative derivative of the state at time (may be null)
     */
    public ComplexODEStateAndDerivative(double time, Complex[] primaryState, Complex[] primaryDerivative,
                                        Complex[][] secondaryState, Complex[][] secondaryDerivative) {
        super(time, primaryState, secondaryState);
        this.primaryDerivative   = primaryDerivative.clone();
        this.secondaryDerivative = copy(secondaryDerivative);
    }

    /** Get derivative of the primary state at time.
     * @return derivative of the primary state at time
     * @see #getSecondaryDerivative(int)
     * @see #getCompleteDerivative()
     */
    public Complex[] getPrimaryDerivative() {
        return primaryDerivative.clone();
    }

    /** Get derivative of the secondary state at time.
     * @param index index of the secondary set as returned
     * by {@link ExpandableODE#addSecondaryEquations(SecondaryODE)}
     * (beware index 0 corresponds to primary state, secondary states start at 1)
     * @return derivative of the secondary state at time
     * @see #getPrimaryDerivative()
     * @see #getCompleteDerivative()
     */
    public Complex[] getSecondaryDerivative(final int index) {
        return index == 0 ? primaryDerivative.clone() : secondaryDerivative[index - 1].clone();
    }

    /** Get complete derivative at time.
     * @return complete derivative at time, starting with
     * {@link #getPrimaryDerivative() primary derivative}, followed
     * by all {@link #getSecondaryDerivative(int) secondary derivatives} in
     * increasing index order
     * @see #getPrimaryDerivative()
     * @see #getSecondaryDerivative(int)
     */
    public Complex[] getCompleteDerivative() {
        final Complex[] completeDerivative = new Complex[getCompleteStateDimension()];
        System.arraycopy(primaryDerivative, 0, completeDerivative, 0, primaryDerivative.length);
        int offset = primaryDerivative.length;
        if (secondaryDerivative != null) {
            for (int index = 0; index < secondaryDerivative.length; ++index) {
                System.arraycopy(secondaryDerivative[index], 0,
                                 completeDerivative, offset,
                                 secondaryDerivative[index].length);
                offset += secondaryDerivative[index].length;
            }
        }
        return completeDerivative;
    }

}
