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

import java.io.Serializable;

import org.hipparchus.complex.Complex;

/** Container for time, main and secondary state vectors.

 * @see ComplexOrdinaryDifferentialEquation
 * @see SecondaryODE
 * @see ODEIntegrator
 * @see ODEStateAndDerivative
 * @since 1.4
 */

public class ComplexODEState implements Serializable {

    /** Serializable UID. */
    private static final long serialVersionUID = 20180902;

    /** Time. */
    private final double time;

    /** Primary state at time. */
    private final Complex[] primaryState;

    /** Secondary state at time. */
    private final Complex[][] secondaryState;

    /** Complete dimension. */
    private final int completeDimension;

    /** Simple constructor.
     * <p>Calling this constructor is equivalent to call {@link
     * #ComplexODEState(double, Complex[], Complex[][])
     * ComplexODEState(time, state, null)}.</p>
     * @param time time
     * @param primaryState primary state at time
     */
    public ComplexODEState(double time, Complex[] primaryState) {
        this(time, primaryState, null);
    }

    /** Simple constructor.
     * @param time time
     * @param primaryState state at time
     * @param secondaryState primary state at time (may be null)
     */
    public ComplexODEState(double time, Complex[] primaryState, Complex[][] secondaryState) {

        this.time           = time;
        this.primaryState   = primaryState.clone();
        this.secondaryState = copy(secondaryState);

        // compute once and for all the complete dimension
        int dimension = primaryState.length;
        if (secondaryState != null) {
            for (final Complex[] secondary : secondaryState) {
                dimension += secondary.length;
            }
        }
        this.completeDimension = dimension;

    }

    /** Copy a two-dimensions array.
     * @param original original array (may be null)
     * @return copied array or null if original array was null
     */
    protected Complex[][] copy(final Complex[][] original) {

        // special handling of null arrays
        if (original == null) {
            return null; // NOPMD
        }

        // allocate the array
        final Complex[][] copied = new Complex[original.length][];

        // copy content
        for (int i = 0; i < original.length; ++i) {
            copied[i] = original[i].clone();
        }

        return copied;

    }

    /** Get time.
     * @return time
     */
    public double getTime() {
        return time;
    }

    /** Get primary state dimension.
     * @return primary state dimension
     * @see #getSecondaryStateDimension(int)
     * @see #getCompleteStateDimension()
     */
    public int getPrimaryStateDimension() {
        return primaryState.length;
    }

    /** Get primary state at time.
     * @return primary state at time
     * @see #getSecondaryState(int)
     * @see #getCompleteState()
     */
    public Complex[] getPrimaryState() {
        return primaryState.clone();
    }

    /** Get the number of secondary states.
     * @return number of secondary states.
     */
    public int getNumberOfSecondaryStates() {
        return secondaryState == null ? 0 : secondaryState.length;
    }

    /** Get secondary state dimension.
     * @param index index of the secondary set as returned
     * by {@link ExpandableODE#addSecondaryEquations(SecondaryODE)}
     * (beware index 0 corresponds to primary state, secondary states start at 1)
     * @return secondary state dimension
     * @see #getPrimaryStateDimension()
     * @see #getCompleteStateDimension()
     */
    public int getSecondaryStateDimension(final int index) {
        return index == 0 ? primaryState.length : secondaryState[index - 1].length;
    }

    /** Get secondary state at time.
     * @param index index of the secondary set as returned
     * by {@link ExpandableODE#addSecondaryEquations(SecondaryODE)}
     * (beware index 0 corresponds to primary state, secondary states start at 1)
     * @return secondary state at time
     * @see #getPrimaryState()
     * @see #getCompleteState()
     */
    public Complex[] getSecondaryState(final int index) {
        return index == 0 ? primaryState.clone() : secondaryState[index - 1].clone();
    }

    /** Return the dimension of the complete set of equations.
     * <p>
     * The complete set of equations correspond to the primary set plus all secondary sets.
     * </p>
     * @return dimension of the complete set of equations
     * @see #getPrimaryStateDimension()
     * @see #getSecondaryStateDimension(int)
     */
    public int getCompleteStateDimension() {
        return completeDimension;
    }

    /** Get complete state at time.
     * @return complete state at time, starting with
     * {@link #getPrimaryState() primary state}, followed
     * by all {@link #getSecondaryState(int) secondary states} in
     * increasing index order
     * @see #getPrimaryState()
     * @see #getSecondaryState(int)
     */
    public Complex[] getCompleteState() {
        final Complex[] completeState = new Complex[getCompleteStateDimension()];
        System.arraycopy(primaryState, 0, completeState, 0, primaryState.length);
        int offset = primaryState.length;
        if (secondaryState != null) {
            for (int index = 0; index < secondaryState.length; ++index) {
                System.arraycopy(secondaryState[index], 0,
                                 completeState, offset,
                                 secondaryState[index].length);
                offset += secondaryState[index].length;
            }
        }
        return completeState;
    }

}
