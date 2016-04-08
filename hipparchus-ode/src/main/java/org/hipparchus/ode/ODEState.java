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

import java.io.Serializable;

/** Container for time, main and secondary state vectors.

 * @see FirstOrderDifferentialEquations
 * @see SecondaryEquations
 * @see FirstOrderIntegrator
 * @see ODEStateAndDerivative
 */

public class ODEState implements Serializable {

    /** Serializable UID. */
    private static final long serialVersionUID = 20160408L;

    /** Time. */
    private final double time;

    /** Main state at time. */
    private final double[] state;

    /** Secondary state at time. */
    private final double[][] secondaryState;

    /** Simple constructor.
     * <p>Calling this constructor is equivalent to call {@link
     * #ODEState(double, double[], double[][])
     * ODEState(time, state, null)}.</p>
     * @param time time
     * @param state state at time
     */
    public ODEState(double time, double[] state) {
        this(time, state, null);
    }

    /** Simple constructor.
     * @param time time
     * @param state state at time
     * @param secondaryState state at time (may be null)
     */
    public ODEState(double time, double[] state, double[][] secondaryState) {
        this.time           = time;
        this.state          = state.clone();
        this.secondaryState = copy(secondaryState);
    }

    /** Copy a two-dimensions array.
     * @param original original array (may be null)
     * @return copied array or null if original array was null
     */
    protected double[][] copy(final double[][] original) {

        // special handling of null arrays
        if (original == null) {
            return null;
        }

        // allocate the array
        final double[][] copied = new double[original.length][];

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

    /** Get main state dimension.
     * @return main state dimension
     */
    public int getStateDimension() {
        return state.length;
    }

    /** Get main state at time.
     * @return main state at time
     */
    public double[] getState() {
        return state.clone();
    }

    /** Get the number of secondary states.
     * @return number of secondary states.
     */
    public int getNumberOfSecondaryStates() {
        return secondaryState == null ? 0 : secondaryState.length;
    }

    /** Get secondary state dimension.
     * @param index index of the secondary set as returned
     * by {@link ExpandableODE#addSecondaryEquations(SecondaryEquations)}
     * (beware index 0 corresponds to main state, additional states start at 1)
     * @return secondary state dimension
     */
    public int getSecondaryStateDimension(final int index) {
        return index == 0 ? state.length : secondaryState[index - 1].length;
    }

    /** Get secondary state at time.
     * @param index index of the secondary set as returned
     * by {@link ExpandableODE#addSecondaryEquations(SecondaryEquations)}
     * (beware index 0 corresponds to main state, additional states start at 1)
     * @return secondary state at time
     */
    public double[] getSecondaryState(final int index) {
        return index == 0 ? state.clone() : secondaryState[index - 1].clone();
    }

    /** Return the dimension of the complete set of equations.
     * <p>
     * The complete set of equations correspond to the primary set plus all secondary sets.
     * </p>
     * @return dimension of the complete set of equations
     */
    public int getTotalDimension() {
        int dimension = state.length;
        if (secondaryState != null) {
            for (final double[] secondary : secondaryState) {
                dimension += secondary.length;
            }
        }
        return dimension;
    }

}
