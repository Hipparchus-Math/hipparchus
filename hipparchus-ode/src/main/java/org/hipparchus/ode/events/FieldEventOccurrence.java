/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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

package org.hipparchus.ode.events;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.ode.FieldODEState;

/**
 * Class to hold the data related to an event occurrence that is needed to decide how
 * to modify integration.
 * @param <T> the type of the field elements
 */
public class FieldEventOccurrence<T extends CalculusFieldElement<T>> {

    /** User requested action. */
    private final Action action;
    /** New state for a reset action. */
    private final FieldODEState<T> newState;
    /** The time to stop propagation if the action is a stop event. */
    private final T stopTime;

    /**
     * Create a new occurrence of an event.
     *
     * @param action   the user requested action.
     * @param newState for a reset event. Should be the current state unless the
     *                 action is {@link Action#RESET_STATE}.
     * @param stopTime to stop propagation if the action is {@link Action#STOP}. Used
     *                 to move the stop time to just after the root.
     */
    public FieldEventOccurrence(final Action action,
                                final FieldODEState<T> newState,
                                final T stopTime) {
        this.action = action;
        this.newState = newState;
        this.stopTime = stopTime;
    }

    /**
     * Get the user requested action.
     *
     * @return the action.
     */
    public Action getAction() {
        return action;
    }

    /**
     * Get the new state for a reset action.
     *
     * @return the new state.
     */
    public FieldODEState<T> getNewState() {
        return newState;
    }

    /**
     * Get the new time for a stop action.
     *
     * @return when to stop propagation.
     */
    public T getStopTime() {
        return stopTime;
    }

}
