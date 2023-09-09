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

package org.hipparchus.ode.events;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.sampling.FieldODEStateInterpolator;

/** This interface handles the state for either one {@link FieldODEEventHandler
 * event handler} or one {@link FieldODEStepEndHandler step end handler}
 * during integration steps.
 * @since 3.0
 * @param <T> the type of the field elements
 */
public interface FieldEventState<T extends CalculusFieldElement<T>> {

    /** Initialize handler at the start of an integration.
     * <p>
     * This method is called once at the start of the integration. It
     * may be used by the event handler to initialize some internal data
     * if needed.
     * </p>
     * @param s0 initial state
     * @param t target time for the integration
     *
     */
    void init(FieldODEStateAndDerivative<T> s0, T t);

    /** Get the occurrence time of the event triggered in the current step.
     * @return occurrence time of the event triggered in the current
     * step or infinity if no events are triggered
     */
    T getEventTime();

    /** Evaluate the impact of the proposed step on the event handler.
     * @param interpolator step interpolator for the proposed step
     * @return true if the event handler triggers an event before
     * the end of the proposed step
     * @exception MathIllegalStateException if the interpolator throws one because
     * the number of functions evaluations is exceeded
     * @exception MathIllegalArgumentException if the event cannot be bracketed
     */
    boolean evaluateStep(FieldODEStateInterpolator<T> interpolator)
        throws MathIllegalArgumentException, MathIllegalStateException;

    /**
     * Notify the user's listener of the event. The event occurs wholly within this method
     * call including a call to {@link FieldODEEventHandler#resetState(FieldODEEventDetector,
     * FieldODEStateAndDerivative)} if necessary.
     *
     * @param state the state at the time of the event. This must be at the same time as
     *              the current value of {@link #getEventTime()}.
     * @return the user's requested action and the new state if the action is {@link
     * Action#RESET_STATE}. Otherwise the new state is {@code state}. The stop time
     * indicates what time propagation should stop if the action is {@link Action#STOP}.
     * This guarantees the integration will stop on or after the root, so that integration
     * may be restarted safely.
     */
    FieldEventOccurrence<T> doEvent(FieldODEStateAndDerivative<T> state);

}
