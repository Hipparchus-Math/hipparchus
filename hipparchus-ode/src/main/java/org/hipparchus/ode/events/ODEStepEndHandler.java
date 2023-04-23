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

import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;

/** This interface represents a handler for discrete events triggered
 * during ODE integration at each step end.
 * @see org.hipparchus.ode.events
 * @since 3.0
 */
public interface ODEStepEndHandler  {

    /** Initialize step end handler at the start of an ODE integration.
     * <p>
     * This method is called once at the start of the integration. It
     * may be used by the step end handler to initialize some internal data
     * if needed.
     * </p>
     * <p>
     * The default implementation does nothing
     * </p>
     * @param initialState initial time, state vector and derivative
     * @param finalTime target time for the integration
     */
    default void init(ODEStateAndDerivative initialState, double finalTime) {
        // nothing by default
    }

    /** Handle an event and choose what to do next.

     * <p>This method is called when the integrator has accepted a step
     * ending exactly on step end, just <em>after</em>
     * the step handler itself is called (see below for scheduling). It
     * allows the user to update his internal data to acknowledge the fact
     * the event has been handled (for example setting a flag in the {@link
     * org.hipparchus.ode.OrdinaryDifferentialEquation
     * differential equations} to switch the derivatives computation in
     * case of discontinuity), or to direct the integrator to either stop
     * or continue integration, possibly with a reset state or derivatives.</p>
     *
     * <ul>
     *   <li>if {@link Action#STOP} is returned, the integration will be stopped,</li>
     *   <li>if {@link Action#RESET_STATE} is returned, the {@link #resetState
     *   resetState} method will be called once the step handler has
     *   finished its task, and the integrator will also recompute the
     *   derivatives,</li>
     *   <li>if {@link Action#RESET_DERIVATIVES} is returned, the integrator
     *   will recompute the derivatives,
     *   <li>if {@link Action#RESET_EVENTS} is returned, the integrator
     *   will recheck all event handlers,
     *   <li>if {@link Action#CONTINUE} is returned, no specific action will
     *   be taken (apart from having called this method) and integration
     *   will continue.</li>
     * </ul>
     *
     * <p>The scheduling between this method and the {@link
     * org.hipparchus.ode.sampling.ODEStepHandler ODEStepHandler} method {@link
     * org.hipparchus.ode.sampling.ODEStepHandler#handleStep(org.hipparchus.ode.sampling.ODEStateInterpolator)
     * handleStep(interpolator)} is to call {@code handleStep} first and this method afterwards.
     * This scheduling allows user code called by this method and user code called by step
     * handlers to get values of the independent time variable consistent with integration direction.</p>
     *
     * @param state current value of the independent <i>time</i> variable, state vector
     * and derivative at step end
     * @param forward if true, propagation is forward
     * @return indication of what the integrator should do next, this
     * value must be one of {@link Action#STOP}, {@link Action#RESET_STATE},
     * {@link Action#RESET_DERIVATIVES}, {@link Action#RESET_EVENTS}, or
     * {@link Action#CONTINUE}
     */
    Action stepEndOccurred(ODEStateAndDerivative state, boolean forward);

    /** Reset the state prior to continue the integration.
     *
     * <p>This method is called after the step handler has returned and
     * before the next step is started, but only when {@link
     * ODEEventHandler#eventOccurred(ODEStateAndDerivative, ODEEventDetector, boolean)}
     * has itself returned the {@link Action#RESET_STATE}
     * indicator. It allows the user to reset the state vector for the
     * next step, without perturbing the step handler of the finishing
     * step.</p>
     * <p>The default implementation returns its argument.</p>
     * @param state current value of the independent <i>time</i> variable, state vector
     * and derivative at step end
     * @return reset state (note that it does not include the derivatives, they will
     * be added automatically by the integrator afterwards)
     */
    default ODEState resetState(ODEStateAndDerivative state) {
        return state;
    }

}
