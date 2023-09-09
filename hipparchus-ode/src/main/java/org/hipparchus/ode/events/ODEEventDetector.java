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

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.solvers.BracketedUnivariateSolver;
import org.hipparchus.ode.ODEStateAndDerivative;

/** This interface represents a detector for discrete events triggered
 * during ODE integration.
 *
 * <p>Some events can be triggered at discrete times as an ODE problem
 * is solved. This occurs for example when the integration process
 * should be stopped as some state is reached (G-stop facility) when the
 * precise date is unknown a priori, or when the derivatives have
 * discontinuities, or simply when the user wants to monitor some
 * states boundaries crossings.
 * </p>
 *
 * <p>These events are defined as occurring when a <code>g</code>
 * switching function sign changes.</p>
 *
 * <p>Since events are only problem-dependent and are triggered by the
 * independent <i>time</i> variable and the state vector, they can
 * occur at virtually any time, unknown in advance. The integrators will
 * take care to avoid sign changes inside the steps, they will reduce
 * the step size when such an event is detected in order to put this
 * event exactly at the end of the current step. This guarantees that
 * step interpolation (which always has a one step scope) is relevant
 * even in presence of discontinuities. This is independent from the
 * stepsize control provided by integrators that monitor the local
 * error (this event handling feature is available for all integrators,
 * including fixed step ones).</p>
 *
 * <p>
 * Note that prior to Hipparchus 3.0, the methods in this interface were
 * in the {@link ODEEventHandler} interface and the defunct
 * {@code EventHandlerConfiguration} interface. The interfaces have been
 * reorganized to allow different objects to be used in event detection
 * and event handling, hence allowing users to reuse predefined events
 * detectors with custom handlers.
 * </p>
 *
 * @see org.hipparchus.ode.events
 * @since 3.0
 */
public interface ODEEventDetector  {

    /** Get the maximal time interval between events handler checks.
     * @return maximal time interval between events handler checks
     */
    AdaptableInterval getMaxCheckInterval();

    /** Get the upper limit in the iteration count for event localization.
     * @return upper limit in the iteration count for event localization
     */
    int getMaxIterationCount();

    /** Get the root-finding algorithm to use to detect state events.
     * @return root-finding algorithm to use to detect state events
     */
    BracketedUnivariateSolver<UnivariateFunction> getSolver();

    /** Get the underlying event handler.
     * @return underlying event handler
     */
    ODEEventHandler getHandler();

    /** Initialize event handler at the start of an ODE integration.
     * <p>
     * This method is called once at the start of the integration. It
     * may be used by the event handler to initialize some internal data
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

    /** Compute the value of the switching function.

     * <p>The discrete events are generated when the sign of this
     * switching function changes. The integrator will take care to change
     * the stepsize in such a way these events occur exactly at step boundaries.
     * The switching function must be continuous in its roots neighborhood
     * (but not necessarily smooth), as the integrator will need to find its
     * roots to locate precisely the events.</p>
     *
     * <p>Also note that for the integrator to detect an event the sign of the switching
     * function must have opposite signs just before and after the event. If this
     * consistency is not preserved the integrator may not detect any events.
     *
     * <p>This need for consistency is sometimes tricky to achieve. A typical
     * example is using an event to model a ball bouncing on the floor. The first
     * idea to represent this would be to have {@code g(state) = h(state)} where h is the
     * height above the floor at time {@code state.getTime()}. When {@code g(state)}
     * reaches 0, the ball is on the floor, so it should bounce and the typical way to do this is
     * to reverse its vertical velocity. However, this would mean that before the
     * event {@code g(state)} was decreasing from positive values to 0, and after the
     * event {@code g(state)} would be increasing from 0 to positive values again.
     * Consistency is broken here! The solution here is to have {@code g(state) = sign
     * * h(state)}, where sign is a variable with initial value set to {@code +1}. Each
     * time {@link ODEEventHandler#eventOccurred(ODEStateAndDerivative,
     * ODEEventDetector, boolean) eventOccurred} is called,
     * {@code sign} is reset to {@code -sign}. This allows the {@code g(state)}
     * function to remain continuous (and even smooth) even across events, despite
     * {@code h(state)} is not. Basically, the event is used to <em>fold</em> {@code h(state)}
     * at bounce points, and {@code sign} is used to <em>unfold</em> it back, so the
     * solvers sees a {@code g(state)} function which behaves smoothly even across events.</p>
     *
     * <p>This method is idempotent, that is calling this multiple times with the same
     * state will result in the same value, with two exceptions. First, the definition of
     * the g function may change when an {@link ODEEventHandler#eventOccurred(ODEStateAndDerivative,
     * ODEEventDetector, boolean) event occurs} on the handler, as in the above example.
     * Second, the definition of the g function may change when the {@link
     * ODEEventHandler#eventOccurred(ODEStateAndDerivative, ODEEventDetector, boolean) eventOccurred}
     * method of any other event handler in the same integrator returns {@link Action#RESET_EVENTS},
     * {@link Action#RESET_DERIVATIVES}, or {@link Action#RESET_STATE}.
     *
     * @param state current value of the independent <i>time</i> variable, state vector
     * and derivative
     * @return value of the g switching function
     * @see org.hipparchus.ode.events
     */
    double g(ODEStateAndDerivative state);

}
