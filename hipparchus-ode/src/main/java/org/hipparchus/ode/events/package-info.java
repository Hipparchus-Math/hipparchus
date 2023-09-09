/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
/**
 *
 * <h2>Events</h2>
 *
 * <p>
 * This package provides classes to handle discrete events occurring during
 * Ordinary Differential Equations integration.
 * </p>
 *
 * <p>
 * Discrete events detection is based on switching functions. The user provides
 * a simple {@link org.hipparchus.ode.events.ODEEventDetector#g(org.hipparchus.ode.ODEStateAndDerivative)
 * g(state)} function depending on the current time and state in a detector object and a
 * {@link org.hipparchus.ode.events.ODEEventHandler#eventOccurred(org.hipparchus.ode.ODEStateAndDerivative,
 * ODEEventDetector, boolean) handler function} in a handler object. The integrator will monitor
 * the value of the function throughout integration range and will trigger the
 * event when its sign changes. The magnitude of the value is almost irrelevant,
 * it should however be continuous (but not necessarily smooth) for the sake of
 * root finding.
 * </p>
 *
 * <p>
 * Events detection is based on two embedded search algorithms. A top level loop,
 * driven by the maxCheck interval, samples the g function during integration steps
 * to identify sign changes and therefore bracket roots. When a sign change has been
 * detected between two successive samples, then a lower level root finding
 * algorithm is triggered to precisely locate the root, using the provided threshold
 * as a convergence criterion. The maxCheck interval should be large enough to avoid
 * sampling too finely the integration step and evaluating the g function too often,
 * but it should be small enough to still separate successive roots, otherwise some
 * events may be missed. See below for a formal definition of the maxCheck behavior.
 * The threshold on the other hand should generally be very small, it really corresponds
 * to the accuracy at which events must be located.
 * </p>
 *
 * <p>
 * When an event is triggered, several different options are available:
 * </p>
 * <ul>
 *  <li>integration can be stopped (this is called a G-stop facility),</li>
 *  <li>the state vector or the derivatives can be changed,</li>
 *  <li>or integration can simply go on.</li>
 * </ul>
 *
 * <p>
 * The first case, G-stop, is the most common one. A typical use case is when an
 * ODE must be solved up to some target state is reached, with a known value of
 * the state but an unknown occurrence time. As an example, if we want to monitor
 * a chemical reaction up to some predefined concentration for the first substance,
 * we can use the following switching function setting in the detector:</p>
 * <pre>
 *  public double g(final ODEStateAndDerivative state) {
 *    return state.getState()[0] - targetConcentration;
 *  }
 * </pre>
 * <p>and the following setting in the event handler:</p>
 * <pre>
 *  public Action eventOccurred(final ODEStateAndDerivative state, final ODEEventDetector detector, final boolean increasing) {
 *    return STOP;
 *  }
 * </pre>
 *
 * <p>
 * The second case, change state vector or derivatives is encountered when dealing
 * with discontinuous dynamical models. A typical case would be the motion of a
 * spacecraft when thrusters are fired for orbital maneuvers. The acceleration is
 * smooth as long as no maneuver are performed, depending only on gravity, drag,
 * third body attraction, radiation pressure. Firing a thruster introduces a
 * discontinuity that must be handled appropriately by the integrator. In such a case,
 * we would use a switching function setting similar to this:</p>
 * <pre>
 *  public double g(final ODEStateAndDerivative state) {
 *    return (state.getTime() - tManeuverStart) &lowast; (state.getTime() - tManeuverStop);
 *  }
 * </pre>
 * <p>and the following setting in the event handler:</p>
 * <pre>
 *  public Action eventOccurred(final ODEStateAndDerivative state, final ODEEventDetector detector, final boolean increasing) {
 *    return RESET_DERIVATIVES;
 *  }
 * </pre>
 *
 * <p>
 * The third case is useful mainly for monitoring purposes, a simple example is:</p>
 * <pre>
 *  public double g(final ODEStateAndDerivative state) {
 *  final double[] y = state.getState();
 *    return y[0] - y[1];
 *  }
 * </pre>
 * <p>and the following setting in the event handler:</p>
 * <pre>
 *  public Action eventOccurred(final ODEStateAndDerivative state, final ODEEventDetector detector, final boolean increasing) {
 *    logger.log("y0(t) and y1(t) curves cross at t = " + t);
 *    return CONTINUE;
 *  }
 * </pre>
 *
 * <h2>Rules of Event Handling</h2>
 *
 * <p> These rules formalize the concept of event detection and are used to determine when
 * an event must be reported to the user and the order in which events must occur. These
 * rules assume the event handler and g function conform to the documentation on
 * {@link org.hipparchus.ode.events.ODEEventHandler ODEEventHandler} and
 * {@link org.hipparchus.ode.ODEIntegrator ODEIntegrator}.
 *
 * <ol>
 *     <li> An event must be detected if the g function has changed signs for longer than
 *     maxCheck(t, y(t)). Formally, given times t, t_e1, t_e2 such that t &lt; t_e1 &lt; t_e2,
 *     g(t_e1) = 0, g(t_e2) = 0, and g(t_i) != 0 on the intervals {[t, t_e1&gt;, &lt;t_e1, t_e2&gt;}
 *     (i.e. t_e1, t_e2 are the next two event times after t) then t_e1 will be detected if
 *     maxCheck(t, y(t)) &lt; t_e2. (I.e. the max check interval must be less than the time until
 *     the second event.)</li>
 *     <li> MaxCheck(t, y(t)) is evaluated at step start and called again at the end of the
 *     previous maxCheck; this implies that if maxCheck depends a lot on state, care must be
 *     taken to return conservative values, i.e. it is better to return small maxCheck and hence
 *     perform a lot of checks than missing an event because maxCheck was too large.</li>
 *     <li> MaxCheck should always return positive values, even if integration is backward</li>
 *     <li> For a given tolerance, h, and root, r, the event may occur at any point on the
 *     interval [r-h, r+h]. The tolerance is the larger of the {@code convergence}
 *     parameter and the convergence settings of the root finder specified when
 *     {@link org.hipparchus.ode.ODEIntegrator#addEventDetector(ODEEventDetector) adding}
 *     the event detector. </li>
 *     <li> At most one event is triggered per root. </li>
 *     <li> Events from the same event detector must alternate between increasing and
 *     decreasing events. That is, for every pair of increasing events there must exist an
 *     intervening decreasing event and vice-versa. </li>
 *     <li> An event starts occurring when the
 *     {@link org.hipparchus.ode.events.ODEEventHandler#eventOccurred(org.hipparchus.ode.ODEStateAndDerivative,
 *     ODEEventDetector, boolean) eventOccured()}
 *     method is called. An event stops occurring when eventOccurred() returns or when the
 *     handler's
 *     {@link org.hipparchus.ode.events.ODEEventHandler#resetState(ODEEventDetector,
 *     org.hipparchus.ode.ODEStateAndDerivative) resetState()} method returns if eventOccured() returned
 *     {@link org.hipparchus.ode.events.Action#RESET_STATE RESET_STATE}. </li>
 *     <li> If event A happens before event B then the effects of A occurring are visible
 *     to B. (Including resetting the state or derivatives, or stopping) </li>
 *     <li> Events occur in chronological order. If integration is forward and event A
 *     happens before event B then the time of event B is greater than or equal to the
 *     time of event A. </li>
 *     <li> There is a total order on events. That is for two events A and B either A
 *     happens before B or B happens before A. </li>
 * </ol>
 *
 */
package org.hipparchus.ode.events;
