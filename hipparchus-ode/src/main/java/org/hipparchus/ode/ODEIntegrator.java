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

import java.util.List;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.events.ODEEventDetector;
import org.hipparchus.ode.events.ODEStepEndHandler;
import org.hipparchus.ode.sampling.ODEStepHandler;

/** This interface represents a first order integrator for
 * differential equations.

 * <p>The classes which are devoted to solve first order differential
 * equations should implement this interface. The problems which can
 * be handled should implement the {@link
 * OrdinaryDifferentialEquation} interface.</p>
 *
 * @see OrdinaryDifferentialEquation
 * @see org.hipparchus.ode.sampling.ODEStepHandler
 * @see org.hipparchus.ode.events.ODEEventHandler
 */
public interface ODEIntegrator  {

    /** Get the name of the method.
     * @return name of the method
     */
    String getName();

    /** Add a step handler to this integrator.
     * <p>The handler will be called by the integrator for each accepted
     * step.</p>
     * @param handler handler for the accepted steps
     * @see #getStepHandlers()
     * @see #clearStepHandlers()
     */
    void addStepHandler(ODEStepHandler handler);

    /** Get all the step handlers that have been added to the integrator.
     * @return an unmodifiable collection of the added events handlers
     * @see #addStepHandler(ODEStepHandler)
     * @see #clearStepHandlers()
     */
    List<ODEStepHandler> getStepHandlers();

    /** Remove all the step handlers that have been added to the integrator.
     * @see #addStepHandler(ODEStepHandler)
     * @see #getStepHandlers()
     */
    void clearStepHandlers();

    /** Add an event detector to the integrator.
     * @param detector event detector
     * @see #getEventDetectors()
     * @see #clearEventDetectors()
     * @since 3.0
     */
    void addEventDetector(ODEEventDetector detector);

    /** Get all the event detectors that have been added to the integrator.
     * @return an unmodifiable list of the added events detectors
     * @see #addEventDetector(ODEEventDetector)
     * @see #clearEventDetectors()
     * @since 3.0
     */
    List<ODEEventDetector> getEventDetectors();

    /** Remove all the event handlers that have been added to the integrator.
     * @see #addEventDetector(ODEEventDetector)
     * @see #getEventDetectors()
     * @since 3.0
     */
    void clearEventDetectors();

    /** Add a handler for step ends to the integrator.
     * <p>
     * The {@link ODEStepEndHandler#stepEndOccurred(ODEStateAndDerivative, boolean)
     * stepEndOccurred(state, forward)} method of the {@code handler} will be called
     * at each step end.
     * </p>
     * @param handler handler for step ends
     * @see #getStepEndHandlers()
     * @see #clearStepEndHandlers()
     * @since 3.0
     */
    void addStepEndHandler(ODEStepEndHandler handler);

    /** Get all the handlers for step ends that have been added to the integrator.
     * @return an unmodifiable list of the added step end handlers
     * @see #addStepEndHandler(ODEStepEndHandler)
     * @see #clearStepEndHandlers()
     * @since 3.0
     */
    List<ODEStepEndHandler> getStepEndHandlers();

    /** Remove all the handlers for step ends that have been added to the integrator.
     * @see #addStepEndHandler(ODEStepEndHandler)
     * @see #getStepEndHandlers()
     * @since 3.0
     */
    void clearStepEndHandlers();

    /** Get the state at step start time t<sub>i</sub>.
     * <p>This method can be called during integration (typically by
     * the object implementing the {@link OrdinaryDifferentialEquation
     * differential equations} problem) if the value of the current step that
     * is attempted is needed.</p>
     * <p>The result is undefined if the method is called outside of
     * calls to <code>integrate</code>.</p>
     * @return state at step start time t<sub>i</sub>
     */
    ODEStateAndDerivative getStepStart();

    /** Get the current signed value of the integration stepsize.
     * <p>This method can be called during integration (typically by
     * the object implementing the {@link OrdinaryDifferentialEquation
     * differential equations} problem) if the signed value of the current stepsize
     * that is tried is needed.</p>
     * <p>The result is undefined if the method is called outside of
     * calls to <code>integrate</code>.</p>
     * @return current signed value of the stepsize
     */
    double getCurrentSignedStepsize();

    /** Set the maximal number of differential equations function evaluations.
     * <p>The purpose of this method is to avoid infinite loops which can occur
     * for example when stringent error constraints are set or when lots of
     * discrete events are triggered, thus leading to many rejected steps.</p>
     * @param maxEvaluations maximal number of function evaluations (negative
     * values are silently converted to maximal integer value, thus representing
     * almost unlimited evaluations)
     */
    void setMaxEvaluations(int maxEvaluations);

    /** Get the maximal number of functions evaluations.
     * @return maximal number of functions evaluations
     */
    int getMaxEvaluations();

    /** Get the number of evaluations of the differential equations function.
     * <p>
     * The number of evaluations corresponds to the last call to the
     * <code>integrate</code> method. It is 0 if the method has not been called yet.
     * </p>
     * @return number of evaluations of the differential equations function
     */
    int getEvaluations();

    /** Integrate the differential equations up to the given time.
     * <p>This method solves an Initial Value Problem (IVP).</p>
     * <p>Since this method stores some internal state variables made
     * available in its public interface during integration ({@link
     * #getCurrentSignedStepsize()}), it is <em>not</em> thread-safe.</p>
     * @param equations differential equations to integrate
     * @param initialState initial state (time, primary and secondary state vectors)
     * @param finalTime target time for the integration
     * (can be set to a value smaller than {@code t0} for backward integration)
     * @return final state, its time will be the same as {@code finalTime} if
     * integration reached its target, but may be different if some {@link
     * org.hipparchus.ode.events.ODEEventHandler} stops it at some point.
     * @exception MathIllegalArgumentException if integration step is too small
     * @exception MathIllegalStateException if the number of functions evaluations is exceeded
     * @exception MathIllegalArgumentException if the location of an event cannot be bracketed
     */
    ODEStateAndDerivative integrate(ExpandableODE equations,
                                    ODEState initialState, double finalTime)
        throws MathIllegalArgumentException, MathIllegalStateException;

    /** Integrate the differential equations up to the given time.
     * <p>This method solves an Initial Value Problem (IVP).</p>
     * <p>Since this method stores some internal state variables made
     * available in its public interface during integration ({@link
     * #getCurrentSignedStepsize()}), it is <em>not</em> thread-safe.</p>
     * @param equations differential equations to integrate
     * @param initialState initial state (time, primary and secondary state vectors)
     * @param finalTime target time for the integration
     * (can be set to a value smaller than {@code t0} for backward integration)
     * @return final state, its time will be the same as {@code finalTime} if
     * integration reached its target, but may be different if some {@link
     * org.hipparchus.ode.events.ODEEventHandler} stops it at some point.
     * @exception MathIllegalArgumentException if integration step is too small
     * @exception MathIllegalStateException if the number of functions evaluations is exceeded
     * @exception MathIllegalArgumentException if the location of an event cannot be bracketed
     */
    default ODEStateAndDerivative integrate(OrdinaryDifferentialEquation equations,
                                            ODEState initialState, double finalTime)
        throws MathIllegalArgumentException, MathIllegalStateException {
        return integrate(new ExpandableODE(equations), initialState, finalTime);
    }

}
