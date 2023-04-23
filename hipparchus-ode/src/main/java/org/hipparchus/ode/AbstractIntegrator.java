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

package org.hipparchus.ode;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.events.Action;
import org.hipparchus.ode.events.DetectorBasedEventState;
import org.hipparchus.ode.events.EventOccurrence;
import org.hipparchus.ode.events.EventState;
import org.hipparchus.ode.events.ODEEventDetector;
import org.hipparchus.ode.events.ODEStepEndHandler;
import org.hipparchus.ode.events.StepEndEventState;
import org.hipparchus.ode.sampling.AbstractODEStateInterpolator;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Incrementor;

/**
 * Base class managing common boilerplate for all integrators.
 */
public abstract class AbstractIntegrator implements ODEIntegrator {

    /** Step handler. */
    private List<ODEStepHandler> stepHandlers;

    /** Current step start time. */
    private ODEStateAndDerivative stepStart;

    /** Current stepsize. */
    private double stepSize;

    /** Indicator for last step. */
    private boolean isLastStep;

    /** Indicator that a state or derivative reset was triggered by some event. */
    private boolean resetOccurred;

    /** Events states related to event detectors. */
    private List<DetectorBasedEventState> detectorBasedEventsStates;

    /** Events states related to step end. */
    private List<StepEndEventState> stepEndEventsStates;

    /** Initialization indicator of events states. */
    private boolean statesInitialized;

    /** Name of the method. */
    private final String name;

    /** Counter for number of evaluations. */
    private Incrementor evaluations;

    /** Differential equations to integrate. */
    private transient ExpandableODE equations;

    /** Build an instance.
     * @param name name of the method
     */
    protected AbstractIntegrator(final String name) {
        this.name                 = name;
        stepHandlers              = new ArrayList<>();
        stepStart                 = null;
        stepSize                  = Double.NaN;
        detectorBasedEventsStates = new ArrayList<>();
        stepEndEventsStates       = new ArrayList<>();
        statesInitialized         = false;
        evaluations               = new Incrementor();
    }

    /** {@inheritDoc} */
    @Override
    public String getName() {
        return name;
    }

    /** {@inheritDoc} */
    @Override
    public void addStepHandler(final ODEStepHandler handler) {
        stepHandlers.add(handler);
    }

    /** {@inheritDoc} */
    @Override
    public List<ODEStepHandler> getStepHandlers() {
        return Collections.unmodifiableList(stepHandlers);
    }

    /** {@inheritDoc} */
    @Override
    public void clearStepHandlers() {
        stepHandlers.clear();
    }

    /** {@inheritDoc} */
    @Override
    public void addEventDetector(final ODEEventDetector detector) {
        detectorBasedEventsStates.add(new DetectorBasedEventState(detector));
    }

    /** {@inheritDoc} */
    @Override
    public List<ODEEventDetector> getEventDetectors() {
        return detectorBasedEventsStates.stream().map(es -> es.getEventDetector()).collect(Collectors.toList());
    }

    /** {@inheritDoc} */
    @Override
    public void clearEventDetectors() {
        detectorBasedEventsStates.clear();
    }

    /** {@inheritDoc} */
    @Override
    public void addStepEndHandler(ODEStepEndHandler handler) {
        stepEndEventsStates.add(new StepEndEventState(handler));
    }

    /** {@inheritDoc} */
    @Override
    public List<ODEStepEndHandler> getStepEndHandlers() {
        return stepEndEventsStates.stream().map(es -> es.getHandler()).collect(Collectors.toList());
    }

    /** {@inheritDoc} */
    @Override
    public void clearStepEndHandlers() {
        stepEndEventsStates.clear();
    }

    /** {@inheritDoc} */
    @Override
    public double getCurrentSignedStepsize() {
        return stepSize;
    }

    /** {@inheritDoc} */
    @Override
    public void setMaxEvaluations(int maxEvaluations) {
        evaluations = evaluations.withMaximalCount((maxEvaluations < 0) ? Integer.MAX_VALUE : maxEvaluations);
    }

    /** {@inheritDoc} */
    @Override
    public int getMaxEvaluations() {
        return evaluations.getMaximalCount();
    }

    /** {@inheritDoc} */
    @Override
    public int getEvaluations() {
        return evaluations.getCount();
    }

    /**
     * Prepare the start of an integration.
     *
     * @param eqn equations to integrate
     * @param s0  initial state vector
     * @param t   target time for the integration
     * @return Initial state with computed derivatives.
     */
    protected ODEStateAndDerivative initIntegration(final ExpandableODE eqn,
                                                    final ODEState s0, final double t) {

        this.equations = eqn;
        evaluations    = evaluations.withCount(0);

        // initialize ODE
        eqn.init(s0, t);

        // set up derivatives of initial state (including primary and secondary components)
        final double   t0    = s0.getTime();
        final double[] y0    = s0.getCompleteState();
        final double[] y0Dot = computeDerivatives(t0, y0);

        // built the state
        final ODEStateAndDerivative s0WithDerivatives =
                        eqn.getMapper().mapStateAndDerivative(t0, y0, y0Dot);

        // initialize detector based event states (both  and step end based)
        detectorBasedEventsStates.stream().forEach(s -> {
            s.init(s0WithDerivatives, t);
            s.getEventDetector().getHandler().init(s0WithDerivatives, t, s.getEventDetector());
        });

        // initialize step end based event states
        stepEndEventsStates.stream().forEach(s -> {
            s.init(s0WithDerivatives, t);
            s.getHandler().init(s0WithDerivatives, t);
        });

        // initialize step handlers
        for (ODEStepHandler handler : stepHandlers) {
            handler.init(s0WithDerivatives, t);
        }

        setStateInitialized(false);

        return s0WithDerivatives;

    }

    /** Get the differential equations to integrate.
     * @return differential equations to integrate
     */
    protected ExpandableODE getEquations() {
        return equations;
    }

    /** Get the evaluations counter.
     * @return evaluations counter
     */
    protected Incrementor getEvaluationsCounter() {
        return evaluations;
    }

    /** Compute the derivatives and check the number of evaluations.
     * @param t current value of the independent <I>time</I> variable
     * @param y array containing the current value of the state vector
     * @return state completed with derivatives
     * @exception MathIllegalArgumentException if arrays dimensions do not match equations settings
     * @exception MathIllegalStateException if the number of functions evaluations is exceeded
     * @exception NullPointerException if the ODE equations have not been set (i.e. if this method
     * is called outside of a call to {@link #integrate(ExpandableODE, ODEState, double) integrate}
     */
    public double[] computeDerivatives(final double t, final double[] y)
        throws MathIllegalArgumentException, MathIllegalStateException, NullPointerException {
        evaluations.increment();
        return equations.computeDerivatives(t, y);
    }

    /** Set the stateInitialized flag.
     * <p>This method must be called by integrators with the value
     * {@code false} before they start integration, so a proper lazy
     * initialization is done automatically on the first step.</p>
     * @param stateInitialized new value for the flag
     */
    protected void setStateInitialized(final boolean stateInitialized) {
        this.statesInitialized = stateInitialized;
    }

    /** Accept a step, triggering events and step handlers.
     * @param interpolator step interpolator
     * @param tEnd final integration time
     * @return state at end of step
     * @exception MathIllegalStateException if the interpolator throws one because
     * the number of functions evaluations is exceeded
     * @exception MathIllegalArgumentException if the location of an event cannot be bracketed
     * @exception MathIllegalArgumentException if arrays dimensions do not match equations settings
     */
    protected ODEStateAndDerivative acceptStep(final AbstractODEStateInterpolator interpolator,
                                               final double tEnd)
            throws MathIllegalArgumentException, MathIllegalStateException {

        ODEStateAndDerivative previousState = interpolator.getGlobalPreviousState();
        final ODEStateAndDerivative currentState = interpolator.getGlobalCurrentState();
        AbstractODEStateInterpolator restricted = interpolator;


        // initialize the events states if needed
        if (!statesInitialized) {
            // initialize event states
            detectorBasedEventsStates.stream().forEach(s -> s.reinitializeBegin(interpolator));
            statesInitialized = true;
        }

        // set end of step
        stepEndEventsStates.stream().forEach(s -> s.setStepEnd(currentState.getTime()));

        // search for next events that may occur during the step
        final int orderingSign = interpolator.isForward() ? +1 : -1;
        final Queue<EventState> occurringEvents = new PriorityQueue<>(new Comparator<EventState>() {
            /** {@inheritDoc} */
            @Override
            public int compare(final EventState es0, final EventState es1) {
                return orderingSign * Double.compare(es0.getEventTime(), es1.getEventTime());
            }
        });

        resetOccurred = false;
        boolean doneWithStep = false;
        resetEvents:
        do {

            // Evaluate all event detectors and end steps for events
            occurringEvents.clear();
            final AbstractODEStateInterpolator finalRestricted = restricted;
            Stream.concat(detectorBasedEventsStates.stream(), stepEndEventsStates.stream()).
            forEach(s -> { if (s.evaluateStep(finalRestricted)) {
                    // the event occurs during the current step
                    occurringEvents.add(s);
                }
            });

            do {

                eventLoop:
                while (!occurringEvents.isEmpty()) {

                    // handle the chronologically first event
                    final EventState currentEvent = occurringEvents.poll();

                    // get state at event time
                    ODEStateAndDerivative eventState = restricted.getInterpolatedState(currentEvent.getEventTime());

                    // restrict the interpolator to the first part of the step, up to the event
                    restricted = restricted.restrictStep(previousState, eventState);

                    // try to advance all event states related to detectors to current time
                    for (final DetectorBasedEventState state : detectorBasedEventsStates) {
                        if (state != currentEvent && state.tryAdvance(eventState, interpolator)) {
                            // we need to handle another event first
                            // remove event we just updated to prevent heap corruption
                            occurringEvents.remove(state);
                            // add it back to update its position in the heap
                            occurringEvents.add(state);
                            // re-queue the event we were processing
                            occurringEvents.add(currentEvent);
                            continue eventLoop;
                        }
                    }
                    // all event detectors agree we can advance to the current event time

                    // handle the first part of the step, up to the event
                    for (final ODEStepHandler handler : stepHandlers) {
                        handler.handleStep(restricted);
                    }

                    // acknowledge event occurrence
                    final EventOccurrence occurrence = currentEvent.doEvent(eventState);
                    final Action action = occurrence.getAction();
                    isLastStep = action == Action.STOP;

                    if (isLastStep) {

                        // ensure the event is after the root if it is returned STOP
                        // this lets the user integrate to a STOP event and then restart
                        // integration from the same time.
                        final ODEStateAndDerivative savedState = eventState;
                        eventState = interpolator.getInterpolatedState(occurrence.getStopTime());
                        restricted = interpolator.restrictStep(savedState, eventState);

                        // handle the almost zero size last part of the final step, at event time
                        for (final ODEStepHandler handler : stepHandlers) {
                            handler.handleStep(restricted);
                            handler.finish(restricted.getCurrentState());
                        }

                    }

                    if (isLastStep) {
                        // the event asked to stop integration
                        return eventState;
                    }

                    if (action == Action.RESET_DERIVATIVES || action == Action.RESET_STATE) {
                        // some event handler has triggered changes that
                        // invalidate the derivatives, we need to recompute them
                        final ODEState newState = occurrence.getNewState();
                        final double[] y = newState.getCompleteState();
                        final double[] yDot = computeDerivatives(newState.getTime(), y);
                        resetOccurred = true;
                        return equations.getMapper().mapStateAndDerivative(newState.getTime(), y, yDot);
                    }
                    // at this point action == Action.CONTINUE or Action.RESET_EVENTS

                    // prepare handling of the remaining part of the step
                    previousState = eventState;
                    restricted = restricted.restrictStep(eventState, currentState);

                    if (action == Action.RESET_EVENTS) {
                        continue resetEvents;
                    }

                    // at this point action == Action.CONTINUE
                    // check if the same event occurs again in the remaining part of the step
                    if (currentEvent.evaluateStep(restricted)) {
                        // the event occurs during the current step
                        occurringEvents.add(currentEvent);
                    }

                }

                // last part of the step, after the last event. Advance all event
                // detectors to the end of the step. Should never find new events unless
                // a previous event modified the g function of another event detector when
                // it returned Action.CONTINUE. Detecting such events here is unreliable
                // and RESET_EVENTS should be used instead. Other option is to replace
                // tryAdvance(...) with a doAdvance(...) that throws an exception when
                // the g function sign is not as expected.
                for (final DetectorBasedEventState state : detectorBasedEventsStates) {
                    if (state.tryAdvance(currentState, interpolator)) {
                        occurringEvents.add(state);
                    }
                }

            } while (!occurringEvents.isEmpty());

            doneWithStep = true;
        } while (!doneWithStep);

        isLastStep = isLastStep || FastMath.abs(currentState.getTime() - tEnd) < FastMath.ulp(tEnd);

        // handle the remaining part of the step, after all events if any
        for (ODEStepHandler handler : stepHandlers) {
            handler.handleStep(restricted);
            if (isLastStep) {
                handler.finish(restricted.getCurrentState());
            }
        }

        return currentState;

    }

    /** Check the integration span.
     * @param initialState initial state
     * @param t target time for the integration
     * @exception MathIllegalArgumentException if integration span is too small
     * @exception MathIllegalArgumentException if adaptive step size integrators
     * tolerance arrays dimensions are not compatible with equations settings
     */
    protected void sanityChecks(final ODEState initialState, final double t)
        throws MathIllegalArgumentException {

        final double threshold = 1000 * FastMath.ulp(FastMath.max(FastMath.abs(initialState.getTime()),
                                                                  FastMath.abs(t)));
        final double dt = FastMath.abs(initialState.getTime() - t);
        if (dt < threshold) {
            throw new MathIllegalArgumentException(LocalizedODEFormats.TOO_SMALL_INTEGRATION_INTERVAL,
                                                   dt, threshold, false);
        }

    }

    /** Check if a reset occurred while last step was accepted.
     * @return true if a reset occurred while last step was accepted
     */
    protected boolean resetOccurred() {
        return resetOccurred;
    }

    /** Set the current step size.
     * @param stepSize step size to set
     */
    protected void setStepSize(final double stepSize) {
        this.stepSize = stepSize;
    }

    /** Get the current step size.
     * @return current step size
     */
    protected double getStepSize() {
        return stepSize;
    }
    /** Set current step start.
     * @param stepStart step start
     */
    protected void setStepStart(final ODEStateAndDerivative stepStart) {
        this.stepStart = stepStart;
    }

    /**  {@inheritDoc} */
    @Override
    public ODEStateAndDerivative getStepStart() {
        return stepStart;
    }

    /** Set the last state flag.
     * @param isLastStep if true, this step is the last one
     */
    protected void setIsLastStep(final boolean isLastStep) {
        this.isLastStep = isLastStep;
    }

    /** Check if this step is the last one.
     * @return true if this step is the last one
     */
    protected boolean isLastStep() {
        return isLastStep;
    }

}
