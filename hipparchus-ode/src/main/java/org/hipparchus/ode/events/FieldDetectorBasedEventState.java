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

package org.hipparchus.ode.events;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver;
import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver.Interval;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.LocalizedODEFormats;
import org.hipparchus.ode.sampling.FieldODEStateInterpolator;
import org.hipparchus.util.FastMath;

/** This class handles the state for one {@link FieldODEEventHandler
 * event handler} during integration steps.
 *
 * <p>Each time the integrator proposes a step, the event handler
 * switching function should be checked. This class handles the state
 * of one handler during one integration step, with references to the
 * state at the end of the preceding step. This information is used to
 * decide if the handler should trigger an event or not during the
 * proposed step.</p>
 *
 * @param <T> the type of the field elements
 */
public class FieldDetectorBasedEventState<T extends CalculusFieldElement<T>> implements FieldEventState<T> {

    /** Event detector.
     * @since 3.0
     */
    private final FieldODEEventDetector<T> detector;

    /** Event solver.
     * @since 3.0
     */
    private final BracketedRealFieldUnivariateSolver<T> solver;

    /** Event handler. */
    private final FieldODEEventHandler<T> handler;

    /** Time of the previous call to g. */
    private T lastT;

    /** Value from the previous call to g. */
    private T lastG;

    /** Time at the beginning of the step. */
    private T t0;

    /** Value of the events handler at the beginning of the step. */
    private T g0;

    /** Sign of g0. */
    private boolean g0Positive;

    /** Indicator of event expected during the step. */
    private boolean pendingEvent;

    /** Occurrence time of the pending event. */
    private T pendingEventTime;

    /**
     * Time to stop propagation if the event is a stop event. Used to enable stopping at
     * an event and then restarting after that event.
     */
    private T stopTime;

    /** Time after the current event. */
    private T afterEvent;

    /** Value of the g function after the current event. */
    private T afterG;

    /** The earliest time considered for events. */
    private T earliestTimeConsidered;

    /** Integration direction. */
    private boolean forward;

    /** Variation direction around pending event.
     *  (this is considered with respect to the integration direction)
     */
    private boolean increasing;

    /** Simple constructor.
     * @param detector event detector
     * @since 3.0
     */
    public FieldDetectorBasedEventState(final FieldODEEventDetector<T> detector) {

        this.detector     = detector;
        this.solver       = detector.getSolver();
        this.handler      = detector.getHandler();

        // some dummy values ...
        t0                = null;
        g0                = null;
        g0Positive        = true;
        pendingEvent      = false;
        pendingEventTime  = null;
        increasing        = true;
        earliestTimeConsidered = null;
        afterEvent = null;
        afterG = null;

    }

    /** Get the underlying event detector.
     * @return underlying event detector
     * @since 3.0
     */
    public FieldODEEventDetector<T> getEventDetector() {
        return detector;
    }

    /** Initialize event handler at the start of an integration.
     * <p>
     * This method is called once at the start of the integration. It
     * may be used by the event handler to initialize some internal data
     * if needed.
     * </p>
     * @param s0 initial state
     * @param t target time for the integration
     *
     */
    @Override
    public void init(final FieldODEStateAndDerivative<T> s0, final T t) {
        detector.init(s0, t);
        lastT = t.getField().getZero().newInstance(Double.NEGATIVE_INFINITY);
        lastG = null;
    }

    /** Compute the value of the switching function.
     * This function must be continuous (at least in its roots neighborhood),
     * as the integrator will need to find its roots to locate the events.
     * @param s the current state information: date, kinematics, attitude
     * @return value of the switching function
     */
    private T g(final FieldODEStateAndDerivative<T> s) {
        if (!s.getTime().subtract(lastT).isZero()) {
            lastG = detector.g(s);
            lastT = s.getTime();
        }
        return lastG;
    }

    /** Reinitialize the beginning of the step.
     * @param interpolator valid for the current step
     * @exception MathIllegalStateException if the interpolator throws one because
     * the number of functions evaluations is exceeded
     */
    public void reinitializeBegin(final FieldODEStateInterpolator<T> interpolator)
        throws MathIllegalStateException {

        forward = interpolator.isForward();
        final FieldODEStateAndDerivative<T> s0 = interpolator.getPreviousState();
        t0 = s0.getTime();
        g0 = g(s0);
        while (g0.isZero()) {
            // excerpt from MATH-421 issue:
            // If an ODE solver is setup with a FieldODEEventHandler that return STOP
            // when the even is triggered, the integrator stops (which is exactly
            // the expected behavior). If however the user wants to restart the
            // solver from the final state reached at the event with the same
            // configuration (expecting the event to be triggered again at a
            // later time), then the integrator may fail to start. It can get stuck
            // at the previous event. The use case for the bug MATH-421 is fairly
            // general, so events occurring exactly at start in the first step should
            // be ignored.

            // extremely rare case: there is a zero EXACTLY at interval start
            // we will use the sign slightly after step beginning to force ignoring this zero
            T tStart = t0.add(solver.getAbsoluteAccuracy().multiply(forward ? 0.5 : -0.5));
            if (tStart.equals(t0)) {
                tStart = nextAfter(t0);
            }
            t0 = tStart;
            g0 = g(interpolator.getInterpolatedState(tStart));
        }
        g0Positive = g0.getReal() > 0;
        // "last" event was increasing
        increasing = g0Positive;

    }

    /** Evaluate the impact of the proposed step on the event handler.
     * @param interpolator step interpolator for the proposed step
     * @return true if the event handler triggers an event before
     * the end of the proposed step
     * @exception MathIllegalStateException if the interpolator throws one because
     * the number of functions evaluations is exceeded
     * @exception MathIllegalArgumentException if the event cannot be bracketed
     */
    @Override
    public boolean evaluateStep(final FieldODEStateInterpolator<T> interpolator)
        throws MathIllegalArgumentException, MathIllegalStateException {

        forward = interpolator.isForward();
        final FieldODEStateAndDerivative<T> s0 = interpolator.getPreviousState();
        final FieldODEStateAndDerivative<T> s1 = interpolator.getCurrentState();
        final T t1 = s1.getTime();
        final T dt = t1.subtract(t0);
        if (dt.abs().subtract(solver.getAbsoluteAccuracy()).getReal() < 0) {
            // we cannot do anything on such a small step, don't trigger any events
            pendingEvent     = false;
            pendingEventTime = null;
            return false;
        }

        T ta = t0;
        T ga = g0;
        for (FieldODEStateAndDerivative<T>sb = nextCheck(s0, s1, interpolator);
             sb != null;
             sb = nextCheck(sb, s1, interpolator)) {

            // evaluate handler value at the end of the substep
            final T tb = sb.getTime();
            final T gb = g(sb);

            // check events occurrence
            if (gb.getReal() == 0.0 || (g0Positive ^ (gb.getReal() > 0))) {
                // there is a sign change: an event is expected during this step
                if (findRoot(interpolator, ta, ga, tb, gb)) {
                    return true;
                }
            } else {
                // no sign change: there is no event for now
                ta = tb;
                ga = gb;
            }

        }

        // no event during the whole step
        pendingEvent     = false;
        pendingEventTime = null;
        return false;

    }

    /** Estimate next state to check.
     * @param done state already checked
     * @param target target state towards which we are checking
     * @param interpolator step interpolator for the proposed step
     * @return intermediate state to check, or exactly {@code null}
     * if we already have {@code done == target}
     * @since 3.0
     */
    private FieldODEStateAndDerivative<T> nextCheck(final FieldODEStateAndDerivative<T> done, final FieldODEStateAndDerivative<T> target,
                                                    final FieldODEStateInterpolator<T> interpolator) {
        if (done == target) {
            // we have already reached target
            return null;
        } else {
            // we have to select some intermediate state
            // attempting to split the remaining time in an integer number of checks
            final T dt       = target.getTime().subtract(done.getTime());
            final double maxCheck = detector.getMaxCheckInterval().currentInterval(done);
            final int    n        = FastMath.max(1, (int) FastMath.ceil(dt.abs().divide(maxCheck).getReal()));
            return n == 1 ? target : interpolator.getInterpolatedState(done.getTime().add(dt.divide(n)));
        }
    }

    /**
     * Find a root in a bracketing interval.
     *
     * <p> When calling this method one of the following must be true. Either ga == 0, gb
     * == 0, (ga < 0  and gb > 0), or (ga > 0 and gb < 0).
     *
     * @param interpolator that covers the interval.
     * @param ta           earliest possible time for root.
     * @param ga           g(ta).
     * @param tb           latest possible time for root.
     * @param gb           g(tb).
     * @return if a zero crossing was found.
     */
    private boolean findRoot(final FieldODEStateInterpolator<T> interpolator,
                             final T ta,
                             final T ga,
                             final T tb,
                             final T gb) {
        // check there appears to be a root in [ta, tb]
        check(ga.getReal() == 0.0 || gb.getReal() == 0.0 ||
                (ga.getReal() > 0.0 && gb.getReal() < 0.0) ||
                (ga.getReal() < 0.0 && gb.getReal() > 0.0));

        final int maxIterationCount = detector.getMaxIterationCount();
        final CalculusFieldUnivariateFunction<T> f = t -> g(interpolator.getInterpolatedState(t));

        // prepare loop below
        T loopT = ta;
        T loopG = ga;

        // event time, just at or before the actual root.
        T beforeRootT = null;
        T beforeRootG = null;
        // time on the other side of the root.
        // Initialized the the loop below executes once.
        T afterRootT = ta;
        T afterRootG = ta.getField().getZero();

        // check for some conditions that the root finders don't like
        // these conditions cannot not happen in the loop below
        // the ga == 0.0 case is handled by the loop below
        if (ta.getReal() == tb.getReal()) {
            // both non-zero but times are the same. Probably due to reset state
            beforeRootT = ta;
            beforeRootG = ga;
            afterRootT = shiftedBy(beforeRootT, solver.getAbsoluteAccuracy());
            afterRootG = f.value(afterRootT);
        } else if (!ga.isZero() && gb.isZero()) {
            // hard: ga != 0.0 and gb == 0.0
            // look past gb by up to convergence to find next sign
            // throw an exception if g(t) = 0.0 in [tb, tb + convergence]
            beforeRootT = tb;
            beforeRootG = gb;
            afterRootT = shiftedBy(beforeRootT, solver.getAbsoluteAccuracy());
            afterRootG = f.value(afterRootT);
        } else if (!ga.isZero()) {
            final T newGa = f.value(ta);
            if (ga.getReal() > 0 != newGa.getReal() > 0) {
                // both non-zero, step sign change at ta, possibly due to reset state
                final T nextT = minTime(shiftedBy(ta, solver.getAbsoluteAccuracy()), tb);
                final T nextG = f.value(nextT);
                if (nextG.getReal() > 0.0 == g0Positive) {
                    // the sign change between ga and new Ga just moved the root less than one convergence
                    // threshold later, we are still in a regular search for another root before tb,
                    // we just need to fix the bracketing interval
                    // (see issue https://github.com/Hipparchus-Math/hipparchus/issues/184)
                    loopT = nextT;
                    loopG = nextG;
                } else {
                    beforeRootT = ta;
                    beforeRootG = newGa;
                    afterRootT  = nextT;
                    afterRootG  = nextG;
                }
            }
        }

        // loop to skip through "fake" roots, i.e. where g(t) = g'(t) = 0.0
        // executed once if we didn't hit a special case above
        while ((afterRootG.isZero() || afterRootG.getReal() > 0.0 == g0Positive) &&
               strictlyAfter(afterRootT, tb)) {
            if (loopG.isZero()) {
                // ga == 0.0 and gb may or may not be 0.0
                // handle the root at ta first
                beforeRootT = loopT;
                beforeRootG = loopG;
                afterRootT = minTime(shiftedBy(beforeRootT, solver.getAbsoluteAccuracy()), tb);
                afterRootG = f.value(afterRootT);
            } else {
                // both non-zero, the usual case, use a root finder.
                if (forward) {
                    try {
                        final Interval<T> interval =
                                        solver.solveInterval(maxIterationCount, f, loopT, tb);
                        beforeRootT = interval.getLeftAbscissa();
                        beforeRootG = interval.getLeftValue();
                        afterRootT = interval.getRightAbscissa();
                        afterRootG = interval.getRightValue();
                        // CHECKSTYLE: stop IllegalCatch check
                    } catch (RuntimeException e) { // NOPMD
                        // CHECKSTYLE: resume IllegalCatch check
                        throw new MathIllegalStateException(e, LocalizedODEFormats.FIND_ROOT,
                                                            detector, loopT.getReal(), loopG.getReal(),
                                                            tb.getReal(), gb.getReal(),
                                                            lastT.getReal(), lastG.getReal());
                    }
                } else {
                    try {
                        final Interval<T> interval =
                                        solver.solveInterval(maxIterationCount, f, tb, loopT);
                        beforeRootT = interval.getRightAbscissa();
                        beforeRootG = interval.getRightValue();
                        afterRootT = interval.getLeftAbscissa();
                        afterRootG = interval.getLeftValue();
                        // CHECKSTYLE: stop IllegalCatch check
                    } catch (RuntimeException e) { // NOPMD
                        // CHECKSTYLE: resume IllegalCatch check
                        throw new MathIllegalStateException(e, LocalizedODEFormats.FIND_ROOT,
                                                            detector, tb.getReal(), gb.getReal(),
                                                            loopT.getReal(), loopG.getReal(),
                                                            lastT.getReal(), lastG.getReal());
                    }
                }
            }
            // tolerance is set to less than 1 ulp
            // assume tolerance is 1 ulp
            if (beforeRootT == afterRootT) {
                afterRootT = nextAfter(afterRootT);
                afterRootG = f.value(afterRootT);
            }
            // check loop is making some progress
            check((forward && afterRootT.getReal() > beforeRootT.getReal()) ||
                  (!forward && afterRootT.getReal() < beforeRootT.getReal()));
            // setup next iteration
            loopT = afterRootT;
            loopG = afterRootG;
        }

        // figure out the result of root finding, and return accordingly
        if (afterRootG.isZero() || afterRootG.getReal() > 0.0 == g0Positive) {
            // loop gave up and didn't find any crossing within this step
            return false;
        } else {
            // real crossing
            check(beforeRootT != null && beforeRootG != null);
            // variation direction, with respect to the integration direction
            increasing = !g0Positive;
            pendingEventTime = beforeRootT;
            stopTime = beforeRootG.isZero() ? beforeRootT : afterRootT;
            pendingEvent = true;
            afterEvent = afterRootT;
            afterG = afterRootG;

            // check increasing set correctly
            check(afterG.getReal() > 0 == increasing);
            check(increasing == gb.getReal() >= ga.getReal());

            return true;
        }
    }

    /**
     * Try to accept the current history up to the given time.
     *
     * <p> It is not necessary to call this method before calling {@link
     * #doEvent(FieldODEStateAndDerivative)} with the same state. It is necessary to call this
     * method before you call {@link #doEvent(FieldODEStateAndDerivative)} on some other event
     * detector.
     *
     * @param state        to try to accept.
     * @param interpolator to use to find the new root, if any.
     * @return if the event detector has an event it has not detected before that is on or
     * before the same time as {@code state}. In other words {@code false} means continue
     * on while {@code true} means stop and handle my event first.
     */
    public boolean tryAdvance(final FieldODEStateAndDerivative<T> state,
                              final FieldODEStateInterpolator<T> interpolator) {
        final T t = state.getTime();
        // check this is only called before a pending event.
        check(!pendingEvent || !strictlyAfter(pendingEventTime, t));

        final boolean meFirst;

        // just found an event and we know the next time we want to search again
        if (earliestTimeConsidered != null && strictlyAfter(t, earliestTimeConsidered)) {
            meFirst = false;
        } else {
            // check g function to see if there is a new event
            final T g = g(state);
            final boolean positive = g.getReal() > 0;

            if (positive == g0Positive) {
                // g function has expected sign
                g0 = g; // g0Positive is the same
                meFirst = false;
            } else {
                // found a root we didn't expect -> find precise location
                final T oldPendingEventTime = pendingEventTime;
                final boolean foundRoot = findRoot(interpolator, t0, g0, t, g);
                // make sure the new root is not the same as the old root, if one exists
                meFirst = foundRoot && !pendingEventTime.equals(oldPendingEventTime);
            }
        }

        if (!meFirst) {
            // advance t0 to the current time so we can't find events that occur before t
            t0 = t;
        }

        return meFirst;
    }

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
    @Override
    public FieldEventOccurrence<T> doEvent(final FieldODEStateAndDerivative<T> state) {
        // check event is pending and is at the same time
        check(pendingEvent);
        check(state.getTime() == this.pendingEventTime);

        final Action action = handler.eventOccurred(state, detector, increasing == forward);
        final FieldODEState<T> newState;
        if (action == Action.RESET_STATE) {
            newState = handler.resetState(detector, state);
        } else {
            newState = state;
        }
        // clear pending event
        pendingEvent = false;
        pendingEventTime = null;
        // setup for next search
        earliestTimeConsidered = afterEvent;
        t0 = afterEvent;
        g0 = afterG;
        g0Positive = increasing;
        // check g0Positive set correctly
        check(g0.isZero() || g0Positive == (g0.getReal() > 0));
        return new FieldEventOccurrence<>(action, newState, stopTime);
    }

    /**
     * Check the ordering of two times.
     *
     * @param t1 the first time.
     * @param t2 the second time.
     * @return true if {@code t2} is strictly after {@code t1} in the propagation
     * direction.
     */
    private boolean strictlyAfter(final T t1, final T t2) {
        return forward ? t1.getReal() < t2.getReal() : t2.getReal() < t1.getReal();
    }

    /**
     * Get the next number after the given number in the current propagation direction.
     *
     * <p> Assumes T has the same precision as a double.
     *
     * @param t input time
     * @return t +/- 1 ulp depending on the direction.
     */
    private T nextAfter(final T t) {
        // direction
        final int sign = forward ? 1 : -1;
        final double ulp = FastMath.ulp(t.getReal());
        return t.add(sign * ulp);
    }

    /**
     * Same as keyword assert, but throw a {@link MathRuntimeException}.
     *
     * @param condition to check
     * @throws MathRuntimeException if {@code condition} is false.
     */
    private void check(final boolean condition) throws MathRuntimeException {
        if (!condition) {
            throw MathRuntimeException.createInternalError();
        }
    }

    /**
     * Get the time that happens first along the current propagation direction: {@link
     * #forward}.
     *
     * @param a first time
     * @param b second time
     * @return min(a, b) if forward, else max (a, b)
     */
    private T minTime(final T a, final T b) {
        return forward ? FastMath.min(a, b) : FastMath.max(a, b);
    }

    /**
     * Shift a time value along the current integration direction: {@link #forward}.
     *
     * @param t     the time to shift.
     * @param delta the amount to shift.
     * @return t + delta if forward, else t - delta. If the result has to be rounded it
     * will be rounded to be before the true value of t + delta.
     */
    private T shiftedBy(final T t, final T delta) {
        if (forward) {
            final T ret = t.add(delta);
            if (ret.subtract(t).getReal() > delta.getReal()) {
                // nextDown(ret)
                return ret.subtract(FastMath.ulp(ret.getReal()));
            } else {
                return ret;
            }
        } else {
            final T ret = t.subtract(delta);
            if (t.subtract(ret).getReal() > delta.getReal()) {
                // nextUp(ret)
                return ret.add(FastMath.ulp(ret.getReal()));
            } else {
                return ret;
            }
        }
    }

    /** Get the occurrence time of the event triggered in the current step.
     * @return occurrence time of the event triggered in the current
     * step or infinity if no events are triggered
     */
    @Override
    public T getEventTime() {
        return pendingEvent ?
               pendingEventTime :
               t0.getField().getZero().add(forward ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY);
    }

}
