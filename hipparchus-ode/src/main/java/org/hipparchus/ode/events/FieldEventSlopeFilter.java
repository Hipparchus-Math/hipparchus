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

import java.util.Arrays;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.util.MathArrays;

/** Wrapper used to detect only increasing or decreasing events.
 *
 * <p>General {@link FieldODEEventDetector events} are defined implicitly
 * by a {@link FieldODEEventDetector#g(FieldODEStateAndDerivative) g function} crossing
 * zero. This function needs to be continuous in the event neighborhood,
 * and its sign must remain consistent between events. This implies that
 * during an ODE integration, events triggered are alternately events
 * for which the function increases from negative to positive values,
 * and events for which the function decreases from positive to
 * negative values.
 * </p>
 *
 * <p>Sometimes, users are only interested in one type of event (say
 * increasing events for example) and not in the other type. In these
 * cases, looking precisely for all events location and triggering
 * events that will later be ignored is a waste of computing time.</p>
 *
 * <p>Users can wrap a regular {@link FieldODEEventDetector event detector} in
 * an instance of this class and provide this wrapping instance to
 * the {@link org.hipparchus.ode.FieldODEIntegrator ODE solver}
 * in order to avoid wasting time looking for uninteresting events.
 * The wrapper will intercept the calls to the {@link
 * FieldODEEventDetector#g(FieldODEStateAndDerivative) g function} and to the {@link
 * FieldODEEventHandler#eventOccurred(FieldODEStateAndDerivative, FieldODEEventDetector, boolean)
 * eventOccurred} method in order to ignore uninteresting events. The
 * wrapped regular {@link FieldODEEventDetector event detector} will the see only
 * the interesting events, i.e. either only {@code increasing} events or
 * {@code decreasing} events. the number of calls to the {@link
 * FieldODEEventDetector#g(FieldODEStateAndDerivative) g function} will also be reduced.</p>
 *
 * @param <T> type of the event detector
 * @param <E> the type of the field elements
 * @since 3.0
 */

public class FieldEventSlopeFilter<T extends FieldODEEventDetector<E>, E extends CalculusFieldElement<E>>
    extends AbstractFieldODEDetector<FieldEventSlopeFilter<T, E>, E> {

    /** Number of past transformers updates stored. */
    private static final int HISTORY_SIZE = 100;

    /** Wrapped event detector.
     * @since 3.0
     */
    private final T rawDetector;

    /** Filter to use. */
    private final FilterType filter;

    /** Transformers of the g function. */
    private final Transformer[] transformers;

    /** Update time of the transformers. */
    private final E[] updates;

    /** Indicator for forward integration. */
    private boolean forward;

    /** Extreme time encountered so far. */
    private E extremeT;

    /** Wrap a {@link FieldODEEventDetector eve,t detector}.
     * @param field field to which array elements belong
     * @param rawDetector event detector to wrap
     * @param filter filter to use
     */
    public FieldEventSlopeFilter(final Field<E> field, final T rawDetector, final FilterType filter) {
        this(field,
             rawDetector.getMaxCheckInterval(), rawDetector.getMaxIterationCount(),
             rawDetector.getSolver(), new LocalHandler<>(rawDetector.getHandler()),
             rawDetector, filter);
    }

    /** Private constructor with full parameters.
     * <p>
     * This constructor is private as users are expected to use the builder
     * API with the various {@code withXxx()} methods to set up the instance
     * in a readable manner without using a huge amount of parameters.
     * </p>
     * @param field field to which array elements belong
     * @param maxCheck maximum checking interval (s)
     * @param maxIter maximum number of iterations in the event time search
     * @param solver solver to user for locating event
     * @param handler event handler to call at event occurrences
     * @param rawDetector event detector to wrap
     * @param filter filter to use
     */
    private FieldEventSlopeFilter(final Field<E> field,
                                  final FieldAdaptableInterval<E> maxCheck, final int maxIter,
                                  final BracketedRealFieldUnivariateSolver<E> solver,
                                  final FieldODEEventHandler<E> handler,
                                  final T rawDetector, final FilterType filter) {
        super(maxCheck, maxIter, solver, handler);
        this.rawDetector  = rawDetector;
        this.filter       = filter;
        this.transformers = new Transformer[HISTORY_SIZE];
        this.updates      = MathArrays.buildArray(field, HISTORY_SIZE);
    }

    /** {@inheritDoc} */
    @Override
    protected FieldEventSlopeFilter<T, E> create(final FieldAdaptableInterval<E> newMaxCheck, final int newMaxIter,
                                                 final BracketedRealFieldUnivariateSolver<E> newSolver,
                                                 final FieldODEEventHandler<E> newHandler) {
        return new FieldEventSlopeFilter<T, E>(newSolver.getAbsoluteAccuracy().getField(), newMaxCheck, newMaxIter,
                                               newSolver, newHandler, rawDetector, filter);
    }

    /**
     * Get the wrapped raw detector.
     * @return the wrapped raw detector
     */
    public T getDetector() {
        return rawDetector;
    }

    /**  {@inheritDoc} */
    @Override
    public void init(final FieldODEStateAndDerivative<E> initialState, E finalTime) {

        // delegate to raw handler
        rawDetector.init(initialState, finalTime);

        // initialize events triggering logic
        forward  = finalTime.subtract(initialState.getTime()).getReal() >= 0;
        extremeT = finalTime.getField().getZero().newInstance(forward ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY);
        Arrays.fill(transformers, Transformer.UNINITIALIZED);
        Arrays.fill(updates, extremeT);

    }

    /**  {@inheritDoc} */
    @Override
    public E g(final FieldODEStateAndDerivative<E> state) {

        final E rawG = rawDetector.g(state);

        // search which transformer should be applied to g
        if (forward) {
            final int last = transformers.length - 1;
            if (extremeT.subtract(state.getTime()).getReal() < 0) {
                // we are at the forward end of the history

                // check if a new rough root has been crossed
                final Transformer previous = transformers[last];
                final Transformer next     = filter.selectTransformer(previous, rawG.getReal(), forward);
                if (next != previous) {
                    // there is a root somewhere between extremeT and t.
                    // the new transformer is valid for t (this is how we have just computed
                    // it above), but it is in fact valid on both sides of the root, so
                    // it was already valid before t and even up to previous time. We store
                    // the switch at extremeT for safety, to ensure the previous transformer
                    // is not applied too close of the root
                    System.arraycopy(updates,      1, updates,      0, last);
                    System.arraycopy(transformers, 1, transformers, 0, last);
                    updates[last]      = extremeT;
                    transformers[last] = next;
                }

                extremeT = state.getTime();

                // apply the transform
                return next.transformed(rawG);

            } else {
                // we are in the middle of the history

                // select the transformer
                for (int i = last; i > 0; --i) {
                    if (updates[i].subtract(state.getTime()).getReal() <= 0) {
                        // apply the transform
                        return transformers[i].transformed(rawG);
                    }
                }

                return transformers[0].transformed(rawG);

            }
        } else {
            if (state.getTime().subtract(extremeT).getReal() < 0) {
                // we are at the backward end of the history

                // check if a new rough root has been crossed
                final Transformer previous = transformers[0];
                final Transformer next     = filter.selectTransformer(previous, rawG.getReal(), forward);
                if (next != previous) {
                    // there is a root somewhere between extremeT and t.
                    // the new transformer is valid for t (this is how we have just computed
                    // it above), but it is in fact valid on both sides of the root, so
                    // it was already valid before t and even up to previous time. We store
                    // the switch at extremeT for safety, to ensure the previous transformer
                    // is not applied too close of the root
                    System.arraycopy(updates,      0, updates,      1, updates.length - 1);
                    System.arraycopy(transformers, 0, transformers, 1, transformers.length - 1);
                    updates[0]      = extremeT;
                    transformers[0] = next;
                }

                extremeT = state.getTime();

                // apply the transform
                return next.transformed(rawG);

            } else {
                // we are in the middle of the history

                // select the transformer
                for (int i = 0; i < updates.length - 1; ++i) {
                    if (state.getTime().subtract(updates[i]).getReal() <= 0) {
                        // apply the transform
                        return transformers[i].transformed(rawG);
                    }
                }

                return transformers[updates.length - 1].transformed(rawG);

            }
       }

    }

    /** Local handler.
     * @param <T> type of the event detector
     * @param <E> the type of the field elements
     */
    private static class LocalHandler<T extends FieldODEEventDetector<E>, E extends CalculusFieldElement<E>>
        implements FieldODEEventHandler<E> {

        /** Raw handler. */
        private final FieldODEEventHandler<E> rawHandler;

        /** Simple constructor.
         * @param rawHandler raw handler
         */
        LocalHandler(final FieldODEEventHandler<E> rawHandler) {
            this.rawHandler = rawHandler;
        }

        /**  {@inheritDoc} */
        @Override
        public Action eventOccurred(final FieldODEStateAndDerivative<E> state,
                                    final FieldODEEventDetector<E> detector,
                                    final boolean increasing) {
            // delegate to raw handler, fixing increasing status on the fly
            @SuppressWarnings("unchecked")
            final FieldEventSlopeFilter<T, E> esf = (FieldEventSlopeFilter<T, E>) detector;
            return rawHandler.eventOccurred(state, esf, esf.filter.isTriggeredOnIncreasing());
        }

        /**  {@inheritDoc} */
        @Override
        public FieldODEState<E> resetState(final FieldODEEventDetector<E> detector,
                                           final FieldODEStateAndDerivative<E> state) {
            // delegate to raw handler
            @SuppressWarnings("unchecked")
            final FieldEventSlopeFilter<T, E> esf = (FieldEventSlopeFilter<T, E>) detector;
            return rawHandler.resetState(esf, state);
        }

    }

}
