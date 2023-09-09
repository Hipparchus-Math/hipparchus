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
import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver;
import org.hipparchus.analysis.solvers.FieldBracketingNthOrderBrentSolver;
import org.hipparchus.ode.FieldODEStateAndDerivative;

/** Base class for #@link {@link FieldODEEventDetector}.
 * @param <T> type of the detector
 * @param <E> type of the field elements
 * @since 3.0
 */
public abstract class AbstractFieldODEDetector<T extends AbstractFieldODEDetector<T, E>, E extends CalculusFieldElement<E>>
    implements FieldODEEventDetector<E> {

    /** Default maximum checking interval (s). */
    public static final double DEFAULT_MAXCHECK = 600;

    /** Default convergence threshold (s). */
    public static final double DEFAULT_THRESHOLD = 1.e-6;

    /** Default maximum number of iterations in the event time search. */
    public static final int DEFAULT_MAX_ITER = 100;

    /** Max check interval. */
    private final FieldAdaptableInterval<E> maxCheck;

    /** Maximum number of iterations in the event time search. */
    private final int maxIter;

    /** Root-finding algorithm to use to detect state events. */
    private final BracketedRealFieldUnivariateSolver<E> solver;

    /** Default handler for event overrides. */
    private final FieldODEEventHandler<E> handler;

    /** Propagation direction. */
    private boolean forward;

    /** Build a new instance.
     * @param maxCheck maximum checking interval, must be strictly positive (s)
     * @param maxIter maximum number of iterations in the event time search
     * @param solver root-finding algorithm to use to detect state events
     * @param handler event handler to call at event occurrences
     */
    protected AbstractFieldODEDetector(final FieldAdaptableInterval<E> maxCheck, final int maxIter,
                                       final BracketedRealFieldUnivariateSolver<E> solver,
                                       final FieldODEEventHandler<E> handler) {
        this.maxCheck  = maxCheck;
        this.maxIter   = maxIter;
        this.solver    = solver;
        this.handler   = handler;
        this.forward   = true;
    }

    /**
     * {@inheritDoc}
     *
     * <p> This implementation sets the direction of integration and initializes the event
     * handler. If a subclass overrides this method it should call {@code
     * super.init(s0, t)}.
     */
    @Override
    public void init(final FieldODEStateAndDerivative<E> s0, final E t) {
        forward = t.subtract(s0.getTime()).getReal() >= 0;
        getHandler().init(s0, t, this);
    }

    /** {@inheritDoc} */
    @Override
    public abstract E g(FieldODEStateAndDerivative<E> s);

    /** {@inheritDoc} */
    @Override
    public FieldAdaptableInterval<E> getMaxCheckInterval() {
        return maxCheck;
    }

    /** {@inheritDoc} */
    @Override
    public int getMaxIterationCount() {
        return maxIter;
    }

    /** {@inheritDoc} */
    @Override
    public BracketedRealFieldUnivariateSolver<E> getSolver() {
        return solver;
    }

    /**
     * Setup the maximum checking interval.
     * <p>
     * This will override a maximum checking interval if it has been configured previously.
     * </p>
     * @param newMaxCheck maximum checking interval (s)
     * @return a new detector with updated configuration (the instance is not changed)
     */
    public T withMaxCheck(final E newMaxCheck) {
        return withMaxCheck(s -> newMaxCheck.getReal());
    }

    /**
     * Setup the maximum checking interval.
     * <p>
     * This will override a maximum checking interval if it has been configured previously.
     * </p>
     * @param newMaxCheck maximum checking interval (s)
     * @return a new detector with updated configuration (the instance is not changed)
     * @since 3.0
     */
    public T withMaxCheck(final FieldAdaptableInterval<E> newMaxCheck) {
        return create(newMaxCheck, getMaxIterationCount(), getSolver(), getHandler());
    }

    /**
     * Setup the maximum number of iterations in the event time search.
     * <p>
     * This will override a number of iterations if it has been configured previously.
     * </p>
     * @param newMaxIter maximum number of iterations in the event time search
     * @return a new detector with updated configuration (the instance is not changed)
     */
    public T withMaxIter(final int newMaxIter) {
        return create(getMaxCheckInterval(), newMaxIter, getSolver(), getHandler());
    }

    /**
     * Setup the convergence threshold.
     * <p>
     * This is equivalent to call {@code withSolver(new FieldBracketingNthOrderBrentSolver<>(zero,
     * newThreshold, zero, 5)}, so it will override a solver if one has been configured previously.
     * </p>
     * @param newThreshold convergence threshold
     * @return a new detector with updated configuration (the instance is not changed)
     * @see #withSolver(BracketedRealFieldUnivariateSolver)
     */
    public T withThreshold(final E newThreshold) {
        final E zero = newThreshold.getField().getZero();
        return withSolver(new FieldBracketingNthOrderBrentSolver<>(zero, newThreshold, zero, 5));
    }

    /**
     * Setup the root-finding algorithm to use to detect state events.
     * <p>
     * This will override a solver if it has been configured previously.
     * </p>
     * @param newSolver root-finding algorithm to use to detect state events
     * @return a new detector with updated configuration (the instance is not changed)
     * @see #withThreshold(CalculusFieldElement)
     */
    public T withSolver(final BracketedRealFieldUnivariateSolver<E> newSolver) {
        return create(getMaxCheckInterval(), getMaxIterationCount(), newSolver, getHandler());
    }

    /**
     * Setup the event handler to call at event occurrences.
     * <p>
     * This will override a handler if it has been configured previously.
     * </p>
     * @param newHandler event handler to call at event occurrences
     * @return a new detector with updated configuration (the instance is not changed)
     */
    public T withHandler(final FieldODEEventHandler<E> newHandler) {
        return create(getMaxCheckInterval(), getMaxIterationCount(), getSolver(), newHandler);
    }

    /** {@inheritDoc} */
    @Override
    public FieldODEEventHandler<E> getHandler() {
        return handler;
    }

    /** Build a new instance.
     * @param newMaxCheck maximum checking interval
     * @param newMaxIter maximum number of iterations in the event time search
     * @param newSolver root-finding algorithm to use to detect state events
     * @param newHandler event handler to call at event occurrences
     * @return a new instance of the appropriate sub-type
     */
    protected abstract T create(FieldAdaptableInterval<E> newMaxCheck, int newMaxIter,
                                BracketedRealFieldUnivariateSolver<E> newSolver,
                                FieldODEEventHandler<E> newHandler);

    /** Check if the current propagation is forward or backward.
     * @return true if the current propagation is forward
     */
    public boolean isForward() {
        return forward;
    }

}
