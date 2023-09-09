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

import java.lang.reflect.Array;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.ode.events.FieldODEEventDetector;
import org.hipparchus.util.MathArrays;

/**
 * This class is used as the base class of the problems that are
 * integrated during the junit tests for the ODE integrators.
 * @param <T> the type of the field elements
 */
public abstract class TestFieldProblemAbstract<T extends CalculusFieldElement<T>>
    implements FieldOrdinaryDifferentialEquation<T> {

    /** Number of functions calls. */
    private int calls;

    /** Initial state */
    private final FieldODEState<T> s0;

    /** Final time */
    private final T t1;

    /** Error scale */
    private final T[] errorScale;

    /**
     * Simple constructor.
     * @param t0 initial time
     * @param y0 initial state
     * @param t1 final time
     * @param errorScale error scale
     */
    protected TestFieldProblemAbstract(T t0, T[] y0, T t1, T[] errorScale) {
        calls      = 0;
        s0         = new FieldODEState<T>(t0, y0);
        this.t1    = t1;
        this.errorScale = errorScale.clone();
    }

    /** get the filed to which elements belong.
     * @return field to which elements belong
     */
    public Field<T> getField() {
        return s0.getTime().getField();
    }

    /** Get the problem dimension.
     * @return problem dimension
     */
    public int getDimension() {
        return s0.getPrimaryStateDimension();
    }

    /**
     * Get the initial time.
     * @return initial time
     */
    public T getInitialTime() {
        return s0.getTime();
    }

   /**
     * Get the initial state.
     * @return initial state
     */
    public FieldODEState<T> getInitialState() {
        return s0;
    }

    /**
     * Get the final time.
     * @return final time
     */
    public T getFinalTime() {
        return t1;
    }

    /**
     * Get the error scale.
     * @return error scale
     */
    public T[] getErrorScale() {
        return errorScale;
    }

    /**
     * Get the event detectors.
     * @param maxCheck maximum checking interval, must be strictly positive
     * @param threshold convergence threshold (s)
     * @param maxIter maximum number of iterations in the event time search
     * @return events detectors   */
    public FieldODEEventDetector<T>[] getEventDetectors(final double maxCheck, final T threshold, final int maxIter) {
        @SuppressWarnings("unchecked")
        final FieldODEEventDetector<T>[] empty =
                        (FieldODEEventDetector<T>[]) Array.newInstance(FieldODEEventDetector.class, 0);
        return empty;
    }

    /**
     * Get the theoretical events times.
     * @return theoretical events times
     */
    public T[] getTheoreticalEventsTimes() {
        return MathArrays.buildArray(s0.getTime().getField(), 0);
    }

    /**
     * Get the number of calls.
     * @return nuber of calls
     */
    public int getCalls() {
        return calls;
    }

    /** {@inheritDoc} */
    public void init(T t0, T[] y0, T t) {
    }

    /** {@inheritDoc} */
    public T[] computeDerivatives(T t, T[] y) {
        ++calls;
        return doComputeDerivatives(t, y);
    }

    abstract public T[] doComputeDerivatives(T t, T[] y);

    /**
     * Compute the theoretical state at the specified time.
     * @param t time at which the state is required
     * @return state vector at time t
     */
    abstract public T[] computeTheoreticalState(T t);

    /** Convert a double.
     * @param d double to convert
     * @return converted double
     */
    protected static <T extends CalculusFieldElement<T>> T convert(Field<T> field, double d) {
        return field.getZero().add(d);
    }

    /** Convert a one dimension array.
     * @param elements array elements
     * @return converted array
     */
    protected static <T extends CalculusFieldElement<T>> T[] convert(Field<T> field, double ... elements) {
        T[] array = MathArrays.buildArray(field, elements.length);
        for (int i = 0; i < elements.length; ++i) {
            array[i] = convert(field, elements[i]);
        }
        return array;
    }

}
