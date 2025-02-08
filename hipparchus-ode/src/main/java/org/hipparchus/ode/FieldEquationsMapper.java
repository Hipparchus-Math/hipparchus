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

import java.io.Serializable;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;

/**
 * Class mapping the part of a complete state or derivative that pertains
 * to a set of differential equations.
 * <p>
 * Instances of this class are guaranteed to be immutable.
 * </p>
 * @see FieldExpandableODE
 * @param <T> the type of the field elements
 */
public class FieldEquationsMapper<T extends CalculusFieldElement<T>> implements Serializable {

    /** Serializable UID. */
    private static final long serialVersionUID = 20151114L;

    /** Start indices of the components. */
    private final int[] start;

    /** Create a mapper by adding a new equation to another mapper.
     * <p>
     * The new equation will have index {@code mapper.}{@link #getNumberOfEquations()},
     * or 0 if {@code mapper} is null.
     * </p>
     * @param mapper former mapper, with one equation less (null for first equation)
     * @param dimension dimension of the equation state vector
     */
    FieldEquationsMapper(final FieldEquationsMapper<T> mapper, final int dimension) {
        final int index = (mapper == null) ? 0 : mapper.getNumberOfEquations();
        this.start = new int[index + 2];
        if (mapper == null) {
            start[0] = 0;
        } else {
            System.arraycopy(mapper.start, 0, start, 0, index + 1);
        }
        start[index + 1] = start[index] + dimension;
    }

    /** Get the number of equations mapped.
     * @return number of equations mapped
     */
    public int getNumberOfEquations() {
        return start.length - 1;
    }

    /** Return the dimension of the complete set of equations.
     * <p>
     * The complete set of equations correspond to the primary set plus all secondary sets.
     * </p>
     * @return dimension of the complete set of equations
     */
    public int getTotalDimension() {
        return start[start.length - 1];
    }

    /** Map flat arrays to a state and derivative.
     * @param t time
     * @param y state array to map, including primary and secondary components
     * @param yDot state derivative array to map, including primary and secondary components
     * @return mapped state
     * @exception MathIllegalArgumentException if an array does not match total dimension
     */
    public FieldODEStateAndDerivative<T> mapStateAndDerivative(final T t, final T[] y, final T[] yDot)
        throws MathIllegalArgumentException {

        if (y.length != getTotalDimension()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   y.length, getTotalDimension());
        }

        if (yDot.length != getTotalDimension()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   yDot.length, getTotalDimension());
        }

        final int n = getNumberOfEquations();
        final T[] state      = extractEquationData(0, y);
        final T[] derivative = extractEquationData(0, yDot);
        if (n < 2) {
            return new FieldODEStateAndDerivative<>(t, state, derivative);
        } else {
            final T[][] secondaryState      = MathArrays.buildArray(t.getField(), n - 1, -1);
            final T[][] secondaryDerivative = MathArrays.buildArray(t.getField(), n - 1, -1);
            for (int index = 1; index < getNumberOfEquations(); ++index) {
                secondaryState[index - 1]      = extractEquationData(index, y);
                secondaryDerivative[index - 1] = extractEquationData(index, yDot);
            }
            return new FieldODEStateAndDerivative<>(t, state, derivative, secondaryState, secondaryDerivative);
        }
    }

    /** Extract equation data from a complete state or derivative array.
     * @param index index of the equation, must be between 0 included and
     * {@link #getNumberOfEquations()} (excluded)
     * @param complete complete state or derivative array from which
     * equation data should be retrieved
     * @return equation data
     * @exception MathIllegalArgumentException if index is out of range
     * @exception MathIllegalArgumentException if complete state has not enough elements
     */
    public T[] extractEquationData(final int index, final T[] complete)
        throws MathIllegalArgumentException {
        checkIndex(index);
        final int begin     = start[index];
        final int end       = start[index + 1];
        if (complete.length < end) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   complete.length, end);
        }
        final int dimension = end - begin;
        final T[] equationData = MathArrays.buildArray(complete[0].getField(), dimension);
        System.arraycopy(complete, begin, equationData, 0, dimension);
        return equationData;
    }

    /** Insert equation data into a complete state or derivative array.
     * @param index index of the equation, must be between 0 included and
     * {@link #getNumberOfEquations()} (excluded)
     * @param equationData equation data to be inserted into the complete array
     * @param complete placeholder where to put equation data (only the
     * part corresponding to the equation will be overwritten)
     * @exception MathIllegalArgumentException if either array has not enough elements
     */
    public void insertEquationData(final int index, T[] equationData, T[] complete)
        throws MathIllegalArgumentException {
        checkIndex(index);
        final int begin     = start[index];
        final int end       = start[index + 1];
        final int dimension = end - begin;
        if (complete.length < end) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   complete.length, end);
        }
        if (equationData.length != dimension) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   equationData.length, dimension);
        }
        System.arraycopy(equationData, 0, complete, begin, dimension);
    }

    /** Check equation index.
     * @param index index of the equation, must be between 0 included and
     * {@link #getNumberOfEquations()} (excluded)
     * @exception MathIllegalArgumentException if index is out of range
     */
    private void checkIndex(final int index) throws MathIllegalArgumentException {
        MathUtils.checkRangeInclusive(index, 0, start.length - 2);
    }

}
