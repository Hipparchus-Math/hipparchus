/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.analysis.differentiation;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Precision;

/** Interface representing both the value and the differentials of a function.
 * @param <T> the type of the field elements
 * @since 1.7
 */
public interface Derivative<T extends CalculusFieldElement<T>> extends CalculusFieldElement<T> {

    /** Get the number of free parameters.
     * @return number of free parameters
     */
    int getFreeParameters();

    /** Get the derivation order.
     * @return derivation order
     */
    int getOrder();

    /** Get the value part of the function.
     * @return value part of the value of the function
     */
    double getValue();

    /** Get a partial derivative.
     * @param orders derivation orders with respect to each variable (if all orders are 0,
     * the value is returned)
     * @return partial derivative
     * @see #getValue()
     * @exception MathIllegalArgumentException if the numbers of variables does not
     * match the instance
     * @exception MathIllegalArgumentException if sum of derivation orders is larger
     * than the instance limits
     */
    double getPartialDerivative(int ... orders)
        throws MathIllegalArgumentException;

    /** Compute composition of the instance by a univariate function.
     * @param f array of value and derivatives of the function at
     * the current point (i.e. [f({@link #getValue()}),
     * f'({@link #getValue()}), f''({@link #getValue()})...]).
     * @return f(this)
     * @exception MathIllegalArgumentException if the number of derivatives
     * in the array is not equal to {@link #getOrder() order} + 1
     */
    T compose(double... f)
        throws MathIllegalArgumentException;

    /** Check if all derivatives are null.
     * <p>
     * In practice, it checks if the sum of the absolute values of derivatives is lower than a user-defined
     * tolerance to decide if they are "null" or not.
     *
     * @param tolerance tolerance on the sum of absolute values of derivatives
     * @return true if sum of all absolute values of derivatives is lower than tolerance, false otherwise
     */
    boolean hasNullDerivatives(double tolerance);

    /** Check if all derivatives are null.
     * <p>
     * Realization of method {@link #hasNullDerivatives(double)} with a tolerance of {@link Precision#SAFE_MIN}
     * @return true if all derivatives are null, false otherwise
     */
    default boolean hasNullDerivatives() {
        return hasNullDerivatives(Precision.SAFE_MIN);
    }
}
