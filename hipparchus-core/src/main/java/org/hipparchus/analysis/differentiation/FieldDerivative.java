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
package org.hipparchus.analysis.differentiation;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;

/** Interface representing both the value and the differentials of a function.
 * @param <S> the type of the field elements
 * @param <T> the type of the function derivative
 * @see Derivative
 * @since 1.7
 */
public interface FieldDerivative<S extends CalculusFieldElement<S>, T extends FieldDerivative<S, T>>
        extends CalculusFieldElement<T>, DifferentialAlgebra {

    /** Get the value part of the function.
     * @return value part of the value of the function
     */
    S getValue();

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
    S getPartialDerivative(int ... orders)
        throws MathIllegalArgumentException;

    /** {@inheritDoc} */
    @Override
    default double getReal() {
        return getValue().getReal();
    }

    /** Create an instance corresponding to a constant Field value.
     * <p>
     * This default implementation is there so that no API gets broken
     * by the next release, which is not a major one. Custom inheritors
     * should probably overwrite it.
     * </p>
     * @param value constant value
     * @return instance corresponding to a constant Field value
     * @since 3.1
     */
    default T newInstance(S value) {
        return newInstance(value.getReal());
    }

    /** Create a new object with new value (zeroth-order derivative, as passed as input)
     * and same derivatives of order one and above.
     * <p>
     * This default implementation is there so that no API gets broken
     * by the next release, which is not a major one. Custom inheritors
     * should probably overwrite it.
     * </p>
     * @param value zeroth-order derivative of new represented function
     * @return new object with changed value
     * @since 3.1
     */
    default T withValue(S value) {
        return add(newInstance(value.subtract(getValue())));
    }

    /** '+' operator.
     * @param a right hand side parameter of the operator
     * @return this+a
     * @since 3.1
     */
    default T add(S a) {
        return withValue(getValue().add(a));
    }

    /** '-' operator.
     * @param a right hand side parameter of the operator
     * @return this-a
     * @since 3.1
     */
    default T subtract(S a) {
        return withValue(getValue().subtract(a));
    }

    /** {@inheritDoc} */
    @Override
    default T log10() {
        return log().divide(FastMath.log(10.));
    }

    /** {@inheritDoc} */
    @Override
    default T pow(T e) {
        return log().multiply(e).exp();
    }

    /** {@inheritDoc} */
    @Override
    default T cosh() {
        return (exp().add(negate().exp())).half();
    }

    /** {@inheritDoc} */
    @Override
    default T sinh() {
        return (exp().subtract(negate().exp())).half();
    }

    /** {@inheritDoc} */
    @Override
    default T acos() {
        return asin().negate().add(getPi().half());
    }

    /** {@inheritDoc} */
    @Override
    default T ceil() {
        return newInstance(getValue().ceil());
    }

    /** {@inheritDoc} */
    @Override
    default T floor() {
        return newInstance(getValue().floor());
    }

    /** {@inheritDoc} */
    @Override
    default T rint() {
        return newInstance(getValue().rint());
    }

    /** {@inheritDoc} */
    @Override
    default T ulp() {
        return newInstance(getValue().ulp());
    }

    /** {@inheritDoc} */
    @Override
    default T sign() {
        return newInstance(getValue().sign());
    }

    /** {@inheritDoc} */
    @Override
    default int getExponent() {
        return getValue().getExponent();
    }

}
