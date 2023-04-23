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
package org.hipparchus;

import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.exception.NullArgumentException;


/**
 * Interface representing <a href="http://mathworld.wolfram.com/Field.html">field</a> elements.
 * @param <T> the type of the field elements
 * @see Field
 */
public interface FieldElement<T extends FieldElement<T>> {

    /** Get the real value of the number.
     * @return real value
     */
    double getReal();

    /** Compute this + a.
     * @param a element to add
     * @return a new element representing this + a
     * @throws NullArgumentException if {@code a} is {@code null}.
     */
    T add(T a) throws NullArgumentException;

    /** Compute this - a.
     * @param a element to subtract
     * @return a new element representing this - a
     * @throws NullArgumentException if {@code a} is {@code null}.
     */
    T subtract(T a) throws NullArgumentException;

    /**
     * Returns the additive inverse of {@code this} element.
     * @return the opposite of {@code this}.
     */
    T negate();

    /** Compute n &times; this. Multiplication by an integer number is defined
     * as the following sum
     * \[
     * n \times \mathrm{this} = \sum_{i=1}^n \mathrm{this}
     * \]
     * @param n Number of times {@code this} must be added to itself.
     * @return A new element representing n &times; this.
     */
    T multiply(int n);

    /** Compute this &times; a.
     * @param a element to multiply
     * @return a new element representing this &times; a
     * @throws NullArgumentException if {@code a} is {@code null}.
     */
    T multiply(T a) throws NullArgumentException;

    /** Compute this &divide; a.
     * @param a element to divide by
     * @return a new element representing this &divide; a
     * @throws NullArgumentException if {@code a} is {@code null}.
     * @throws MathRuntimeException if {@code a} is zero
     */
    T divide(T a) throws NullArgumentException, MathRuntimeException;

    /**
     * Returns the multiplicative inverse of {@code this} element.
     * @return the inverse of {@code this}.
     * @throws MathRuntimeException if {@code this} is zero
     */
    T reciprocal() throws MathRuntimeException;

    /** Get the {@link Field} to which the instance belongs.
     * @return {@link Field} to which the instance belongs
     */
    Field<T> getField();

    /** Check if an element is semantically equal to zero.
     * <p>
     * The default implementation simply calls {@code equals(getField().getZero())}.
     * However, this may need to be overridden in some cases as due to
     * compatibility with {@code hashCode()} some classes implements
     * {@code equals(Object)} in such a way that -0.0 and +0.0 are different,
     * which may be a problem. It prevents for example identifying a diagonal
     * element is zero and should be avoided when doing partial pivoting in
     * LU decomposition.
     * </p>
     * @return true if the element is semantically equal to zero
     * @since 1.8
     */
    default boolean isZero() {
        return equals(getField().getZero());
    }

}
