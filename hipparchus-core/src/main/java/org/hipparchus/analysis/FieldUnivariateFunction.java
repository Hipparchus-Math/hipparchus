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
package org.hipparchus.analysis;

import org.hipparchus.Field;
import org.hipparchus.RealFieldElement;

/**
 * An interface representing a univariate real function for any field type.
 * <p>
 * This interface is more general than {@link RealFieldUnivariateFunction} because
 * the same instance can accept any field type, not just one.
 * </p>
 * @see UnivariateFunction
 * @see RealFieldUnivariateFunction
 * @since 1.3
 */
public interface FieldUnivariateFunction {

    /** Convert to a {@link RealFieldUnivariateFunction} with a specific type.
     * @param <T> the type of the field elements
     * @param field field for the argument and value
     * @return converted function
     */
    default <T extends RealFieldElement<T>> RealFieldUnivariateFunction<T> toRealFieldUnivariateFunction(Field<T> field) {
        return this::value;
    }

    /**
     * Compute the value of the function.
     *
     * @param <T> the type of the field elements
     * @param x Point at which the function value should be computed.
     * @return the value of the function.
     * @throws IllegalArgumentException when the activated method itself can
     * ascertain that a precondition, specified in the API expressed at the
     * level of the activated method, has been violated.
     * When Hipparchus throws an {@code IllegalArgumentException}, it is
     * usually the consequence of checking the actual parameters passed to
     * the method.
     */
    <T extends RealFieldElement<T>> T value(T x);

}
