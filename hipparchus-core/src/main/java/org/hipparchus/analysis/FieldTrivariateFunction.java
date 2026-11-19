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
package org.hipparchus.analysis;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;

/**
 * An interface representing a trivariate field function.
 * @since 4.1
 */
@FunctionalInterface
public interface FieldTrivariateFunction
{

    /** Convert to a {@link CalculusFieldTrivariateFunction} with a specific type.
     * @param <T> the type of the field elements
     * @param field field for the argument and value
     * @return converted function
     */
    default <T extends CalculusFieldElement<T>> CalculusFieldTrivariateFunction<T> toCalculusFieldTrivariateFunction(Field<T> field) {
        return this::value;
    }

    /**
     * Compute the value for the function.
     *
     * @param x x-coordinate of the interpolation point.
     * @param y y-coordinate of the interpolation point.
     * @param z z-coordinate of the interpolation point.
     * @param <T> type of the field elements
     * @return the value.
     */
    <T extends CalculusFieldElement<T>> T value(T x, T y, T z);

}
