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
package org.hipparchus.analysis.interpolation;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * Interface representing a univariate field interpolating function.
 * @since 1.5
 */
public interface FieldUnivariateInterpolator {
    /**
     * Compute an interpolating function for the dataset.
     *
     * @param xval Arguments for the interpolation points.
     * @param yval Values for the interpolation points.
     * @param <T> the type of the field elements
     * @return a function which interpolates the dataset.
     * @throws MathIllegalArgumentException
     * if the arguments violate assumptions made by the interpolation
     * algorithm.
     * @throws MathIllegalArgumentException if arrays lengthes do not match
     */
    <T extends CalculusFieldElement<T>> CalculusFieldUnivariateFunction<T> interpolate(T[] xval, T[] yval)
        throws MathIllegalArgumentException;
}
