/*
 * Licensed to the Hipparchus project under one or more
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
package org.hipparchus.analysis.interpolation;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.analysis.CalculusFieldTrivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * Interface representing a trivariate field interpolating function where the
 * sample points must be specified on a regular grid.
 * @param <T> type of the field elements
 * @since 4.1
 */
@FunctionalInterface
public interface FieldTrivariateGridInterpolator<T extends CalculusFieldElement<T>> {
    /**
     * Compute an interpolating function for the dataset.
     *
     * @param xval All the x-coordinates of the interpolation points, sorted
     * in increasing order.
     * @param yval All the y-coordinates of the interpolation points, sorted
     * in increasing order.
     * @param zval All the z-coordinates of the interpolation points, sorted
     * in increasing order.
     * @param fval The values of the interpolation points on all the grid knots:
     * {@code fval[i][j][k] = f(xval[i], yval[j], zval[k])}.
     * @return a function which interpolates the dataset.
     * @throws MathIllegalArgumentException if any of the arrays has zero length.
     * @throws MathIllegalArgumentException if the array lengths are inconsistent.
     * @throws MathIllegalArgumentException if the array is not sorted.
     * @throws MathIllegalArgumentException if the number of points is too small for
     * the order of the interpolation
     */
    CalculusFieldTrivariateFunction<T> interpolate(T[] xval, T[] yval, T[] zval, T[][][] fval)
        throws MathIllegalArgumentException;
}
