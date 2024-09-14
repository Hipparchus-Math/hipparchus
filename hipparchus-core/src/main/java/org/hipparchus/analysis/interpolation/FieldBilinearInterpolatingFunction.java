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
import org.hipparchus.Field;
import org.hipparchus.analysis.CalculusFieldBivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.MathArrays;

/**
 * Interpolate grid data using bi-linear interpolation.
 * <p>
 * This interpolator is thread-safe.
 * </p>
 * @param <T> Type of the field elements.
 * @since 4.0
 */
public class FieldBilinearInterpolatingFunction<T extends CalculusFieldElement<T>>
    implements CalculusFieldBivariateFunction<T> {

    /** Grid along the x axis. */
    private final FieldGridAxis<T> xGrid;

    /** Grid along the y axis. */
    private final FieldGridAxis<T> yGrid;

    /** Grid size along the y axis. */
    private final int ySize;

    /** Values of the interpolation points on all the grid knots (in a flatten array). */
    private final T[] fVal;

    /** Simple constructor.
     * @param xVal All the x-coordinates of the interpolation points, sorted
     * in increasing order.
     * @param yVal All the y-coordinates of the interpolation points, sorted
     * in increasing order.
     * @param fVal The values of the interpolation points on all the grid knots:
     * {@code fVal[i][j] = f(xVal[i], yVal[j])}.
     * @exception MathIllegalArgumentException if grid size is smaller than 2
     * or if the grid is not sorted in strict increasing order
     */
    public FieldBilinearInterpolatingFunction(final T[] xVal, final T[] yVal, final T[][] fVal)
        throws MathIllegalArgumentException {
        final Field<T> field = fVal[0][0].getField();
        this.xGrid = new FieldGridAxis<>(xVal, 2);
        this.yGrid = new FieldGridAxis<>(yVal, 2);
        this.ySize = yVal.length;
        this.fVal  = MathArrays.buildArray(field, xVal.length * ySize);
        int k = 0;
        for (int i = 0; i < xVal.length; ++i) {
            final T[] fi = fVal[i];
            for (int j = 0; j < ySize; ++j) {
                this.fVal[k++] = fi[j];
            }
        }
    }

    /** Get the lowest grid x coordinate.
     * @return lowest grid x coordinate
     */
    public T getXInf() {
        return xGrid.node(0);
    }

    /** Get the highest grid x coordinate.
     * @return highest grid x coordinate
     */
    public T getXSup() {
        return xGrid.node(xGrid.size() - 1);
    }

    /** Get the lowest grid y coordinate.
     * @return lowest grid y coordinate
     */
    public T getYInf() {
        return yGrid.node(0);
    }

    /** Get the highest grid y coordinate.
     * @return highest grid y coordinate
     */
    public T getYSup() {
        return yGrid.node(yGrid.size() - 1);
    }

    /** {@inheritDoc} */
    @Override
    public T value(T x, T y) {

        // get the interpolation nodes
        final int i   = xGrid.interpolationIndex(x);
        final int j   = yGrid.interpolationIndex(y);
        final T   x0  = xGrid.node(i);
        final T   x1  = xGrid.node(i + 1);
        final T   y0  = yGrid.node(j);
        final T   y1  = yGrid.node(j + 1);

        // get the function values at interpolation nodes
        final int k0  = i * ySize + j;
        final int k1  = k0 + ySize;
        final T   z00 = fVal[k0];
        final T   z01 = fVal[k0 + 1];
        final T   z10 = fVal[k1];
        final T   z11 = fVal[k1 + 1];

        // interpolate
        final T   dx0  = x.subtract(x0);
        final T   mdx1 = x.subtract(x1);
        final T   dx10 = x1.subtract(x0);
        final T   dy0  = y.subtract(y0);
        final T   mdy1 = y.subtract(y1);
        final T   dy10 = y1.subtract(y0);

        return          dy0.multiply(z11).subtract(mdy1.multiply(z10)).multiply(dx0).
               subtract(dy0.multiply(z01).subtract(mdy1.multiply(z00)).multiply(mdx1)).
               divide(dx10.multiply(dy10));

    }

}
