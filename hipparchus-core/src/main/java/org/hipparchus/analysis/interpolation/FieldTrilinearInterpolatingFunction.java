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
import org.hipparchus.analysis.CalculusFieldTrivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.MathArrays;

/**
 * Interpolate grid data using tri-linear interpolation.
 * <p>
 * This interpolator is thread-safe.
 * </p>
 * @param <T> Type of the field elements.
 * @since 4.1
 */
public class FieldTrilinearInterpolatingFunction<T extends CalculusFieldElement<T>>
    implements CalculusFieldTrivariateFunction<T>
{
    /** Grid along the x axis. */
    private final FieldGridAxis<T> xGrid;

    /** Grid along the y axis. */
    private final FieldGridAxis<T> yGrid;

    /** Grid along the z axis. */
    private final FieldGridAxis<T> zGrid;

    /** Values of the interpolation points on all the grid knots */
    private final T[][][] fVal;

    /**
     * Simple constructor.
     * @param xVal All the x-coordinates of the interpolation points, sorted
     * in increasing order.
     * @param yVal All the y-coordinates of the interpolation points, sorted
     * in increasing order.
     * @param zVal All the z-coordinates of the interpolation points, sorted
     * in increasing order.
     * @param fVal The values of the interpolation points on all the grid knots:
     * {@code fVal[i][j][k] = f(xVal[i], yVal[j], zVal[k])}.
     * @exception MathIllegalArgumentException if grid size is smaller than 2
     * or if the grid is not sorted in strict increasing order
     */
    public FieldTrilinearInterpolatingFunction(final T[] xVal, final T[] yVal, final T[] zVal, final T[][][] fVal)
        throws MathIllegalArgumentException {
        final Field<T> field = fVal[0][0][0].getField();
        this.xGrid = new FieldGridAxis<>(xVal, 2);
        this.yGrid = new FieldGridAxis<>(yVal, 2);
        this.zGrid = new FieldGridAxis<>(zVal, 2);

        this.fVal = MathArrays.buildArray(field, xVal.length, yVal.length, zVal.length);
        for (int i = 0; i < xVal.length; i++) {
            for (int j = 0; j < yVal.length; j++) {
                System.arraycopy(fVal[i][j], 0, this.fVal[i][j], 0, zVal.length);
            }
        }
    }

    /**
     * Get the lowest grid x coordinate.
     *
     * @return lowest grid x coordinate
     */
    public T getXInf() {
        return xGrid.node(0);
    }

    /**
     * Get the highest grid x coordinate.
     *
     * @return highest grid x coordinate
     */
    public T getXSup() {
        return xGrid.node(xGrid.size() - 1);
    }

    /**
     * Get the lowest grid y coordinate.
     *
     * @return lowest grid y coordinate
     */
    public T getYInf() {
        return yGrid.node(0);
    }

    /**
     * Get the highest grid y coordinate.
     *
     * @return highest grid y coordinate
     */
    public T getYSup() {
        return yGrid.node(yGrid.size() - 1);
    }

    /**
     * Get the lowest grid z coordinate.
     *
     * @return lowest grid z coordinate
     */
    public T getZInf() {
        return zGrid.node(0);
    }

    /**
     * Get the highest grid z coordinate.
     *
     * @return highest grid z coordinate
     */
    public T getZSup() {
        return zGrid.node(zGrid.size() - 1);
    }

    /** {@inheritDoc} */
    @Override
    public T value(T x, T y, T z) {
        // get the interpolation nodes
        final int i = xGrid.interpolationIndex(x);
        final int j = yGrid.interpolationIndex(y);
        final int k = zGrid.interpolationIndex(z);
        final T x0  = xGrid.node(i);
        final T x1  = xGrid.node(i + 1);
        final T y0  = yGrid.node(j);
        final T y1  = yGrid.node(j + 1);
        final T z0  = zGrid.node(k);
        final T z1  = zGrid.node(k + 1);

        // get the function values at interpolation nodes
        final T c000 = fVal[i][j][k];
        final T c100 = fVal[i + 1][j][k];
        final T c010 = fVal[i][j + 1][k];
        final T c110 = fVal[i + 1][j + 1][k];
        final T c001 = fVal[i][j][k + 1];
        final T c101 = fVal[i + 1][j][k + 1];
        final T c011 = fVal[i][j + 1][k + 1];
        final T c111 = fVal[i + 1][j + 1][k + 1];

        // interpolate
        final T dx0      = x.subtract(x0);
        final T dx1      = x1.subtract(x);
        final T dx10     = x1.subtract(x0);
        final T dy0      = y.subtract(y0);
        final T dy1      = y1.subtract(y);
        final T dy10     = y1.subtract(y0);
        final T dx10dy10 = dx10.multiply(dy10);
        final T c0 = dy0.multiply(c110).add(dy1.multiply(c100)).multiply(dx0).
                add(dy0.multiply(c010).add(dy1.multiply(c000)).multiply(dx1)).
                divide(dx10dy10);
        final T c1 = dy0.multiply(c111).add(dy1.multiply(c101)).multiply(dx0).
                add(dy0.multiply(c011).add(dy1.multiply(c001)).multiply(dx1)).
                divide(dx10dy10);

        // interpolate along z
        final T dz0  = z.subtract(z0);
        final T dz10 = z1.subtract(z0);
        final T dc10 = c1.subtract(c0);

        return c0.add(dz0.divide(dz10).multiply(dc10));
    }
}
