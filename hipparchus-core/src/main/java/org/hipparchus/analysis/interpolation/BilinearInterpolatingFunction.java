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
package org.hipparchus.analysis.interpolation;

import java.io.Serializable;

import org.hipparchus.RealFieldElement;
import org.hipparchus.analysis.BivariateFunction;
import org.hipparchus.analysis.FieldBivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * Interpolate grid data using bi-linear interpolation.
 * <p>
 * This interpolator is thread-safe.
 * </p>
 * @since 1.4
 */
public class BilinearInterpolatingFunction implements BivariateFunction, FieldBivariateFunction, Serializable {

    /** Serializable UID. */
    private static final long serialVersionUID = 20180926L;

    /** Grid along the x axis. */
    private final GridAxis xGrid;

    /** Grid along the y axis. */
    private final GridAxis yGrid;

    /** Grid size along the y axis. */
    private final int ySize;

    /** Values of the interpolation points on all the grid knots (in a flatten array). */
    private final double[] fVal;

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
    public BilinearInterpolatingFunction(final double[] xVal, final double[] yVal,
                                         final double[][] fVal)
        throws MathIllegalArgumentException {
        this.xGrid = new GridAxis(xVal, 2);
        this.yGrid = new GridAxis(yVal, 2);
        this.ySize = yVal.length;
        this.fVal  = new double[xVal.length * ySize];
        int k = 0;
        for (int i = 0; i < xVal.length; ++i) {
            final double[] fi = fVal[i];
            for (int j = 0; j < ySize; ++j) {
                this.fVal[k++] = fi[j];
            }
        }
    }

    /** Get the lowest grid x coordinate.
     * @return lowest grid x coordinate
     */
    public double getXInf() {
        return xGrid.node(0);
    }

    /** Get the highest grid x coordinate.
     * @return highest grid x coordinate
     */
    public double getXSup() {
        return xGrid.node(xGrid.size() - 1);
    }

    /** Get the lowest grid y coordinate.
     * @return lowest grid y coordinate
     */
    public double getYInf() {
        return yGrid.node(0);
    }

    /** Get the highest grid y coordinate.
     * @return highest grid y coordinate
     */
    public double getYSup() {
        return yGrid.node(yGrid.size() - 1);
    }

    /** {@inheritDoc} */
    @Override
    public double value(final double x, final double y) {

        // get the interpolation nodes
        final int    i   = xGrid.interpolationIndex(x);
        final int    j   = yGrid.interpolationIndex(y);
        final double x0  = xGrid.node(i);
        final double x1  = xGrid.node(i + 1);
        final double y0  = yGrid.node(j);
        final double y1  = yGrid.node(j + 1);

        // get the function values at interpolation nodes
        final int    k0  = i * ySize + j;
        final int    k1  = k0 + ySize;
        final double z00 = fVal[k0];
        final double z01 = fVal[k0 + 1];
        final double z10 = fVal[k1];
        final double z11 = fVal[k1 + 1];

        // interpolate
        final double dx0  = x  - x0;
        final double dx1  = x1 - x;
        final double dx10 = x1 - x0;
        final double dy0  = y  - y0;
        final double dy1  = y1 - y;
        final double dy10 = y1 - y0;
        return (dx0 * (dy0 * z11 + dy1 * z10) + dx1 * (dy0 * z01 + dy1 * z00)) /
               (dx10 * dy10);

    }

    /** {@inheritDoc}
     * @since 1.5
     */
    @Override
    public <T extends RealFieldElement<T>> T value(T x, T y) {

        // get the interpolation nodes
        final int    i   = xGrid.interpolationIndex(x.getReal());
        final int    j   = yGrid.interpolationIndex(y.getReal());
        final double x0  = xGrid.node(i);
        final double x1  = xGrid.node(i + 1);
        final double y0  = yGrid.node(j);
        final double y1  = yGrid.node(j + 1);

        // get the function values at interpolation nodes
        final int    k0  = i * ySize + j;
        final int    k1  = k0 + ySize;
        final double z00 = fVal[k0];
        final double z01 = fVal[k0 + 1];
        final double z10 = fVal[k1];
        final double z11 = fVal[k1 + 1];

        // interpolate
        final T      dx0   = x.subtract(x0);
        final T      mdx1  = x.subtract(x1);
        final double dx10  = x1 - x0;
        final T      dy0   = y.subtract(y0);
        final T      mdy1  = y.subtract(y1);
        final double dy10  = y1 - y0;
        return          dy0.multiply(z11).subtract(mdy1.multiply(z10)).multiply(dx0).
               subtract(dy0.multiply(z01).subtract(mdy1.multiply(z00)).multiply(mdx1)).
               divide(dx10 * dy10);

    }

}
