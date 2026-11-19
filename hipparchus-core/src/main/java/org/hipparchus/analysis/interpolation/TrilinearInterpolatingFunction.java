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
import org.hipparchus.analysis.FieldTrivariateFunction;
import org.hipparchus.analysis.TrivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;

import java.io.Serializable;

/**
 * Interpolate grid data using tri-linear interpolation.
 * <p>
 * This interpolator is thread-safe.
 * </p>
 * @since 4.1
 */
public class TrilinearInterpolatingFunction implements TrivariateFunction, FieldTrivariateFunction, Serializable {

    /**
     * Serializable UID.
     */
    private static final long serialVersionUID = 20180926L;

    /**
     * Grid along the x axis.
     */
    private final GridAxis xGrid;

    /**
     * Grid along the y axis.
     */
    private final GridAxis yGrid;

    /**
     * Grid along the z axis.
     */
    private final GridAxis zGrid;

    /**
     * Values of the interpolation points on all the grid knots
     */
    private final double[][][] fVal;

    /**
     * Simple constructor.
     *
     * @param xVal All the x-coordinates of the interpolation points, sorted
     *             in increasing order.
     * @param yVal All the y-coordinates of the interpolation points, sorted
     *             in increasing order.
     * @param zVal All the z-coordinates of the interpolation points, sorted
     *             in increasing order.
     * @param fVal The values of the interpolation points on all the grid knots:
     *             {@code fVal[i][j][k] = f(xVal[i], yVal[j], zVal[k])}.
     * @throws MathIllegalArgumentException if grid size is smaller than 2
     *                                      or if the grid is not sorted in strict increasing order
     */
    public TrilinearInterpolatingFunction(final double[] xVal, final double[] yVal, final double[] zVal,
                                          final double[][][] fVal)
            throws MathIllegalArgumentException {
        this.xGrid = new GridAxis(xVal, 2);
        this.yGrid = new GridAxis(yVal, 2);
        this.zGrid = new GridAxis(zVal, 2);

        // deep copy of the array
        this.fVal = new double[xVal.length][yVal.length][zVal.length];
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
    public double getXInf() {
        return xGrid.node(0);
    }

    /**
     * Get the highest grid x coordinate.
     *
     * @return highest grid x coordinate
     */
    public double getXSup() {
        return xGrid.node(xGrid.size() - 1);
    }

    /**
     * Get the lowest grid y coordinate.
     *
     * @return lowest grid y coordinate
     */
    public double getYInf() {
        return yGrid.node(0);
    }

    /**
     * Get the highest grid y coordinate.
     *
     * @return highest grid y coordinate
     */
    public double getYSup() {
        return yGrid.node(yGrid.size() - 1);
    }

    /**
     * Get the lowest grid z coordinate.
     *
     * @return lowest grid z coordinate
     */
    public double getZInf() {
        return zGrid.node(0);
    }

    /**
     * Get the highest grid z coordinate.
     *
     * @return highest grid z coordinate
     */
    public double getZSup() {
        return zGrid.node(zGrid.size() - 1);
    }


    /** {@inheritDoc} */
    @Override
    public double value(final double x, final double y, final double z) {

        // get the interpolation nodes
        final int i = xGrid.interpolationIndex(x);
        final int j = yGrid.interpolationIndex(y);
        final int k = zGrid.interpolationIndex(z);
        final double x0 = xGrid.node(i);
        final double x1 = xGrid.node(i + 1);
        final double y0 = yGrid.node(j);
        final double y1 = yGrid.node(j + 1);
        final double z0 = zGrid.node(k);
        final double z1 = zGrid.node(k + 1);

        // get the function values at interpolation nodes
        final double c000 = fVal[i][j][k];
        final double c100 = fVal[i + 1][j][k];
        final double c010 = fVal[i][j + 1][k];
        final double c110 = fVal[i + 1][j + 1][k];
        final double c001 = fVal[i][j][k + 1];
        final double c101 = fVal[i + 1][j][k + 1];
        final double c011 = fVal[i][j + 1][k + 1];
        final double c111 = fVal[i + 1][j + 1][k + 1];

        // bilinear interpolations on (x, y)
        final double dx0  = x  - x0;
        final double dx1  = x1 - x;
        final double dx10 = x1 - x0;
        final double dy0  = y  - y0;
        final double dy1  = y1 - y;
        final double dy10 = y1 - y0;
        final double c0 = (dx0 * (dy0 * c110 + dy1 * c100) + dx1 * (dy0 * c010 + dy1 * c000)) / (dx10 * dy10);
        final double c1 = (dx0 * (dy0 * c111 + dy1 * c101) + dx1 * (dy0 * c011 + dy1 * c001)) / (dx10 * dy10);

        // interpolate along z
        final double t = (z  - z0) / (z1 - z0);
        return c0 + t * (c1 - c0);
    }

    /**
     * {@inheritDoc}
     *
     * @since 4.1
     */
    @Override
    public <T extends CalculusFieldElement<T>> T value(T x, T y, T z) {

        // get the interpolation nodes
        final int i = xGrid.interpolationIndex(x.getReal());
        final int j = yGrid.interpolationIndex(y.getReal());
        final int k = zGrid.interpolationIndex(z.getReal());
        final double x0 = xGrid.node(i);
        final double x1 = xGrid.node(i + 1);
        final double y0 = yGrid.node(j);
        final double y1 = yGrid.node(j + 1);
        final double z0 = zGrid.node(k);
        final double z1 = zGrid.node(k + 1);

        // get the function values at interpolation nodes
        final double c000 = fVal[i][j][k];
        final double c100 = fVal[i + 1][j][k];
        final double c010 = fVal[i][j + 1][k];
        final double c110 = fVal[i + 1][j + 1][k];
        final double c001 = fVal[i][j][k + 1];
        final double c101 = fVal[i + 1][j][k + 1];
        final double c011 = fVal[i][j + 1][k + 1];
        final double c111 = fVal[i + 1][j + 1][k + 1];

        // interpolate
        final T      dx0  = x.subtract(x0);
        final T      mdx1 = x.subtract(x1);
        final double dx10 = x1 - x0;
        final T      dy0  = y.subtract(y0);
        final T      mdy1 = y.subtract(y1);
        final double dy10 = y1 - y0;
        final T c0 = dy0.multiply(c110).subtract(mdy1.multiply(c100)).multiply(dx0).
                subtract(dy0.multiply(c010).subtract(mdy1.multiply(c000)).multiply(mdx1)).
                divide(dx10 * dy10);
        final T c1 = dy0.multiply(c111).subtract(mdy1.multiply(c101)).multiply(dx0).
                subtract(dy0.multiply(c011).subtract(mdy1.multiply(c001)).multiply(mdx1)).
                divide(dx10 * dy10);

        // interpolate along z
        final T t = z.subtract(z0).divide(z1 - z0);
        return c0.add(t.multiply(c1.subtract(c0)));
    }
}
