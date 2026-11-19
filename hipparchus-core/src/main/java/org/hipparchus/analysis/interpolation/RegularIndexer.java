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

import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/**
 * Helper for sampled data along a regular grid.
 * @see AxisChecker
 * @since 4.1
 */
public class RegularIndexer {

    /** Minimum value. */
    private final double min;

    /** Maximum value. */
    private final double max;

    /** Number of points (including min and max). */
    private final int n;

    /** Step size. */
    private final double step;

    /** Simple constructor.
     * @param min minimum value
     * @param max maximum value
     * @param n number of points (including min and max)
     */
    public RegularIndexer(final double min, final double max, final int n) {
        this.min  = min;
        this.max  = max;
        this.n    = n;
        this.step = (max - min) / (n - 1);
    }

    /** Get the minimum coordinate.
     * @return minimum coordinate
     */
    public double getMin() {
        return min;
    }

    /** Get the maximum coordinate.
     * @return maximum coordinate
     */
    public double getMax() {
        return max;
    }

    /** Get the number of different coordinates.
     * @return number of different coordinates (including {@link #getMin()} and {@link #getMax()})
     */
    public int getN() {
        return n;
    }

    /** Get the step size between two successive coordinates.
     * @return step size between two successive coordinates
     */
    public double getStep() {
        return step;
    }

    /** Get index corresponding to coordinate.
     * <p>
     * The method rounds to the closest index, so it may be called with
     * {@code coordinate} slightly overshooting nominal range. If called
     * with {@code coordinate} between {@link #getMin()} - {@link #getStep()}/2
     * and {@link #getMin()} + {@link #getStep()}/2, it will return 0. If called
     * with {@code coordinate} between {@link #getMax()} - {@link #getStep()}/2
     * and {@link #getMax()} + {@link #getStep()}/2, it will return {@link #getN()} - 1.
     * If called outside of this extended range, it will throw an exception.
     * </p>
     * @param coordinate coordinate
     * @return index corresponding to coordinate (i.e. {@code 0} for {@code min}, {@code n-1} for {@code max})
     */
    public int index(final double coordinate) {
        final int index = (int) FastMath.rint((coordinate - min) / step);
        MathUtils.checkRangeInclusive(index, 0, n - 1);
        return index;
    }

    /** Get coordinates corresponding to index.
     * @param index index
     * @return coordinate (i.e. {@code min} for {@code 0}, {@code max} for {@code n-1})
     */
    public double coordinate(final int index) {
        MathUtils.checkRangeInclusive(index, 0, n - 1);
        // we count from the closest end to avoid numerical noise
        return index < n / 2 ?
               min + index * step :
               max - (n - 1 - index) * step;
    }

}
