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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;

import java.util.function.ToDoubleFunction;

/**
 * Checker for sampled data along a regular grid.
 * <p>
 * This class is intended to check one axis of n-dimensional grid data.
 * </p>
 * <p>
 * A typical use case is loading data from a file with entries of the form:
 * </p>
 * <pre>
 *   x₁ y₁ valueA₁₁ valueB₁₁
 *   x₁ y₂ valueA₁₂ valueB₁₂
 *   ⋮
 *   x₁ yₙ valueA₁ₙ valueB₁ₙ
 *   x₂ y₁ valueA₂₁ valueB₂₁
 *   x₂ y₂ valueA₂₂ valueB₂₂
 *   ⋮
 *   x₂ yₙ valueA₂ₙ valueB₂ₙ
 *   x₃ yₙ valueA₃ₙ valueB₃ₙ
 *   ⋮
 *   xₘ yₙ valueAₘₙ valueBₘₙ
 * </pre>
 * <p>
 * In this case, the grid is two-dimensional, with the x-axis having m steps
 *  ranging from x₁ to xₘ and the y-axis having n steps ranging from y₁ to yₙ.
 * </p>
 * <p>
 * In this case, users could store the entries in some custom container
 * class {@code Entry} with getters {@code getX()}, {@code getY()}, {@code getValueA()}, and
 * {@code getValueB()}. Then the checks for both x and y axes could be built as follows:
 * </p>
 * <pre>{@code
 *   List<Entry>    data     = readData();
 *   RegularIndexer xIndexer = new AxisChecker<>(Entry::getX(), toleranceX).checkGridData(data);
 *   RegularIndexer yIndexer = new AxisChecker<>(Entry::getY(), toleranceY).checkGridData(data);
 * }</pre>
 * @param <T> grid data type
 * @since 4.1
 */
public class AxisChecker<T> {

    /** Coordinate extractor. */
    private final ToDoubleFunction<T> extractor;

    /** Tolerance below which extracted coordinates are considered equal. */
    private final double tolerance;

    /** Simple constructor.
     * @param extractor coordinate extractor
     * @param tolerance tolerance below which extracted coordinates are considered equal
     */
    public AxisChecker(final ToDoubleFunction<T> extractor, final double tolerance) {
        this.extractor = extractor;
        this.tolerance = tolerance;
    }

    /** Check grid data for regularity (i.e., regular sampling, no missing points,
     * same number of points at all coordinates…).
     * <p>
     * The check for one axis performs two passes over the full grid data.
     * </p>
     * @param gridData grid data to check
     * @return regular sampling for the axis configured at construction
     */
    public RegularIndexer checkGridData(final Iterable<T> gridData) {

        // first pass, extract statistics
        double currentMin = Double.POSITIVE_INFINITY;
        double currentMax = Double.NEGATIVE_INFINITY;
        int    minCount   = 0;
        int    size       = 0;
        for (final T t : gridData) {
            final double coordinate = extractor.applyAsDouble(t);
            if (Double.isNaN(coordinate)) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.NAN_ELEMENT_AT_INDEX, size);
            }
            ++size;
            if (coordinate < currentMin - tolerance) {
                // this is a new minimum, we reset counter
                currentMin = coordinate;
                minCount   = 1;
            } else if (coordinate - currentMin < tolerance) {
                // we found again an already known minimum, update the count
                ++minCount;
            }
            currentMax = FastMath.max(currentMax, coordinate);
        }

        // store global indexing data
        final RegularIndexer regularIndexer = new RegularIndexer(currentMin, currentMax, size / minCount);

        // second pass, check grid data is regular
        final int[] count = new int[regularIndexer.getN()];
        for (final T t : gridData) {
            final double coordinate = extractor.applyAsDouble(t);
            final int    index      = regularIndexer.index(coordinate);
            final double expected   = regularIndexer.coordinate(index);
            if (FastMath.abs(coordinate - expected) > tolerance) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.MISALIGNED_GRID_POINT,
                                                       index, expected, coordinate - expected);
            }
            count[index]++;
        }
        for (int i = 0; i < regularIndexer.getN(); i++) {
            if (count[i] != minCount) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.IRREGULAR_GRID, minCount, i, count[i]);
            }
        }

        return regularIndexer;

    }

}
