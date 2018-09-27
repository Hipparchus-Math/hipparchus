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
import java.util.concurrent.atomic.AtomicInteger;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;

/**
 * Helper for finding interpolation nodes along one axis of grid data.
 * <p>
 * This class is intended to be used for interpolating inside grids.
 * It works on any sorted data without duplication and size at least
 * {@code n} where {@code n} is the number of points required for
 * interpolation (i.e. 2 for linear interpolation, 3 for quadratic...)
 * <p>
 * </p>
 * The method uses linear interpolation to select the nodes indices.
 * It should be O(1) for sufficiently regular data, therefore much faster
 * than bisection. It also features caching, which improves speed when
 * interpolating several points in raw in the close locations, i.e. when
 * successive calls have a high probability to return the same interpolation
 * nodes. This occurs for example when scanning with small steps a loose
 * grid. The method also works on non-regular grids, but may be slower in
 * this case.
 * </p>
 * <p>
 * This class is thread-safe.
 * </p>
 * @since 1.4
 */
public class GridAxis implements Serializable {

    /** Serializable UID. */
    private static final long serialVersionUID = 20180926L;

    /** All the coordinates of the interpolation points, sorted in increasing order. */
    private final double[] grid;

    /** Number of points required for interpolation. */
    private int n;

    /** Cached value of last x index. */
    private final AtomicInteger cache;

    /** Simple constructor.
     * @param grid coordinates of the interpolation points, sorted in increasing order
     * @param n number of points required for interpolation, i.e. 2 for linear, 3
     * for quadratic...
     * @exception MathIllegalArgumentException if grid size is smaller than {@code n}
     * or if the grid is not sorted in strict increasing order
     */
    public GridAxis(final double[] grid, final int n)
        throws MathIllegalArgumentException {

        // safety checks
        if (grid.length < n) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.INSUFFICIENT_DIMENSION,
                                                   grid.length, n);
        }
        MathArrays.checkOrder(grid);

        this.grid  = grid.clone();
        this.n     = n;
        this.cache = new AtomicInteger(0);

    }

    /** Get the number of points of the grid.
     * @return number of points of the grid
     */
    public int size() {
        return grid.length;
    }

    /** Get the number of points required for interpolation.
     * @return number of points required for interpolation
     */
    public int getN() {
        return n;
    }

    /** Get the interpolation node at specified index.
     * @param index node index
     * @return coordinate of the node at specified index
     */
    public double node(final int index) {
        return grid[index];
    }

    /** Get the index of the first interpolation node for some coordinate along the grid.
     * <p>
     * The index return is the one for the lowest interpolation node suitable for
     * {@code t}. This means that if {@code i} is returned the nodes to use for
     * interpolation at coordinate {@code t} are at indices {@code i}, {@code i+1},
     * ..., {@code i+n-1}, where {@code n} is the number of points required for
     * interpolation passed at construction.
     * </p>
     * <p>
     * The index is selected in order to have the subset of nodes from {@code i} to
     * {@code i+n-1} as balanced as possible around {@code t}:
     * </p>
     * <ul>
     *   <li>
     *     if {@code t} is inside the grid and sufficiently far from the endpoints
     *     <ul>
     *       <li>
     *         if {@code n} is even, the returned nodes will be perfectly balanced:
     *         there will be {@code n/2} nodes smaller than {@code t} and {@code n/2}
     *         nodes larger than {@code t}
     *       </li>
     *       <li>
     *         if {@code n} is odd, the returned nodes will be slightly unbalanced by
     *         one point: there will be {@code (n+1)/2} nodes smaller than {@code t}
     *         and {@code (n-1)/2} nodes larger than {@code t}
     *       </li>
     *     </ul>
     *   </li>
     *   <li>
     *     if {@code t} is inside the grid and close to endpoints, the returned nodes
     *     will be unbalanced: there will be less nodes on the endpoints side and
     *     more nodes on the interior side
     *   </li>
     *   <li>
     *     if {@code t} is outside of the grid, the returned nodes will completely
     *     off balance: all nodes will be on the same side with respect to {@code t}
     *   </li>
     * </ul>
     * <p>
     * It is <em>not</em> an error to call this method with {@code t} outside of the grid,
     * it simply implies that the interpolation will become an extrapolation and accuracy
     * will decrease as {@code t} goes farther from the grid points. This is intended so
     * interpolation does not fail near the end of the grid.
     * </p>
     * @param t coordinate of the point to interpolate
     * @return index {@code i} such {@link #node(int) node(i)}, {@link #node(int) node(i+1)},
     * ... {@link #node(int) node(i+n-1)} can be used for interpolating a value at
     * coordinate {@code t}
     * @since 1.4
     */
    public int interpolationIndex(final double t) {

        final int middleOffset = (n - 1) / 2;
        int iInf = middleOffset;
        int iSup = grid.length - (n - 1) + middleOffset;

        // first try to simply reuse the cached index,
        // for faster return in a common case
        final int    cached = cache.get();
        final int    middle = cached + middleOffset;
        final double aMid0  = grid[middle];
        final double aMid1  = grid[middle + 1];
        if (t < aMid0) {
            if (middle == iInf) {
                // we are in the unbalanced low area
                return cached;
            }
        } else if (t < aMid1) {
            // we are in the balanced middle area
            return cached;
        } else {
            if (middle == iSup - 1) {
                // we are in the unbalanced high area
                return cached;
            }
        }

        // we need to find a new index
        double aInf = grid[iInf];
        double aSup = grid[iSup];
        while (iSup - iInf > 1) {
            final int iInterp = (int) ((iInf * (aSup - t) + iSup * (t - aInf)) / (aSup - aInf));
            final int iMed    = FastMath.max(iInf + 1, FastMath.min(iInterp, iSup - 1));
            if (t < grid[iMed]) {
                // keeps looking in the lower part of the grid
                iSup = iMed;
                aSup = grid[iSup];
            } else {
                // keeps looking in the upper part of the grid
                iInf = iMed;
                aInf = grid[iInf];
            }
        }

       final int newCached = iInf - middleOffset;
       cache.compareAndSet(cached, newCached);
       return newCached;

    }

}
