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
package org.hipparchus.geometry.partitioning;

import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.geometry.LocalizedGeometryFormats;
import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.Space;

/**
 * This class find points that are inside a convex cell.
 * @param <S> Type of the space.
 * @param <P> Type of the points in space.
 * @param <H> Type of the hyperplane.
 * @param <I> Type of the sub-hyperplane.
 * @since 4.0
 */
class InsideCellFinder<S extends Space,
                       P extends Point<S, P>,
                       H extends Hyperplane<S, P, H, I>,
                       I extends SubHyperplane<S, P, H, I>> {

    /** Maximum number of target offset stages for inside points search. */
    private static final int MAX_STAGES = 500;

    /** Factor at each stage for inside points search. */
    private static final int FACTOR = 20;

    /** Cell where point should be searched. */
    private final BSPTree<S, P, H, I> cell;

    /** Cell depth in the tree. */
    private final int depth;

    /** Simple constructor.
     * @param cell cell where point should be searched
     */
    InsideCellFinder(final BSPTree<S, P, H, I> cell) {
        this.cell  = cell;
        this.depth = cell.getDepth();
    }

    /** Find an inside point, starting from a specified point.
     * @param start starting point for the search
     * @return inside point
     */
    public P findInsidePoint(final P start) {

        // special handling for tree root
        if (cell.getParent() == null) {
            return start;
        }

        P point = start;

        // for each stage (i.e. multiplication factor wrt tolerance), we allow as many
        // attempts as there are parent nodes to move the test point. The goal is
        // to avoid infinite loops (same point/hyperplane used again and again)
        // when we exhaust the number of attempts, we change the stage factor (reducing it)
        for (int stage = 0; stage < MAX_STAGES; ++stage) {

            // the last stage factor (the smallest one), should be 2
            final int k           = MAX_STAGES - 1 - stage;
            final int stageFactor = 2 + FACTOR * k * k;

            for (int attempt = 0; attempt < depth; ++attempt) {
                // look for the two hyperplanes that give maximum offsets for current point
                // i.e. the worst cell boundaries for current point
                Selection<S, P, H, I> max  = null;
                Selection<S, P, H, I> next = null;
                for (BSPTree<S, P, H, I> node = cell; node.getParent() != null; node = node.getParent()) {
                    if (node.getParent().getCut() != null) {
                        final Selection<S, P, H, I> current = new Selection<>(node, point);
                        if (max == null || current.offset > max.offset) {
                            next = max;
                            max  = current;
                        } else if (next == null || current.offset > next.offset) {
                            next = current;
                        }
                    }
                }

                if (max == null || max.pointIsInside()) {
                    // either there are no parent (i.e. cell is the whole space)
                    // or the current point is known to be inside the cell
                    // in both cases, current point is what we were looking for
                    return point;
                }

                // move the point in an attempt to reduce maximum offset
                if (next == null || next.pointIsInside()) {
                    // we have only one bad hyperplane
                    // we target a point on the interior side of the cell wrt this bad hyperplane
                    point = max.movePoint(point, stageFactor);
                } else {
                    // we have two bad hyperplanes
                    // we target a point on the bisector
                    // this is an attempt to work around cases where two hyperplanes
                    // form a very thin wedge and their tolerances create an overlapping
                    // zone where we can't find any interior points;
                    // using a single hyperplane in this case (as per the alternative above)
                    // would result in ping-pong test points
                    final P point1 = max.movePoint(point, stageFactor);
                    final P point2 = next.movePoint(point, stageFactor);
                    point = point1.moveTowards(point2, 0.5);
                }

            }
        }

        // we were not able to find an inside point after many iterations
        throw new MathIllegalStateException(LocalizedGeometryFormats.CANNOT_FIND_INSIDE_POINT, MAX_STAGES);

    }

    /** Container for hyperplanes metadata using during search.
     * @param <S> Type of the space.
     * @param <P> Type of the points in space.
     * @param <H> Type of the hyperplane.
     * @param <I> Type of the sub-hyperplane.
     */
    private static class Selection<S extends Space,
                                   P extends Point<S, P>,
                                   H extends Hyperplane<S, P, H, I>,
                                   I extends SubHyperplane<S, P, H, I>> {

        /** Hyperplane used as a reference for moving search point. */
        private final H hyperplane;

        /** Sign correction to apply to offset. */
        private final double sign;

        /** Sign-corrected offset of current point. */
        private final double offset;

        /** Simple constructor.
         * @param node node containing the reference hyperplane (must have a non-null parent)
         * @param point current point
         */
        private Selection(final BSPTree<S, P, H, I> node, final P point) {

            hyperplane = node.getParent().getCut().getHyperplane();

            // we want negative offsets within the cell, positive offsets outside the cell
            sign   = node == node.getParent().getMinus() ? 1 : -1;
            offset = sign * hyperplane.getOffset(point);

        }

        /** Check if offset is negative enough to ensure the current point is inside the cell.
         * @return true if offset is negative enough to ensure the current point is inside the cell
         */
        private boolean pointIsInside() {
            return offset < -hyperplane.getTolerance();
        }

        /** Move point, attempting to get a negative offset.
         * @param point point to move
         * @param stageFactor multiplicative factor wrt. tolerance
         * @return moved point
         */
        private P movePoint(final P point, final double stageFactor) {
            final double tolerance = sign * hyperplane.getTolerance();
            return hyperplane.moveToOffset(point, -stageFactor * tolerance);
        }

    }

}
