/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.geometry.partitioning;

import java.util.ArrayList;
import java.util.List;

import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.Space;
import org.hipparchus.geometry.partitioning.Region.Location;
import org.hipparchus.util.FastMath;

/** Local tree visitor to compute projection on boundary.
 * @param <S> Type of the space.
 * @param <P> Type of the points in space.
 * @param <H> Type of the hyperplane.
 * @param <I> Type of the sub-hyperplane.
 * @param <T> Type of the sub-space.
 * @param <Q> Type of the points in sub-space.
 * @param <F> Type of the hyperplane in the destination sub-space.
 * @param <J> Type of the sub-hyperplane in the destination sub-space.
 */
class BoundaryProjector<S extends Space,
                        P extends Point<S, P>,
                        H extends Hyperplane<S, P, H, I>,
                        I extends SubHyperplane<S, P, H, I>,
                        T extends Space,
                        Q extends Point<T, Q>,
                        F extends Hyperplane<T, Q, F, J>,
                        J extends SubHyperplane<T, Q, F, J>>
    implements BSPTreeVisitor<S, P, H, I> {

    /** Original point. */
    private final P original;

    /** Current best projected point. */
    private P projected;

    /** Leaf node closest to the test point. */
    private BSPTree<S, P, H, I> leaf;

    /** Current offset. */
    private double offset;

    /** Simple constructor.
     * @param original original point
     */
    BoundaryProjector(final P original) {
        this.original  = original;
        this.projected = null;
        this.leaf      = null;
        this.offset    = Double.POSITIVE_INFINITY;
    }

    /** {@inheritDoc} */
    @Override
    public Order visitOrder(final BSPTree<S, P, H, I> node) {
        // we want to visit the tree so that the first encountered
        // leaf is the one closest to the test point
        if (node.getCut().getHyperplane().getOffset(original) <= 0) {
            return Order.MINUS_SUB_PLUS;
        } else {
            return Order.PLUS_SUB_MINUS;
        }
    }

    /** {@inheritDoc} */
    @Override
    public void visitInternalNode(final BSPTree<S, P, H, I> node) {

        // project the point on the cut sub-hyperplane
        final H hyperplane = node.getCut().getHyperplane();
        final double signedOffset = hyperplane.getOffset(original);
        if (FastMath.abs(signedOffset) < offset) {

            // project point
            final P regular = hyperplane.project(original);

            // get boundary parts
            final List<Region<T, Q, F, J>> boundaryParts = boundaryRegions(node);

            // check if regular projection really belongs to the boundary
            boolean regularFound = false;
            for (final Region<T, Q, F, J> part : boundaryParts) {
                if (!regularFound && belongsToPart(regular, hyperplane, part)) {
                    // the projected point lies in the boundary
                    projected    = regular;
                    offset       = FastMath.abs(signedOffset);
                    regularFound = true;
                }
            }

            if (!regularFound) {
                // the regular projected point is not on boundary,
                // so we have to check further if a singular point
                // (i.e. a vertex in 2D case) is a possible projection
                for (final Region<T, Q, F, J> part : boundaryParts) {
                    final P spI = singularProjection(regular, hyperplane, part);
                    if (spI != null) {
                        final double distance = original.distance(spI);
                        if (distance < offset) {
                            projected = spI;
                            offset    = distance;
                        }
                    }
                }

            }

        }

    }

    /** {@inheritDoc} */
    @Override
    public void visitLeafNode(final BSPTree<S, P, H, I> node) {
        if (leaf == null) {
            // this is the first leaf we visit,
            // it is the closest one to the original point
            leaf = node;
        }
    }

    /** Get the projection.
     * @return projection
     */
    public BoundaryProjection<S, P> getProjection() {

        // fix offset sign
        offset = FastMath.copySign(offset, (Boolean) leaf.getAttribute() ? -1 : +1);

        return new BoundaryProjection<>(original, projected, offset);

    }

    /** Extract the regions of the boundary on an internal node.
     * @param node internal node
     * @return regions in the node sub-hyperplane
     */
    private List<Region<T, Q, F, J>> boundaryRegions(final BSPTree<S, P, H, I> node) {

        final List<Region<T, Q, F, J>> regions = new ArrayList<>(2);

        @SuppressWarnings("unchecked")
        final BoundaryAttribute<S, P, H, I> ba = (BoundaryAttribute<S, P, H, I>) node.getAttribute();
        addRegion(ba.getPlusInside(),  regions);
        addRegion(ba.getPlusOutside(), regions);

        return regions;

    }

    /** Add a boundary region to a list.
     * @param sub sub-hyperplane defining the region
     * @param list to fill up
     */
    private void addRegion(final SubHyperplane<S, P, H, I> sub, final List<Region<T, Q, F, J>> list) {
        if (sub != null) {
            final Region<T, Q, F, J> region = ((AbstractSubHyperplane<S, P, H, I, T, Q, F, J>) sub).getRemainingRegion();
            if (region != null) {
                list.add(region);
            }
        }
    }

    /** Check if a projected point lies on a boundary part.
     * @param point projected point to check
     * @param hyperplane hyperplane into which the point was projected
     * @param part boundary part
     * @return true if point lies on the boundary part
     */
    private boolean belongsToPart(final P point, final Hyperplane<S, P, H, I> hyperplane,
                                  final Region<T, Q, F, J> part) {

        // there is a non-null sub-space, we can dive into smaller dimensions
        @SuppressWarnings("unchecked")
        final Embedding<S, P, T, Q> embedding = (Embedding<S, P, T, Q>) hyperplane;
        return part.checkPoint(embedding.toSubSpace(point)) != Location.OUTSIDE;

    }

    /** Get the projection to the closest boundary singular point.
     * @param point projected point to check
     * @param hyperplane hyperplane into which the point was projected
     * @param part boundary part
     * @return projection to a singular point of boundary part (may be null)
     */
    private P singularProjection(final P point, final Hyperplane<S, P, H, I> hyperplane,
                                        final Region<T, Q, F, J> part) {

        // there is a non-null sub-space, we can dive into smaller dimensions
        @SuppressWarnings("unchecked")
        final Embedding<S, P, T, Q> embedding = (Embedding<S, P, T, Q>) hyperplane;
        final BoundaryProjection<T, Q> bp = part.projectToBoundary(embedding.toSubSpace(point));

        // back to initial dimension
        return (bp.getProjected() == null) ? null : embedding.toSpace(bp.getProjected());

    }

}
