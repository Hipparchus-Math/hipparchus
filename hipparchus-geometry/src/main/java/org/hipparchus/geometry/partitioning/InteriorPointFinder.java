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

import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.Space;

/**
 * Finder for interior points.
 * @param <S> Type of the space.
 * @param <P> Type of the points in space.
 * @param <H> Type of the hyperplane.
 * @param <I> Type of the sub-hyperplane.
 * @since 4.0
 */
public class InteriorPointFinder<S extends Space,
                                 P extends Point<S, P>,
                                 H extends Hyperplane<S, P, H, I>,
                                 I extends SubHyperplane<S, P, H, I>>
    implements BSPTreeVisitor<S, P, H, I> {

    /** Default point to use for whole space. */
    private final P defaultPoint;

    /** Selected point. */
    private BSPTree.InteriorPoint<S, P> selected;

    /**
     * Simple constructor.
     * @param defaultPoint default point to use for whole space
     */
    public InteriorPointFinder(final P defaultPoint) {
        this.defaultPoint = defaultPoint;
        this.selected     = null;
    }

    /** {@inheritDoc} */
    @Override
    public Order visitOrder(final BSPTree<S, P, H, I> node) {
        return Order.MINUS_PLUS_SUB;
    }

    /** {@inheritDoc} */
    @Override
    public void visitInternalNode(final BSPTree<S, P, H, I> node) {
        // nothing to do
    }

    /** {@inheritDoc} */
    @Override
    public void visitLeafNode(final BSPTree<S, P, H, I> node) {
        if ((Boolean) node.getAttribute()) {
            // this is an inside cell, look for the barycenter of edges/facets interior points
            final BSPTree.InteriorPoint<S, P> interior = node.getInteriorPoint(defaultPoint);
            if (selected == null || interior.getDistance() > selected.getDistance()) {
                // this new point is farther away to edges/facets than the selected one, change selection
                selected = interior;
            }
        }
    }

    /**
     * Get the point found.
     * @return found point (null if tree was empty)
     */
    public BSPTree.InteriorPoint<S, P> getPoint() {
        return selected;
    }

}
