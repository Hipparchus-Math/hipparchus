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

import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.Space;

/** Visitor building boundary shell tree.
 * <p>
 * The boundary shell is represented as {@link BoundaryAttribute boundary attributes}
 * at each internal node.
 * </p>
 * @param <S> Type of the space.
 * @param <P> Type of the points in space.
 * @param <H> Type of the hyperplane.
 * @param <I> Type of the sub-hyperplane.
 */
class BoundaryBuilder<S extends Space, P extends Point<S, P>, H extends Hyperplane<S, P, H, I>, I extends SubHyperplane<S, P, H, I>>
    implements BSPTreeVisitor<S, P, H, I> {

    /** {@inheritDoc} */
    @Override
    public Order visitOrder(BSPTree<S, P, H, I> node) {
        return Order.PLUS_MINUS_SUB;
    }

    /** {@inheritDoc} */
    @Override
    public void visitInternalNode(BSPTree<S, P, H, I> node) {

        I                    plusOutside = null;
        I                    plusInside  = null;
        NodesSet<S, P, H, I> splitters   = null;

        // characterize the cut sub-hyperplane,
        // first with respect to the plus sub-tree
        final Characterization<S, P, H, I> plusChar = new Characterization<>(node.getPlus(), node.getCut().copySelf());

        if (plusChar.touchOutside()) {
            // plusChar.outsideTouching() corresponds to a subset of the cut sub-hyperplane
            // known to have outside cells on its plus side, we want to check if parts
            // of this subset do have inside cells on their minus side
            final Characterization<S, P, H, I> minusChar = new Characterization<>(node.getMinus(), plusChar.outsideTouching());
            if (minusChar.touchInside()) {
                // this part belongs to the boundary,
                // it has the outside on its plus side and the inside on its minus side
                plusOutside = minusChar.insideTouching();
                splitters = new NodesSet<>();
                splitters.addAll(minusChar.getInsideSplitters());
                splitters.addAll(plusChar.getOutsideSplitters());
            }
        }

        if (plusChar.touchInside()) {
            // plusChar.insideTouching() corresponds to a subset of the cut sub-hyperplane
            // known to have inside cells on its plus side, we want to check if parts
            // of this subset do have outside cells on their minus side
            final Characterization<S, P, H, I> minusChar = new Characterization<>(node.getMinus(), plusChar.insideTouching());
            if (minusChar.touchOutside()) {
                // this part belongs to the boundary,
                // it has the inside on its plus side and the outside on its minus side
                plusInside = minusChar.outsideTouching();
                if (splitters == null) {
                    splitters = new NodesSet<>();
                }
                splitters.addAll(minusChar.getOutsideSplitters());
                splitters.addAll(plusChar.getInsideSplitters());
            }
        }

        if (splitters != null) {
            // the parent nodes are natural splitters for boundary sub-hyperplanes
            for (BSPTree<S, P, H, I> up = node.getParent(); up != null; up = up.getParent()) {
                splitters.add(up);
            }
        }

        // set the boundary attribute at non-leaf nodes
        node.setAttribute(new BoundaryAttribute<>(plusOutside, plusInside, splitters));

    }

    /** {@inheritDoc} */
    @Override
    public void visitLeafNode(BSPTree<S, P, H, I> node) {
    }

}
