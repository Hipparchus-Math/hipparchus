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

/** Visitor computing the boundary size.
 * @param <S> Type of the space.
 * @param <P> Type of the points in space.
 * @param <H> Type of the hyperplane.
 * @param <I> Type of the sub-hyperplane.
 */
class BoundarySizeVisitor<S extends Space, P extends Point<S, P>, H extends Hyperplane<S, P, H, I>, I extends SubHyperplane<S, P, H, I>>
    implements BSPTreeVisitor<S, P, H, I> {

    /** Size of the boundary. */
    private double boundarySize;

    /** Simple constructor.
     */
    BoundarySizeVisitor() {
        boundarySize = 0;
    }

    /** {@inheritDoc}*/
    @Override
    public Order visitOrder(final BSPTree<S, P, H, I> node) {
        return Order.MINUS_SUB_PLUS;
    }

    /** {@inheritDoc}*/
    @Override
    public void visitInternalNode(final BSPTree<S, P, H, I> node) {
        @SuppressWarnings("unchecked")
        final BoundaryAttribute<S, P, H, I> attribute =
            (BoundaryAttribute<S, P, H, I>) node.getAttribute();
        if (attribute.getPlusOutside() != null) {
            boundarySize += attribute.getPlusOutside().getSize();
        }
        if (attribute.getPlusInside() != null) {
            boundarySize += attribute.getPlusInside().getSize();
        }
    }

    /** {@inheritDoc}*/
    @Override
    public void visitLeafNode(final BSPTree<S, P, H, I> node) {
    }

    /** Get the size of the boundary.
     * @return size of the boundary
     */
    public double getSize() {
        return boundarySize;
    }

}
