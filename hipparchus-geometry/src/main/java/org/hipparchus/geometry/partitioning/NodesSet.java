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
import java.util.Iterator;
import java.util.List;

import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.Space;

/** Set of {@link BSPTree BSP tree} nodes.
 * @see BoundaryAttribute
 * @param <S> Type of the space.
 * @param <P> Type of the points in space.
 * @param <H> Type of the hyperplane.
 * @param <I> Type of the sub-hyperplane.
 */
public class NodesSet<S extends Space, P extends Point<S, P>, H extends Hyperplane<S, P, H, I>, I extends SubHyperplane<S, P, H, I>>
    implements Iterable<BSPTree<S, P, H, I>> {

    /** List of sub-hyperplanes. */
    private final List<BSPTree<S, P, H, I>> list;

    /** Simple constructor.
     */
    public NodesSet() {
        list = new ArrayList<>();
    }

    /** Add a node if not already known.
     * @param node node to add
     */
    public void add(final BSPTree<S, P, H, I> node) {

        for (final BSPTree<S, P, H, I> existing : list) {
            if (node == existing) {
                // the node is already known, don't add it
                return;
            }
        }

        // the node was not known, add it
        list.add(node);

    }

    /** Add nodes if they are not already known.
     * @param iterator nodes iterator
     */
    public void addAll(final Iterable<BSPTree<S, P, H, I>> iterator) {
        for (final BSPTree<S, P, H, I> node : iterator) {
            add(node);
        }
    }

    /** {@inheritDoc} */
    @Override
    public Iterator<BSPTree<S, P, H, I>> iterator() {
        return list.iterator();
    }

}
