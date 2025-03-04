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

import java.util.HashMap;
import java.util.Map;

import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.Space;

/** This class implements the dimension-independent parts of {@link SubHyperplane}.

 * <p>sub-hyperplanes are obtained when parts of an {@link
 * Hyperplane hyperplane} are chopped off by other hyperplanes that
 * intersect it. The remaining part is a convex region. Such objects
 * appear in {@link BSPTree BSP trees} as the intersection of a cut
 * hyperplane with the convex region which it splits, the chopping
 * hyperplanes are the cut hyperplanes closer to the tree root.</p>

 * @param <S> Type of the space.
 * @param <P> Type of the points in space.
 * @param <H> Type of the hyperplane.
 * @param <I> Type of the sub-hyperplane.
 * @param <T> Type of the sub-space.
 * @param <Q> Type of the points in sub-space.
 * @param <F> Type of the hyperplane.
 * @param <J> Type of the sub-hyperplane.

 */
public abstract class AbstractSubHyperplane<S extends Space,
                                            P extends Point<S, P>,
                                            H extends Hyperplane<S, P, H, I>,
                                            I extends SubHyperplane<S, P, H, I>,
                                            T extends Space, Q extends Point<T, Q>,
                                            F extends Hyperplane<T, Q, F, J>,
                                            J extends SubHyperplane<T, Q, F, J>>
    implements SubHyperplane<S, P, H, I> {

    /** Underlying hyperplane. */
    private final H hyperplane;

    /** Remaining region of the hyperplane. */
    private final Region<T, Q, F, J> remainingRegion;

    /** Build a sub-hyperplane from an hyperplane and a region.
     * @param hyperplane underlying hyperplane
     * @param remainingRegion remaining region of the hyperplane
     */
    protected AbstractSubHyperplane(final H hyperplane,
                                    final Region<T, Q, F, J> remainingRegion) {
        this.hyperplane      = hyperplane;
        this.remainingRegion = remainingRegion;
    }

    /** Build a sub-hyperplane from an hyperplane and a region.
     * @param hyper underlying hyperplane
     * @param remaining remaining region of the hyperplane
     * @return a new sub-hyperplane
     */
    protected abstract I buildNew(H hyper, Region<T, Q, F, J> remaining);

    /** {@inheritDoc} */
    @Override
    public I copySelf() {
        return buildNew(hyperplane.copySelf(), remainingRegion);
    }

    /** Get the underlying hyperplane.
     * @return underlying hyperplane
     */
    @Override
    public H getHyperplane() {
        return hyperplane;
    }

    /** Get the remaining region of the hyperplane.
     * <p>The returned region is expressed in the canonical hyperplane
     * frame and has the hyperplane dimension. For example a chopped
     * hyperplane in the 3D euclidean is a 2D plane and the
     * corresponding region is a convex 2D polygon.</p>
     * @return remaining region of the hyperplane
     */
    public Region<T, Q, F, J> getRemainingRegion() {
        return remainingRegion;
    }

    /** {@inheritDoc} */
    @Override
    public double getSize() {
        return remainingRegion.getSize();
    }

    /** {@inheritDoc} */
    @Override
    public I reunite(final I other) {
        @SuppressWarnings("unchecked")
        AbstractSubHyperplane<S, P, H, I, T, Q, F, J> o = (AbstractSubHyperplane<S, P, H, I, T, Q, F, J>) other;
        return buildNew(hyperplane,
                        new RegionFactory<T, Q, F, J>().union(remainingRegion, o.remainingRegion));
    }

    /** Apply a transform to the instance.
     * <p>The instance must be a (D-1)-dimension sub-hyperplane with
     * respect to the transform <em>not</em> a (D-2)-dimension
     * sub-hyperplane the transform knows how to transform by
     * itself. The transform will consist in transforming first the
     * hyperplane and then the all region using the various methods
     * provided by the transform.</p>
     * @param transform D-dimension transform to apply
     * @return the transformed instance
     */
    public I applyTransform(final Transform<S, P, H, I, T, Q, F, J> transform) {
        final H tHyperplane = transform.apply(hyperplane);

        // transform the tree, except for boundary attribute splitters
        final Map<BSPTree<T, Q, F, J>, BSPTree<T, Q, F, J>> map = new HashMap<>();
        final BSPTree<T, Q, F, J> tTree =
            recurseTransform(remainingRegion.getTree(false), tHyperplane, transform, map);

        // set up the boundary attributes splitters
        for (final Map.Entry<BSPTree<T, Q, F, J>, BSPTree<T, Q, F, J>> entry : map.entrySet()) {
            if (entry.getKey().getCut() != null) {
                @SuppressWarnings("unchecked")
                BoundaryAttribute<T, Q, F, J> original = (BoundaryAttribute<T, Q, F, J>) entry.getKey().getAttribute();
                if (original != null) {
                    @SuppressWarnings("unchecked")
                    BoundaryAttribute<T, Q, F, J> transformed = (BoundaryAttribute<T, Q, F, J>) entry.getValue().getAttribute();
                    for (final BSPTree<T, Q, F, J> splitter : original.getSplitters()) {
                        transformed.getSplitters().add(map.get(splitter));
                    }
                }
            }
        }

        return buildNew(tHyperplane, remainingRegion.buildNew(tTree));

    }

    /** Recursively transform a BSP-tree from a sub-hyperplane.
     * @param node current BSP tree node
     * @param transformed image of the instance hyperplane by the transform
     * @param transform transform to apply
     * @param map transformed nodes map
     * @return a new tree
     */
    private BSPTree<T, Q, F, J> recurseTransform(final BSPTree<T, Q, F, J> node,
                                                 final H transformed,
                                                 final Transform<S, P, H, I, T, Q, F, J> transform,
                                                 final Map<BSPTree<T, Q, F, J>, BSPTree<T, Q, F, J>> map) {

        final BSPTree<T, Q, F, J> transformedNode;
        if (node.getCut() == null) {
            transformedNode = new BSPTree<>(node.getAttribute());
        } else {

            @SuppressWarnings("unchecked")
            BoundaryAttribute<T, Q, F, J> attribute = (BoundaryAttribute<T, Q, F, J>) node.getAttribute();
            if (attribute != null) {
                final J tPO = (attribute.getPlusOutside() == null) ?
                    null : transform.apply(attribute.getPlusOutside(), hyperplane, transformed);
                final J tPI = (attribute.getPlusInside() == null) ?
                    null : transform.apply(attribute.getPlusInside(), hyperplane, transformed);
                // we start with an empty list of splitters, it will be filled in out of recursion
                attribute = new BoundaryAttribute<>(tPO, tPI, new NodesSet<>());
            }

            transformedNode = new BSPTree<>(transform.apply(node.getCut(), hyperplane, transformed),
                                            recurseTransform(node.getPlus(),  transformed, transform, map),
                                            recurseTransform(node.getMinus(), transformed, transform, map),
                                            attribute);
        }

        map.put(node, transformedNode);
        return transformedNode;

    }

    /** {@inheritDoc} */
    @Override
    public abstract SplitSubHyperplane<S, P, H, I> split(H hyper);

    /** {@inheritDoc} */
    @Override
    public boolean isEmpty() {
        return remainingRegion.isEmpty();
    }

}
