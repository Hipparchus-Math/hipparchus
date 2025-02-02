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
 * @param <T> Type of the sub-space.
 * @param <Q> Type of the points in sub-space.

 */
public abstract class AbstractSubHyperplane<S extends Space, P extends Point<S>, T extends Space, Q extends Point<T>>
    implements SubHyperplane<S, P> {

    /** Underlying hyperplane. */
    private final Hyperplane<S, P> hyperplane;

    /** Remaining region of the hyperplane. */
    private final Region<T, Q> remainingRegion;

    /** Build a sub-hyperplane from an hyperplane and a region.
     * @param hyperplane underlying hyperplane
     * @param remainingRegion remaining region of the hyperplane
     */
    protected AbstractSubHyperplane(final Hyperplane<S, P> hyperplane,
                                    final Region<T, Q> remainingRegion) {
        this.hyperplane      = hyperplane;
        this.remainingRegion = remainingRegion;
    }

    /** Build a sub-hyperplane from an hyperplane and a region.
     * @param hyper underlying hyperplane
     * @param remaining remaining region of the hyperplane
     * @return a new sub-hyperplane
     */
    protected abstract AbstractSubHyperplane<S, P, T, Q> buildNew(Hyperplane<S, P> hyper,
                                                                  Region<T, Q> remaining);

    /** {@inheritDoc} */
    @Override
    public AbstractSubHyperplane<S, P, T, Q> copySelf() {
        return buildNew(hyperplane.copySelf(), remainingRegion);
    }

    /** Get the underlying hyperplane.
     * @return underlying hyperplane
     */
    @Override
    public Hyperplane<S, P> getHyperplane() {
        return hyperplane;
    }

    /** Get the remaining region of the hyperplane.
     * <p>The returned region is expressed in the canonical hyperplane
     * frame and has the hyperplane dimension. For example a chopped
     * hyperplane in the 3D euclidean is a 2D plane and the
     * corresponding region is a convex 2D polygon.</p>
     * @return remaining region of the hyperplane
     */
    public Region<T, Q> getRemainingRegion() {
        return remainingRegion;
    }

    /** {@inheritDoc} */
    @Override
    public double getSize() {
        return remainingRegion.getSize();
    }

    /** {@inheritDoc} */
    @Override
    public AbstractSubHyperplane<S, P, T, Q> reunite(final SubHyperplane<S, P> other) {
        AbstractSubHyperplane<S, P, T, Q> o = (AbstractSubHyperplane<S, P, T, Q>) other;
        return buildNew(hyperplane,
                        new RegionFactory<T, Q>().union(remainingRegion, o.remainingRegion));
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
    public AbstractSubHyperplane<S, P, T, Q> applyTransform(final Transform<S, P, T, Q> transform) {
        final Hyperplane<S, P> tHyperplane = transform.apply(hyperplane);

        // transform the tree, except for boundary attribute splitters
        final Map<BSPTree<T, Q>, BSPTree<T, Q>> map = new HashMap<>();
        final BSPTree<T, Q> tTree =
            recurseTransform(remainingRegion.getTree(false), tHyperplane, transform, map);

        // set up the boundary attributes splitters
        for (final Map.Entry<BSPTree<T, Q>, BSPTree<T, Q>> entry : map.entrySet()) {
            if (entry.getKey().getCut() != null) {
                @SuppressWarnings("unchecked")
                BoundaryAttribute<T, Q> original = (BoundaryAttribute<T, Q>) entry.getKey().getAttribute();
                if (original != null) {
                    @SuppressWarnings("unchecked")
                    BoundaryAttribute<T, Q> transformed = (BoundaryAttribute<T, Q>) entry.getValue().getAttribute();
                    for (final BSPTree<T, Q> splitter : original.getSplitters()) {
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
    private BSPTree<T, Q> recurseTransform(final BSPTree<T, Q> node,
                                           final Hyperplane<S, P> transformed,
                                           final Transform<S, P, T, Q> transform,
                                           final Map<BSPTree<T, Q>, BSPTree<T, Q>> map) {

        final BSPTree<T, Q> transformedNode;
        if (node.getCut() == null) {
            transformedNode = new BSPTree<>(node.getAttribute());
        } else {

            @SuppressWarnings("unchecked")
            BoundaryAttribute<T, Q> attribute = (BoundaryAttribute<T, Q>) node.getAttribute();
            if (attribute != null) {
                final SubHyperplane<T, Q> tPO = (attribute.getPlusOutside() == null) ?
                    null : transform.apply(attribute.getPlusOutside(), hyperplane, transformed);
                final SubHyperplane<T, Q> tPI = (attribute.getPlusInside() == null) ?
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
    public abstract SplitSubHyperplane<S, P> split(Hyperplane<S, P> hyper);

    /** {@inheritDoc} */
    @Override
    public boolean isEmpty() {
        return remainingRegion.isEmpty();
    }

}
