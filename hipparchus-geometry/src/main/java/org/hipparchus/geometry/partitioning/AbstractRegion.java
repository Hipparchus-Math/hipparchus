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
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.Space;

/** Abstract class for all regions, independently of geometry type or dimension.

 * @param <S> Type of the space.
 * @param <P> Type of the points in space.
 * @param <H> Type of the hyperplane.
 * @param <I> Type of the sub-hyperplane.
 * @param <T> Type of the sub-space.
 * @param <Q> Type of the points in sub-space.
 * @param <F> Type of the hyperplane.
 * @param <J> Type of the sub-hyperplane.

 */
public abstract class AbstractRegion<S extends Space,
                                     P extends Point<S, P>,
                                     H extends Hyperplane<S, P, H, I>,
                                     I extends SubHyperplane<S, P, H, I>,
                                     T extends Space, Q extends Point<T, Q>,
                                     F extends Hyperplane<T, Q, F, J>,
                                     J extends SubHyperplane<T, Q, F, J>>
    implements Region<S, P, H, I> {

    /** Inside/Outside BSP tree. */
    private final BSPTree<S, P, H, I> tree;

    /** Tolerance below which points are considered to belong to hyperplanes. */
    private final double tolerance;

    /** Size of the instance. */
    private double size;

    /** Barycenter. */
    private P barycenter;

    /** Build a region representing the whole space.
     * @param tolerance tolerance below which points are considered identical.
     */
    protected AbstractRegion(final double tolerance) {
        this.tree      = new BSPTree<>(Boolean.TRUE);
        this.tolerance = tolerance;
    }

    /** Build a region from an inside/outside BSP tree.
     * <p>The leaf nodes of the BSP tree <em>must</em> have a
     * {@code Boolean} attribute representing the inside status of
     * the corresponding cell (true for inside cells, false for outside
     * cells). In order to avoid building too many small objects, it is
     * recommended to use the predefined constants
     * {@code Boolean.TRUE} and {@code Boolean.FALSE}. The
     * tree also <em>must</em> have either null internal nodes or
     * internal nodes representing the boundary as specified in the
     * {@link #getTree getTree} method).</p>
     * @param tree inside/outside BSP tree representing the region
     * @param tolerance tolerance below which points are considered identical.
     */
    protected AbstractRegion(final BSPTree<S, P, H, I> tree, final double tolerance) {
        this.tree      = tree;
        this.tolerance = tolerance;
    }

    /** Build a Region from a Boundary REPresentation (B-rep).
     * <p>The boundary is provided as a collection of {@link
     * SubHyperplane sub-hyperplanes}. Each sub-hyperplane has the
     * interior part of the region on its minus side and the exterior on
     * its plus side.</p>
     * <p>The boundary elements can be in any order, and can form
     * several non-connected sets (like for example polygons with holes
     * or a set of disjoints polyhedrons considered as a whole). In
     * fact, the elements do not even need to be connected together
     * (their topological connections are not used here). However, if the
     * boundary does not really separate an inside open from an outside
     * open (open having here its topological meaning), then subsequent
     * calls to the {@link #checkPoint(Point) checkPoint} method will not be
     * meaningful anymore.</p>
     * <p>If the boundary is empty, the region will represent the whole
     * space.</p>
     * @param boundary collection of boundary elements, as a
     * collection of {@link SubHyperplane SubHyperplane} objects
     * @param tolerance tolerance below which points are considered identical.
     */
    protected AbstractRegion(final Collection<I> boundary, final double tolerance) {

        this.tolerance = tolerance;

        if (boundary.isEmpty()) {

            // the tree represents the whole space
            tree = new BSPTree<>(Boolean.TRUE);

        } else {

            // sort the boundary elements in decreasing size order
            // (we don't want equal size elements to be removed, so
            // we use a trick to fool the TreeSet)
            final TreeSet<I> ordered = new TreeSet<>(new Comparator<I>() {
                /** {@inheritDoc} */
                @Override
                public int compare(final I o1, final I o2) {
                    final double size1 = o1.getSize();
                    final double size2 = o2.getSize();
                    return (size2 < size1) ? -1 : ((o1 == o2) ? 0 : +1);
                }
            });
            ordered.addAll(boundary);

            // build the tree top-down
            tree = new BSPTree<>();
            insertCuts(tree, ordered);

            // set up the inside/outside flags
            tree.visit(new BSPTreeVisitor<S, P, H, I>() {

                /** {@inheritDoc} */
                @Override
                public Order visitOrder(final BSPTree<S, P, H, I> node) {
                    return Order.PLUS_SUB_MINUS;
                }

                /** {@inheritDoc} */
                @Override
                public void visitInternalNode(final BSPTree<S, P, H, I> node) {
                }

                /** {@inheritDoc} */
                @Override
                public void visitLeafNode(final BSPTree<S, P, H, I> node) {
                    if (node.getParent() == null || node == node.getParent().getMinus()) {
                        node.setAttribute(Boolean.TRUE);
                    } else {
                        node.setAttribute(Boolean.FALSE);
                    }
                }
            });

        }

    }

    /** Build a convex region from an array of bounding hyperplanes.
     * @param hyperplanes array of bounding hyperplanes (if null, an
     * empty region will be built)
     * @param tolerance tolerance below which points are considered identical.
     */
    public AbstractRegion(final H[] hyperplanes, final double tolerance) {
        this.tolerance = tolerance;
        if ((hyperplanes == null) || (hyperplanes.length == 0)) {
            tree = new BSPTree<>(Boolean.FALSE);
        } else {

            // use the first hyperplane to build the right class
            tree = hyperplanes[0].wholeSpace().getTree(false);

            // chop off parts of the space
            BSPTree<S, P, H, I> node = tree;
            node.setAttribute(Boolean.TRUE);
            for (final H hyperplane : hyperplanes) {
                if (node.insertCut(hyperplane)) {
                    node.setAttribute(null);
                    node.getPlus().setAttribute(Boolean.FALSE);
                    node = node.getMinus();
                    node.setAttribute(Boolean.TRUE);
                }
            }

        }

    }

    /** {@inheritDoc} */
    @Override
    public abstract AbstractRegion<S, P, H, I, T, Q, F, J> buildNew(BSPTree<S, P, H, I> newTree);

    /** Get the tolerance below which points are considered to belong to hyperplanes.
     * @return tolerance below which points are considered to belong to hyperplanes
     */
    public double getTolerance() {
        return tolerance;
    }

    /** Recursively build a tree by inserting cut sub-hyperplanes.
     * @param node current tree node (it is a leaf node at the beginning
     * of the call)
     * @param boundary collection of edges belonging to the cell defined
     * by the node
     */
    private void insertCuts(final BSPTree<S, P, H, I> node, final Collection<I> boundary) {

        final Iterator<I> iterator = boundary.iterator();

        // build the current level
        H inserted = null;
        while ((inserted == null) && iterator.hasNext()) {
            inserted = iterator.next().getHyperplane();
            if (!node.insertCut(inserted.copySelf())) {
                inserted = null;
            }
        }

        if (!iterator.hasNext()) {
            return;
        }

        // distribute the remaining edges in the two sub-trees
        final ArrayList<I> plusList  = new ArrayList<>();
        final ArrayList<I> minusList = new ArrayList<>();
        while (iterator.hasNext()) {
            final I other = iterator.next();
            final SubHyperplane.SplitSubHyperplane<S, P, H, I> split = other.split(inserted);
            switch (split.getSide()) {
            case PLUS:
                plusList.add(other);
                break;
            case MINUS:
                minusList.add(other);
                break;
            case BOTH:
                plusList.add(split.getPlus());
                minusList.add(split.getMinus());
                break;
            default:
                // ignore the sub-hyperplanes belonging to the cut hyperplane
            }
        }

        // recurse through lower levels
        insertCuts(node.getPlus(),  plusList);
        insertCuts(node.getMinus(), minusList);

    }

    /** {@inheritDoc} */
    @Override
    public AbstractRegion<S, P, H, I, T, Q, F, J> copySelf() {
        return buildNew(tree.copySelf());
    }

    /** {@inheritDoc} */
    @Override
    public boolean isEmpty() {
        return isEmpty(tree);
    }

    /** {@inheritDoc} */
    @Override
    public boolean isEmpty(final BSPTree<S, P, H, I> node) {

        // we use a recursive function rather than the BSPTreeVisitor
        // interface because we can stop visiting the tree as soon as we
        // have found an inside cell

        if (node.getCut() == null) {
            // if we find an inside node, the region is not empty
            return !((Boolean) node.getAttribute());
        }

        // check both sides of the sub-tree
        return isEmpty(node.getMinus()) && isEmpty(node.getPlus());

    }

    /** {@inheritDoc} */
    @Override
    public boolean isFull() {
        return isFull(tree);
    }

    /** {@inheritDoc} */
    @Override
    public boolean isFull(final BSPTree<S, P, H, I> node) {

        // we use a recursive function rather than the BSPTreeVisitor
        // interface because we can stop visiting the tree as soon as we
        // have found an outside cell

        if (node.getCut() == null) {
            // if we find an outside node, the region does not cover full space
            return (Boolean) node.getAttribute();
        }

        // check both sides of the sub-tree
        return isFull(node.getMinus()) && isFull(node.getPlus());

    }

    /** {@inheritDoc} */
    @Override
    public boolean contains(final Region<S, P, H, I> region) {
        return new RegionFactory<S, P, H, I>().difference(region, this).isEmpty();
    }

    /** {@inheritDoc}
     */
    @Override
    public BoundaryProjection<S, P> projectToBoundary(final P point) {
        final BoundaryProjector<S, P, H, I, T, Q, F, J> projector = new BoundaryProjector<>(point);
        getTree(true).visit(projector);
        return projector.getProjection();
    }

    /** {@inheritDoc} */
    @Override
    public Location checkPoint(final P point) {
        return checkPoint(tree, point);
    }

    /** Check a point with respect to the region starting at a given node.
     * @param node root node of the region
     * @param point point to check
     * @return a code representing the point status: either {@link
     * Region.Location#INSIDE INSIDE}, {@link Region.Location#OUTSIDE
     * OUTSIDE} or {@link Region.Location#BOUNDARY BOUNDARY}
     */
    protected Location checkPoint(final BSPTree<S, P, H, I> node, final P point) {
        final BSPTree<S, P, H, I> cell = node.getCell(point, tolerance);
        if (cell.getCut() == null) {
            // the point is in the interior of a cell, just check the attribute
            return ((Boolean) cell.getAttribute()) ? Location.INSIDE : Location.OUTSIDE;
        }

        // the point is on a cut-sub-hyperplane, is it on a boundary ?
        final Location minusCode = checkPoint(cell.getMinus(), point);
        final Location plusCode  = checkPoint(cell.getPlus(),  point);
        return (minusCode == plusCode) ? minusCode : Location.BOUNDARY;

    }

    /** {@inheritDoc} */
    @Override
    public BSPTree<S, P, H, I> getTree(final boolean includeBoundaryAttributes) {
        if (includeBoundaryAttributes && (tree.getCut() != null) && (tree.getAttribute() == null)) {
            // compute the boundary attributes
            tree.visit(new BoundaryBuilder<>());
        }
        return tree;
    }

    /** {@inheritDoc} */
    @Override
    public double getBoundarySize() {
        final BoundarySizeVisitor<S, P, H, I> visitor = new BoundarySizeVisitor<>();
        getTree(true).visit(visitor);
        return visitor.getSize();
    }

    /** {@inheritDoc} */
    @Override
    public double getSize() {
        if (barycenter == null) {
            computeGeometricalProperties();
        }
        return size;
    }

    /** Set the size of the instance.
     * @param size size of the instance
     */
    protected void setSize(final double size) {
        this.size = size;
    }

    /** {@inheritDoc} */
    @Override
    public P getBarycenter() {
        if (barycenter == null) {
            computeGeometricalProperties();
        }
        return barycenter;
    }

    /** Set the barycenter of the instance.
     * @param barycenter barycenter of the instance
     */
    protected void setBarycenter(final P barycenter) {
        this.barycenter = barycenter;
    }

    /** Compute some geometrical properties.
     * <p>The properties to compute are the barycenter and the size.</p>
     */
    protected abstract void computeGeometricalProperties();

    /** {@inheritDoc} */
    @Override
    public I intersection(final I sub) {
        return recurseIntersection(tree, sub);
    }

    /** Recursively compute the parts of a sub-hyperplane that are
     * contained in the region.
     * @param node current BSP tree node
     * @param sub sub-hyperplane traversing the region
     * @return filtered sub-hyperplane
     */
    private I recurseIntersection(final BSPTree<S, P, H, I> node, final I sub) {

        if (node.getCut() == null) {
            return (Boolean) node.getAttribute() ? sub.copySelf() : null;
        }

        final H hyperplane = node.getCut().getHyperplane();
        final SubHyperplane.SplitSubHyperplane<S, P, H, I> split = sub.split(hyperplane);
        if (split.getPlus() != null) {
            if (split.getMinus() != null) {
                // both sides
                final I plus  = recurseIntersection(node.getPlus(),  split.getPlus());
                final I minus = recurseIntersection(node.getMinus(), split.getMinus());
                if (plus == null) {
                    return minus;
                } else if (minus == null) {
                    return plus;
                } else {
                    return plus.reunite(minus);
                }
            } else {
                // only on plus side
                return recurseIntersection(node.getPlus(), sub);
            }
        } else if (split.getMinus() != null) {
            // only on minus side
            return recurseIntersection(node.getMinus(), sub);
        } else {
            // on hyperplane
            return recurseIntersection(node.getPlus(),
                                       recurseIntersection(node.getMinus(), sub));
        }

    }

    /** Transform a region.
     * <p>Applying a transform to a region consist in applying the
     * transform to all the hyperplanes of the underlying BSP tree and
     * of the boundary (and also to the sub-hyperplanes embedded in
     * these hyperplanes) and to the barycenter. The instance is not
     * modified, a new instance is built.</p>
     * @param transform transform to apply
     * @return a new region, resulting from the application of the
     * transform to the instance
     */
    public AbstractRegion<S, P, H, I, T, Q, F, J> applyTransform(final Transform<S, P, H, I, T, Q, F, J> transform) {

        // transform the tree, except for boundary attribute splitters
        final Map<BSPTree<S, P, H, I>, BSPTree<S, P, H, I>> map = new HashMap<>();
        final BSPTree<S, P, H, I> transformedTree = recurseTransform(getTree(false), transform, map);

        // set up the boundary attributes splitters
        for (final Map.Entry<BSPTree<S, P, H, I>, BSPTree<S, P, H, I>> entry : map.entrySet()) {
            if (entry.getKey().getCut() != null) {
                @SuppressWarnings("unchecked")
                BoundaryAttribute<S, P, H, I> original = (BoundaryAttribute<S, P, H, I>) entry.getKey().getAttribute();
                if (original != null) {
                    @SuppressWarnings("unchecked")
                    BoundaryAttribute<S, P, H, I> transformed = (BoundaryAttribute<S, P, H, I>) entry.getValue().getAttribute();
                    for (final BSPTree<S, P, H, I> splitter : original.getSplitters()) {
                        transformed.getSplitters().add(map.get(splitter));
                    }
                }
            }
        }

        return buildNew(transformedTree);

    }

    /** Recursively transform an inside/outside BSP-tree.
     * @param node current BSP tree node
     * @param transform transform to apply
     * @param map transformed nodes map
     * @return a new tree
     */
    @SuppressWarnings("unchecked")
    private BSPTree<S, P, H, I> recurseTransform(final BSPTree<S, P, H, I> node, final Transform<S, P, H, I, T, Q, F, J> transform,
                                        final Map<BSPTree<S, P, H, I>, BSPTree<S, P, H, I>> map) {

        final BSPTree<S, P, H, I> transformedNode;
        if (node.getCut() == null) {
            transformedNode = new BSPTree<>(node.getAttribute());
        } else {

            final I  sub = node.getCut();
            final I tSub = ((AbstractSubHyperplane<S, P, H, I, T, Q, F, J>) sub).applyTransform(transform);
            BoundaryAttribute<S, P, H, I> attribute = (BoundaryAttribute<S, P, H, I>) node.getAttribute();
            if (attribute != null) {
                final I tPO = (attribute.getPlusOutside() == null) ?
                    null : ((AbstractSubHyperplane<S, P, H, I, T, Q, F, J>) attribute.getPlusOutside()).applyTransform(transform);
                final I tPI = (attribute.getPlusInside()  == null) ?
                    null  : ((AbstractSubHyperplane<S, P, H, I, T, Q, F, J>) attribute.getPlusInside()).applyTransform(transform);
                // we start with an empty list of splitters, it will be filled in out of recursion
                attribute = new BoundaryAttribute<>(tPO, tPI, new NodesSet<>());
            }

            transformedNode = new BSPTree<>(tSub,
                                            recurseTransform(node.getPlus(),  transform, map),
                                            recurseTransform(node.getMinus(), transform, map),
                                            attribute);
        }

        map.put(node, transformedNode);
        return transformedNode;

    }

}
