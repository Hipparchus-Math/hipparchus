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
package org.hipparchus.geometry.spherical.twod;

import java.util.ArrayList;
import java.util.List;

import org.hipparchus.geometry.partitioning.AbstractSubHyperplane;
import org.hipparchus.geometry.partitioning.BSPTree;
import org.hipparchus.geometry.partitioning.BSPTreeVisitor;
import org.hipparchus.geometry.partitioning.BoundaryAttribute;
import org.hipparchus.geometry.partitioning.SubHyperplane;
import org.hipparchus.geometry.spherical.oned.Arc;
import org.hipparchus.geometry.spherical.oned.ArcsSet;
import org.hipparchus.geometry.spherical.oned.S1Point;
import org.hipparchus.geometry.spherical.oned.Sphere1D;
import org.hipparchus.util.FastMath;

/** Visitor building edges.
 * @since 1.4
 */
class EdgesWithNodeInfoBuilder implements BSPTreeVisitor<Sphere2D> {

    /** Tolerance for close nodes connection. */
    private final double tolerance;

    /** Built segments. */
    private final List<EdgeWithNodeInfo> edges;

    /** Simple constructor.
     * @param tolerance below which points are consider to be identical
     */
    EdgesWithNodeInfoBuilder(final double tolerance) {
        this.tolerance = tolerance;
        this.edges     = new ArrayList<>();
    }

    /** {@inheritDoc} */
    @Override
    public Order visitOrder(final BSPTree<Sphere2D> node) {
        return Order.MINUS_SUB_PLUS;
    }

    /** {@inheritDoc} */
    @Override
    public void visitInternalNode(final BSPTree<Sphere2D> node) {
        @SuppressWarnings("unchecked")
        final BoundaryAttribute<Sphere2D> attribute = (BoundaryAttribute<Sphere2D>) node.getAttribute();
        final Iterable<BSPTree<Sphere2D>> splitters = attribute.getSplitters();
        if (attribute.getPlusOutside() != null) {
            addContribution(attribute.getPlusOutside(), node, splitters, false);
        }
        if (attribute.getPlusInside() != null) {
            addContribution(attribute.getPlusInside(), node, splitters, true);
        }
    }

    /** {@inheritDoc} */
    @Override
    public void visitLeafNode(final BSPTree<Sphere2D> node) {
    }

    /** Add the contribution of a boundary edge.
     * @param sub boundary facet
     * @param node node to which the edge belongs
     * @param splitters splitters for the boundary facet
     * @param reversed if true, the facet has the inside on its plus side
     */
    private void addContribution(final SubHyperplane<Sphere2D> sub, final BSPTree<Sphere2D> node,
                                 final Iterable<BSPTree<Sphere2D>> splitters,
                                 final boolean reversed) {
        final AbstractSubHyperplane<Sphere2D, Sphere1D> absSub =
                        (AbstractSubHyperplane<Sphere2D, Sphere1D>) sub;
        final Circle circle  = (Circle) sub.getHyperplane();
        final List<Arc> arcs = ((ArcsSet) absSub.getRemainingRegion()).asList();
        for (final Arc a : arcs) {

            // find the 2D points
            final Vertex startS = new Vertex(circle.toSpace(new S1Point(a.getInf())));
            final Vertex endS   = new Vertex(circle.toSpace(new S1Point(a.getSup())));

            // recover the connectivity information
            final BSPTree<Sphere2D> startN = selectClosest(startS.getLocation(), splitters);
            final BSPTree<Sphere2D> endN   = selectClosest(endS.getLocation(), splitters);

            if (reversed) {
                edges.add(new EdgeWithNodeInfo(endS, startS, a.getSize(), circle.getReverse(),
                                              node, endN, startN));
            } else {
                edges.add(new EdgeWithNodeInfo(startS, endS, a.getSize(), circle,
                                              node, startN, endN));
            }

        }
    }

    /** Select the node whose cut sub-hyperplane is closest to specified point.
     * @param point reference point
     * @param candidates candidate nodes
     * @return node closest to point, or null if no node is closer than tolerance
     */
    private BSPTree<Sphere2D> selectClosest(final S2Point point, final Iterable<BSPTree<Sphere2D>> candidates) {

        if (point == null) {
            return null;
        }

        BSPTree<Sphere2D> selected = null;

        double min = Double.POSITIVE_INFINITY;
        for (final BSPTree<Sphere2D> node : candidates) {
            final double distance = FastMath.abs(node.getCut().getHyperplane().getOffset(point));
            if (distance < min) {
                selected = node;
                min      = distance;
            }
        }

        return min <= tolerance ? selected : null;

    }

    /** Get the edges.
     * @return built edges
     */
    public List<EdgeWithNodeInfo> getEdges() {
        return edges;
    }

}
