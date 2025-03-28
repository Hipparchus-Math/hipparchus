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
package org.hipparchus.geometry.euclidean.threed;

import java.util.ArrayList;

import org.hipparchus.geometry.euclidean.twod.Euclidean2D;
import org.hipparchus.geometry.euclidean.twod.Line;
import org.hipparchus.geometry.euclidean.twod.PolygonsSet;
import org.hipparchus.geometry.euclidean.twod.SubLine;
import org.hipparchus.geometry.euclidean.twod.Vector2D;
import org.hipparchus.geometry.partitioning.BSPTree;
import org.hipparchus.geometry.partitioning.BSPTreeVisitor;
import org.hipparchus.geometry.partitioning.BoundaryAttribute;
import org.hipparchus.geometry.partitioning.RegionFactory;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/** Extractor for {@link PolygonsSet polyhedrons sets} outlines.
 * <p>This class extracts the 2D outlines from {{@link PolygonsSet
 * polyhedrons sets} in a specified projection plane.</p>
 */
public class OutlineExtractor {

    /** Abscissa axis of the projection plane. */
    private final Vector3D u;

    /** Ordinate axis of the projection plane. */
    private final Vector3D v;

    /** Normal of the projection plane (viewing direction). */
    private final Vector3D w;

    /** Build an extractor for a specific projection plane.
     * @param u abscissa axis of the projection point
     * @param v ordinate axis of the projection point
     */
    public OutlineExtractor(final Vector3D u, final Vector3D v) {
        this.u = u;
        this.v = v;
        w = Vector3D.crossProduct(u, v);
    }

    /** Extract the outline of a polyhedrons set.
     * @param polyhedronsSet polyhedrons set whose outline must be extracted
     * @return an outline, as an array of loops.
     */
    public Vector2D[][] getOutline(final PolyhedronsSet polyhedronsSet) {

        // project all boundary facets into one polygons set
        final BoundaryProjector projector = new BoundaryProjector(polyhedronsSet.getTolerance());
        polyhedronsSet.getTree(true).visit(projector);
        final PolygonsSet projected = projector.getProjected();

        // Remove the spurious intermediate vertices from the outline
        final Vector2D[][] outline = projected.getVertices();
        for (int i = 0; i < outline.length; ++i) {
            final Vector2D[] rawLoop = outline[i];
            int end = rawLoop.length;
            int j = 0;
            while (j < end) {
                if (pointIsBetween(rawLoop, end, j)) {
                    // the point should be removed
                    for (int k = j; k < (end - 1); ++k) {
                        rawLoop[k] = rawLoop[k + 1];
                    }
                    --end;
                } else {
                    // the point remains in the loop
                    ++j;
                }
            }
            if (end != rawLoop.length) {
                // resize the array
                outline[i] = new Vector2D[end];
                System.arraycopy(rawLoop, 0, outline[i], 0, end);
            }
        }

        return outline;

    }

    /** Check if a point is geometrically between its neighbor in an array.
     * <p>The neighbors are computed considering the array is a loop
     * (i.e. point at index (n-1) is before point at index 0)</p>
     * @param loop points array
     * @param n number of points to consider in the array
     * @param i index of the point to check (must be between 0 and n-1)
     * @return true if the point is exactly between its neighbors
     */
    private boolean pointIsBetween(final Vector2D[] loop, final int n, final int i) {
        final Vector2D previous = loop[(i + n - 1) % n];
        final Vector2D current  = loop[i];
        final Vector2D next     = loop[(i + 1) % n];
        final double dx1       = current.getX() - previous.getX();
        final double dy1       = current.getY() - previous.getY();
        final double dx2       = next.getX()    - current.getX();
        final double dy2       = next.getY()    - current.getY();
        final double cross     = dx1 * dy2 - dx2 * dy1;
        final double dot       = dx1 * dx2 + dy1 * dy2;
        final double d1d2      = FastMath.sqrt((dx1 * dx1 + dy1 * dy1) * (dx2 * dx2 + dy2 * dy2));
        return (FastMath.abs(cross) <= (1.0e-6 * d1d2)) && (dot >= 0.0);
    }

    /** Visitor projecting the boundary facets on a plane. */
    private class BoundaryProjector implements BSPTreeVisitor<Euclidean3D, Vector3D, Plane, SubPlane> {

        /** Projection of the polyhedrons set on the plane. */
        private PolygonsSet projected;

        /** Tolerance below which points are considered identical. */
        private final double tolerance;

        /** Simple constructor.
         * @param tolerance tolerance below which points are considered identical
         */
        BoundaryProjector(final double tolerance) {
            this.projected = new PolygonsSet(new BSPTree<>(Boolean.FALSE), tolerance);
            this.tolerance = tolerance;
        }

        /** {@inheritDoc} */
        @Override
        public Order visitOrder(final BSPTree<Euclidean3D, Vector3D, Plane, SubPlane> node) {
            return Order.MINUS_SUB_PLUS;
        }

        /** {@inheritDoc} */
        @Override
        public void visitInternalNode(final BSPTree<Euclidean3D, Vector3D, Plane, SubPlane> node) {
            @SuppressWarnings("unchecked")
            final BoundaryAttribute<Euclidean3D, Vector3D, Plane, SubPlane> attribute =
                (BoundaryAttribute<Euclidean3D, Vector3D, Plane, SubPlane>) node.getAttribute();
            if (attribute.getPlusOutside() != null) {
                addContribution(attribute.getPlusOutside());
            }
            if (attribute.getPlusInside() != null) {
                addContribution(attribute.getPlusInside());
            }
        }

        /** {@inheritDoc} */
        @Override
        public void visitLeafNode(final BSPTree<Euclidean3D, Vector3D, Plane, SubPlane> node) {
        }

        /** Add he contribution of a boundary facet.
         * @param facet boundary facet
         */
        private void addContribution(final SubPlane facet) {

            final double scal = facet.getHyperplane().getNormal().dotProduct(w);
            if (FastMath.abs(scal) > 1.0e-3) {
                Vector2D[][] vertices =
                    ((PolygonsSet) facet.getRemainingRegion()).getVertices();

                if (scal < 0) {
                    // the facet is seen from the back of the plane,
                    // we need to invert its boundary orientation
                    final Vector2D[][] newVertices = new Vector2D[vertices.length][];
                    for (int i = 0; i < vertices.length; ++i) {
                        final Vector2D[] loop = vertices[i];
                        final Vector2D[] newLoop = new Vector2D[loop.length];
                        if (loop[0] == null) {
                            newLoop[0] = null;
                            for (int j = 1; j < loop.length; ++j) {
                                newLoop[j] = loop[loop.length - j];
                            }
                        } else {
                            for (int j = 0; j < loop.length; ++j) {
                                newLoop[j] = loop[loop.length - (j + 1)];
                            }
                        }
                        newVertices[i] = newLoop;
                    }

                    // use the reverted vertices
                    vertices = newVertices;

                }

                // compute the projection of the facet in the outline plane
                final ArrayList<SubLine> edges = new ArrayList<>();
                for (Vector2D[] loop : vertices) {
                    final boolean closed = loop[0] != null;
                    int previous         = closed ? (loop.length - 1) : 1;
                    final Vector3D previous3D = facet.getHyperplane().toSpace(loop[previous]);
                    int current          = (previous + 1) % loop.length;
                    Vector2D pPoint      = new Vector2D(previous3D.dotProduct(u), previous3D.dotProduct(v));
                    while (current < loop.length) {

                        final Vector3D current3D = facet.getHyperplane().toSpace(loop[current]);
                        final Vector2D  cPoint    = new Vector2D(current3D.dotProduct(u),
                                                                 current3D.dotProduct(v));
                        final Line line = new Line(pPoint, cPoint, tolerance);
                        SubLine edge = line.wholeHyperplane();

                        if (closed || (previous != 1)) {
                            // the previous point is a real vertex
                            // it defines one bounding point of the edge
                            final double angle = line.getAngle() + MathUtils.SEMI_PI;
                            final Line l = new Line(pPoint, angle, tolerance);
                            edge = edge.split(l).getPlus();
                        }

                        if (closed || (current != (loop.length - 1))) {
                            // the current point is a real vertex
                            // it defines one bounding point of the edge
                            final double angle = line.getAngle() + MathUtils.SEMI_PI;
                            final Line l = new Line(cPoint, angle, tolerance);
                            edge = edge.split(l).getMinus();
                        }

                        edges.add(edge);

                        previous   = current++;
                        pPoint     = cPoint;

                    }
                }
                final PolygonsSet projectedFacet = new PolygonsSet(edges, tolerance);

                // add the contribution of the facet to the global outline
                final RegionFactory<Euclidean2D, Vector2D, Line, SubLine> factory = new RegionFactory<>();
                projected = (PolygonsSet) factory.union(projected, projectedFacet);

            }
        }

        /** Get the projection of the polyhedrons set on the plane.
         * @return projection of the polyhedrons set on the plane
         */
        public PolygonsSet getProjected() {
            return projected;
        }

    }

}
