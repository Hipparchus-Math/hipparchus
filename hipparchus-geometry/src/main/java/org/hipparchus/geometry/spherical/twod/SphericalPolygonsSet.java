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
package org.hipparchus.geometry.spherical.twod;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.function.IntPredicate;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.geometry.LocalizedGeometryFormats;
import org.hipparchus.geometry.enclosing.EnclosingBall;
import org.hipparchus.geometry.enclosing.WelzlEncloser;
import org.hipparchus.geometry.euclidean.threed.Euclidean3D;
import org.hipparchus.geometry.euclidean.threed.Rotation;
import org.hipparchus.geometry.euclidean.threed.RotationConvention;
import org.hipparchus.geometry.euclidean.threed.SphereGenerator;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.geometry.partitioning.AbstractRegion;
import org.hipparchus.geometry.partitioning.BSPTree;
import org.hipparchus.geometry.partitioning.BoundaryProjection;
import org.hipparchus.geometry.partitioning.RegionFactory;
import org.hipparchus.geometry.partitioning.SubHyperplane;
import org.hipparchus.geometry.spherical.oned.Arc;
import org.hipparchus.geometry.spherical.oned.Sphere1D;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Precision;

/** This class represents a region on the 2-sphere: a set of spherical polygons.
 */
public class SphericalPolygonsSet extends AbstractRegion<Sphere2D, Sphere1D> {

    /** Boundary defined as an array of closed loops start vertices. */
    private List<Vertex> loops;

    /** Build a polygons set representing the whole real 2-sphere.
     * @param tolerance below which points are consider to be identical
     * @exception MathIllegalArgumentException if tolerance is smaller than {@link Sphere1D#SMALLEST_TOLERANCE}
     */
    public SphericalPolygonsSet(final double tolerance) throws MathIllegalArgumentException {
        super(tolerance);
        Sphere2D.checkTolerance(tolerance);
    }

    /** Build a polygons set representing a hemisphere.
     * @param pole pole of the hemisphere (the pole is in the inside half)
     * @param tolerance below which points are consider to be identical
     * @exception MathIllegalArgumentException if tolerance is smaller than {@link Sphere1D#SMALLEST_TOLERANCE}
     */
    public SphericalPolygonsSet(final Vector3D pole, final double tolerance)
        throws MathIllegalArgumentException {
        super(new BSPTree<Sphere2D>(new Circle(pole, tolerance).wholeHyperplane(),
                                    new BSPTree<Sphere2D>(Boolean.FALSE),
                                    new BSPTree<Sphere2D>(Boolean.TRUE),
                                    null),
              tolerance);
        Sphere2D.checkTolerance(tolerance);
    }

    /** Build a polygons set representing a regular polygon.
     * @param center center of the polygon (the center is in the inside half)
     * @param meridian point defining the reference meridian for first polygon vertex
     * @param outsideRadius distance of the vertices to the center
     * @param n number of sides of the polygon
     * @param tolerance below which points are consider to be identical
     * @exception MathIllegalArgumentException if tolerance is smaller than {@link Sphere1D#SMALLEST_TOLERANCE}
     */
    public SphericalPolygonsSet(final Vector3D center, final Vector3D meridian,
                                final double outsideRadius, final int n,
                                final double tolerance)
        throws MathIllegalArgumentException {
        this(tolerance, createRegularPolygonVertices(center, meridian, outsideRadius, n));
    }

    /** Build a polygons set from a BSP tree.
     * <p>The leaf nodes of the BSP tree <em>must</em> have a
     * {@code Boolean} attribute representing the inside status of
     * the corresponding cell (true for inside cells, false for outside
     * cells). In order to avoid building too many small objects, it is
     * recommended to use the predefined constants
     * {@code Boolean.TRUE} and {@code Boolean.FALSE}</p>
     * @param tree inside/outside BSP tree representing the region
     * @param tolerance below which points are consider to be identical
     * @exception MathIllegalArgumentException if tolerance is smaller than {@link Sphere1D#SMALLEST_TOLERANCE}
     */
    public SphericalPolygonsSet(final BSPTree<Sphere2D> tree, final double tolerance)
        throws MathIllegalArgumentException {
        super(tree, tolerance);
        Sphere2D.checkTolerance(tolerance);
    }

    /** Build a polygons set from a Boundary REPresentation (B-rep).
     * <p>The boundary is provided as a collection of {@link
     * SubHyperplane sub-hyperplanes}. Each sub-hyperplane has the
     * interior part of the region on its minus side and the exterior on
     * its plus side.</p>
     * <p>The boundary elements can be in any order, and can form
     * several non-connected sets (like for example polygons with holes
     * or a set of disjoint polygons considered as a whole). In
     * fact, the elements do not even need to be connected together
     * (their topological connections are not used here). However, if the
     * boundary does not really separate an inside open from an outside
     * open (open having here its topological meaning), then subsequent
     * calls to the {@link
     * org.hipparchus.geometry.partitioning.Region#checkPoint(org.hipparchus.geometry.Point)
     * checkPoint} method will not be meaningful anymore.</p>
     * <p>If the boundary is empty, the region will represent the whole
     * space.</p>
     * @param boundary collection of boundary elements, as a
     * collection of {@link SubHyperplane SubHyperplane} objects
     * @param tolerance below which points are consider to be identical
     * @exception MathIllegalArgumentException if tolerance is smaller than {@link Sphere1D#SMALLEST_TOLERANCE}
     */
    public SphericalPolygonsSet(final Collection<SubHyperplane<Sphere2D>> boundary, final double tolerance)
        throws MathIllegalArgumentException {
        super(boundary, tolerance);
        Sphere2D.checkTolerance(tolerance);
    }

    /** Build a polygon from a simple list of vertices.
     * <p>The boundary is provided as a list of points considering to
     * represent the vertices of a simple loop. The interior part of the
     * region is on the left side of this path and the exterior is on its
     * right side.</p>
     * <p>This constructor does not handle polygons with a boundary
     * forming several disconnected paths (such as polygons with holes).</p>
     * <p>For cases where this simple constructor applies, it is expected to
     * be numerically more robust than the {@link #SphericalPolygonsSet(Collection,
     * double) general constructor} using {@link SubHyperplane subhyperplanes}.</p>
     * <p>If the list is empty, the region will represent the whole
     * space.</p>
     * <p>This constructor assumes that edges between {@code vertices}, including the edge
     * between the last and the first vertex, are shorter than pi. If edges longer than pi
     * are used it may produce unintuitive results, such as reversing the direction of the
     * edge. This implies using a {@code vertices} array of length 1 or 2 in this
     * constructor produces an ill-defined region. Use one of the other constructors or
     * {@link RegionFactory} instead.</p>
     * <p>The list of {@code vertices} is reduced by selecting a sub-set of vertices
     * before creating the boundary set. Every point in {@code vertices} will be on the
     * {@link #checkPoint(org.hipparchus.geometry.Point) boundary} of the constructed polygon set, but not
     * necessarily the center-line of the boundary.</p>
     * <p>
     * Polygons with thin pikes or dents are inherently difficult to handle because
     * they involve circles with almost opposite directions at some vertices. Polygons
     * whose vertices come from some physical measurement with noise are also
     * difficult because an edge that should be straight may be broken in lots of
     * different pieces with almost equal directions. In both cases, computing the
     * circles intersections is not numerically robust due to the almost 0 or almost
     * &pi; angle. Such cases need to carefully adjust the {@code hyperplaneThickness}
     * parameter. A too small value would often lead to completely wrong polygons
     * with large area wrongly identified as inside or outside. Large values are
     * often much safer. As a rule of thumb, a value slightly below the size of the
     * most accurate detail needed is a good value for the {@code hyperplaneThickness}
     * parameter.
     * </p>
     * @param hyperplaneThickness tolerance below which points are considered to
     * belong to the hyperplane (which is therefore more a slab). Should be greater than
     * {@code FastMath.ulp(4 * FastMath.PI)} for meaningful results.
     * @param vertices vertices of the simple loop boundary
     * @exception MathIllegalArgumentException if tolerance is smaller than {@link Sphere1D#SMALLEST_TOLERANCE}
     * @exception org.hipparchus.exception.MathRuntimeException if {@code vertices}
     * contains only a single vertex or repeated vertices.
     */
    public SphericalPolygonsSet(final double hyperplaneThickness, final S2Point ... vertices)
        throws MathIllegalArgumentException {
        super(verticesToTree(hyperplaneThickness, vertices), hyperplaneThickness);
        Sphere2D.checkTolerance(hyperplaneThickness);
    }

    /** Build the vertices representing a regular polygon.
     * @param center center of the polygon (the center is in the inside half)
     * @param meridian point defining the reference meridian for first polygon vertex
     * @param outsideRadius distance of the vertices to the center
     * @param n number of sides of the polygon
     * @return vertices array
     */
    private static S2Point[] createRegularPolygonVertices(final Vector3D center, final Vector3D meridian,
                                                          final double outsideRadius, final int n) {
        final S2Point[] array = new S2Point[n];
        final Rotation r0 = new Rotation(Vector3D.crossProduct(center, meridian),
                                         outsideRadius, RotationConvention.VECTOR_OPERATOR);
        array[0] = new S2Point(r0.applyTo(center));

        final Rotation r = new Rotation(center, MathUtils.TWO_PI / n, RotationConvention.VECTOR_OPERATOR);
        for (int i = 1; i < n; ++i) {
            array[i] = new S2Point(r.applyTo(array[i - 1].getVector()));
        }

        return array;
    }

    /** Build the BSP tree of a polygons set from a simple list of vertices.
     * <p>The boundary is provided as a list of points considering to
     * represent the vertices of a simple loop. The interior part of the
     * region is on the left side of this path and the exterior is on its
     * right side.</p>
     * <p>This constructor does not handle polygons with a boundary
     * forming several disconnected paths (such as polygons with holes).</p>
     * <p>This constructor handles only polygons with edges strictly shorter
     * than \( \pi \). If longer edges are needed, they need to be broken up
     * in smaller sub-edges so this constraint holds.</p>
     * <p>For cases where this simple constructor applies, it is expected to
     * be numerically more robust than the {@link #PolygonsSet(Collection) general
     * constructor} using {@link SubHyperplane subhyperplanes}.</p>
     * @param hyperplaneThickness tolerance below which points are consider to
     * belong to the hyperplane (which is therefore more a slab)
     * @param vertices vertices of the simple loop boundary
     * @return the BSP tree of the input vertices
     */
    private static BSPTree<Sphere2D> verticesToTree(final double hyperplaneThickness,
                                                    S2Point ... vertices) {
        // thin vertices to those that define distinct circles
        vertices = reduce(hyperplaneThickness, vertices).toArray(new S2Point[0]);
        final int n = vertices.length;
        if (n == 0) {
            // the tree represents the whole space
            return new BSPTree<Sphere2D>(Boolean.TRUE);
        }

        // build the vertices
        final Vertex[] vArray = new Vertex[n];
        for (int i = 0; i < n; ++i) {
            vArray[i] = new Vertex(vertices[i]);
        }

        // build the edges
        final List<Edge> edges = new ArrayList<>(n);
        Vertex end = vArray[n - 1];
        for (int i = 0; i < n; ++i) {

            // get the endpoints of the edge
            final Vertex start = end;
            end = vArray[i];

            // get the circle supporting the edge
            final Circle circle = new Circle(start.getLocation(), end.getLocation(), hyperplaneThickness);

            // create the edge and store it
            edges.add(new Edge(start, end,
                               Vector3D.angle(start.getLocation().getVector(),
                                              end.getLocation().getVector()),
                               circle));

        }

        // build the tree top-down
        final BSPTree<Sphere2D> tree = new BSPTree<>();
        insertEdges(hyperplaneThickness, tree, edges);

        return tree;

    }

    /**
     * Compute a subset of vertices that define the boundary to within the given
     * tolerance. This method partitions {@code vertices} into segments that all lie same
     * arc to within {@code hyperplaneThickness}, and then returns the end points of the
     * arcs. Combined arcs are limited to length of pi. If the input vertices has arcs
     * longer than pi these will be preserved in the returned data.
     *
     * @param hyperplaneThickness of each circle in radians.
     * @param vertices            to decimate.
     * @return a subset of {@code vertices}.
     */
    private static List<S2Point> reduce(final double hyperplaneThickness,
                                        final S2Point[] vertices) {
        final int n = vertices.length;
        if (n <= 3) {
            // can't reduce to fewer than three points
            return Arrays.asList(vertices.clone());
        }
        final List<S2Point> points = new ArrayList<>();
        /* Use a simple greedy search to add points to a circle s.t. all intermediate
         * points are within the thickness. Running time is O(n lg n) worst case.
         * Since the first vertex may be the middle of a straight edge, look backward
         * and forward to establish the first edge.
         * Uses the property that any two points define a circle, so don't check
         * circles that just span two points.
         */
        // first look backward
        final IntPredicate onCircleBackward = j -> {
            final int i = n - 2 - j;
            // circle spanning considered points
            final Circle circle = new Circle(vertices[0], vertices[i], hyperplaneThickness);
            final Arc arc = circle.getArc(vertices[0], vertices[i]);
            if (arc.getSize() >= FastMath.PI) {
                return false;
            }
            for (int k = i + 1; k < n; k++) {
                final S2Point vertex = vertices[k];
                if (FastMath.abs(circle.getOffset(vertex)) > hyperplaneThickness ||
                        arc.getOffset(circle.toSubSpace(vertex)) > 0) {
                    // point is not within the thickness or arc, start new edge
                    return false;
                }
            }
            return true;
        };
        // last index in vertices of last entry added to points
        int last = n - 2 - searchHelper(onCircleBackward, 0, n - 2);
        if (last > 1) {
            points.add(vertices[last]);
        } else {
            // all points lie on one semi-circle, distance from 0 to 1 is > pi
            // ill-defined case, just return three points from the list
            return Arrays.asList(Arrays.copyOfRange(vertices, 0, 3));
        }
        final int first = last;
        // then build edges forward
        for (int j = 1; ; j += 2) {
            final int lastFinal = last;
            final IntPredicate onCircle = i -> {
                // circle spanning considered points
                final Circle circle = new Circle(vertices[lastFinal], vertices[i], hyperplaneThickness);
                final Arc arc = circle.getArc(vertices[lastFinal], vertices[i]);
                if (arc.getSize() >= FastMath.PI) {
                    return false;
                }
                final int end = lastFinal < i ? i : i + n;
                for (int k = lastFinal + 1; k < end; k++) {
                    final S2Point vertex = vertices[k % n];
                    if (FastMath.abs(circle.getOffset(vertex)) > hyperplaneThickness ||
                            arc.getOffset(circle.toSubSpace(vertex)) > 0) {
                        // point is not within the thickness or arc, start new edge
                        return false;
                    }
                }
                return true;
            };
            j = searchHelper(onCircle, j, first + 1);
            if (j >= first) {
                break;
            }
            last = j;
            points.add(vertices[last]);
        }
        // put first point last
        final S2Point swap = points.remove(0);
        points.add(swap);
        return points;
    }

    /**
     * Search {@code items} for the first item where {@code predicate} is false between
     * {@code a} and {@code b}. Assumes that predicate switches from true to false at
     * exactly one location in [a, b]. Similar to {@link Arrays#binarySearch(int[], int,
     * int, int)} except that 1. it operates on indices, not elements, 2. there is not a
     * shortcut for equality, and 3. it is optimized for cases where the return value is
     * close to a.
     *
     * <p> This method achieves O(lg n) performance in the worst case, where n = b - a.
     * Performance improves to O(lg(i-a)) when i is close to a, where i is the return
     * value.
     *
     * @param predicate to apply.
     * @param a         start, inclusive.
     * @param b         end, exclusive.
     * @return a if a==b, a-1 if predicate.test(a) == false, b - 1 if predicate.test(b-1),
     * otherwise i s.t. predicate.test(i) == true && predicate.test(i + 1) == false.
     * @throws MathIllegalArgumentException if a > b.
     */
    private static int searchHelper(final IntPredicate predicate,
                                    final int a,
                                    final int b) {
        if (a > b) {
            throw new MathIllegalArgumentException(
                    LocalizedCoreFormats.LOWER_ENDPOINT_ABOVE_UPPER_ENDPOINT, a, b);
        }
        // Argument checks and special cases
        if (a == b) {
            return a;
        }
        if (!predicate.test(a)) {
            return a - 1;
        }

        // start with exponential search
        int start = a;
        int end = b;
        for (int i = 2; a + i < b; i *= 2) {
            if (predicate.test(a + i)) {
                // update lower bound of search
                start = a + i;
            } else {
                // found upper bound of search
                end = a + i;
                break;
            }
        }

        // next binary search
        // copied from Arrays.binarySearch() and modified to work on indices alone
        int low = start;
        int high = end - 1;
        while (low <= high) {
            final int mid = (low + high) >>> 1;
            if (predicate.test(mid)) {
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        // low is now insertion point, according to Arrays.binarySearch()
        return low - 1;
    }

    /** Recursively build a tree by inserting cut sub-hyperplanes.
     * @param hyperplaneThickness tolerance below which points are considered to
     * belong to the hyperplane (which is therefore more a slab)
     * @param node current tree node (it is a leaf node at the beginning
     * of the call)
     * @param edges list of edges to insert in the cell defined by this node
     * (excluding edges not belonging to the cell defined by this node)
     */
    private static void insertEdges(final double hyperplaneThickness,
                                    final BSPTree<Sphere2D> node,
                                    final List<Edge> edges) {

        // find an edge with an hyperplane that can be inserted in the node
        int index = 0;
        Edge inserted = null;
        while (inserted == null && index < edges.size()) {
            inserted = edges.get(index++);
            if (!node.insertCut(inserted.getCircle())) {
                inserted = null;
            }
        }

        if (inserted == null) {
            // no suitable edge was found, the node remains a leaf node
            // we need to set its inside/outside boolean indicator
            final BSPTree<Sphere2D> parent = node.getParent();
            if (parent == null || node == parent.getMinus()) {
                node.setAttribute(Boolean.TRUE);
            } else {
                node.setAttribute(Boolean.FALSE);
            }
            return;
        }

        // we have split the node by inserting an edge as a cut sub-hyperplane
        // distribute the remaining edges in the two sub-trees
        final List<Edge> outsideList = new ArrayList<>();
        final List<Edge> insideList  = new ArrayList<>();
        for (final Edge edge : edges) {
            if (edge != inserted) {
                edge.split(inserted.getCircle(), outsideList, insideList);
            }
        }

        // recurse through lower levels
        if (!outsideList.isEmpty()) {
            insertEdges(hyperplaneThickness, node.getPlus(), outsideList);
        } else {
            node.getPlus().setAttribute(Boolean.FALSE);
        }
        if (!insideList.isEmpty()) {
            insertEdges(hyperplaneThickness, node.getMinus(),  insideList);
        } else {
            node.getMinus().setAttribute(Boolean.TRUE);
        }

    }

    /** {@inheritDoc} */
    @Override
    public SphericalPolygonsSet buildNew(final BSPTree<Sphere2D> tree) {
        return new SphericalPolygonsSet(tree, getTolerance());
    }

    /** {@inheritDoc}
     * @exception MathIllegalStateException if the tolerance setting does not allow to build
     * a clean non-ambiguous boundary
     */
    @Override
    protected void computeGeometricalProperties() throws MathIllegalStateException {

        final BSPTree<Sphere2D> tree = getTree(true);

        if (tree.getCut() == null) {

            // the instance has a single cell without any boundaries

            if (tree.getCut() == null && (Boolean) tree.getAttribute()) {
                // the instance covers the whole space
                setSize(4 * FastMath.PI);
                setBarycenter(new S2Point(0, 0));
            } else {
                setSize(0);
                setBarycenter(S2Point.NaN);
            }

        } else {

            // the instance has a boundary
            final PropertiesComputer pc = new PropertiesComputer(getTolerance());
            tree.visit(pc);
            setSize(pc.getArea());
            setBarycenter(pc.getBarycenter());

        }

    }

    /** Get the boundary loops of the polygon.
     * <p>The polygon boundary can be represented as a list of closed loops,
     * each loop being given by exactly one of its vertices. From each loop
     * start vertex, one can follow the loop by finding the outgoing edge,
     * then the end vertex, then the next outgoing edge ... until the start
     * vertex of the loop (exactly the same instance) is found again once
     * the full loop has been visited.</p>
     * <p>If the polygon has no boundary at all, a zero length loop
     * array will be returned.</p>
     * <p>If the polygon is a simple one-piece polygon, then the returned
     * array will contain a single vertex.
     * </p>
     * <p>All edges in the various loops have the inside of the region on
     * their left side (i.e. toward their pole) and the outside on their
     * right side (i.e. away from their pole) when moving in the underlying
     * circle direction. This means that the closed loops obey the direct
     * trigonometric orientation.</p>
     * @return boundary of the polygon, organized as an unmodifiable list of loops start vertices.
     * @exception MathIllegalStateException if the tolerance setting does not allow to build
     * a clean non-ambiguous boundary
     * @see Vertex
     * @see Edge
     */
    public List<Vertex> getBoundaryLoops() throws MathIllegalStateException {

        if (loops == null) {
            if (getTree(false).getCut() == null) {
                loops = Collections.emptyList();
            } else {

                // sort the arcs according to their start point
                final EdgesWithNodeInfoBuilder visitor = new EdgesWithNodeInfoBuilder(getTolerance());
                getTree(true).visit(visitor);
                final List<EdgeWithNodeInfo> edges = visitor.getEdges();

                // connect all edges, using topological criteria first
                // and using Euclidean distance only as a last resort
                int pending = edges.size();
                pending -= naturalFollowerConnections(edges);
                if (pending > 0) {
                    pending -= splitEdgeConnections(edges);
                }
                if (pending > 0) {
                    closeVerticesConnections(edges);
                }

                // extract the edges loops
                loops = new ArrayList<>();
                for (EdgeWithNodeInfo s = getUnprocessed(edges); s != null; s = getUnprocessed(edges)) {
                    loops.add(s.getStart());
                    followLoop(s);
                }

            }
        }

        return Collections.unmodifiableList(loops);

    }

    /** Connect the edges using only natural follower information.
     * @param edges edges complete edges list
     * @return number of connections performed
     */
    private int naturalFollowerConnections(final List<EdgeWithNodeInfo> edges) {
        int connected = 0;
        for (final EdgeWithNodeInfo edge : edges) {
            if (edge.getEnd().getOutgoing() == null) {
                for (final EdgeWithNodeInfo candidateNext : edges) {
                    if (EdgeWithNodeInfo.areNaturalFollowers(edge, candidateNext)) {
                        // connect the two edges
                        edge.setNextEdge(candidateNext);
                        ++connected;
                        break;
                    }
                }
            }
        }
        return connected;
    }

    /** Connect the edges resulting from a circle splitting a circular edge.
     * @param edges edges complete edges list
     * @return number of connections performed
     */
    private int splitEdgeConnections(final List<EdgeWithNodeInfo> edges) {
        int connected = 0;
        for (final EdgeWithNodeInfo edge : edges) {
            if (edge.getEnd().getOutgoing() == null) {
                for (final EdgeWithNodeInfo candidateNext : edges) {
                    if (EdgeWithNodeInfo.resultFromASplit(edge, candidateNext)) {
                        // connect the two edges
                        edge.setNextEdge(candidateNext);
                        ++connected;
                        break;
                    }
                }
            }
        }
        return connected;
    }

    /** Connect the edges using spherical distance.
     * <p>
     * This connection heuristic should be used last, as it relies
     * only on a fuzzy distance criterion.
     * </p>
     * @param edges edges complete edges list
     * @return number of connections performed
     */
    private int closeVerticesConnections(final List<EdgeWithNodeInfo> edges) {
        int connected = 0;
        for (final EdgeWithNodeInfo edge : edges) {
            if (edge.getEnd().getOutgoing() == null && edge.getEnd() != null) {
                final Vector3D end = edge.getEnd().getLocation().getVector();
                EdgeWithNodeInfo selectedNext = null;
                double min = Double.POSITIVE_INFINITY;
                for (final EdgeWithNodeInfo candidateNext : edges) {
                    if (candidateNext.getStart().getIncoming() == null) {
                        final double distance = Vector3D.distance(end, candidateNext.getStart().getLocation().getVector());
                        if (distance < min) {
                            selectedNext = candidateNext;
                            min          = distance;
                        }
                    }
                }
                if (min <= getTolerance()) {
                    // connect the two edges
                    edge.setNextEdge(selectedNext);
                    ++connected;
                }
            }
        }
        return connected;
    }

    /** Get first unprocessed edge from a list.
     * @param edges edges list
     * @return first edge that has not been processed yet
     * or null if all edges have been processed
     */
    private EdgeWithNodeInfo getUnprocessed(final List<EdgeWithNodeInfo> edges) {
        for (final EdgeWithNodeInfo edge : edges) {
            if (!edge.isProcessed()) {
                return edge;
            }
        }
        return null;
    }

    /** Build the loop containing a edge.
     * <p>
     * All edges put in the loop will be marked as processed.
     * </p>
     * @param defining edge used to define the loop
     */
    private void followLoop(final EdgeWithNodeInfo defining) {

        defining.setProcessed(true);

        // process edges in connection order
        EdgeWithNodeInfo previous = defining;
        EdgeWithNodeInfo next     = (EdgeWithNodeInfo) defining.getEnd().getOutgoing();
        while (next != defining) {
            if (next == null) {
                // this should not happen
                throw new MathIllegalStateException(LocalizedGeometryFormats.OUTLINE_BOUNDARY_LOOP_OPEN);
            }
            next.setProcessed(true);

            // filter out spurious vertices
            if (Vector3D.angle(previous.getCircle().getPole(), next.getCircle().getPole()) <= Precision.EPSILON) {
                // the vertex between the two edges is a spurious one
                // replace the two edges by a single one
                previous.setNextEdge(next.getEnd().getOutgoing());
                previous.setLength(previous.getLength() + next.getLength());
            }

            previous = next;
            next     = (EdgeWithNodeInfo) next.getEnd().getOutgoing();

        }

    }

    /** Get a spherical cap enclosing the polygon.
     * <p>
     * This method is intended as a first test to quickly identify points
     * that are guaranteed to be outside of the region, hence performing a full
     * {@link #checkPoint(org.hipparchus.geometry.Vector) checkPoint}
     * only if the point status remains undecided after the quick check. It is
     * is therefore mostly useful to speed up computation for small polygons with
     * complex shapes (say a country boundary on Earth), as the spherical cap will
     * be small and hence will reliably identify a large part of the sphere as outside,
     * whereas the full check can be more computing intensive. A typical use case is
     * therefore:
     * </p>
     * <pre>
     *   // compute region, plus an enclosing spherical cap
     *   SphericalPolygonsSet complexShape = ...;
     *   EnclosingBall&lt;Sphere2D, S2Point&gt; cap = complexShape.getEnclosingCap();
     *
     *   // check lots of points
     *   for (Vector3D p : points) {
     *
     *     final Location l;
     *     if (cap.contains(p)) {
     *       // we cannot be sure where the point is
     *       // we need to perform the full computation
     *       l = complexShape.checkPoint(v);
     *     } else {
     *       // no need to do further computation,
     *       // we already know the point is outside
     *       l = Location.OUTSIDE;
     *     }
     *
     *     // use l ...
     *
     *   }
     * </pre>
     * <p>
     * In the special cases of empty or whole sphere polygons, special
     * spherical caps are returned, with angular radius set to negative
     * or positive infinity so the {@link
     * EnclosingBall#contains(org.hipparchus.geometry.Point) ball.contains(point)}
     * method return always false or true.
     * </p>
     * <p>
     * This method is <em>not</em> guaranteed to return the smallest enclosing cap.
     * </p>
     * @return a spherical cap enclosing the polygon
     */
    public EnclosingBall<Sphere2D, S2Point> getEnclosingCap() {

        // handle special cases first
        if (isEmpty()) {
            return new EnclosingBall<Sphere2D, S2Point>(S2Point.PLUS_K, Double.NEGATIVE_INFINITY);
        }
        if (isFull()) {
            return new EnclosingBall<Sphere2D, S2Point>(S2Point.PLUS_K, Double.POSITIVE_INFINITY);
        }

        // as the polygons is neither empty nor full, it has some boundaries and cut hyperplanes
        final BSPTree<Sphere2D> root = getTree(false);
        if (isEmpty(root.getMinus()) && isFull(root.getPlus())) {
            // the polygon covers an hemisphere, and its boundary is one 2π long edge
            final Circle circle = (Circle) root.getCut().getHyperplane();
            return new EnclosingBall<Sphere2D, S2Point>(new S2Point(circle.getPole()).negate(),
                                                        0.5 * FastMath.PI);
        }
        if (isFull(root.getMinus()) && isEmpty(root.getPlus())) {
            // the polygon covers an hemisphere, and its boundary is one 2π long edge
            final Circle circle = (Circle) root.getCut().getHyperplane();
            return new EnclosingBall<Sphere2D, S2Point>(new S2Point(circle.getPole()),
                                                        0.5 * FastMath.PI);
        }

        // gather some inside points, to be used by the encloser
        final List<Vector3D> points = getInsidePoints();

        // extract points from the boundary loops, to be used by the encloser as well
        final List<Vertex> boundary = getBoundaryLoops();
        for (final Vertex loopStart : boundary) {
            int count = 0;
            for (Vertex v = loopStart; count == 0 || v != loopStart; v = v.getOutgoing().getEnd()) {
                ++count;
                points.add(v.getLocation().getVector());
            }
        }

        // find the smallest enclosing 3D sphere
        final SphereGenerator generator = new SphereGenerator();
        final WelzlEncloser<Euclidean3D, Vector3D> encloser =
                new WelzlEncloser<>(getTolerance(), generator);
        EnclosingBall<Euclidean3D, Vector3D> enclosing3D = encloser.enclose(points);
        final Vector3D[] support3D = enclosing3D.getSupport();

        // convert to 3D sphere to spherical cap
        final double r = enclosing3D.getRadius();
        final double h = enclosing3D.getCenter().getNorm();
        if (h < getTolerance()) {
            // the 3D sphere is centered on the unit sphere and covers it
            // fall back to a crude approximation, based only on outside convex cells
            EnclosingBall<Sphere2D, S2Point> enclosingS2 =
                    new EnclosingBall<>(S2Point.PLUS_K, Double.POSITIVE_INFINITY);
            for (Vector3D outsidePoint : getOutsidePoints()) {
                final S2Point outsideS2 = new S2Point(outsidePoint);
                final BoundaryProjection<Sphere2D> projection = projectToBoundary(outsideS2);
                if (FastMath.PI - projection.getOffset() < enclosingS2.getRadius()) {
                    enclosingS2 = new EnclosingBall<>(outsideS2.negate(),
                                                      FastMath.PI - projection.getOffset(),
                                                      (S2Point) projection.getProjected());
                }
            }
            return enclosingS2;
        }
        final S2Point[] support = new S2Point[support3D.length];
        for (int i = 0; i < support3D.length; ++i) {
            support[i] = new S2Point(support3D[i]);
        }

        return new EnclosingBall<>(new S2Point(enclosing3D.getCenter()),
                                   FastMath.acos((1 + h * h - r * r) / (2 * h)),
                                   support);


    }

    /** Gather some inside points.
     * @return list of points known to be strictly in all inside convex cells
     */
    private List<Vector3D> getInsidePoints() {
        final PropertiesComputer pc = new PropertiesComputer(getTolerance());
        getTree(true).visit(pc);
        return pc.getConvexCellsInsidePoints();
    }

    /** Gather some outside points.
     * @return list of points known to be strictly in all outside convex cells
     */
    private List<Vector3D> getOutsidePoints() {
        final SphericalPolygonsSet complement =
                (SphericalPolygonsSet) new RegionFactory<Sphere2D>().getComplement(this);
        final PropertiesComputer pc = new PropertiesComputer(getTolerance());
        complement.getTree(true).visit(pc);
        return pc.getConvexCellsInsidePoints();
    }

}
