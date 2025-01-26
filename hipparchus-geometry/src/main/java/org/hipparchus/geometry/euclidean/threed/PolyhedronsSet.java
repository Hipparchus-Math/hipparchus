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
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.geometry.LocalizedGeometryFormats;
import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.euclidean.oned.Euclidean1D;
import org.hipparchus.geometry.euclidean.oned.Vector1D;
import org.hipparchus.geometry.euclidean.twod.Euclidean2D;
import org.hipparchus.geometry.euclidean.twod.PolygonsSet;
import org.hipparchus.geometry.euclidean.twod.SubLine;
import org.hipparchus.geometry.euclidean.twod.Vector2D;
import org.hipparchus.geometry.partitioning.AbstractRegion;
import org.hipparchus.geometry.partitioning.BSPTree;
import org.hipparchus.geometry.partitioning.BSPTreeVisitor;
import org.hipparchus.geometry.partitioning.BoundaryAttribute;
import org.hipparchus.geometry.partitioning.Hyperplane;
import org.hipparchus.geometry.partitioning.Region;
import org.hipparchus.geometry.partitioning.RegionFactory;
import org.hipparchus.geometry.partitioning.SubHyperplane;
import org.hipparchus.geometry.partitioning.Transform;
import org.hipparchus.util.FastMath;

/** This class represents a 3D region: a set of polyhedrons.
 */
public class PolyhedronsSet extends AbstractRegion<Euclidean3D, Vector3D, Euclidean2D, Vector2D> {

    /** Build a polyhedrons set representing the whole real line.
     * @param tolerance tolerance below which points are considered identical
     */
    public PolyhedronsSet(final double tolerance) {
        super(tolerance);
    }

    /** Build a polyhedrons set from a BSP tree.
     * <p>The leaf nodes of the BSP tree <em>must</em> have a
     * {@code Boolean} attribute representing the inside status of
     * the corresponding cell (true for inside cells, false for outside
     * cells). In order to avoid building too many small objects, it is
     * recommended to use the predefined constants
     * {@code Boolean.TRUE} and {@code Boolean.FALSE}</p>
     * <p>
     * This constructor is aimed at expert use, as building the tree may
     * be a difficult task. It is not intended for general use and for
     * performances reasons does not check thoroughly its input, as this would
     * require walking the full tree each time. Failing to provide a tree with
     * the proper attributes, <em>will</em> therefore generate problems like
     * {@link NullPointerException} or {@link ClassCastException} only later on.
     * This limitation is known and explains why this constructor is for expert
     * use only. The caller does have the responsibility to provided correct arguments.
     * </p>
     * @param tree inside/outside BSP tree representing the region
     * @param tolerance tolerance below which points are considered identical
     */
    public PolyhedronsSet(final BSPTree<Euclidean3D, Vector3D> tree, final double tolerance) {
        super(tree, tolerance);
    }

    /** Build a polyhedrons set from a Boundary REPresentation (B-rep) specified by sub-hyperplanes.
     * <p>The boundary is provided as a collection of {@link
     * SubHyperplane sub-hyperplanes}. Each sub-hyperplane has the
     * interior part of the region on its minus side and the exterior on
     * its plus side.</p>
     * <p>The boundary elements can be in any order, and can form
     * several non-connected sets (like for example polyhedrons with holes
     * or a set of disjoint polyhedrons considered as a whole). In
     * fact, the elements do not even need to be connected together
     * (their topological connections are not used here). However, if the
     * boundary does not really separate an inside open from an outside
     * open (open having here its topological meaning), then subsequent
     * calls to the {@link Region#checkPoint(Point) checkPoint} method will
     * not be meaningful anymore.</p>
     * <p>If the boundary is empty, the region will represent the whole
     * space.</p>
     * @param boundary collection of boundary elements, as a
     * collection of {@link SubHyperplane SubHyperplane} objects
     * @param tolerance tolerance below which points are considered identical
     */
    public PolyhedronsSet(final Collection<SubHyperplane<Euclidean3D, Vector3D>> boundary,
                          final double tolerance) {
        super(boundary, tolerance);
    }

    /** Build a polyhedrons set from a Boundary REPresentation (B-rep) specified by connected vertices.
     * <p>
     * The boundary is provided as a list of vertices and a list of facets.
     * Each facet is specified as an integer array containing the arrays vertices
     * indices in the vertices list. Each facet normal is oriented by right hand
     * rule to the facet vertices list.
     * </p>
     * <p>
     * Some basic sanity checks are performed but not everything is thoroughly
     * assessed, so it remains under caller responsibility to ensure the vertices
     * and facets are consistent and properly define a polyhedrons set.
     * </p>
     * @param vertices list of polyhedrons set vertices
     * @param facets list of facets, as vertices indices in the vertices list
     * @param tolerance tolerance below which points are considered identical
     * @exception MathIllegalArgumentException if some basic sanity checks fail
     */
    public PolyhedronsSet(final List<Vector3D> vertices, final List<int[]> facets,
                          final double tolerance) {
        super(buildBoundary(vertices, facets, tolerance), tolerance);
    }

    /** Build a polyhedrons set from a Boundary REPresentation (B-rep) specified by connected vertices.
     * <p>
     * Some basic sanity checks are performed but not everything is thoroughly
     * assessed, so it remains under caller responsibility to ensure the vertices
     * and facets are consistent and properly define a polyhedrons set.
     * </p>
     * @param brep Boundary REPresentation of the polyhedron to build
     * @param tolerance tolerance below which points are considered identical
     * @exception MathIllegalArgumentException if some basic sanity checks fail
     * @since 1.2
     */
    public PolyhedronsSet(final BRep brep, final double tolerance) {
        super(buildBoundary(brep.getVertices(), brep.getFacets(), tolerance), tolerance);
    }

    /** Build a parallellepipedic box.
     * @param xMin low bound along the x direction
     * @param xMax high bound along the x direction
     * @param yMin low bound along the y direction
     * @param yMax high bound along the y direction
     * @param zMin low bound along the z direction
     * @param zMax high bound along the z direction
     * @param tolerance tolerance below which points are considered identical
     */
    public PolyhedronsSet(final double xMin, final double xMax,
                          final double yMin, final double yMax,
                          final double zMin, final double zMax,
                          final double tolerance) {
        super(buildBoundary(xMin, xMax, yMin, yMax, zMin, zMax, tolerance), tolerance);
    }

    /** Build a parallellepipedic box boundary.
     * @param xMin low bound along the x direction
     * @param xMax high bound along the x direction
     * @param yMin low bound along the y direction
     * @param yMax high bound along the y direction
     * @param zMin low bound along the z direction
     * @param zMax high bound along the z direction
     * @param tolerance tolerance below which points are considered identical
     * @return boundary tree
     */
    private static BSPTree<Euclidean3D, Vector3D> buildBoundary(final double xMin, final double xMax,
                                                                final double yMin, final double yMax,
                                                                final double zMin, final double zMax,
                                                                final double tolerance) {
        if ((xMin >= xMax - tolerance) || (yMin >= yMax - tolerance) || (zMin >= zMax - tolerance)) {
            // too thin box, build an empty polygons set
            return new BSPTree<>(Boolean.FALSE);
        }
        final Plane pxMin = new Plane(new Vector3D(xMin, 0,    0),   Vector3D.MINUS_I, tolerance);
        final Plane pxMax = new Plane(new Vector3D(xMax, 0,    0),   Vector3D.PLUS_I,  tolerance);
        final Plane pyMin = new Plane(new Vector3D(0,    yMin, 0),   Vector3D.MINUS_J, tolerance);
        final Plane pyMax = new Plane(new Vector3D(0,    yMax, 0),   Vector3D.PLUS_J,  tolerance);
        final Plane pzMin = new Plane(new Vector3D(0,    0,   zMin), Vector3D.MINUS_K, tolerance);
        final Plane pzMax = new Plane(new Vector3D(0,    0,   zMax), Vector3D.PLUS_K,  tolerance);
        final Region<Euclidean3D, Vector3D> boundary =
        new RegionFactory<Euclidean3D, Vector3D>().buildConvex(pxMin, pxMax, pyMin, pyMax, pzMin, pzMax);
        return boundary.getTree(false);
    }

    /** Build boundary from vertices and facets.
     * @param vertices list of polyhedrons set vertices
     * @param facets list of facets, as vertices indices in the vertices list
     * @param tolerance tolerance below which points are considered identical
     * @return boundary as a list of sub-hyperplanes
     * @exception MathIllegalArgumentException if some basic sanity checks fail
     */
    private static List<SubHyperplane<Euclidean3D, Vector3D>> buildBoundary(final List<Vector3D> vertices,
                                                                            final List<int[]> facets,
                                                                            final double tolerance) {

        // check vertices distances
        for (int i = 0; i < vertices.size() - 1; ++i) {
            final Vector3D vi = vertices.get(i);
            for (int j = i + 1; j < vertices.size(); ++j) {
                if (Vector3D.distance(vi, vertices.get(j)) <= tolerance) {
                    throw new MathIllegalArgumentException(LocalizedGeometryFormats.CLOSE_VERTICES,
                                                           vi.getX(), vi.getY(), vi.getZ());
                }
            }
        }

        // find how vertices are referenced by facets
        final int[][] references = findReferences(vertices, facets);

        // find how vertices are linked together by edges along the facets they belong to
        final int[][] successors = successors(vertices, facets, references);

        // check edges orientations
        for (int vA = 0; vA < vertices.size(); ++vA) {
            for (final int vB : successors[vA]) {

                if (vB >= 0) {
                    // when facets are properly oriented, if vB is the successor of vA on facet f1,
                    // then there must be an adjacent facet f2 where vA is the successor of vB
                    boolean found = false;
                    for (final int v : successors[vB]) {
                        found = found || (v == vA);
                    }
                    if (!found) {
                        final Vector3D start = vertices.get(vA);
                        final Vector3D end   = vertices.get(vB);
                        throw new MathIllegalArgumentException(LocalizedGeometryFormats.EDGE_CONNECTED_TO_ONE_FACET,
                                                               start.getX(), start.getY(), start.getZ(),
                                                               end.getX(),   end.getY(),   end.getZ());
                    }
                }
            }
        }

        final List<SubHyperplane<Euclidean3D, Vector3D>> boundary = new ArrayList<>();

        for (final int[] facet : facets) {

            // define facet plane from the first 3 points
            Plane plane = new Plane(vertices.get(facet[0]), vertices.get(facet[1]), vertices.get(facet[2]),
                                    tolerance);

            // check all points are in the plane
            final Vector2D[] two2Points = new Vector2D[facet.length];
            for (int i = 0 ; i < facet.length; ++i) {
                final Vector3D v = vertices.get(facet[i]);
                if (!plane.contains(v)) {
                    throw new MathIllegalArgumentException(LocalizedGeometryFormats.OUT_OF_PLANE,
                                                           v.getX(), v.getY(), v.getZ());
                }
                two2Points[i] = plane.toSubSpace(v);
            }

            // create the polygonal facet
            boundary.add(new SubPlane(plane, new PolygonsSet(tolerance, two2Points)));

        }

        return boundary;

    }

    /** Find the facets that reference each edges.
     * @param vertices list of polyhedrons set vertices
     * @param facets list of facets, as vertices indices in the vertices list
     * @return references array such that r[v][k] = f for some k if facet f contains vertex v
     * @exception MathIllegalArgumentException if some facets have fewer than 3 vertices
     */
    private static int[][] findReferences(final List<Vector3D> vertices, final List<int[]> facets) {

        // find the maximum number of facets a vertex belongs to
        final int[] nbFacets = new int[vertices.size()];
        int maxFacets  = 0;
        for (final int[] facet : facets) {
            if (facet.length < 3) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.WRONG_NUMBER_OF_POINTS,
                                                    3, facet.length, true);
            }
            for (final int index : facet) {
                maxFacets = FastMath.max(maxFacets, ++nbFacets[index]);
            }
        }

        // set up the references array
        final int[][] references = new int[vertices.size()][maxFacets];
        for (int[] r : references) {
            Arrays.fill(r, -1);
        }
        for (int f = 0; f < facets.size(); ++f) {
            for (final int v : facets.get(f)) {
                // vertex v is referenced by facet f
                int k = 0;
                while (k < maxFacets && references[v][k] >= 0) {
                    ++k;
                }
                references[v][k] = f;
            }
        }

        return references;

    }

    /** Find the successors of all vertices among all facets they belong to.
     * @param vertices list of polyhedrons set vertices
     * @param facets list of facets, as vertices indices in the vertices list
     * @param references facets references array
     * @return indices of vertices that follow vertex v in some facet (the array
     * may contain extra entries at the end, set to negative indices)
     * @exception MathIllegalArgumentException if the same vertex appears more than
     * once in the successors list (which means one facet orientation is wrong)
     */
    private static int[][] successors(final List<Vector3D> vertices, final List<int[]> facets,
                                      final int[][] references) {

        // create an array large enough
        final int[][] successors = new int[vertices.size()][references[0].length];
        for (final int[] s : successors) {
            Arrays.fill(s, -1);
        }

        for (int v = 0; v < vertices.size(); ++v) {
            for (int k = 0; k < successors[v].length && references[v][k] >= 0; ++k) {

                // look for vertex v
                final int[] facet = facets.get(references[v][k]);
                int i = 0;
                while (i < facet.length && facet[i] != v) {
                    ++i;
                }

                // we have found vertex v, we deduce its successor on current facet
                successors[v][k] = facet[(i + 1) % facet.length];
                for (int l = 0; l < k; ++l) {
                    if (successors[v][l] == successors[v][k]) {
                        final Vector3D start = vertices.get(v);
                        final Vector3D end   = vertices.get(successors[v][k]);
                        throw new MathIllegalArgumentException(LocalizedGeometryFormats.FACET_ORIENTATION_MISMATCH,
                                                               start.getX(), start.getY(), start.getZ(),
                                                               end.getX(),   end.getY(),   end.getZ());
                    }
                }

            }
        }

        return successors;

    }

    /** {@inheritDoc} */
    @Override
    public PolyhedronsSet buildNew(final BSPTree<Euclidean3D, Vector3D> tree) {
        return new PolyhedronsSet(tree, getTolerance());
    }

    /** Get the boundary representation of the instance.
     * <p>
     * The boundary representation can be extracted <em>only</em> from
     * bounded polyhedrons sets. If the polyhedrons set is unbounded,
     * a {@link MathRuntimeException} will be thrown.
     * </p>
     * <p>
     * The boundary representation extracted is not minimal, as for
     * example canonical facets may be split into several smaller
     * independent sub-facets sharing the same plane and connected by
     * their edges.
     * </p>
     * <p>
     * As the {@link BRep B-Rep} representation does not support
     * facets with several boundary loops (for example facets with
     * holes), an exception is triggered when attempting to extract
     * B-Rep from such complex polyhedrons sets.
     * </p>
     * @return boundary representation of the instance
     * @exception MathRuntimeException if polyhedrons is unbounded
     * @since 1.2
     */
    public BRep getBRep() throws MathRuntimeException {
        BRepExtractor extractor = new BRepExtractor(getTolerance());
        getTree(true).visit(extractor);
        return extractor.getBRep();
    }

    /** Visitor extracting BRep. */
    private static class BRepExtractor implements BSPTreeVisitor<Euclidean3D, Vector3D> {

        /** Tolerance for vertices identification. */
        private final double tolerance;

        /** Extracted vertices. */
        private final List<Vector3D> vertices;

        /** Extracted facets. */
        private final List<int[]> facets;

        /** Simple constructor.
         * @param tolerance tolerance for vertices identification
         */
        BRepExtractor(final double tolerance) {
            this.tolerance = tolerance;
            this.vertices  = new ArrayList<>();
            this.facets    = new ArrayList<>();
        }

        /** Get the BRep.
         * @return extracted BRep
         */
        public BRep getBRep() {
            return new BRep(vertices, facets);
        }

        /** {@inheritDoc} */
        @Override
        public Order visitOrder(final BSPTree<Euclidean3D, Vector3D> node) {
            return Order.MINUS_SUB_PLUS;
        }

        /** {@inheritDoc} */
        @Override
        public void visitInternalNode(final BSPTree<Euclidean3D, Vector3D> node) {
            @SuppressWarnings("unchecked")
            final BoundaryAttribute<Euclidean3D, Vector3D> attribute =
                (BoundaryAttribute<Euclidean3D, Vector3D>) node.getAttribute();
            if (attribute.getPlusOutside() != null) {
                addContribution(attribute.getPlusOutside(), false);
            }
            if (attribute.getPlusInside() != null) {
                addContribution(attribute.getPlusInside(), true);
            }
        }

        /** {@inheritDoc} */
        @Override
        public void visitLeafNode(final BSPTree<Euclidean3D, Vector3D> node) {
        }

        /** Add the contribution of a boundary facet.
         * @param facet boundary facet
         * @param reversed if true, the facet has the inside on its plus side
         * @exception MathRuntimeException if facet is unbounded
         */
        private void addContribution(final SubHyperplane<Euclidean3D, Vector3D> facet, final boolean reversed)
            throws MathRuntimeException {

            final Plane plane = (Plane) facet.getHyperplane();
            final PolygonsSet polygon = (PolygonsSet) ((SubPlane) facet).getRemainingRegion();
            final Vector2D[][] loops2D = polygon.getVertices();
            if (loops2D.length == 0) {
                throw new MathRuntimeException(LocalizedGeometryFormats.OUTLINE_BOUNDARY_LOOP_OPEN);
            } else if (loops2D.length > 1) {
                throw new MathRuntimeException(LocalizedGeometryFormats.FACET_WITH_SEVERAL_BOUNDARY_LOOPS);
            } else {
                for (final Vector2D[] loop2D : polygon.getVertices()) {
                    final int[] loop3D = new int[loop2D.length];
                    for (int i = 0; i < loop2D.length ; ++i) {
                        if (loop2D[i] == null) {
                            throw new MathRuntimeException(LocalizedGeometryFormats.OUTLINE_BOUNDARY_LOOP_OPEN);
                        }
                        loop3D[reversed ? loop2D.length - 1 - i : i] = getVertexIndex(plane.toSpace(loop2D[i]));
                    }
                    facets.add(loop3D);
                }
            }

        }

        /** Get the index of a vertex.
         * @param vertex vertex as a 3D point
         * @return index of the vertex
         */
        private int getVertexIndex(final Vector3D vertex) {

            for (int i = 0; i < vertices.size(); ++i) {
                if (Vector3D.distance(vertex, vertices.get(i)) <= tolerance) {
                    // the vertex is already known
                    return i;
                }
            }

            // the vertex is a new one, add it
            vertices.add(vertex);
            return vertices.size() - 1;

        }

    }

    /** {@inheritDoc} */
    @Override
    protected void computeGeometricalProperties() {

        // compute the contribution of all boundary facets
        getTree(true).visit(new FacetsContributionVisitor());

        if (getSize() < 0) {
            // the polyhedrons set as a finite outside
            // surrounded by an infinite inside
            setSize(Double.POSITIVE_INFINITY);
            setBarycenter(Vector3D.NaN);
        } else {
            // the polyhedrons set is finite, apply the remaining scaling factors
            setSize(getSize() / 3.0);
            setBarycenter(new Vector3D(1.0 / (4 * getSize()), getBarycenter()));
        }

    }

    /** Visitor computing geometrical properties. */
    private class FacetsContributionVisitor implements BSPTreeVisitor<Euclidean3D, Vector3D> {

        /** Simple constructor. */
        FacetsContributionVisitor() {
            setSize(0);
            setBarycenter(new Vector3D(0, 0, 0));
        }

        /** {@inheritDoc} */
        @Override
        public Order visitOrder(final BSPTree<Euclidean3D, Vector3D> node) {
            return Order.MINUS_SUB_PLUS;
        }

        /** {@inheritDoc} */
        @Override
        public void visitInternalNode(final BSPTree<Euclidean3D, Vector3D> node) {
            @SuppressWarnings("unchecked")
            final BoundaryAttribute<Euclidean3D, Vector3D> attribute =
                (BoundaryAttribute<Euclidean3D, Vector3D>) node.getAttribute();
            if (attribute.getPlusOutside() != null) {
                addContribution(attribute.getPlusOutside(), false);
            }
            if (attribute.getPlusInside() != null) {
                addContribution(attribute.getPlusInside(), true);
            }
        }

        /** {@inheritDoc} */
        @Override
        public void visitLeafNode(final BSPTree<Euclidean3D, Vector3D> node) {
        }

        /** Add he contribution of a boundary facet.
         * @param facet boundary facet
         * @param reversed if true, the facet has the inside on its plus side
         */
        private void addContribution(final SubHyperplane<Euclidean3D, Vector3D> facet, final boolean reversed) {

            final Region<Euclidean2D, Vector2D> polygon = ((SubPlane) facet).getRemainingRegion();
            final double area    = polygon.getSize();

            if (Double.isInfinite(area)) {
                setSize(Double.POSITIVE_INFINITY);
                setBarycenter(Vector3D.NaN);
            } else {

                final Plane    plane  = (Plane) facet.getHyperplane();
                final Vector3D facetB = plane.toSpace(polygon.getBarycenter());
                double   scaled = area * facetB.dotProduct(plane.getNormal());
                if (reversed) {
                    scaled = -scaled;
                }

                setSize(getSize() + scaled);
                setBarycenter(new Vector3D(1.0, getBarycenter(), scaled, facetB));

            }

        }

    }

    /** Get the first sub-hyperplane crossed by a semi-infinite line.
     * @param point start point of the part of the line considered
     * @param line line to consider (contains point)
     * @return the first sub-hyperplane crossed by the line after the
     * given point, or null if the line does not intersect any
     * sub-hyperplane
     */
    public SubHyperplane<Euclidean3D, Vector3D> firstIntersection(final Vector3D point, final Line line) {
        return recurseFirstIntersection(getTree(true), point, line);
    }

    /** Get the first sub-hyperplane crossed by a semi-infinite line.
     * @param node current node
     * @param point start point of the part of the line considered
     * @param line line to consider (contains point)
     * @return the first sub-hyperplane crossed by the line after the
     * given point, or null if the line does not intersect any
     * sub-hyperplane
     */
    private SubHyperplane<Euclidean3D, Vector3D> recurseFirstIntersection(final BSPTree<Euclidean3D, Vector3D> node,
                                                                          final Vector3D point, final Line line) {

        final SubHyperplane<Euclidean3D, Vector3D> cut = node.getCut();
        if (cut == null) {
            return null;
        }
        final BSPTree<Euclidean3D, Vector3D> minus = node.getMinus();
        final BSPTree<Euclidean3D, Vector3D> plus  = node.getPlus();
        final Plane                plane = (Plane) cut.getHyperplane();

        // establish search order
        final double offset = plane.getOffset(point);
        final boolean in    = FastMath.abs(offset) < getTolerance();
        final BSPTree<Euclidean3D, Vector3D> near;
        final BSPTree<Euclidean3D, Vector3D> far;
        if (offset < 0) {
            near = minus;
            far  = plus;
        } else {
            near = plus;
            far  = minus;
        }

        if (in) {
            // search in the cut hyperplane
            final SubHyperplane<Euclidean3D, Vector3D> facet = boundaryFacet(point, node);
            if (facet != null) {
                return facet;
            }
        }

        // search in the near branch
        final SubHyperplane<Euclidean3D, Vector3D> crossed = recurseFirstIntersection(near, point, line);
        if (crossed != null) {
            return crossed;
        }

        if (!in) {
            // search in the cut hyperplane
            final Vector3D hit3D = plane.intersection(line);
            if (hit3D != null && line.getAbscissa(hit3D) > line.getAbscissa(point)) {
                final SubHyperplane<Euclidean3D, Vector3D> facet = boundaryFacet(hit3D, node);
                if (facet != null) {
                    return facet;
                }
            }
        }

        // search in the far branch
        return recurseFirstIntersection(far, point, line);

    }

    /** Check if a point belongs to the boundary part of a node.
     * @param point point to check
     * @param node node containing the boundary facet to check
     * @return the boundary facet this points belongs to (or null if it
     * does not belong to any boundary facet)
     */
    private SubHyperplane<Euclidean3D, Vector3D> boundaryFacet(final Vector3D point,
                                                               final BSPTree<Euclidean3D, Vector3D> node) {
        final Vector2D point2D = ((Plane) node.getCut().getHyperplane()).toSubSpace(point);
        @SuppressWarnings("unchecked")
        final BoundaryAttribute<Euclidean3D, Vector3D> attribute =
            (BoundaryAttribute<Euclidean3D, Vector3D>) node.getAttribute();
        if ((attribute.getPlusOutside() != null) &&
            (((SubPlane) attribute.getPlusOutside()).getRemainingRegion().checkPoint(point2D) != Location.OUTSIDE)) {
            return attribute.getPlusOutside();
        }
        if ((attribute.getPlusInside() != null) &&
            (((SubPlane) attribute.getPlusInside()).getRemainingRegion().checkPoint(point2D) != Location.OUTSIDE)) {
            return attribute.getPlusInside();
        }
        return null;
    }

    /** Rotate the region around the specified point.
     * <p>The instance is not modified, a new instance is created.</p>
     * @param center rotation center
     * @param rotation vectorial rotation operator
     * @return a new instance representing the rotated region
     */
    public PolyhedronsSet rotate(final Vector3D center, final Rotation rotation) {
        return (PolyhedronsSet) applyTransform(new RotationTransform(center, rotation));
    }

    /** 3D rotation as a Transform. */
    private static class RotationTransform
        implements Transform<Euclidean3D, Vector3D, Euclidean2D, Vector2D> {

        /** Center point of the rotation. */
        private final Vector3D   center;

        /** Vectorial rotation. */
        private final Rotation   rotation;

        /** Cached original hyperplane. */
        private Plane cachedOriginal;

        /** Cached 2D transform valid inside the cached original hyperplane. */
        private Transform<Euclidean2D, Vector2D, Euclidean1D, Vector1D> cachedTransform;

        /** Build a rotation transform.
         * @param center center point of the rotation
         * @param rotation vectorial rotation
         */
        RotationTransform(final Vector3D center, final Rotation rotation) {
            this.center   = center;
            this.rotation = rotation;
        }

        /** {@inheritDoc} */
        @Override
        public Vector3D apply(final Vector3D point) {
            final Vector3D delta = point.subtract(center);
            return new Vector3D(1.0, center, 1.0, rotation.applyTo(delta));
        }

        /** {@inheritDoc} */
        @Override
        public Plane apply(final Hyperplane<Euclidean3D, Vector3D> hyperplane) {
            return ((Plane) hyperplane).rotate(center, rotation);
        }

        /** {@inheritDoc} */
        @Override
        public SubHyperplane<Euclidean2D, Vector2D> apply(final SubHyperplane<Euclidean2D, Vector2D> sub,
                                                          final Hyperplane<Euclidean3D, Vector3D> original,
                                                          final Hyperplane<Euclidean3D, Vector3D> transformed) {
            if (original != cachedOriginal) {
                // we have changed hyperplane, reset the in-hyperplane transform

                final Plane    oPlane = (Plane) original;
                final Plane    tPlane = (Plane) transformed;
                final Vector3D p00    = oPlane.getOrigin();
                final Vector3D p10    = oPlane.toSpace(new Vector2D(1.0, 0.0));
                final Vector3D p01    = oPlane.toSpace(new Vector2D(0.0, 1.0));
                final Vector2D tP00   = tPlane.toSubSpace(apply(p00));
                final Vector2D tP10   = tPlane.toSubSpace(apply(p10));
                final Vector2D tP01   = tPlane.toSubSpace(apply(p01));

                cachedOriginal  = (Plane) original;
                cachedTransform =
                        org.hipparchus.geometry.euclidean.twod.Line.getTransform(tP10.getX() - tP00.getX(),
                                                                                 tP10.getY() - tP00.getY(),
                                                                                 tP01.getX() - tP00.getX(),
                                                                                 tP01.getY() - tP00.getY(),
                                                                                 tP00.getX(),
                                                                                 tP00.getY());

            }
            return ((SubLine) sub).applyTransform(cachedTransform);
        }

    }

    /** Translate the region by the specified amount.
     * <p>The instance is not modified, a new instance is created.</p>
     * @param translation translation to apply
     * @return a new instance representing the translated region
     */
    public PolyhedronsSet translate(final Vector3D translation) {
        return (PolyhedronsSet) applyTransform(new TranslationTransform(translation));
    }

    /** 3D translation as a transform. */
    private static class TranslationTransform
        implements Transform<Euclidean3D, Vector3D, Euclidean2D, Vector2D> {

        /** Translation vector. */
        private final Vector3D   translation;

        /** Cached original hyperplane. */
        private Plane cachedOriginal;

        /** Cached 2D transform valid inside the cached original hyperplane. */
        private Transform<Euclidean2D, Vector2D, Euclidean1D, Vector1D>  cachedTransform;

        /** Build a translation transform.
         * @param translation translation vector
         */
        TranslationTransform(final Vector3D translation) {
            this.translation = translation;
        }

        /** {@inheritDoc} */
        @Override
        public Vector3D apply(final Vector3D point) {
            return new Vector3D(1.0, point, 1.0, translation);
        }

        /** {@inheritDoc} */
        @Override
        public Plane apply(final Hyperplane<Euclidean3D, Vector3D> hyperplane) {
            return ((Plane) hyperplane).translate(translation);
        }

        /** {@inheritDoc} */
        @Override
        public SubHyperplane<Euclidean2D, Vector2D> apply(final SubHyperplane<Euclidean2D, Vector2D> sub,
                                                          final Hyperplane<Euclidean3D, Vector3D> original,
                                                          final Hyperplane<Euclidean3D, Vector3D> transformed) {
            if (original != cachedOriginal) {
                // we have changed hyperplane, reset the in-hyperplane transform

                final Plane   oPlane = (Plane) original;
                final Plane   tPlane = (Plane) transformed;
                final Vector2D shift  = tPlane.toSubSpace(apply(oPlane.getOrigin()));

                cachedOriginal  = (Plane) original;
                cachedTransform =
                        org.hipparchus.geometry.euclidean.twod.Line.getTransform(1, 0, 0, 1,
                                                                                           shift.getX(),
                                                                                           shift.getY());

            }

            return ((SubLine) sub).applyTransform(cachedTransform);

        }

    }

    /** Container for Boundary REPresentation (B-Rep).
     * <p>
     * The boundary is provided as a list of vertices and a list of facets.
     * Each facet is specified as an integer array containing the arrays vertices
     * indices in the vertices list. Each facet normal is oriented by right hand
     * rule to the facet vertices list.
     * </p>
     * @see PolyhedronsSet#PolyhedronsSet(BSPTree, double)
     * @see PolyhedronsSet#getBRep()
     * @since 1.2
     */
    public static class BRep {

        /** List of polyhedrons set vertices. */
        private final List<Vector3D> vertices;

        /** List of facets, as vertices indices in the vertices list. */
        private final List<int[]> facets;

        /** Simple constructor.
         * @param vertices list of polyhedrons set vertices
         * @param facets list of facets, as vertices indices in the vertices list
         */
        public BRep(final List<Vector3D> vertices, final List<int[]> facets) {
            this.vertices = vertices;
            this.facets   = facets;
        }

        /** Get the extracted vertices.
         * @return extracted vertices
         */
        public List<Vector3D> getVertices() {
            return vertices;
        }

        /** Get the extracted facets.
         * @return extracted facets
         */
        public List<int[]> getFacets() {
            return facets;
        }

    }

}
