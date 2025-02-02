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

import org.hipparchus.exception.Localizable;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.geometry.LocalizedGeometryFormats;
import org.hipparchus.geometry.euclidean.twod.Euclidean2D;
import org.hipparchus.geometry.euclidean.twod.PolygonsSet;
import org.hipparchus.geometry.euclidean.twod.SubLine;
import org.hipparchus.geometry.euclidean.twod.Vector2D;
import org.hipparchus.geometry.partitioning.BSPTree;
import org.hipparchus.geometry.partitioning.BSPTreeVisitor;
import org.hipparchus.geometry.partitioning.BoundaryAttribute;
import org.hipparchus.geometry.partitioning.Region;
import org.hipparchus.geometry.partitioning.RegionDumper;
import org.hipparchus.geometry.partitioning.RegionFactory;
import org.hipparchus.geometry.partitioning.RegionParser;
import org.hipparchus.geometry.partitioning.SubHyperplane;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

class PolyhedronsSetTest {

    @Test
    void testBox() {
        PolyhedronsSet tree = new PolyhedronsSet(0, 1, 0, 1, 0, 1, 1.0e-10);
        assertEquals(1.0, tree.getSize(), 1.0e-10);
        assertEquals(6.0, tree.getBoundarySize(), 1.0e-10);
        Vector3D barycenter = tree.getBarycenter();
        assertEquals(0.5, barycenter.getX(), 1.0e-10);
        assertEquals(0.5, barycenter.getY(), 1.0e-10);
        assertEquals(0.5, barycenter.getZ(), 1.0e-10);
        for (double x = -0.25; x < 1.25; x += 0.1) {
            boolean xOK = (x >= 0.0) && (x <= 1.0);
            for (double y = -0.25; y < 1.25; y += 0.1) {
                boolean yOK = (y >= 0.0) && (y <= 1.0);
                for (double z = -0.25; z < 1.25; z += 0.1) {
                    boolean zOK = (z >= 0.0) && (z <= 1.0);
                    Region.Location expected =
                        (xOK && yOK && zOK) ? Region.Location.INSIDE : Region.Location.OUTSIDE;
                    assertEquals(expected, tree.checkPoint(new Vector3D(x, y, z)));
                }
            }
        }
        checkPoints(Region.Location.BOUNDARY, tree, new Vector3D[] {
            new Vector3D(0.0, 0.5, 0.5),
            new Vector3D(1.0, 0.5, 0.5),
            new Vector3D(0.5, 0.0, 0.5),
            new Vector3D(0.5, 1.0, 0.5),
            new Vector3D(0.5, 0.5, 0.0),
            new Vector3D(0.5, 0.5, 1.0)
        });
        checkPoints(Region.Location.OUTSIDE, tree, new Vector3D[] {
            new Vector3D(0.0, 1.2, 1.2),
            new Vector3D(1.0, 1.2, 1.2),
            new Vector3D(1.2, 0.0, 1.2),
            new Vector3D(1.2, 1.0, 1.2),
            new Vector3D(1.2, 1.2, 0.0),
            new Vector3D(1.2, 1.2, 1.0)
        });
    }

    @Test
    void testBRepExtractor() {
        double x = 1.0;
        double y = 2.0;
        double z = 3.0;
        double w = 0.1;
        double l = 1.0;
        PolyhedronsSet polyhedron =
            new PolyhedronsSet(x - l, x + l, y - w, y + w, z - w, z + w, 1.0e-10);
        PolyhedronsSet.BRep brep = polyhedron.getBRep();
        assertEquals(6, brep.getFacets().size());
        assertEquals(8, brep.getVertices().size());
    }

    @Test
    void testEmptyBRepIfEmpty() {
        PolyhedronsSet empty = (PolyhedronsSet) new RegionFactory<Euclidean3D, Vector3D, Plane, SubPlane>().
                getComplement(new PolyhedronsSet(1.0e-10));
        assertTrue(empty.isEmpty());
        assertEquals(0.0, empty.getSize(), 1.0e-10);
        PolyhedronsSet.BRep brep = empty.getBRep();
        assertEquals(0, brep.getFacets().size());
        assertEquals(0, brep.getVertices().size());
    }

    @Test
    void testNoBRepHalfSpace() {
        BSPTree<Euclidean3D, Vector3D, Plane, SubPlane> bsp = new BSPTree<>();
        bsp.insertCut(new Plane(Vector3D.PLUS_K, 1.0e-10));
        bsp.getPlus().setAttribute(Boolean.FALSE);
        bsp.getMinus().setAttribute(Boolean.TRUE);
        PolyhedronsSet polyhedron = new PolyhedronsSet(bsp, 1.0e-10);
        assertEquals(Double.POSITIVE_INFINITY, polyhedron.getSize(), 1.0e-10);
        try {
            polyhedron.getBRep();
            fail("an exception should have been thrown");
        } catch (MathRuntimeException mre) {
            assertEquals(LocalizedGeometryFormats.OUTLINE_BOUNDARY_LOOP_OPEN, mre.getSpecifier());
        }
    }

    @Test
    void testNoBRepUnboundedOctant() {
        BSPTree<Euclidean3D, Vector3D, Plane, SubPlane> bsp = new BSPTree<>();
        bsp.insertCut(new Plane(Vector3D.PLUS_K, 1.0e-10));
        bsp.getPlus().setAttribute(Boolean.FALSE);
        bsp.getMinus().insertCut(new Plane(Vector3D.PLUS_I, 1.0e-10));
        bsp.getMinus().getPlus().setAttribute(Boolean.FALSE);
        bsp.getMinus().getMinus().insertCut(new Plane(Vector3D.PLUS_J, 1.0e-10));
        bsp.getMinus().getMinus().getPlus().setAttribute(Boolean.FALSE);
        bsp.getMinus().getMinus().getMinus().setAttribute(Boolean.TRUE);
        PolyhedronsSet polyhedron = new PolyhedronsSet(bsp, 1.0e-10);
        assertEquals(Double.POSITIVE_INFINITY, polyhedron.getSize(), 1.0e-10);
        try {
            polyhedron.getBRep();
            fail("an exception should have been thrown");
        } catch (MathRuntimeException mre) {
            assertEquals(LocalizedGeometryFormats.OUTLINE_BOUNDARY_LOOP_OPEN, mre.getSpecifier());
        }
    }

    @Test
    void testNoBRepHolesInFacet() {
        double tolerance = 1.0e-10;
        PolyhedronsSet cube       = new PolyhedronsSet(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, tolerance);
        PolyhedronsSet tubeAlongX = new PolyhedronsSet(-2.0, 2.0, -0.5, 0.5, -0.5, 0.5, tolerance);
        PolyhedronsSet tubeAlongY = new PolyhedronsSet(-0.5, 0.5, -2.0, 2.0, -0.5, 0.5, tolerance);
        PolyhedronsSet tubeAlongZ = new PolyhedronsSet(-0.5, 0.5, -0.5, 0.5, -2.0, 2.0, tolerance);
        RegionFactory<Euclidean3D, Vector3D, Plane, SubPlane> factory = new RegionFactory<>();
        PolyhedronsSet cubeWithHoles = (PolyhedronsSet) factory.difference(cube,
                                                                           factory.union(tubeAlongX,
                                                                                         factory.union(tubeAlongY, tubeAlongZ)));
        assertEquals(4.0, cubeWithHoles.getSize(), 1.0e-10);
        try {
            cubeWithHoles.getBRep();
            fail("an exception should have been thrown");
        } catch (MathRuntimeException mre) {
            assertEquals(LocalizedGeometryFormats.FACET_WITH_SEVERAL_BOUNDARY_LOOPS, mre.getSpecifier());
        }
    }

    @Test
    void testTetrahedron() throws MathRuntimeException {
        Vector3D vertex1 = new Vector3D(1, 2, 3);
        Vector3D vertex2 = new Vector3D(2, 2, 4);
        Vector3D vertex3 = new Vector3D(2, 3, 3);
        Vector3D vertex4 = new Vector3D(1, 3, 4);
        PolyhedronsSet tree =
            (PolyhedronsSet) new RegionFactory<Euclidean3D, Vector3D, Plane, SubPlane>().
                    buildConvex(new Plane(vertex3, vertex2, vertex1, 1.0e-10),
                                new Plane(vertex2, vertex3, vertex4, 1.0e-10),
                                new Plane(vertex4, vertex3, vertex1, 1.0e-10),
                                new Plane(vertex1, vertex2, vertex4, 1.0e-10));
        assertEquals(1.0 / 3.0, tree.getSize(), 1.0e-10);
        assertEquals(2.0 * FastMath.sqrt(3.0), tree.getBoundarySize(), 1.0e-10);
        Vector3D barycenter = tree.getBarycenter();
        assertEquals(1.5, barycenter.getX(), 1.0e-10);
        assertEquals(2.5, barycenter.getY(), 1.0e-10);
        assertEquals(3.5, barycenter.getZ(), 1.0e-10);
        double third = 1.0 / 3.0;
        checkPoints(Region.Location.BOUNDARY, tree, new Vector3D[] {
            vertex1, vertex2, vertex3, vertex4,
            new Vector3D(third, vertex1, third, vertex2, third, vertex3),
            new Vector3D(third, vertex2, third, vertex3, third, vertex4),
            new Vector3D(third, vertex3, third, vertex4, third, vertex1),
            new Vector3D(third, vertex4, third, vertex1, third, vertex2)
        });
        checkPoints(Region.Location.OUTSIDE, tree, new Vector3D[] {
            new Vector3D(1, 2, 4),
            new Vector3D(2, 2, 3),
            new Vector3D(2, 3, 4),
            new Vector3D(1, 3, 3)
        });
    }

    @Test
    void testIsometry() throws MathRuntimeException {
        Vector3D vertex1 = new Vector3D(1.1, 2.2, 3.3);
        Vector3D vertex2 = new Vector3D(2.0, 2.4, 4.2);
        Vector3D vertex3 = new Vector3D(2.8, 3.3, 3.7);
        Vector3D vertex4 = new Vector3D(1.0, 3.6, 4.5);
        PolyhedronsSet tree =
            (PolyhedronsSet) new RegionFactory<Euclidean3D, Vector3D, Plane, SubPlane>().
                    buildConvex(new Plane(vertex3, vertex2, vertex1, 1.0e-10),
                                new Plane(vertex2, vertex3, vertex4, 1.0e-10),
                                new Plane(vertex4, vertex3, vertex1, 1.0e-10),
                                new Plane(vertex1, vertex2, vertex4, 1.0e-10));
        Vector3D barycenter = tree.getBarycenter();
        Vector3D s = new Vector3D(10.2, 4.3, -6.7);
        Vector3D c = new Vector3D(-0.2, 2.1, -3.2);
        Rotation r = new Rotation(new Vector3D(6.2, -4.4, 2.1), 0.12, RotationConvention.VECTOR_OPERATOR);

        tree = tree.rotate(c, r).translate(s);

        Vector3D newB = new Vector3D(1.0, s,
                                     1.0, c,
                                     1.0, r.applyTo(barycenter.subtract(c)));
        assertEquals(0.0, newB.subtract(tree.getBarycenter()).getNorm(), 1.0e-10);

        final Vector3D[] expectedV = new Vector3D[] {
            new Vector3D(1.0, s,
                         1.0, c,
                         1.0, r.applyTo(vertex1.subtract(c))),
                         new Vector3D(1.0, s,
                                      1.0, c,
                                      1.0, r.applyTo(vertex2.subtract(c))),
                                      new Vector3D(1.0, s,
                                                   1.0, c,
                                                   1.0, r.applyTo(vertex3.subtract(c))),
                                                   new Vector3D(1.0, s,
                                                                1.0, c,
                                                                1.0, r.applyTo(vertex4.subtract(c)))
        };
        tree.getTree(true).visit(new BSPTreeVisitor<Euclidean3D, Vector3D, Plane, SubPlane>() {

            public Order visitOrder(BSPTree<Euclidean3D, Vector3D, Plane, SubPlane> node) {
                return Order.MINUS_SUB_PLUS;
            }

            public void visitInternalNode(BSPTree<Euclidean3D, Vector3D, Plane, SubPlane> node) {
                @SuppressWarnings("unchecked")
                BoundaryAttribute<Euclidean3D, Vector3D, Plane, SubPlane> attribute =
                    (BoundaryAttribute<Euclidean3D, Vector3D, Plane, SubPlane>) node.getAttribute();
                if (attribute.getPlusOutside() != null) {
                    checkFacet(attribute.getPlusOutside());
                }
                if (attribute.getPlusInside() != null) {
                    checkFacet(attribute.getPlusInside());
                }
            }

            public void visitLeafNode(BSPTree<Euclidean3D, Vector3D, Plane, SubPlane> node) {
            }

            private void checkFacet(SubPlane facet) {
                Plane plane = facet.getHyperplane();
                Vector2D[][] vertices =
                    ((PolygonsSet) facet.getRemainingRegion()).getVertices();
                assertEquals(1, vertices.length);
                for (int i = 0; i < vertices[0].length; ++i) {
                    Vector3D v = plane.toSpace(vertices[0][i]);
                    double d = Double.POSITIVE_INFINITY;
                    for (final Vector3D u : expectedV) {
                        d = FastMath.min(d, v.subtract(u).getNorm());
                    }
                    assertEquals(0, d, 1.0e-10);
                }
            }

        });

    }

    @Test
    void testBuildBox() {
        double x = 1.0;
        double y = 2.0;
        double z = 3.0;
        double w = 0.1;
        double l = 1.0;
        PolyhedronsSet tree =
            new PolyhedronsSet(x - l, x + l, y - w, y + w, z - w, z + w, 1.0e-10);
        Vector3D barycenter = tree.getBarycenter();
        assertEquals(x, barycenter.getX(), 1.0e-10);
        assertEquals(y, barycenter.getY(), 1.0e-10);
        assertEquals(z, barycenter.getZ(), 1.0e-10);
        assertEquals(8 * l * w * w, tree.getSize(), 1.0e-10);
        assertEquals(8 * w * (2 * l + w), tree.getBoundarySize(), 1.0e-10);
    }

    @Test
    void testCross() {

        double x = 1.0;
        double y = 2.0;
        double z = 3.0;
        double w = 0.1;
        double l = 1.0;
        PolyhedronsSet xBeam =
            new PolyhedronsSet(x - l, x + l, y - w, y + w, z - w, z + w, 1.0e-10);
        PolyhedronsSet yBeam =
            new PolyhedronsSet(x - w, x + w, y - l, y + l, z - w, z + w, 1.0e-10);
        PolyhedronsSet zBeam =
            new PolyhedronsSet(x - w, x + w, y - w, y + w, z - l, z + l, 1.0e-10);
        RegionFactory<Euclidean3D, Vector3D, Plane, SubPlane> factory = new RegionFactory<>();
        PolyhedronsSet tree = (PolyhedronsSet) factory.union(xBeam, factory.union(yBeam, zBeam));
        Vector3D barycenter = tree.getBarycenter();

        assertEquals(x, barycenter.getX(), 1.0e-10);
        assertEquals(y, barycenter.getY(), 1.0e-10);
        assertEquals(z, barycenter.getZ(), 1.0e-10);
        assertEquals(8 * w * w * (3 * l - 2 * w), tree.getSize(), 1.0e-10);
        assertEquals(24 * w * (2 * l - w), tree.getBoundarySize(), 1.0e-10);

    }

    @Test
    void testIssue780() throws MathRuntimeException {
        float[] coords = {
            1.000000f, -1.000000f, -1.000000f,
            1.000000f, -1.000000f, 1.000000f,
            -1.000000f, -1.000000f, 1.000000f,
            -1.000000f, -1.000000f, -1.000000f,
            1.000000f, 1.000000f, -1f,
            0.999999f, 1.000000f, 1.000000f,   // 1.000000f, 1.000000f, 1.000000f,
            -1.000000f, 1.000000f, 1.000000f,
            -1.000000f, 1.000000f, -1.000000f};
        int[] indices = {
            0, 1, 2, 0, 2, 3,
            4, 7, 6, 4, 6, 5,
            0, 4, 5, 0, 5, 1,
            1, 5, 6, 1, 6, 2,
            2, 6, 7, 2, 7, 3,
            4, 0, 3, 4, 3, 7};
        ArrayList<SubPlane> subHyperplaneList = new ArrayList<>();
        for (int idx = 0; idx < indices.length; idx += 3) {
            int idxA = indices[idx] * 3;
            int idxB = indices[idx + 1] * 3;
            int idxC = indices[idx + 2] * 3;
            Vector3D v_1 = new Vector3D(coords[idxA], coords[idxA + 1], coords[idxA + 2]);
            Vector3D v_2 = new Vector3D(coords[idxB], coords[idxB + 1], coords[idxB + 2]);
            Vector3D v_3 = new Vector3D(coords[idxC], coords[idxC + 1], coords[idxC + 2]);
            Vector3D[] vertices = {v_1, v_2, v_3};
            Plane polyPlane = new Plane(v_1, v_2, v_3, 1.0e-10);
            ArrayList<SubLine> lines = new ArrayList<>();

            Vector2D[] projPts = new Vector2D[vertices.length];
            for (int ptIdx = 0; ptIdx < projPts.length; ptIdx++) {
                projPts[ptIdx] = polyPlane.toSubSpace(vertices[ptIdx]);
            }

            for (int ptIdx = 0; ptIdx < projPts.length; ptIdx++) {
                lines.add(new SubLine(projPts[ptIdx], projPts[(ptIdx + 1) % projPts.length], 1.0e-10));
            }
            Region<Euclidean2D, Vector2D, org.hipparchus.geometry.euclidean.twod.Line, SubLine> polyRegion =
                    new PolygonsSet(lines, 1.0e-10);
            SubPlane polygon    = new SubPlane(polyPlane, polyRegion);
            subHyperplaneList.add(polygon);
        }
        PolyhedronsSet polyhedronsSet = new PolyhedronsSet(subHyperplaneList, 1.0e-10);
        assertEquals( 8.0, polyhedronsSet.getSize(), 3.0e-6);
        assertEquals(24.0, polyhedronsSet.getBoundarySize(), 5.0e-6);
    }

    @Test
    void testTooThinBox() {
        assertEquals(0.0,
                            new PolyhedronsSet(0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0e-10).getSize(),
                            1.0e-10);
    }

    @Test
    void testWrongUsage() {
        // the following is a wrong usage of the constructor.
        // as explained in the javadoc, the failure is NOT detected at construction
        // time but occurs later on
        PolyhedronsSet ps = new PolyhedronsSet(new BSPTree<>(), 1.0e-10);
        assertNotNull(ps);
        try {
            ps.checkPoint(Vector3D.ZERO);
            fail("an exception should have been thrown");
        } catch (NullPointerException npe) {
            // this is expected
        }
    }

    @Test
    void testDumpParse() throws IOException, ParseException {
        double tol=1e-8;

            Vector3D[] verts=new Vector3D[8];
            double xmin=-1,xmax=1;
            double ymin=-1,ymax=1;
            double zmin=-1,zmax=1;
            verts[0]=new Vector3D(xmin,ymin,zmin);
            verts[1]=new Vector3D(xmax,ymin,zmin);
            verts[2]=new Vector3D(xmax,ymax,zmin);
            verts[3]=new Vector3D(xmin,ymax,zmin);
            verts[4]=new Vector3D(xmin,ymin,zmax);
            verts[5]=new Vector3D(xmax,ymin,zmax);
            verts[6]=new Vector3D(xmax,ymax,zmax);
            verts[7]=new Vector3D(xmin,ymax,zmax);
            //
            int[][] faces=new int[12][];
            faces[0]=new int[]{3,1,0};  // bottom (-z)
            faces[1]=new int[]{1,3,2};  // bottom (-z)
            faces[2]=new int[]{5,7,4};  // top (+z)
            faces[3]=new int[]{7,5,6};  // top (+z)
            faces[4]=new int[]{2,5,1};  // right (+x)
            faces[5]=new int[]{5,2,6};  // right (+x)
            faces[6]=new int[]{4,3,0};  // left (-x)
            faces[7]=new int[]{3,4,7};  // left (-x)
            faces[8]=new int[]{4,1,5};  // front (-y)
            faces[9]=new int[]{1,4,0};  // front (-y)
            faces[10]=new int[]{3,6,2}; // back (+y)
            faces[11]=new int[]{6,3,7}; // back (+y)
            PolyhedronsSet polyset = new PolyhedronsSet(Arrays.asList(verts), Arrays.asList(faces), tol);
            assertEquals(8.0, polyset.getSize(), 1.0e-10);
            assertEquals(24.0, polyset.getBoundarySize(), 1.0e-10);
            String dump = RegionDumper.dump(polyset);
            PolyhedronsSet parsed = RegionParser.parsePolyhedronsSet(dump);
            assertEquals(8.0, parsed.getSize(), 1.0e-10);
            assertEquals(24.0, parsed.getBoundarySize(), 1.0e-10);
            assertTrue(new RegionFactory<Euclidean3D, Vector3D, Plane, SubPlane>().difference(polyset, parsed).isEmpty());
    }

    @Test
    void testConnectedFacets() throws IOException, ParseException {
        InputStream stream = getClass().getResourceAsStream("pentomino-N.ply");
        PLYParser   parser = new PLYParser(stream);
        stream.close();
        PolyhedronsSet polyhedron = new PolyhedronsSet(parser.getVertices(), parser.getFaces(), 1.0e-10);
        assertEquals( 5.0, polyhedron.getSize(), 1.0e-10);
        assertEquals(22.0, polyhedron.getBoundarySize(), 1.0e-10);
    }

    @Test
    void testTooClose() {
        checkError("pentomino-N-too-close.ply", LocalizedGeometryFormats.CLOSE_VERTICES);
    }

    @Test
    void testHole() {
        checkError("pentomino-N-hole.ply", LocalizedGeometryFormats.EDGE_CONNECTED_TO_ONE_FACET);
    }

    @Test
    void testNonPlanar() {
        checkError("pentomino-N-out-of-plane.ply", LocalizedGeometryFormats.OUT_OF_PLANE);
    }

    @Test
    void testOrientation() {
        checkError("pentomino-N-bad-orientation.ply", LocalizedGeometryFormats.FACET_ORIENTATION_MISMATCH);
    }

    @Test
    void testFacet2Vertices() {
        checkError(Arrays.asList(Vector3D.ZERO, Vector3D.PLUS_I, Vector3D.PLUS_J, Vector3D.PLUS_K),
                   Arrays.asList(new int[] { 0, 1, 2 }, new int[] {2, 3}),
                   LocalizedCoreFormats.WRONG_NUMBER_OF_POINTS);
    }

    private void checkError(final String resourceName, final Localizable expected) {
        try {
            InputStream stream = getClass().getResourceAsStream(resourceName);
            PLYParser   parser = new PLYParser(stream);
            stream.close();
            checkError(parser.getVertices(), parser.getFaces(), expected);
        } catch (IOException | ParseException e) {
            fail(e.getLocalizedMessage());
        }
    }

    private void checkError(final List<Vector3D> vertices, final List<int[]> facets,
                            final Localizable expected) {
        try {
            new PolyhedronsSet(vertices, facets, 1.0e-10);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(expected, miae.getSpecifier());
        }
    }

    // issue GEOMETRY-38
    @Test
    void testFirstIntersectionLinesPassThroughBoundaries() {
        // arrange
        Vector3D lowerCorner = Vector3D.ZERO;
        Vector3D upperCorner = new Vector3D(1, 1, 1);
        Vector3D center = new Vector3D(0.5, lowerCorner, 0.5, upperCorner);

        PolyhedronsSet polySet = new PolyhedronsSet(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0e-15);

        Line upDiagonal = new Line(lowerCorner, upperCorner, 1.0e-15);
        Line downDiagonal = upDiagonal.revert();

        // act/assert
        SubPlane upFromOutsideResult = (SubPlane) polySet.firstIntersection(new Vector3D(-1, -1, -1), upDiagonal);
        assertNotNull(upFromOutsideResult);
        assertEquals(0.0,
                            Vector3D.distance(lowerCorner,
                                              upFromOutsideResult.getHyperplane().intersection(upDiagonal)),
                            1.0e-15);

        SubPlane upFromCenterResult = (SubPlane) polySet.firstIntersection(center, upDiagonal);
        assertNotNull(upFromCenterResult);
        assertEquals(0.0,
                            Vector3D.distance(upperCorner,
                                              upFromCenterResult.getHyperplane().intersection(upDiagonal)),
                            1.0e-15);

        SubPlane downFromOutsideResult = (SubPlane) polySet.firstIntersection(new Vector3D(2, 2, 2), downDiagonal);
        assertNotNull(downFromOutsideResult);
        assertEquals(0.0,
                            Vector3D.distance(upperCorner,
                            downFromOutsideResult.getHyperplane().intersection(downDiagonal)),
                            1.0e-15);

        SubPlane downFromCenterResult = (SubPlane) polySet.firstIntersection(center, downDiagonal);
        assertNotNull(downFromCenterResult);
        assertEquals(0.0,
                            Vector3D.distance(lowerCorner,
                                              downFromCenterResult.getHyperplane().intersection(downDiagonal)),
                            1.0e-15);
    }

    @Test
    void testIssue1211() throws IOException, ParseException {

        PolyhedronsSet polyset = RegionParser.parsePolyhedronsSet(loadTestData("issue-1211.bsp"));
        RandomGenerator random = new Well1024a(0xb97c9d1ade21e40aL);
        int nrays = 1000;
        for (int i = 0; i < nrays; i++) {
            Vector3D origin    = Vector3D.ZERO;
            Vector3D direction = new Vector3D(2 * random.nextDouble() - 1,
                                              2 * random.nextDouble() - 1,
                                              2 * random.nextDouble() - 1).normalize();
            Line line = new Line(origin, origin.add(direction), polyset.getTolerance());
            SubHyperplane<Euclidean3D, Vector3D, Plane, SubPlane> plane = polyset.firstIntersection(origin, line);
            if (plane != null) {
                Vector3D intersectionPoint = plane.getHyperplane().intersection(line);
                double dotProduct = direction.dotProduct(intersectionPoint.subtract(origin));
                assertTrue(dotProduct > 0);
            }
        }
    }

    private String loadTestData(final String resourceName)
            throws IOException {
            InputStream stream = getClass().getResourceAsStream(resourceName);
            Reader reader = new InputStreamReader(stream, StandardCharsets.UTF_8);
            StringBuilder builder = new StringBuilder();
            for (int c = reader.read(); c >= 0; c = reader.read()) {
                builder.append((char) c);
            }
            return builder.toString();
        }

    private void checkPoints(Region.Location expected, PolyhedronsSet tree, Vector3D[] points) {
        for (final Vector3D point : points) {
            assertEquals(expected, tree.checkPoint(point));
        }
    }

}
