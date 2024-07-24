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

import org.hipparchus.geometry.euclidean.twod.Vector2D;
import org.hipparchus.geometry.partitioning.RegionFactory;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class OutlineExtractorTest {

    @Test
    public void testBox() {

        PolyhedronsSet tree = new PolyhedronsSet(0, 1, 0, 1, 0, 1, 1.0e-10);
        Assertions.assertEquals(1.0, tree.getSize(), 1.0e-10);
        Assertions.assertEquals(6.0, tree.getBoundarySize(), 1.0e-10);
        Vector2D[][] outline = new OutlineExtractor(Vector3D.PLUS_I, Vector3D.PLUS_J).getOutline(tree);
        Assertions.assertEquals(1, outline.length);
        Assertions.assertEquals(4, outline[0].length);

        Vector2D[] expected = new Vector2D[] {
            new Vector2D(0.0, 0.0),
            new Vector2D(0.0, 1.0),
            new Vector2D(1.0, 1.0),
            new Vector2D(1.0, 0.0)
        };
        for (final Vector2D vertex : outline[0]) {
            for (int j = 0; j < expected.length; ++j) {
                if (expected[j] != null && Vector2D.distance(vertex, expected[j]) < 1.0e-10) {
                    expected[j] = null;
                }
            }
        }
        for (final Vector2D e : expected) {
            Assertions.assertNull(e);
        }

        // observed from a diagonal, the cube appears as an hexagon
        Rotation r = new Rotation(Vector3D.PLUS_I, Vector3D.PLUS_J,
                                  new Vector3D(1, 1, 0), new Vector3D(0, 1, 1));
        Vector2D[][] outlineSkew = new OutlineExtractor(r.applyTo(Vector3D.PLUS_I),
                                                        r.applyTo(Vector3D.PLUS_J)).
                                   getOutline(tree);
        Assertions.assertEquals(1, outlineSkew.length);
        int n = outlineSkew[0].length;
        Assertions.assertEquals(6, n);
        for (int i = 0; i < n; ++i) {
            Vector2D v1 = outlineSkew[0][i];
            Vector2D v2 = outlineSkew[0][(i + n - 1) % n];
            Assertions.assertEquals(FastMath.sqrt(2.0 / 3.0), Vector2D.distance(v1, v2), 1.0e-10);
        }

    }

    @Test
    public void testHolesInFacet() {
        double tolerance = 1.0e-10;
        PolyhedronsSet cube       = new PolyhedronsSet(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, tolerance);
        PolyhedronsSet tubeAlongX = new PolyhedronsSet(-2.0, 2.0, -0.5, 0.5, -0.5, 0.5, tolerance);
        PolyhedronsSet tubeAlongY = new PolyhedronsSet(-0.5, 0.5, -2.0, 2.0, -0.5, 0.5, tolerance);
        PolyhedronsSet tubeAlongZ = new PolyhedronsSet(-0.5, 0.5, -0.5, 0.5, -2.0, 2.0, tolerance);
        RegionFactory<Euclidean3D> factory = new RegionFactory<>();
        PolyhedronsSet cubeWithHoles = (PolyhedronsSet) factory.difference(cube,
                                                                           factory.union(tubeAlongX,
                                                                                         factory.union(tubeAlongY, tubeAlongZ)));
        Assertions.assertEquals(4.0, cubeWithHoles.getSize(), 1.0e-10);
        Vector2D[][] outline = new OutlineExtractor(Vector3D.PLUS_I, Vector3D.PLUS_J).getOutline(cubeWithHoles);
        Assertions.assertEquals(2, outline.length);
        Assertions.assertEquals(4, outline[0].length);
        Assertions.assertEquals(4, outline[1].length);
    }

}
