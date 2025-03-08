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
package org.hipparchus.geometry;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.geometry.euclidean.twod.Vector2D;
import org.hipparchus.geometry.spherical.twod.S2Point;
import org.hipparchus.util.CombinatoricsUtils;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class GeometryTest {

    @Test
    public void testBarycenterException() {
        try {
            final List<S2Point> points = Collections.emptyList();
            Geometry.barycenter(points);
            Assertions.fail("Should have thrown an exception");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedCoreFormats.NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE, miae.getSpecifier());
            Assertions.assertEquals(0, miae.getParts()[0]);
        }
    }

    @Test
    public void testBarycenterIdentity() {
        final Vector3D p = new Vector3D(1.0, 2.0, 3.0);
        Assertions.assertSame(p, Geometry.barycenter(Collections.singletonList(p)));
    }

    @Test
    void testBarycenter2D() {
        final List<Vector2D> points = Arrays.asList(new Vector2D(740.0, 1070.0),
                                                    new Vector2D(1010.0, 1330.0),
                                                    new Vector2D(1240.0, 1160.0),
                                                    new Vector2D(1250.0,  800.0),
                                                    new Vector2D( 890.0,  720.0));
        final Vector2D naive = barycenterNaive2D(points);
        CombinatoricsUtils.
                permutations(points).
                map(Geometry::barycenter).
                forEach(barycenter -> Assertions.assertEquals(0.0, naive.distance(barycenter), 2.0e-13));

    }

    @Test
    void testBarycenter3D() {
        final List<Vector3D> points = Arrays.asList(new Vector3D( 740.0, 1070.0, -12.0),
                                                    new Vector3D(1010.0, 1330.0,   0.0),
                                                    new Vector3D(1240.0, 1160.0,   6.0),
                                                    new Vector3D(1250.0,  800.0,   6.0));
        final Vector3D naive = barycenterNaive3D(points);
        CombinatoricsUtils.
                permutations(points).
                map(Geometry::barycenter).
                forEach(barycenter -> Assertions.assertEquals(0.0, naive.distance(barycenter), 1.0e-15));

    }

    private Vector2D barycenterNaive2D(final List<Vector2D> points) {
        double sumX = 0;
        double sumY = 0;
        for (final Vector2D point : points) {
            sumX += point.getX();
            sumY += point.getY();
        }
        return new Vector2D(sumX / points.size(), sumY / points.size());
    }

    private Vector3D barycenterNaive3D(final List<Vector3D> points) {
        double sumX = 0;
        double sumY = 0;
        double sumZ = 0;
        for (final Vector3D point : points) {
            sumX += point.getX();
            sumY += point.getY();
            sumZ += point.getZ();
        }
        return new Vector3D(sumX / points.size(), sumY / points.size(), sumZ / points.size());
    }

}
