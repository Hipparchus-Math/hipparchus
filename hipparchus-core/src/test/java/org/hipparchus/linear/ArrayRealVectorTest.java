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
package org.hipparchus.linear;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the {@link ArrayRealVector} class.
 *
 */
public class ArrayRealVectorTest extends RealVectorAbstractTest {


    @Override
    public RealVector create(final double[] data) {
        return new ArrayRealVector(data, true);
    }

    @Test
    void testConstructors() {
        final double[] vec1 = {1d, 2d, 3d};
        final double[] vec3 = {7d, 8d, 9d};
        final double[] vec4 = {1d, 2d, 3d, 4d, 5d, 6d, 7d, 8d, 9d};
        final Double[] dvec1 = {1d, 2d, 3d, 4d, 5d, 6d, 7d, 8d, 9d};

        ArrayRealVector v0 = new ArrayRealVector();
        assertEquals(0, v0.getDimension(), "testData len");

        ArrayRealVector v1 = new ArrayRealVector(7);
        assertEquals(7, v1.getDimension(), "testData len");
        assertEquals(0.0, v1.getEntry(6), 0, "testData is 0.0 ");

        ArrayRealVector v2 = new ArrayRealVector(5, 1.23);
        assertEquals(5, v2.getDimension(), "testData len");
        assertEquals(1.23, v2.getEntry(4), 0, "testData is 1.23 ");

        ArrayRealVector v3 = new ArrayRealVector(vec1);
        assertEquals(3, v3.getDimension(), "testData len");
        assertEquals(2.0, v3.getEntry(1), 0, "testData is 2.0 ");

        ArrayRealVector v3_bis = new ArrayRealVector(vec1, true);
        assertEquals(3, v3_bis.getDimension(), "testData len");
        assertEquals(2.0, v3_bis.getEntry(1), 0, "testData is 2.0 ");
        assertNotSame(v3_bis.getDataRef(), vec1);
        assertNotSame(v3_bis.toArray(), vec1);

        ArrayRealVector v3_ter = new ArrayRealVector(vec1, false);
        assertEquals(3, v3_ter.getDimension(), "testData len");
        assertEquals(2.0, v3_ter.getEntry(1), 0, "testData is 2.0 ");
        assertSame(v3_ter.getDataRef(), vec1);
        assertNotSame(v3_ter.toArray(), vec1);

        ArrayRealVector v4 = new ArrayRealVector(vec4, 3, 2);
        assertEquals(2, v4.getDimension(), "testData len");
        assertEquals(4.0, v4.getEntry(0), 0, "testData is 4.0 ");
        try {
            new ArrayRealVector(vec4, 8, 3);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected behavior
        }

        RealVector v5_i = new ArrayRealVector(dvec1);
        assertEquals(9, v5_i.getDimension(), "testData len");
        assertEquals(9.0, v5_i.getEntry(8), 0, "testData is 9.0 ");

        ArrayRealVector v5 = new ArrayRealVector(dvec1);
        assertEquals(9, v5.getDimension(), "testData len");
        assertEquals(9.0, v5.getEntry(8), 0, "testData is 9.0 ");

        ArrayRealVector v6 = new ArrayRealVector(dvec1, 3, 2);
        assertEquals(2, v6.getDimension(), "testData len");
        assertEquals(4.0, v6.getEntry(0), 0, "testData is 4.0 ");
        try {
            new ArrayRealVector(dvec1, 8, 3);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected behavior
        }

        ArrayRealVector v7 = new ArrayRealVector(v1);
        assertEquals(7, v7.getDimension(), "testData len");
        assertEquals(0.0, v7.getEntry(6), 0, "testData is 0.0 ");

        RealVectorTestImpl v7_i = new RealVectorTestImpl(vec1);

        ArrayRealVector v7_2 = new ArrayRealVector(v7_i);
        assertEquals(3, v7_2.getDimension(), "testData len");
        assertEquals(2.0d, v7_2.getEntry(1), 0, "testData is 0.0 ");

        ArrayRealVector v8 = new ArrayRealVector(v1, true);
        assertEquals(7, v8.getDimension(), "testData len");
        assertEquals(0.0, v8.getEntry(6), 0, "testData is 0.0 ");
        assertNotSame(v1.getDataRef(), v8.getDataRef(), "testData not same object ");

        ArrayRealVector v8_2 = new ArrayRealVector(v1, false);
        assertEquals(7, v8_2.getDimension(), "testData len");
        assertEquals(0.0, v8_2.getEntry(6), 0, "testData is 0.0 ");
        assertEquals(v1.getDataRef(), v8_2.getDataRef(), "testData same object ");

        ArrayRealVector v9 = new ArrayRealVector(v1, v3);
        assertEquals(10, v9.getDimension(), "testData len");
        assertEquals(1.0, v9.getEntry(7), 0, "testData is 1.0 ");

        ArrayRealVector v10 = new ArrayRealVector(v2, new RealVectorTestImpl(vec3));
        assertEquals(8, v10.getDimension(), "testData len");
        assertEquals(1.23, v10.getEntry(4), 0, "testData is 1.23 ");
        assertEquals(7.0, v10.getEntry(5), 0, "testData is 7.0 ");

        ArrayRealVector v11 = new ArrayRealVector(new RealVectorTestImpl(vec3), v2);
        assertEquals(8, v11.getDimension(), "testData len");
        assertEquals(9.0, v11.getEntry(2), 0, "testData is 9.0 ");
        assertEquals(1.23, v11.getEntry(3), 0, "testData is 1.23 ");

        ArrayRealVector v12 = new ArrayRealVector(v2, vec3);
        assertEquals(8, v12.getDimension(), "testData len");
        assertEquals(1.23, v12.getEntry(4), 0, "testData is 1.23 ");
        assertEquals(7.0, v12.getEntry(5), 0, "testData is 7.0 ");

        ArrayRealVector v13 = new ArrayRealVector(vec3, v2);
        assertEquals(8, v13.getDimension(), "testData len");
        assertEquals(9.0, v13.getEntry(2), 0, "testData is 9.0 ");
        assertEquals(1.23, v13.getEntry(3), 0, "testData is 1.23 ");

        ArrayRealVector v14 = new ArrayRealVector(vec3, vec4);
        assertEquals(12, v14.getDimension(), "testData len");
        assertEquals(9.0, v14.getEntry(2), 0, "testData is 9.0 ");
        assertEquals(1.0, v14.getEntry(3), 0, "testData is 1.0 ");

    }

    @Test
    void testGetDataRef() {
        final double[] data = {1d, 2d, 3d, 4d};
        final ArrayRealVector v = new ArrayRealVector(data);
        v.getDataRef()[0] = 0d;
        assertEquals(0d, v.getEntry(0), 0, "");
    }

    @Test
    void testPredicates() {

        assertEquals(create(new double[] { Double.NaN, 1, 2 }).hashCode(),
                     create(new double[] { 0, Double.NaN, 2 }).hashCode());

        assertTrue(create(new double[] { Double.NaN, 1, 2 }).hashCode() !=
                   create(new double[] { 0, 1, 2 }).hashCode());
    }

    @Test
    void testZeroVectors() {
        assertEquals(0, new ArrayRealVector(new double[0]).getDimension());
        assertEquals(0, new ArrayRealVector(new double[0], true).getDimension());
        assertEquals(0, new ArrayRealVector(new double[0], false).getDimension());
    }
}
