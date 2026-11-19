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
package org.hipparchus.analysis.interpolation;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class RegularIndexerTest {

    @Test
    public void testGetters() {
        final RegularIndexer sampling = new RegularIndexer(3.0, 5.0, 5);
        Assertions.assertEquals(3.0, sampling.getMin(), 1.0e-15);
        Assertions.assertEquals(5.0, sampling.getMax(), 1.0e-15);
        Assertions.assertEquals(5,   sampling.getN());
        Assertions.assertEquals(0.5, sampling.getStep(), 1.0e-15);
    }

    @Test
    public void testIndexAtNodes() {
        final RegularIndexer sampling = new RegularIndexer(3.0, 5.0, 5);
        Assertions.assertEquals(0, sampling.index(3.0));
        Assertions.assertEquals(1, sampling.index(3.5));
        Assertions.assertEquals(2, sampling.index(4.0));
        Assertions.assertEquals(3, sampling.index(4.5));
        Assertions.assertEquals(4, sampling.index(5.0));
    }

    @Test
    public void testIndexAwayFromNodes() {
        final RegularIndexer sampling = new RegularIndexer(3.0, 5.0, 5);
        Assertions.assertEquals(0, sampling.index(2.750000001));
        Assertions.assertEquals(0, sampling.index(3.249999999));
        Assertions.assertEquals(1, sampling.index(3.250000001));
        Assertions.assertEquals(1, sampling.index(3.749999999));
        Assertions.assertEquals(2, sampling.index(3.750000001));
        Assertions.assertEquals(2, sampling.index(4.249999999));
        Assertions.assertEquals(3, sampling.index(4.250000001));
        Assertions.assertEquals(3, sampling.index(4.749999999));
        Assertions.assertEquals(4, sampling.index(4.750000001));
        Assertions.assertEquals(4, sampling.index(5.249999999));
    }

    @Test
    public void testIndexTooLow() {
        final RegularIndexer sampling = new RegularIndexer(300.0, 500.0, 5);
        try {
            sampling.index(274.9999999);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, miae.getSpecifier());
            Assertions.assertEquals(-1L, miae.getParts()[0]);
            Assertions.assertEquals( 0L, miae.getParts()[1]);
            Assertions.assertEquals( 4L, miae.getParts()[2]);
        }
    }

    @Test
    public void testIndexTooHigh() {
        final RegularIndexer sampling = new RegularIndexer(300.0, 500.0, 5);
        try {
            sampling.index(525.0000001);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, miae.getSpecifier());
            Assertions.assertEquals( 5L, miae.getParts()[0]);
            Assertions.assertEquals( 0L, miae.getParts()[1]);
            Assertions.assertEquals( 4L, miae.getParts()[2]);
        }
    }

    @Test
    public void testCoordinateInRange() {
        final RegularIndexer sampling = new RegularIndexer(300.0, 500.0, 5);
        Assertions.assertEquals(300.0, sampling.coordinate(0), 1.0e-15);
        Assertions.assertEquals(350.0, sampling.coordinate(1), 1.0e-15);
        Assertions.assertEquals(400.0, sampling.coordinate(2), 1.0e-15);
        Assertions.assertEquals(450.0, sampling.coordinate(3), 1.0e-15);
        Assertions.assertEquals(500.0, sampling.coordinate(4), 1.0e-15);
    }

    @Test
    public void testCoordinateTooLow() {
        final RegularIndexer sampling = new RegularIndexer(300.0, 500.0, 5);
        try {
            sampling.coordinate(-1);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, miae.getSpecifier());
            Assertions.assertEquals(-1L, miae.getParts()[0]);
            Assertions.assertEquals( 0L, miae.getParts()[1]);
            Assertions.assertEquals( 4L, miae.getParts()[2]);
        }
    }

    @Test
    public void testCoordinateTooHigh() {
        final RegularIndexer sampling = new RegularIndexer(300.0, 500.0, 5);
        try {
            sampling.coordinate(5);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, miae.getSpecifier());
            Assertions.assertEquals( 5L, miae.getParts()[0]);
            Assertions.assertEquals( 0L, miae.getParts()[1]);
            Assertions.assertEquals( 4L, miae.getParts()[2]);
        }
    }

}
