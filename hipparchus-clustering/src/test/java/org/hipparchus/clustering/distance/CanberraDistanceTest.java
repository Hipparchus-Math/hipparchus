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
package org.hipparchus.clustering.distance;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Tests for {@link CanberraDistance} class.
 */
class CanberraDistanceTest {
    final DistanceMeasure distance = new CanberraDistance();

    @Test
    void testZero() {
        final double[] a = { 0, 1, -2, 3.4, 5, -6.7, 89 };
        assertEquals(0, distance.compute(a, a), 0d);
    }

    @Test
    void testZero2() {
        final double[] a = { 0, 0 };
        assertEquals(0, distance.compute(a, a), 0d);
    }

    @Test
    void test() {
        final double[] a = { 1, 2, 3, 4, 9 };
        final double[] b = { -5, -6, 7, 4, 3 };
        final double expected = 2.9;
        assertEquals(expected, distance.compute(a, b), 0d);
        assertEquals(expected, distance.compute(b, a), 0d);
    }
}
