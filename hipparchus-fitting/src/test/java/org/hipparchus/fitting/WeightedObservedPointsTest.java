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
package org.hipparchus.fitting;

import org.hipparchus.util.Precision;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests {@link WeightedObservedPoints}.
 *
 */
class WeightedObservedPointsTest {
    @Test
    void testAdd1() {
        final WeightedObservedPoints store = new WeightedObservedPoints();

        final double x = 1.2;
        final double y = 34.56;
        final double w = 0.789;

        store.add(w, x, y);

        assertTrue(lastElementIsSame(store, new WeightedObservedPoint(w, x, y)));
    }

    @Test
    void testAdd2() {
        final WeightedObservedPoints store = new WeightedObservedPoints();

        final double x = 1.2;
        final double y = 34.56;
        final double w = 0.789;

        store.add(new WeightedObservedPoint(w, x, y));

        assertTrue(lastElementIsSame(store, new WeightedObservedPoint(w, x, y)));
    }

    @Test
    void testAdd3() {
        final WeightedObservedPoints store = new WeightedObservedPoints();

        final double x = 1.2;
        final double y = 34.56;

        store.add(x, y);

        assertTrue(lastElementIsSame(store, new WeightedObservedPoint(1, x, y)));
    }

    @Test
    void testClear() {
        final WeightedObservedPoints store = new WeightedObservedPoints();

        store.add(new WeightedObservedPoint(1, 2, 3));
        store.add(new WeightedObservedPoint(2, -1, -2));
        assertEquals(2, store.toList().size());

        store.clear();
        assertEquals(0, store.toList().size());
    }

    // Ensure that an instance returned by "toList()" is independent from
    // the original container.
    @Test
    void testToListCopy() {
        final WeightedObservedPoints store = new WeightedObservedPoints();

        store.add(new WeightedObservedPoint(1, 2, 3));
        store.add(new WeightedObservedPoint(2, -3, -4));

        final List<WeightedObservedPoint> list = store.toList();
        assertEquals(2, list.size());

        // Adding an element to "list" has no impact on "store".
        list.add(new WeightedObservedPoint(1.2, 3.4, 5.6));
        assertNotEquals(list.size(), store.toList().size());

        // Clearing "store" has no impact on "list".
        store.clear();
        assertNotEquals(0, list.size());
    }

    /**
     * Checks that the contents of the last element is equal to the
     * contents of {@code p}.
     *
     * @param store Container.
     * @param point Observation.
     * @return {@code true} if both elements have the same contents.
     */
    private boolean lastElementIsSame(WeightedObservedPoints store,
                                      WeightedObservedPoint point) {
        final List<WeightedObservedPoint> list = store.toList();
        final WeightedObservedPoint lastPoint = list.get(list.size() - 1);

        if (!Precision.equals(lastPoint.getX(), point.getX())) {
            return false;
        }
        if (!Precision.equals(lastPoint.getY(), point.getY())) {
            return false;
        }
        if (!Precision.equals(lastPoint.getWeight(), point.getWeight())) {
            return false;
        }

        return true;
    }
}
