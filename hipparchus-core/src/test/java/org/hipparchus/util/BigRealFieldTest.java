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
package org.hipparchus.util;

import org.hipparchus.Field;
import org.hipparchus.UnitTestUtils;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertSame;

class BigRealFieldTest {

    @Test
    void testZero() {
        assertEquals(BigReal.ZERO, BigRealField.getInstance().getZero());
    }

    @Test
    void testOne() {
        assertEquals(BigReal.ONE, BigRealField.getInstance().getOne());
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    void testMap() {
        Map<Field<?>, Integer> map = new HashMap<>();
        for (int i = 1; i < 100; ++i) {
            map.put(new BigReal(i).getField(), 0);
        }
        // there should be only one field for all values
        assertEquals(1, map.size());
        assertEquals(BigRealField.getInstance(), map.entrySet().iterator().next().getKey());
        assertNotEquals(BigRealField.getInstance(), Binary64Field.getInstance());
    }

    @Test
    void testRunTImeClass() {
        assertEquals(BigReal.class, BigRealField.getInstance().getRuntimeClass());
    }

    @Test
    void testSerial() {
        // deserializing the singleton should give the singleton itself back
        BigRealField field = BigRealField.getInstance();
        assertSame(field, UnitTestUtils.serializeAndRecover(field));
    }

}
