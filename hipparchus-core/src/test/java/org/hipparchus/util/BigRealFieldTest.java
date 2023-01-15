/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
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


import java.util.HashMap;
import java.util.Map;

import org.hipparchus.Field;
import org.hipparchus.UnitTestUtils;
import org.junit.Assert;
import org.junit.Test;

public class BigRealFieldTest {

    @Test
    public void testZero() {
        Assert.assertEquals(BigReal.ZERO, BigRealField.getInstance().getZero());
    }

    @Test
    public void testOne() {
        Assert.assertEquals(BigReal.ONE, BigRealField.getInstance().getOne());
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    public void testMap() {
        Map<Field<?>, Integer> map = new HashMap<>();
        for (int i = 1; i < 100; ++i) {
            map.put(new BigReal(i).getField(), 0);
        }
        // there should be only one field for all values
        Assert.assertEquals(1, map.size());
        Assert.assertTrue(BigRealField.getInstance().equals(map.entrySet().iterator().next().getKey()));
        Assert.assertFalse(BigRealField.getInstance().equals(Binary64Field.getInstance()));
    }

    @Test
    public void testRunTImeClass() {
        Assert.assertEquals(BigReal.class, BigRealField.getInstance().getRuntimeClass());
    }

    @Test
    public void testSerial() {
        // deserializing the singleton should give the singleton itself back
        BigRealField field = BigRealField.getInstance();
        Assert.assertTrue(field == UnitTestUtils.serializeAndRecover(field));
    }

}
