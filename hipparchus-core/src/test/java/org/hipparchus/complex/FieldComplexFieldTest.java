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
package org.hipparchus.complex;

import org.hipparchus.Field;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.Map;

public class FieldComplexFieldTest {

    @Test
    public void testZero() {
        Assertions.assertEquals(new FieldComplex<>(Binary64.ZERO), FieldComplexField.getField(Binary64Field.getInstance()).getZero());
    }

    @Test
    public void testOne() {
        Assertions.assertEquals(new FieldComplex<>(Binary64.ONE), FieldComplexField.getField(Binary64Field.getInstance()).getOne());
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    public void testMap() {
        Map<Field<?>, Integer> map = new HashMap<>();
        for (int i = 1; i < 100; ++i) {
            map.put(new FieldComplex<>(new Binary64(i)).getField(), 0);
        }
        // there should be only one field for all values
        FieldComplexField<Binary64> field = FieldComplexField.getField(Binary64Field.getInstance());
        Assertions.assertEquals(1, map.size());
        Assertions.assertEquals(field, map.entrySet().iterator().next().getKey());
        Assertions.assertNotEquals(field, Binary64Field.getInstance());
    }

    @Test
    public void testRunTimeClass() {
        Assertions.assertEquals(Complex.class, ComplexField.getInstance().getRuntimeClass());
    }

}
