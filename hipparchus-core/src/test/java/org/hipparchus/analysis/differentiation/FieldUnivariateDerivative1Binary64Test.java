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

package org.hipparchus.analysis.differentiation;

import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test for class {@link FieldUnivariateDerivative1} on {@link Binary64}.
 */
class FieldUnivariateDerivative1Binary64Test extends FieldUnivariateDerivative1AbstractTest<Binary64> {

    @Test
    void test2Norm() {
        final FieldUnivariateDerivative1<Binary64> fud2 = new FieldUnivariateDerivative1<>(Binary64.ONE, Binary64.PI);
        assertEquals(1 + FastMath.PI, fud2.norm());
    }

    @Override
    protected Binary64Field getValueField() {
        return Binary64Field.getInstance();
    }

    @Test
    void testHashcode() {
        assertEquals(2108686789, build(2, 1).hashCode());
    }

    @Override
    @Test
    public void testLinearCombinationReference() {
        doTestLinearCombinationReference(this::build, 5.0e-9, 1.0);
    }

    @Test
    public void testIsSmall() {
        Assertions.assertTrue(build(1.0e-13, 2.5e-13).isSmall(build(1.0, 1.0), 1.0e-12));
        Assertions.assertFalse(build(1.0e-11, 2.5e-13).isSmall(build(1.0, 1.0), 1.0e-12));
        Assertions.assertFalse(build(1.0e-13, 2.5e-11).isSmall(build(1.0, 1.0), 1.0e-12));
    }

}
