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

import org.hipparchus.Field;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Test for class {@link FieldDerivativeStructure} on {@link Binary64}.
 */
class FieldDerivativeStructureBinary64Test extends FieldDerivativeStructureAbstractTest<Binary64> {

    @Override
    protected Field<Binary64> getField() {
        return Binary64Field.getInstance();
    }

    @Test
    void testNorm() {
        final FDSFactory<Binary64> factory = buildFactory(1, 1);
        final FieldDerivativeStructure<Binary64> fds = new FieldDerivativeStructure<>(factory,
                new Binary64[]{ Binary64.ONE, Binary64.PI});
        final double norm = fds.norm();
        Assertions.assertEquals(1 + FastMath.PI, norm);
    }

    @Override
    @Test
    public void testComposeField() {
        doTestComposeField(new double[] { 1.0e-100, 5.0e-14, 2.0e-13, 3.0e-13, 2.0e-13, 1.0e-100 });
    }

    @Override
    @Test
    public void testComposePrimitive() {
        doTestComposePrimitive(new double[] { 1.0e-100, 5.0e-14, 2.0e-13, 3.0e-13, 2.0e-13, 1.0e-100 });
    }

    @Override
    @Test
    public void testHypotNoOverflow() {
        doTestHypotNoOverflow(250);
    }

    @Test
    public void testIsSmallDimensionMismatch() {
        try {
            buildFactory(3, 2).constant(0.0).isSmall(buildFactory(3, 4).constant(0.0), 1.0e-15);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            Assertions.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, e.getSpecifier());
        }
    }

    @Test
    public void testIsSmall() {
        final FDSFactory<Binary64> factory = buildFactory(2, 2);
        Assertions.assertTrue(factory.variable(1, 2.5).multiply(1.0e-13).isSmall(factory.variable(1, 1.0), 1.0e-12));
        Assertions.assertFalse(factory.variable(1, 2.5).multiply(1.0e-11).isSmall(factory.variable(1, 1.0), 1.0e-12));
    }

}
