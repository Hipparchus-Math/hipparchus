/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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

package org.hipparchus.analysis.differentiation;

import org.hipparchus.Field;
import org.hipparchus.dfp.Dfp;
import org.hipparchus.dfp.DfpField;
import org.junit.Test;

/**
 * Test for class {@link FieldDerivativeStructure} on {@link Dfp}.
 */
public class FieldDerivativeStructureDfpTest extends FieldDerivativeStructureAbstractTest<Dfp> {

    private static final DfpField FIELD = new DfpField(25);

    @Override
    protected Field<Dfp> getField() {
        return FIELD;
    }

    @Override
    @Test
    public void testComposeField() {
        doTestComposeField(new double[] { 5.0e-24, 5.0e-24, 7.0e-24, 2.0e-23, 3.0e-23, 1.0e-100 });
    }

    @Override
    @Test
    public void testComposePrimitive() {
        doTestComposePrimitive(new double[] { 2.0e-14, 3.0e-14, 8.0e-14, 2.0e-13, 9.0e-14, 1.0e-100 });
    }

    @Override
    @Test
    public void testHypotNoOverflow() {
        doTestHypotNoOverflow(65600);
    }

}
