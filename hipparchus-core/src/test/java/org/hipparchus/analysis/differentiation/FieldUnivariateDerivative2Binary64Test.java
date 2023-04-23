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
import org.junit.Assert;
import org.junit.Test;

/**
 * Test for class {@link FieldUnivariateDerivative2} on {@link Binary64}.
 */
public class FieldUnivariateDerivative2Binary64Test extends FieldUnivariateDerivative2AbstractTest<Binary64> {

    @Override
    protected Binary64Field getValueField() {
        return Binary64Field.getInstance();
    }

    @Test
    public void testHashcode() {
        Assert.assertEquals(905969981, build(2, 1, 4).hashCode());
    }

    @Override
    @Test
    public void testLinearCombinationReference() {
        doTestLinearCombinationReference(x -> build(x), 5.0e-16, 1.0);
    }
}
