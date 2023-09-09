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

import org.hipparchus.dfp.Dfp;
import org.hipparchus.dfp.DfpField;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well19937a;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test for class {@link FieldUnivariateDerivative2} on {@link Dfp}.
 */
public class FieldUnivariateDerivative2DfpTest extends FieldUnivariateDerivative2AbstractTest<Dfp> {

    private static final DfpField FIELD = new DfpField(25);

    @Override
    protected DfpField getValueField() {
        return FIELD;
    }

    @Test
    public void testHashcode() {
        Assert.assertEquals(-1300667743, build(2, 1, 4).hashCode());
    }

    @Override
    @Test
    public void testLinearCombinationReference() {
        doTestLinearCombinationReference(x -> build(x), 5.0e-9, 4.212e-9);
    }

    @Override
    @Test
    public void testUlp() {
        final RandomGenerator random = new Well19937a(0x36d4f8862421e0e4l);
        for (int i = -300; i < 300; ++i) {
            final double x = FastMath.scalb(2.0 * random.nextDouble() - 1.0, i);
            Assert.assertTrue(FastMath.ulp(x) >= build(x).ulp().getReal());
        }
    }

    @Override
    @Test
    public void testUlpVsDS() {
        // skipped as Dfp is much higher accuracy than double
    }

}
