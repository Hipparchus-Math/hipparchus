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
package org.hipparchus.analysis.integration.gauss;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test of the {@link FieldLegendreRuleFactory}.
 *
 */
public class FieldLegendreTest {
    private static final FieldGaussIntegratorFactory<Binary64> factory = new FieldGaussIntegratorFactory<>(Binary64Field.getInstance());

    @Test
    public void testTooLArgeNumberOfPoints() {
        try {
            factory.legendre(10000, new Binary64(0), new Binary64(Math.PI / 2));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.NUMBER_TOO_LARGE, miae.getSpecifier());
            Assert.assertEquals(10000, ((Integer) miae.getParts()[0]).intValue());
            Assert.assertEquals(1000,  ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    public void testCos() {
        final FieldGaussIntegrator<Binary64> integrator = factory.legendre(7, new Binary64(0), new Binary64(Math.PI / 2));
        final double s = integrator.integrate(x -> FastMath.cos(x)).getReal();
        // System.out.println("s=" + s + " e=" + 1);
        Assert.assertEquals(1, s, Math.ulp(1d));
    }


    @Test
    public void testInverse() {

        final Binary64 lo = new Binary64(12.34);
        final Binary64 hi = new Binary64(456.78);

        final FieldGaussIntegrator<Binary64> integrator = factory.legendre(60, lo, hi);
        final double s = integrator.integrate(x -> x.reciprocal()).getReal();
        final double expected = FastMath.log(hi).subtract(FastMath.log(lo)).getReal();
        // System.out.println("s=" + s + " e=" + expected);
        Assert.assertEquals(expected, s, 1e-14);
    }
}
