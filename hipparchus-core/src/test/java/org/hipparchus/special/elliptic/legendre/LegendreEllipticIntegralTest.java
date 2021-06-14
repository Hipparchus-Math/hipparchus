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
package org.hipparchus.special.elliptic.legendre;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class LegendreEllipticIntegralTest {

    @Test
    public void testNoConvergence() {
        try {
            LegendreEllipticIntegral.bigK(Double.NaN);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedCoreFormats.CONVERGENCE_FAILED, mise.getSpecifier());
        }
    }

    @Test
    public void testComplementary() {
        for (double k = 0.01; k < 1; k += 0.01) {
            double k1 = LegendreEllipticIntegral.bigK(k);
            double k2 = LegendreEllipticIntegral.bigKPrime(FastMath.sqrt(1 - k * k));
            Assert.assertEquals(k1, k2, FastMath.ulp(k1));
        }
    }

    @Test
    public void testAbramowitzStegunExample3() {
        Assert.assertEquals(3.591545001,
                            LegendreEllipticIntegral.bigK(FastMath.sqrt(80.0 / 81.0)),
                            2.0e-9);
    }

}
