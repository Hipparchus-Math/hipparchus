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
package org.hipparchus.special.elliptic;

import org.hipparchus.complex.Complex;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well19937a;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class CarlsonEllipticIntegralTest {

    @Test
    public void testNoConvergenceRf() {
        try {
            CarlsonEllipticIntegral.rF(new Complex(1), new Complex(2), Complex.NaN);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedCoreFormats.CONVERGENCE_FAILED, mise.getSpecifier());
        }
    }

    @Test
    public void testDlmfRf() {
        Complex rf = CarlsonEllipticIntegral.rF(new Complex(1), new Complex(2), new Complex(4));
        Assert.assertEquals(0.6850858166, rf.getRealPart(),      1.0e-10);
        Assert.assertEquals(0.0,          rf.getImaginaryPart(), 1.0e-10);
    }

    @Test
    public void testCarlson1995rF() {

        Complex rf1 = CarlsonEllipticIntegral.rF(new Complex(1), new Complex(2), new Complex(0));
        Assert.assertEquals( 1.3110287771461,  rf1.getRealPart(),      1.0e-13);
        Assert.assertEquals( 0.0,              rf1.getImaginaryPart(), 1.0e-13);

        Complex rf2 = CarlsonEllipticIntegral.rF(new Complex(0.5), new Complex(1), new Complex(0));
        Assert.assertEquals( 1.8540746773014,  rf2.getRealPart(),      1.0e-13);
        Assert.assertEquals( 0.0,              rf2.getImaginaryPart(), 1.0e-13);

        Complex rf3 = CarlsonEllipticIntegral.rF(new Complex(-1, 1), new Complex(0, 1), new Complex(0));
        Assert.assertEquals( 0.79612586584234, rf3.getRealPart(),      1.0e-13);
        Assert.assertEquals(-1.2138566698365,  rf3.getImaginaryPart(), 1.0e-13);

        Complex rf4 = CarlsonEllipticIntegral.rF(new Complex(2), new Complex(3), new Complex(4));
        Assert.assertEquals( 0.58408284167715, rf4.getRealPart(),      1.0e-13);
        Assert.assertEquals( 0.0,              rf4.getImaginaryPart(), 1.0e-13);

        Complex rf5 = CarlsonEllipticIntegral.rF(new Complex(0, 1), new Complex(0, -1), new Complex(2));
        Assert.assertEquals( 1.0441445654064,  rf5.getRealPart(),      1.0e-13);
        Assert.assertEquals( 0.0,              rf5.getImaginaryPart(), 1.0e-13);

        Complex rf6 = CarlsonEllipticIntegral.rF(new Complex(-1, 1), new Complex(0, 1), new Complex(1, -1));
        Assert.assertEquals( 0.93912050218619, rf6.getRealPart(),      1.0e-13);
        Assert.assertEquals(-0.53296252018635, rf6.getImaginaryPart(), 1.0e-13);

    }

    @Test
    public void testNoConvergenceRc() {
        try {
            CarlsonEllipticIntegral.rC(new Complex(1), Complex.NaN);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedCoreFormats.CONVERGENCE_FAILED, mise.getSpecifier());
        }
    }

    @Test
    public void testCarlson1995rC() {

        Complex rc1 = CarlsonEllipticIntegral.rC(new Complex(0), new Complex(0.25));
        Assert.assertEquals(FastMath.PI,       rc1.getRealPart(),      1.0e-15);
        Assert.assertEquals( 0.0,              rc1.getImaginaryPart(), 1.0e-15);

        Complex rc2 = CarlsonEllipticIntegral.rC(new Complex(2.25), new Complex(2));
        Assert.assertEquals(FastMath.log(2),   rc2.getRealPart(),      1.0e-15);
        Assert.assertEquals( 0.0,              rc2.getImaginaryPart(), 1.0e-15);

        Complex rc3 = CarlsonEllipticIntegral.rC(new Complex(0), new Complex(0, 1));
        Assert.assertEquals( 1.1107207345396,  rc3.getRealPart(),      1.0e-13);
        Assert.assertEquals(-1.1107207345396,  rc3.getImaginaryPart(), 1.0e-13);

        Complex rc4 = CarlsonEllipticIntegral.rC(new Complex(0, -1), new Complex(0, 1));
        Assert.assertEquals( 1.2260849569072,  rc4.getRealPart(),      1.0e-13);
        Assert.assertEquals(-0.34471136988768, rc4.getImaginaryPart(), 1.0e-13);

        Complex rc5 = CarlsonEllipticIntegral.rC(new Complex(0.25), new Complex(-2));
        Assert.assertEquals(FastMath.log(2) / 3.0,  rc5.getRealPart(), 1.0e-15);
        // the 1995 paper does not show any imaginary part for this case,
        // but it seems to be non-zero, which is consistent with R_F(x, y, y)

    }

    @Test
    public void testRfRc() {
        RandomGenerator random = new Well19937a(0x7e8041334a8c20edl);
        for (int i = 0; i < 100000; ++i) {
            final Complex x = new Complex(6 * random.nextDouble() - 3,
                                          6 * random.nextDouble() - 3);
            final Complex y = new Complex(6 * random.nextDouble() - 3,
                                          6 * random.nextDouble() - 3);
            final Complex rf = CarlsonEllipticIntegral.rF(x, y, y);
            final Complex rc = CarlsonEllipticIntegral.rC(x, y);
            Assert.assertEquals(0.0, rf.subtract(rc).norm(), 3.0e-14 * rf.norm());
        }
    }

}
