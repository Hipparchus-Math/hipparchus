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
package org.hipparchus;

import java.util.function.DoubleFunction;

import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.random.Well19937a;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.FieldSinhCosh;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.SinCos;
import org.junit.Assert;
import org.junit.Test;

public abstract class CalculusFieldElementAbstractTest<T extends CalculusFieldElement<T>> {

    protected abstract T build(double x);

    @Test
    public void testAddField() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                checkRelative(x + y, build(x).add(build(y)));
            }
        }
    }

    @Test
    public void testAddDouble() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                checkRelative(x + y, build(x).add(y));
            }
        }
    }

    @Test
    public void testSubtractField() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                checkRelative(x - y, build(x).subtract(build(y)));
            }
        }
    }

    @Test
    public void testSubtractDouble() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                checkRelative(x - y, build(x).subtract(y));
            }
        }
    }

    @Test
    public void testMultiplyField() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                checkRelative(x * y, build(x).multiply(build(y)));
            }
        }
    }

    @Test
    public void testMultiplyDouble() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                checkRelative(x * y, build(x).multiply(y));
            }
        }
    }

    @Test
    public void testMultiplyInt() {
        for (double x = -3; x < 3; x += 0.2) {
            for (int y = -10; y < 10; y += 1) {
                checkRelative(x * y, build(x).multiply(y));
            }
        }
    }

    @Test
    public void testDivideField() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                checkRelative(x / y, build(x).divide(build(y)));
            }
        }
    }

    @Test
    public void testDivideDouble() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                    checkRelative(x / y, build(x).divide(y));
            }
        }
    }

    @Test
    public void testToDegrees() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.toDegrees(x), build(x).toDegrees());
        }
    }

    @Test
    public void testToRadians() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.toRadians(x), build(x).toRadians());
        }
    }

    @Test
    public void testCos() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.cos(x), build(x).cos());
        }
    }

    @Test
    public void testAcos() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.acos(x), build(x).acos());
        }
    }

    @Test
    public void testSin() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.sin(x), build(x).sin());
        }
    }

    @Test
    public void testAsin() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.asin(x), build(x).asin());
        }
    }

    @Test
    public void testSinCos() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            FieldSinCos<T> sinCos = build(x).sinCos();
            checkRelative(FastMath.sin(x), sinCos.sin());
            checkRelative(FastMath.cos(x), sinCos.cos());
        }
    }

    @Test
    public void testSinCosNaN() {
        FieldSinCos<T> sinCos = build(Double.NaN).sinCos();
        Assert.assertTrue(sinCos.sin().isNaN());
        Assert.assertTrue(sinCos.cos().isNaN());
    }

    @Test
    public void testSinCosSum() {
        final RandomGenerator random = new Well19937a(0x4aab62a42c9eb940l);
        for (int i = 0; i < 10000; ++i) {
            final double alpha  = 10.0 * (2.0 * random.nextDouble() - 1.0);
            final double beta   = 10.0 * (2.0 * random.nextDouble() - 1.0);
            final T      alphaT = build(alpha);
            final T      betaT  = build(beta);
            final SinCos scSum = SinCos.sum(FastMath.sinCos(alpha), FastMath.sinCos(beta));
            final FieldSinCos<T> scSumT = FieldSinCos.sum(alphaT.sinCos(), betaT.sinCos());
            checkRelative(scSum.sin(), scSumT.sin());
            checkRelative(scSum.cos(), scSumT.cos());
        }
    }

    @Test
    public void testSinCosdifference() {
        final RandomGenerator random = new Well19937a(0x589aaf49471b03d5l);
        for (int i = 0; i < 10000; ++i) {
            final double alpha  = 10.0 * (2.0 * random.nextDouble() - 1.0);
            final double beta   = 10.0 * (2.0 * random.nextDouble() - 1.0);
            final T      alphaT = build(alpha);
            final T      betaT  = build(beta);
            final SinCos scDifference = SinCos.difference(FastMath.sinCos(alpha), FastMath.sinCos(beta));
            final FieldSinCos<T> scDifferenceT = FieldSinCos.difference(alphaT.sinCos(), betaT.sinCos());
            checkRelative(scDifference.sin(), scDifferenceT.sin());
            checkRelative(scDifference.cos(), scDifferenceT.cos());
        }
    }

    @Test
    public void testTan() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.tan(x), build(x).tan());
        }
    }

    @Test
    public void testAtan() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.atan(x), build(x).atan());
        }
    }

    @Test
    public void testAtan2() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                checkRelative(FastMath.atan2(y, x), build(y).atan2(build(x)));
            }
        }
    }

    @Test
    public void testAtan2SpecialCases() {
        checkRelative(FastMath.atan2(+0.0, +0.0), build(+0.0).atan2(build(+0.0)));
        checkRelative(FastMath.atan2(-0.0, +0.0), build(-0.0).atan2(build(+0.0)));
        checkRelative(FastMath.atan2(+0.0, -0.0), build(+0.0).atan2(build(-0.0)));
        checkRelative(FastMath.atan2(-0.0, -0.0), build(-0.0).atan2(build(-0.0)));
    }

    @Test
    public void testCosh() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.cosh(x), build(x).cosh());
        }
    }

    @Test
    public void testAcosh() {
        for (double x = 1.1; x < 5.0; x += 0.05) {
            checkRelative(FastMath.acosh(x), build(x).acosh());
        }
    }

    @Test
    public void testSinh() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.sinh(x), build(x).sinh());
        }
    }

    @Test
    public void testAsinh() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.asinh(x), build(x).asinh());
        }
    }

    @Test
    public void testSinhCosh() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            FieldSinhCosh<T> sinhCosh = build(x).sinhCosh();
            checkRelative(FastMath.sinh(x), sinhCosh.sinh());
            checkRelative(FastMath.cosh(x), sinhCosh.cosh());
        }
    }

    @Test
    public void testSinhCoshNaN() {
        FieldSinhCosh<T> sinhCosh = build(Double.NaN).sinhCosh();
        Assert.assertTrue(sinhCosh.sinh().isNaN());
        Assert.assertTrue(sinhCosh.cosh().isNaN());
    }

    @Test
    public void testTanh() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.tanh(x), build(x).tanh());
        }
    }

    @Test
    public void testAtanh() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.atanh(x), build(x).atanh());
        }
    }

    @Test
    public void testSqrt() {
        for (double x = 0.01; x < 0.9; x += 0.05) {
            checkRelative(FastMath.sqrt(x), build(x).sqrt());
        }
    }

    @Test
    public void testCbrt() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.cbrt(x), build(x).cbrt());
        }
    }

    @Test
    public void testHypot() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                checkRelative(FastMath.hypot(x, y), build(x).hypot(build(y)));
            }
        }
    }

    @Test
    public void testHypotSpecialCases() {
        Assert.assertTrue(Double.isNaN(build(Double.NaN).hypot(build(0)).getReal()));
        Assert.assertTrue(Double.isNaN(build(0).hypot(build(Double.NaN)).getReal()));
        Assert.assertEquals(Double.POSITIVE_INFINITY, build(Double.POSITIVE_INFINITY).hypot(build(0)).getReal(), 1.0);
        Assert.assertEquals(Double.POSITIVE_INFINITY, build(Double.NEGATIVE_INFINITY).hypot(build(0)).getReal(), 1.0);
        Assert.assertEquals(Double.POSITIVE_INFINITY, build(Double.POSITIVE_INFINITY).hypot(build(Double.NaN)).getReal(), 1.0);
        Assert.assertEquals(Double.POSITIVE_INFINITY, build(Double.NEGATIVE_INFINITY).hypot(build(Double.NaN)).getReal(), 1.0);
        Assert.assertEquals(Double.POSITIVE_INFINITY, build(0).hypot(build(Double.POSITIVE_INFINITY)).getReal(), 1.0);
        Assert.assertEquals(Double.POSITIVE_INFINITY, build(0).hypot(build(Double.NEGATIVE_INFINITY)).getReal(), 1.0);
        Assert.assertEquals(Double.POSITIVE_INFINITY, build(Double.NaN).hypot(build(Double.POSITIVE_INFINITY)).getReal(), 1.0);
        Assert.assertEquals(Double.POSITIVE_INFINITY, build(Double.NaN).hypot(build(Double.NEGATIVE_INFINITY)).getReal(), 1.0);
    }

    @Test
    public void testRootN() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            for (int n = 1; n < 5; ++n) {
                if (x < 0) {
                    if (n % 2 == 1) {
                        checkRelative(-FastMath.pow(-x, 1.0 / n), build(x).rootN(n));
                    }
                } else {
                    checkRelative(FastMath.pow(x, 1.0 / n), build(x).rootN(n));
                }
            }
        }
    }

    @Test
    public void testPowField() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            for (double y = 0.1; y < 4; y += 0.2) {
                checkRelative(FastMath.pow(x, y), build(x).pow(build(y)));
            }
        }
    }

    @Test
    public void testPowDouble() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            for (double y = 0.1; y < 4; y += 0.2) {
                checkRelative(FastMath.pow(x, y), build(x).pow(y));
            }
            checkRelative(FastMath.pow(x, 0.0), build(x).pow(0.0));
        }
    }

    @Test
    public void testPowInt() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            for (int n = 0; n < 5; ++n) {
                checkRelative(FastMath.pow(x, n), build(x).pow(n));
            }
        }
    }

    @Test
    public void testExp() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.exp(x), build(x).exp());
        }
    }

    @Test
    public void testExpm1() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.expm1(x), build(x).expm1());
        }
    }

    @Test
    public void testLog() {
        for (double x = 0.01; x < 0.9; x += 0.05) {
            checkRelative(FastMath.log(x), build(x).log());
        }
    }

    @Test
    public void testLog1p() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.log1p(x), build(x).log1p());
        }
    }

    @Test
    public void testLog10() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.log10(x), build(x).log10());
        }
    }

    @Test
    public void testScalb() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            for (int n = -100; n < 100; ++n) {
                checkRelative(FastMath.scalb(x, n), build(x).scalb(n));
            }
        }
    }

    @Test
    public void testUlp() {
        final RandomGenerator random = new Well19937a(0x36d4f8862421e0e4l);
        for (int i = -300; i < 300; ++i) {
            final double x = FastMath.scalb(2.0 * random.nextDouble() - 1.0, i);
            Assert.assertTrue(FastMath.ulp(x) >= build(x).ulp().getReal());
        }
    }

    @Test
    public void testCeil() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.ceil(x), build(x).ceil());
        }
    }

    @Test
    public void testFloor() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.floor(x), build(x).floor());
        }
    }

    @Test
    public void testRint() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.rint(x), build(x).rint());
        }
    }

    @Test
    public void testRemainderField() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                checkRelative(FastMath.IEEEremainder(x, y), build(x).remainder(build(y)));
            }
        }
    }

    @Test
    public void testRemainderDouble() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3.2; y < 3.2; y += 0.25) {
                checkRelative(FastMath.IEEEremainder(x, y), build(x).remainder(y));
            }
        }
    }

    @Test
    public void testCopySignField() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                checkRelative(FastMath.copySign(x, y), build(x).copySign(build(y)));
            }
        }
    }

    @Test
    public void testCopySignDouble() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                checkRelative(FastMath.copySign(x, y), build(x).copySign(y));
            }
        }
    }

    @Test
    public void testCopySignSpecialField() {

        Assert.assertEquals(-2.0, build(-2.0).copySign(build(-5.0)).getReal(),                     1.0e-10);
        Assert.assertEquals(-2.0, build(+2.0).copySign(build(-5.0)).getReal(),                     1.0e-10);
        Assert.assertEquals(+2.0, build(-2.0).copySign(build(+5.0)).getReal(),                     1.0e-10);
        Assert.assertEquals(+2.0, build(+2.0).copySign(build(+5.0)).getReal(),                     1.0e-10);
        Assert.assertEquals(-2.0, build(-2.0).copySign(build(Double.NEGATIVE_INFINITY)).getReal(), 1.0e-10);
        Assert.assertEquals(-2.0, build(+2.0).copySign(build(Double.NEGATIVE_INFINITY)).getReal(), 1.0e-10);
        Assert.assertEquals(+2.0, build(-2.0).copySign(build(Double.POSITIVE_INFINITY)).getReal(), 1.0e-10);
        Assert.assertEquals(+2.0, build(+2.0).copySign(build(Double.POSITIVE_INFINITY)).getReal(), 1.0e-10);
        Assert.assertEquals(+2.0, build(-2.0).copySign(build(Double.NaN)).getReal(),               1.0e-10);
        Assert.assertEquals(+2.0, build(-2.0).copySign(build(Double.NaN)).getReal(),               1.0e-10);
        Assert.assertEquals(+2.0, build(-2.0).copySign(build(-Double.NaN)).getReal(),              1.0e-10);
        Assert.assertEquals(+2.0, build(-2.0).copySign(build(-Double.NaN)).getReal(),              1.0e-10);
        Assert.assertEquals(-2.0, build(-2.0).copySign(build(-0.0)).getReal(),                     1.0e-10);
        Assert.assertEquals(-2.0, build(+2.0).copySign(build(-0.0)).getReal(),                     1.0e-10);
        Assert.assertEquals(+2.0, build(-2.0).copySign(build(+0.0)).getReal(),                     1.0e-10);
        Assert.assertEquals(+2.0, build(+2.0).copySign(build(+0.0)).getReal(),                     1.0e-10);

        Assert.assertEquals(-3.0, build(+3.0).copySign(build(-0.0).copySign(build(-5.0))).getReal(),                     1.0e-10);
        Assert.assertEquals(-3.0, build(+3.0).copySign(build(+0.0).copySign(build(-5.0))).getReal(),                     1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(-0.0).copySign(build(+5.0))).getReal(),                     1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(+0.0).copySign(build(+5.0))).getReal(),                     1.0e-10);
        Assert.assertEquals(-3.0, build(+3.0).copySign(build(-0.0).copySign(build(Double.NEGATIVE_INFINITY))).getReal(), 1.0e-10);
        Assert.assertEquals(-3.0, build(+3.0).copySign(build(+0.0).copySign(build(Double.NEGATIVE_INFINITY))).getReal(), 1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(-0.0).copySign(build(Double.POSITIVE_INFINITY))).getReal(), 1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(+0.0).copySign(build(Double.POSITIVE_INFINITY))).getReal(), 1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(-0.0).copySign(build(Double.NaN))).getReal(),               1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(-0.0).copySign(build(Double.NaN))).getReal(),               1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(-0.0).copySign(build(-Double.NaN))).getReal(),              1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(-0.0).copySign(build(-Double.NaN))).getReal(),              1.0e-10);
        Assert.assertEquals(-3.0, build(+3.0).copySign(build(-0.0).copySign(build(-0.0))).getReal(),                     1.0e-10);
        Assert.assertEquals(-3.0, build(+3.0).copySign(build(+0.0).copySign(build(-0.0))).getReal(),                     1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(-0.0).copySign(build(+0.0))).getReal(),                     1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(+0.0).copySign(build(+0.0))).getReal(),                     1.0e-10);

    }

    @Test
    public void testCopySignSpecialDouble() {

        Assert.assertEquals(-2.0, build(-2.0).copySign(-5.0).getReal(),                     1.0e-10);
        Assert.assertEquals(-2.0, build(+2.0).copySign(-5.0).getReal(),                     1.0e-10);
        Assert.assertEquals(+2.0, build(-2.0).copySign(+5.0).getReal(),                     1.0e-10);
        Assert.assertEquals(+2.0, build(+2.0).copySign(+5.0).getReal(),                     1.0e-10);
        Assert.assertEquals(-2.0, build(-2.0).copySign(Double.NEGATIVE_INFINITY).getReal(), 1.0e-10);
        Assert.assertEquals(-2.0, build(+2.0).copySign(Double.NEGATIVE_INFINITY).getReal(), 1.0e-10);
        Assert.assertEquals(+2.0, build(-2.0).copySign(Double.POSITIVE_INFINITY).getReal(), 1.0e-10);
        Assert.assertEquals(+2.0, build(+2.0).copySign(Double.POSITIVE_INFINITY).getReal(), 1.0e-10);
        Assert.assertEquals(+2.0, build(-2.0).copySign(Double.NaN).getReal(),               1.0e-10);
        Assert.assertEquals(+2.0, build(-2.0).copySign(Double.NaN).getReal(),               1.0e-10);
        Assert.assertEquals(+2.0, build(-2.0).copySign(-Double.NaN).getReal(),              1.0e-10);
        Assert.assertEquals(+2.0, build(-2.0).copySign(-Double.NaN).getReal(),              1.0e-10);
        Assert.assertEquals(-2.0, build(-2.0).copySign(-0.0).getReal(),                     1.0e-10);
        Assert.assertEquals(-2.0, build(+2.0).copySign(-0.0).getReal(),                     1.0e-10);
        Assert.assertEquals(+2.0, build(-2.0).copySign(+0.0).getReal(),                     1.0e-10);
        Assert.assertEquals(+2.0, build(+2.0).copySign(+0.0).getReal(),                     1.0e-10);

        Assert.assertEquals(-3.0, build(+3.0).copySign(build(-0.0).copySign(-5.0)).getReal(),                     1.0e-10);
        Assert.assertEquals(-3.0, build(+3.0).copySign(build(+0.0).copySign(-5.0)).getReal(),                     1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(-0.0).copySign(+5.0)).getReal(),                     1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(+0.0).copySign(+5.0)).getReal(),                     1.0e-10);
        Assert.assertEquals(-3.0, build(+3.0).copySign(build(-0.0).copySign(Double.NEGATIVE_INFINITY)).getReal(), 1.0e-10);
        Assert.assertEquals(-3.0, build(+3.0).copySign(build(+0.0).copySign(Double.NEGATIVE_INFINITY)).getReal(), 1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(-0.0).copySign(Double.POSITIVE_INFINITY)).getReal(), 1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(+0.0).copySign(Double.POSITIVE_INFINITY)).getReal(), 1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(-0.0).copySign(Double.NaN)).getReal(),               1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(-0.0).copySign(Double.NaN)).getReal(),               1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(-0.0).copySign(-Double.NaN)).getReal(),              1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(-0.0).copySign(-Double.NaN)).getReal(),              1.0e-10);
        Assert.assertEquals(-3.0, build(+3.0).copySign(build(-0.0).copySign(-0.0)).getReal(),                     1.0e-10);
        Assert.assertEquals(-3.0, build(+3.0).copySign(build(+0.0).copySign(-0.0)).getReal(),                     1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(-0.0).copySign(+0.0)).getReal(),                     1.0e-10);
        Assert.assertEquals(+3.0, build(+3.0).copySign(build(+0.0).copySign(+0.0)).getReal(),                     1.0e-10);

    }

    @Test
    public void testSign() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.signum(x), build(x).sign());
        }
    }

    @Test
    public void testLinearCombinationReference() {
        doTestLinearCombinationReference(x -> build(x), 5.0e-16, 1.0);
    }

    protected void doTestLinearCombinationReference(final DoubleFunction<T> builder,
                                                    final double toleranceLinearCombination,
                                                    final double relativeErrorNaiveImplementation) {

        final T[] a = MathArrays.buildArray(build(0).getField(), 3);
        a[0] = builder.apply(-1321008684645961.0);
        a[1] = builder.apply(-5774608829631843.0);
        a[2] = builder.apply(-7645843051051357.0 / 32.0);
        final T[] b = MathArrays.buildArray(build(0).getField(), 3);
        b[0] = builder.apply(-5712344449280879.0 / 16.0);
        b[1] = builder.apply(-4550117129121957.0 / 16.0);
        b[2] = builder.apply(8846951984510141.0);

        final T abSumInline = a[0].linearCombination(a[0], b[0], a[1], b[1], a[2], b[2]);
        final T abSumArray  = a[0].linearCombination(a, b);
        final T abNaive     = a[0].multiply(b[0]).add(a[1].multiply(b[1])).add(a[2].multiply(b[2]));

        Assert.assertEquals(abSumInline.getReal(), abSumArray.getReal(), 0);
        final double reference = -65271563724949.90625;
        Assert.assertEquals(reference, abSumInline.getReal(),
                            toleranceLinearCombination * FastMath.abs(reference));
        Assert.assertEquals(relativeErrorNaiveImplementation * FastMath.abs(reference),
                            FastMath.abs(abNaive.subtract(reference).getReal()),
                            1.0e-3 * relativeErrorNaiveImplementation * FastMath.abs(reference));

    }

    @Test
    public void testLinearCombinationFaFa() {
        RandomGenerator r = new Well1024a(0xfafal);
        for (int i = 0; i < 50; ++i) {
            double[] aD = generateDouble(r, 10);
            double[] bD = generateDouble(r, 10);
            T[] aF      = toFieldArray(aD);
            T[] bF      = toFieldArray(bD);
            checkRelative(MathArrays.linearCombination(aD, bD),
                          aF[0].linearCombination(aF, bF));
        }
    }

    @Test
    public void testLinearCombinationDaFa() {
        RandomGenerator r = new Well1024a(0xdafal);
        for (int i = 0; i < 50; ++i) {
            double[] aD = generateDouble(r, 10);
            double[] bD = generateDouble(r, 10);
            T[] bF      = toFieldArray(bD);
            checkRelative(MathArrays.linearCombination(aD, bD),
                          bF[0].linearCombination(aD, bF));
        }
    }

    @Test
    public void testLinearCombinationFF2() {
        RandomGenerator r = new Well1024a(0xff2l);
        for (int i = 0; i < 50; ++i) {
            double[] aD = generateDouble(r, 2);
            double[] bD = generateDouble(r, 2);
            T[] aF      = toFieldArray(aD);
            T[] bF      = toFieldArray(bD);
            checkRelative(MathArrays.linearCombination(aD[0], bD[0], aD[1], bD[1]),
                          aF[0].linearCombination(aF[0], bF[0], aF[1], bF[1]));
        }
    }

    @Test
    public void testLinearCombinationDF2() {
        RandomGenerator r = new Well1024a(0xdf2l);
        for (int i = 0; i < 50; ++i) {
            double[] aD = generateDouble(r, 2);
            double[] bD = generateDouble(r, 2);
            T[] bF      = toFieldArray(bD);
            checkRelative(MathArrays.linearCombination(aD[0], bD[0], aD[1], bD[1]),
                          bF[0].linearCombination(aD[0], bF[0], aD[1], bF[1]));
        }
    }

    @Test
    public void testLinearCombinationFF3() {
        RandomGenerator r = new Well1024a(0xff3l);
        for (int i = 0; i < 50; ++i) {
            double[] aD = generateDouble(r, 3);
            double[] bD = generateDouble(r, 3);
            T[] aF      = toFieldArray(aD);
            T[] bF      = toFieldArray(bD);
            checkRelative(MathArrays.linearCombination(aD[0], bD[0], aD[1], bD[1], aD[2], bD[2]),
                          aF[0].linearCombination(aF[0], bF[0], aF[1], bF[1], aF[2], bF[2]));
        }
    }

    @Test
    public void testLinearCombinationDF3() {
        RandomGenerator r = new Well1024a(0xdf3l);
        for (int i = 0; i < 50; ++i) {
            double[] aD = generateDouble(r, 3);
            double[] bD = generateDouble(r, 3);
            T[] bF      = toFieldArray(bD);
            checkRelative(MathArrays.linearCombination(aD[0], bD[0], aD[1], bD[1], aD[2], bD[2]),
                          bF[0].linearCombination(aD[0], bF[0], aD[1], bF[1], aD[2], bF[2]));
        }
    }

    @Test
    public void testLinearCombinationFF4() {
        RandomGenerator r = new Well1024a(0xff4l);
        for (int i = 0; i < 50; ++i) {
            double[] aD = generateDouble(r, 4);
            double[] bD = generateDouble(r, 4);
            T[] aF      = toFieldArray(aD);
            T[] bF      = toFieldArray(bD);
            checkRelative(MathArrays.linearCombination(aD[0], bD[0], aD[1], bD[1], aD[2], bD[2], aD[3], bD[3]),
                          aF[0].linearCombination(aF[0], bF[0], aF[1], bF[1], aF[2], bF[2], aF[3], bF[3]));
        }
    }

    @Test
    public void testLinearCombinationDF4() {
        RandomGenerator r = new Well1024a(0xdf4l);
        for (int i = 0; i < 50; ++i) {
            double[] aD = generateDouble(r, 4);
            double[] bD = generateDouble(r, 4);
            T[] bF      = toFieldArray(bD);
            checkRelative(MathArrays.linearCombination(aD[0], bD[0], aD[1], bD[1], aD[2], bD[2], aD[3], bD[3]),
                          bF[0].linearCombination(aD[0], bF[0], aD[1], bF[1], aD[2], bF[2], aD[3], bF[3]));
        }
    }

    @Test
    public void testGetPi() {
        checkRelative(FastMath.PI, build(-10).getPi());
    }

    @Test
    public void testGetField() {
        checkRelative(1.0, build(-10).getField().getOne());
        checkRelative(0.0, build(-10).getField().getZero());
    }

    @Test
    public void testAbs() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            checkRelative(FastMath.abs(x), build(x).abs());
        }
    }

    @Test
    public void testRound() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            Assert.assertEquals(FastMath.round(x), build(x).round());
        }
    }
    protected void checkRelative(double expected, T obtained) {
        Assert.assertEquals(expected, obtained.getReal(), 1.0e-15 * (1 + FastMath.abs(expected)));
    }

    private double[] generateDouble (final RandomGenerator r, int n) {
        double[] a = new double[n];
        for (int i = 0; i < n; ++i) {
            a[i] = r.nextDouble();
        }
        return a;
    }

    private T[] toFieldArray (double[] a) {
        T[] f = MathArrays.buildArray(build(0).getField(), a.length);
        for (int i = 0; i < a.length; ++i) {
            f[i] = build(a[i]);
        }
        return f;
    }

}
