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
import org.hipparchus.util.MathUtils;
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

    @Test
    public void testAbramowitzStegunExample4() {
        Assert.assertEquals(1.019106060,
                            LegendreEllipticIntegral.bigE(FastMath.sqrt(80.0 / 81.0)),
                            2.0e-8);
    }

    @Test
    public void testAbramowitzStegunExample8() {
        final double k    = FastMath.sqrt(1.0 / 5.0);
        final double phi1 = FastMath.acos(FastMath.sqrt(2) / 3.0);
        final double phi2 = FastMath.acos(FastMath.sqrt(2) / 2.0);
        Assert.assertEquals(1.115921, LegendreEllipticIntegral.bigF(phi1, k), 1.0e-6);
        Assert.assertEquals(0.800380, LegendreEllipticIntegral.bigF(phi2, k), 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample9() {
        final double k    = FastMath.sqrt(1.0 / 2.0);
        final double phi1 = MathUtils.SEMI_PI;
        final double phi2 = FastMath.PI / 6.0;
        Assert.assertEquals(1.854075, LegendreEllipticIntegral.bigF(phi1, k), 1.0e-6);
        Assert.assertEquals(0.535623, LegendreEllipticIntegral.bigF(phi2, k), 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample10() {
        final double k    = FastMath.sqrt(4.0 / 5.0);
        final double phi  = FastMath.PI / 6.0;
        Assert.assertEquals(0.543604, LegendreEllipticIntegral.bigF(phi, k), 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample14() {
        final double k    = 3.0 / 5.0;
        final double phi1 = FastMath.asin(FastMath.sqrt(5.0) / 3.0);
        final double phi2 = FastMath.asin(5.0 / (3.0 * FastMath.sqrt(17.0)));
        Assert.assertEquals(0.80904, LegendreEllipticIntegral.bigE(phi1, k), 1.0e-5);
        Assert.assertEquals(0.41192, LegendreEllipticIntegral.bigE(phi2, k), 1.0e-5);
    }

    @Test
    public void testAbramowitzStegunTable175() {
        Assert.assertEquals(0.26263487,
                            LegendreEllipticIntegral.bigF(FastMath.toRadians(15),
                                                          FastMath.sin(FastMath.toRadians(32))),
                            1.0e-8);
        Assert.assertEquals(1.61923762,
                            LegendreEllipticIntegral.bigF(FastMath.toRadians(80),
                                                          FastMath.sin(FastMath.toRadians(46))),
                            1.0e-8);
    }

    @Test
    public void testAbramowitzStegunTable176() {
        Assert.assertEquals(0.42531712,
                            LegendreEllipticIntegral.bigE(FastMath.toRadians(25),
                                                          FastMath.sin(FastMath.toRadians(64))),
                            1.0e-8);
        Assert.assertEquals(0.96208074,
                            LegendreEllipticIntegral.bigE(FastMath.toRadians(70),
                                                          FastMath.sin(FastMath.toRadians(76))),
                            1.0e-8);
    }

    @Test
    public void testAbramowitzStegunTable179() {
        Assert.assertEquals(1.62298,
                            LegendreEllipticIntegral.bigPi(FastMath.toRadians(75),
                                                           0.4,
                                                           FastMath.sin(FastMath.toRadians(15))),
                            1.0e-5);
        Assert.assertEquals(1.03076,
                            LegendreEllipticIntegral.bigPi(FastMath.toRadians(45),
                                                           0.8,
                                                           FastMath.sin(FastMath.toRadians(60))),
                            1.0e-5);
        Assert.assertEquals(2.79990,
                            LegendreEllipticIntegral.bigPi(FastMath.toRadians(75),
                                                           0.9,
                                                           FastMath.sin(FastMath.toRadians(15))),
                            1.0e-5);
    }

    @Test
    public void testCompleteVsIncompleteF() {
        for (double k = 0.01; k < 1; k += 0.01) {
            double complete   = LegendreEllipticIntegral.bigK(k);
            double incomplete = LegendreEllipticIntegral.bigF(MathUtils.SEMI_PI, k);
            Assert.assertEquals(complete, incomplete, FastMath.ulp(complete));
        }
    }

    @Test
    public void testCompleteVsIncompleteE() {
        for (double k = 0.01; k < 1; k += 0.01) {
            double complete   = LegendreEllipticIntegral.bigE(k);
            double incomplete = LegendreEllipticIntegral.bigE(MathUtils.SEMI_PI, k);
            Assert.assertEquals(complete, incomplete, 3 * FastMath.ulp(complete));
        }
    }

    @Test
    public void testCompleteVsIncompleteD() {
        for (double k = 0.01; k < 1; k += 0.01) {
            double complete   = LegendreEllipticIntegral.bigD(k);
            double incomplete = LegendreEllipticIntegral.bigD(MathUtils.SEMI_PI, k);
            Assert.assertEquals(complete, incomplete, FastMath.ulp(complete));
        }
    }

    @Test
    public void testCompleteVsIncompletePi() {
        for (double alpha2 = 0.01; alpha2 < 1; alpha2 += 0.01) {
            for (double k = 0.01; k < 1; k += 0.01) {
                double complete   = LegendreEllipticIntegral.bigPi(alpha2, k);
                double incomplete = LegendreEllipticIntegral.bigPi(MathUtils.SEMI_PI, alpha2, k);
                Assert.assertEquals(complete, incomplete, FastMath.ulp(complete));
            }
        }
    }

}
