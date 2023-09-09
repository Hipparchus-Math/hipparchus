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
package org.hipparchus.special.elliptic.legendre;

import org.hipparchus.special.elliptic.carlson.CarlsonEllipticIntegral;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.junit.Assert;
import org.junit.Test;

public class LegendreEllipticIntegralTest {

    @Test
    public void testNoConvergence() {
        Assert.assertTrue(Double.isNaN(LegendreEllipticIntegral.bigK(Double.NaN)));
    }

    @Test
    public void testComplementary() {
        for (double m = 0.01; m < 1; m += 0.01) {
            double k1 = LegendreEllipticIntegral.bigK(m);
            double k2 = LegendreEllipticIntegral.bigKPrime(1 - m);
            Assert.assertEquals(k1, k2, FastMath.ulp(k1));
        }
    }

    @Test
    public void testAbramowitzStegunExample3() {
        Assert.assertEquals(3.591545001,
                            LegendreEllipticIntegral.bigK(80.0 / 81.0),
                            2.0e-9);
    }

    @Test
    public void testAbramowitzStegunExample4() {
        Assert.assertEquals(1.019106060,
                            LegendreEllipticIntegral.bigE(80.0 / 81.0),
                            2.0e-8);
    }

    @Test
    public void testAbramowitzStegunExample8() {
        final double m    = 1.0 / 5.0;
        final double phi1 = FastMath.acos(FastMath.sqrt(2) / 3.0);
        final double phi2 = FastMath.acos(FastMath.sqrt(2) / 2.0);
        Assert.assertEquals(1.115921, LegendreEllipticIntegral.bigF(phi1, m), 1.0e-6);
        Assert.assertEquals(0.800380, LegendreEllipticIntegral.bigF(phi2, m), 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample9() {
        final double m    = 1.0 / 2.0;
        final double phi1 = MathUtils.SEMI_PI;
        final double phi2 = FastMath.PI / 6.0;
        Assert.assertEquals(1.854075, LegendreEllipticIntegral.bigF(phi1, m), 1.0e-6);
        Assert.assertEquals(0.535623, LegendreEllipticIntegral.bigF(phi2, m), 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample10() {
        final double m    = 4.0 / 5.0;
        final double phi  = FastMath.PI / 6.0;
        Assert.assertEquals(0.543604, LegendreEllipticIntegral.bigF(phi, m), 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample14() {
        final double k    = 3.0 / 5.0;
        final double phi1 = FastMath.asin(FastMath.sqrt(5.0) / 3.0);
        final double phi2 = FastMath.asin(5.0 / (3.0 * FastMath.sqrt(17.0)));
        Assert.assertEquals(0.80904, LegendreEllipticIntegral.bigE(phi1, k * k), 1.0e-5);
        Assert.assertEquals(0.41192, LegendreEllipticIntegral.bigE(phi2, k * k), 1.0e-5);
    }

    @Test
    public void testAbramowitzStegunTable175() {
        final double sinAlpha1 = FastMath.sin(FastMath.toRadians(32));
        Assert.assertEquals(0.26263487,
                            LegendreEllipticIntegral.bigF(FastMath.toRadians(15), sinAlpha1 * sinAlpha1),
                            1.0e-8);
        final double sinAlpha2 = FastMath.sin(FastMath.toRadians(46));
        Assert.assertEquals(1.61923762,
                            LegendreEllipticIntegral.bigF(FastMath.toRadians(80), sinAlpha2 * sinAlpha2),
                            1.0e-8);
    }

    @Test
    public void testAbramowitzStegunTable176() {
        final double sinAlpha1 = FastMath.sin(FastMath.toRadians(64));
        Assert.assertEquals(0.42531712,
                            LegendreEllipticIntegral.bigE(FastMath.toRadians(25), sinAlpha1 * sinAlpha1),
                            1.0e-8);
        final double sinAlpha2 = FastMath.sin(FastMath.toRadians(76));
        Assert.assertEquals(0.96208074,
                            LegendreEllipticIntegral.bigE(FastMath.toRadians(70), sinAlpha2 * sinAlpha2),
                            1.0e-8);
    }

    @Test
    public void testAbramowitzStegunTable179() {
        final double sinAlpha1 = FastMath.sin(FastMath.toRadians(15));
        Assert.assertEquals(1.62298,
                            LegendreEllipticIntegral.bigPi(0.4, FastMath.toRadians(75), sinAlpha1 * sinAlpha1),
                            1.0e-5);
        final double sinAlpha2 = FastMath.sin(FastMath.toRadians(60));
        Assert.assertEquals(1.03076,
                            LegendreEllipticIntegral.bigPi(0.8, FastMath.toRadians(45), sinAlpha2 * sinAlpha2),
                            1.0e-5);
        final double sinAlpha3 = FastMath.sin(FastMath.toRadians(15));
        Assert.assertEquals(2.79990,
                            LegendreEllipticIntegral.bigPi(0.9, FastMath.toRadians(75), sinAlpha3 * sinAlpha3),
                            1.0e-5);
    }

    @Test
    public void testCompleteVsIncompleteF() {
        for (double m = 0.01; m < 1; m += 0.01) {
            double complete   = LegendreEllipticIntegral.bigK(m);
            double incomplete = LegendreEllipticIntegral.bigF(MathUtils.SEMI_PI, m);
            Assert.assertEquals(complete, incomplete, FastMath.ulp(complete));
        }
    }

    @Test
    public void testCompleteVsIncompleteE() {
        for (double m = 0.01; m < 1; m += 0.01) {
            double complete   = LegendreEllipticIntegral.bigE(m);
            double incomplete = LegendreEllipticIntegral.bigE(MathUtils.SEMI_PI, m);
            Assert.assertEquals(complete, incomplete, 4 * FastMath.ulp(complete));
        }
    }

    @Test
    public void testCompleteVsIncompleteD() {
        for (double m = 0.01; m < 1; m += 0.01) {
            double complete   = LegendreEllipticIntegral.bigD(m);
            double incomplete = LegendreEllipticIntegral.bigD(MathUtils.SEMI_PI, m);
            Assert.assertEquals(complete, incomplete, FastMath.ulp(complete));
        }
    }

    @Test
    public void testCompleteVsIncompletePi() {
        for (double alpha2 = 0.01; alpha2 < 1; alpha2 += 0.01) {
            for (double m = 0.01; m < 1; m += 0.01) {
                double complete   = LegendreEllipticIntegral.bigPi(alpha2, m);
                double incomplete = LegendreEllipticIntegral.bigPi(alpha2, MathUtils.SEMI_PI, m);
                Assert.assertEquals(complete, incomplete, FastMath.ulp(complete));
            }
        }
    }

    @Test
    public void testNomeMediumParameter() {
        Assert.assertEquals(0.0857957337021947665168, LegendreEllipticIntegral.nome(0.75), 1.0e-15);
    }

    @Test
    public void testNomeSmallParameter() {
        Assert.assertEquals(5.9375e-18, LegendreEllipticIntegral.nome(0.95e-16), 1.0e-22);
    }

    @Test
    public void testIntegralsSmallParameter() {
        Assert.assertEquals(7.8539816428e-10,
                            LegendreEllipticIntegral.bigK(2.0e-9) - MathUtils.SEMI_PI,
                            1.0e-15);
    }

    @Test
    public void testPrecomputedDelta() {

        double n   = 0.7;
        double m   = 0.2;
        double phi = 1.2;
        double ref = 1.8264362537906997;
        Assert.assertEquals(ref, LegendreEllipticIntegral.bigPi(n, phi, m), 1.0e-15);

        // no argument reduction and no precomputed delta
        final double csc     = 1.0 / FastMath.sin(phi);
        final double csc2    = csc * csc;
        final double cM1     = csc2 - 1;
        final double cMm     = csc2 - m;
        final double cMn     = csc2 - n;
        final double pinphim = CarlsonEllipticIntegral.rF(cM1, cMm, csc2) +
                               CarlsonEllipticIntegral.rJ(cM1, cMm, csc2, cMn) * n / 3;
        Assert.assertEquals(ref, pinphim, 1.0e-15);

    }

}
