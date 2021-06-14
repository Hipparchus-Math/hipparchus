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

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.Decimal64Field;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.junit.Assert;
import org.junit.Test;

public class FieldLegendreEllipticIntegralTest {

    @Test
    public void testNoConvergence() {
        doTestNoConvergence(Decimal64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestNoConvergence(final Field<T> field) {
        try {
            LegendreEllipticIntegral.bigK(field.getZero().newInstance(Double.NaN));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedCoreFormats.CONVERGENCE_FAILED, mise.getSpecifier());
        }
    }

    @Test
    public void testComplementary() {
        doTestComplementary(Decimal64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestComplementary(final Field<T> field) {
        for (double k = 0.01; k < 1; k += 0.01) {
            T k1 = LegendreEllipticIntegral.bigK(field.getZero().newInstance(k));
            T k2 = LegendreEllipticIntegral.bigKPrime(field.getZero().newInstance(FastMath.sqrt(1 - k * k)));
            Assert.assertEquals(k1.getReal(), k2.getReal(), FastMath.ulp(k1).getReal());
        }
    }

    @Test
    public void testAbramowitzStegunExample3() {
        doTestAbramowitzStegunExample3(Decimal64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestAbramowitzStegunExample3(final Field<T> field) {
        T k = LegendreEllipticIntegral.bigK(field.getZero().newInstance(FastMath.sqrt(80.0 / 81.0)));
        Assert.assertEquals(3.591545001, k.getReal(), 2.0e-9);
    }

    public void testAbramowitzStegunExample4() {
        doTestBigE(Decimal64Field.getInstance(), FastMath.sqrt(80.0 / 81.0), 1.019106060, 2.0e-8);
    }

    @Test
    public void testAbramowitzStegunExample8() {
        final double k    = FastMath.sqrt(1.0 / 5.0);
        doTestBigF(Decimal64Field.getInstance(), FastMath.acos(FastMath.sqrt(2) / 3.0), k, 1.115921, 1.0e-6);
        doTestBigF(Decimal64Field.getInstance(), FastMath.acos(FastMath.sqrt(2) / 2.0), k, 0.800380, 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample9() {
        final double k    = FastMath.sqrt(1.0 / 2.0);
        doTestBigF(Decimal64Field.getInstance(), MathUtils.SEMI_PI, k, 1.854075, 1.0e-6);
        doTestBigF(Decimal64Field.getInstance(), FastMath.PI / 6.0, k, 0.535623, 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample10() {
        final double k    = FastMath.sqrt(4.0 / 5.0);
        doTestBigF(Decimal64Field.getInstance(), FastMath.PI / 6.0, k, 0.543604, 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample14() {
        final double k    = 3.0 / 5.0;
        doTestBigE(Decimal64Field.getInstance(), FastMath.asin(FastMath.sqrt(5.0) / 3.0),          k, 0.80904, 1.0e-5);
        doTestBigE(Decimal64Field.getInstance(), FastMath.asin(5.0 / (3.0 * FastMath.sqrt(17.0))), k, 0.41192, 1.0e-5);
    }

    @Test
    public void testAbramowitzStegunTable175() {
        doTestBigF(Decimal64Field.getInstance(), FastMath.toRadians(15), FastMath.sin(FastMath.toRadians(32)), 0.26263487, 1.0e-8);
        doTestBigF(Decimal64Field.getInstance(), FastMath.toRadians(80), FastMath.sin(FastMath.toRadians(46)), 1.61923762, 1.0e-8);
    }

    @Test
    public void testAbramowitzStegunTable176() {
        doTestBigE(Decimal64Field.getInstance(), FastMath.toRadians(25), FastMath.sin(FastMath.toRadians(64)), 0.42531712, 1.0e-8);
        doTestBigE(Decimal64Field.getInstance(), FastMath.toRadians(70), FastMath.sin(FastMath.toRadians(76)), 0.96208074, 1.0e-8);
    }

    @Test
    public void testAbramowitzStegunTable179() {
        doTestBigPi(Decimal64Field.getInstance(), FastMath.toRadians(75), 0.4, FastMath.sin(FastMath.toRadians(15)), 1.62298, 1.0e-5);
        doTestBigPi(Decimal64Field.getInstance(), FastMath.toRadians(45), 0.8, FastMath.sin(FastMath.toRadians(60)), 1.03076, 1.0e-5);
        doTestBigPi(Decimal64Field.getInstance(), FastMath.toRadians(75), 0.9, FastMath.sin(FastMath.toRadians(15)), 2.79990, 1.0e-5);
    }

    @Test
    public void testCompleteVsIncompleteF() {
        doTestCompleteVsIncompleteF(Decimal64Field.getInstance());
    }

    @Test
    public void testCompleteVsIncompleteE() {
        doTestCompleteVsIncompleteE(Decimal64Field.getInstance());
    }

    @Test
    public void testCompleteVsIncompleteD() {
        doTestCompleteVsIncompleteD(Decimal64Field.getInstance());
    }

    @Test
    public void testCompleteVsIncompletePi() {
        doTestCompleteVsIncompletePi(Decimal64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestBigE(final Field<T> field, final double k,
                                                                final double expected, final double tol) {
        Assert.assertEquals(expected,
                            LegendreEllipticIntegral.bigE(field.getZero().newInstance(k)).getReal(),
                            tol);
    }

    private <T extends CalculusFieldElement<T>> void doTestBigE(final Field<T> field,
                                                                final double phi, final double k,
                                                                final double expected, final double tol) {
        Assert.assertEquals(expected,
                            LegendreEllipticIntegral.bigE(field.getZero().newInstance(phi),
                                                          field.getZero().newInstance(k)).getReal(),
                            tol);
    }

    private <T extends CalculusFieldElement<T>> void doTestBigF(final Field<T> field,
                                                                final double phi, final double k,
                                                                final double expected, final double tol) {
        Assert.assertEquals(expected,
                            LegendreEllipticIntegral.bigF(field.getZero().newInstance(phi),
                                                          field.getZero().newInstance(k)).getReal(),
                            tol);
    }

    private <T extends CalculusFieldElement<T>> void doTestBigPi(final Field<T> field,
                                                                 final double phi, final double alpha2, final double k,
                                                                 final double expected, final double tol) {
        Assert.assertEquals(expected,
                            LegendreEllipticIntegral.bigPi(field.getZero().newInstance(phi),
                                                           field.getZero().newInstance(alpha2),
                                                           field.getZero().newInstance(k)).getReal(),
                            tol);
    }

    private <T extends CalculusFieldElement<T>> void doTestCompleteVsIncompleteF(final Field<T> field) {
        for (double k = 0.01; k < 1; k += 0.01) {
            double complete   = LegendreEllipticIntegral.bigK(field.getZero().newInstance(k)).getReal();
            double incomplete = LegendreEllipticIntegral.bigF(field.getZero().newInstance(MathUtils.SEMI_PI),
                                                              field.getZero().newInstance(k)).getReal();
            Assert.assertEquals(complete, incomplete, FastMath.ulp(complete));
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestCompleteVsIncompleteE(final Field<T> field) {
        for (double k = 0.01; k < 1; k += 0.01) {
            double complete   = LegendreEllipticIntegral.bigE(field.getZero().newInstance(k)).getReal();
            double incomplete = LegendreEllipticIntegral.bigE(field.getZero().newInstance(MathUtils.SEMI_PI),
                                                              field.getZero().newInstance(k)).getReal();
            Assert.assertEquals(complete, incomplete, 3 * FastMath.ulp(complete));
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestCompleteVsIncompleteD(final Field<T> field) {
        for (double k = 0.01; k < 1; k += 0.01) {
            double complete   = LegendreEllipticIntegral.bigD(field.getZero().newInstance(k)).getReal();
            double incomplete = LegendreEllipticIntegral.bigD(field.getZero().newInstance(MathUtils.SEMI_PI),
                                                              field.getZero().newInstance(k)).getReal();
            Assert.assertEquals(complete, incomplete, FastMath.ulp(complete));
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestCompleteVsIncompletePi(final Field<T> field) {
        for (double alpha2 = 0.01; alpha2 < 1; alpha2 += 0.01) {
            for (double k = 0.01; k < 1; k += 0.01) {
                double complete   = LegendreEllipticIntegral.bigPi(field.getZero().newInstance(alpha2),
                                                                   field.getZero().newInstance(k)).getReal();
                double incomplete = LegendreEllipticIntegral.bigPi(field.getZero().newInstance(MathUtils.SEMI_PI),
                                                                   field.getZero().newInstance(alpha2),
                                                                   field.getZero().newInstance(k)).getReal();
                Assert.assertEquals(complete, incomplete, FastMath.ulp(complete));
            }
        }
    }

}
