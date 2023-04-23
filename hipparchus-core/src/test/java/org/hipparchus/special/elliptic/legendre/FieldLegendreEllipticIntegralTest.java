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

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.special.elliptic.carlson.CarlsonEllipticIntegral;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.junit.Assert;
import org.junit.Test;

public class FieldLegendreEllipticIntegralTest {

    @Test
    public void testNoConvergence() {
        doTestNoConvergence(Binary64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestNoConvergence(final Field<T> field) {
        Assert.assertTrue(LegendreEllipticIntegral.bigK(field.getZero().newInstance(Double.NaN)).isNaN());
    }

    @Test
    public void testComplementary() {
        doTestComplementary(Binary64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestComplementary(final Field<T> field) {
        for (double m = 0.01; m < 1; m += 0.01) {
            T k1 = LegendreEllipticIntegral.bigK(field.getZero().newInstance(m));
            T k2 = LegendreEllipticIntegral.bigKPrime(field.getZero().newInstance(1 - m));
            Assert.assertEquals(k1.getReal(), k2.getReal(), FastMath.ulp(k1).getReal());
        }
    }

    @Test
    public void testAbramowitzStegunExample3() {
        doTestAbramowitzStegunExample3(Binary64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestAbramowitzStegunExample3(final Field<T> field) {
        T k = LegendreEllipticIntegral.bigK(field.getZero().newInstance(80.0 / 81.0));
        Assert.assertEquals(3.591545001, k.getReal(), 2.0e-9);
    }

    public void testAbramowitzStegunExample4() {
        doTestBigE(Binary64Field.getInstance(), 80.0 / 81.0, 1.019106060, 2.0e-8);
    }

    @Test
    public void testAbramowitzStegunExample8() {
        final double m    = 1.0 / 5.0;
        doTestBigF(Binary64Field.getInstance(), FastMath.acos(FastMath.sqrt(2) / 3.0), m, 1.115921, 1.0e-6);
        doTestBigF(Binary64Field.getInstance(), FastMath.acos(FastMath.sqrt(2) / 2.0), m, 0.800380, 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample9() {
        final double m    = 1.0 / 2.0;
        doTestBigF(Binary64Field.getInstance(), MathUtils.SEMI_PI, m, 1.854075, 1.0e-6);
        doTestBigF(Binary64Field.getInstance(), FastMath.PI / 6.0, m, 0.535623, 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample10() {
        final double m    = 4.0 / 5.0;
        doTestBigF(Binary64Field.getInstance(), FastMath.PI / 6.0, m, 0.543604, 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample14() {
        final double k    = 3.0 / 5.0;
        doTestBigE(Binary64Field.getInstance(), FastMath.asin(FastMath.sqrt(5.0) / 3.0),          k * k, 0.80904, 1.0e-5);
        doTestBigE(Binary64Field.getInstance(), FastMath.asin(5.0 / (3.0 * FastMath.sqrt(17.0))), k * k, 0.41192, 1.0e-5);
    }

    @Test
    public void testAbramowitzStegunTable175() {
        final double sinAlpha1 = FastMath.sin(FastMath.toRadians(32));
        doTestBigF(Binary64Field.getInstance(), FastMath.toRadians(15), sinAlpha1 * sinAlpha1, 0.26263487, 1.0e-8);
        final double sinAlpha2 = FastMath.sin(FastMath.toRadians(46));
        doTestBigF(Binary64Field.getInstance(), FastMath.toRadians(80), sinAlpha2 * sinAlpha2, 1.61923762, 1.0e-8);
    }

    @Test
    public void testAbramowitzStegunTable176() {
        final double sinAlpha1 = FastMath.sin(FastMath.toRadians(64));
        doTestBigE(Binary64Field.getInstance(), FastMath.toRadians(25), sinAlpha1 * sinAlpha1, 0.42531712, 1.0e-8);
        final double sinAlpha2 = FastMath.sin(FastMath.toRadians(76));
        doTestBigE(Binary64Field.getInstance(), FastMath.toRadians(70), sinAlpha2 * sinAlpha2, 0.96208074, 1.0e-8);
    }

    @Test
    public void testAbramowitzStegunTable179() {
        final double sinAlpha1 = FastMath.sin(FastMath.toRadians(15));
        doTestBigPi(Binary64Field.getInstance(), FastMath.toRadians(75), 0.4, sinAlpha1 * sinAlpha1, 1.62298, 1.0e-5);
        final double sinAlpha2 = FastMath.sin(FastMath.toRadians(60));
        doTestBigPi(Binary64Field.getInstance(), FastMath.toRadians(45), 0.8, sinAlpha2 * sinAlpha2, 1.03076, 1.0e-5);
        final double sinAlpha3 = FastMath.sin(FastMath.toRadians(15));
        doTestBigPi(Binary64Field.getInstance(), FastMath.toRadians(75), 0.9, sinAlpha3 * sinAlpha3, 2.79990, 1.0e-5);
    }

    @Test
    public void testCompleteVsIncompleteF() {
        doTestCompleteVsIncompleteF(Binary64Field.getInstance());
    }

    @Test
    public void testCompleteVsIncompleteE() {
        doTestCompleteVsIncompleteE(Binary64Field.getInstance());
    }

    @Test
    public void testCompleteVsIncompleteD() {
        doTestCompleteVsIncompleteD(Binary64Field.getInstance());
    }

    @Test
    public void testCompleteVsIncompletePi() {
        doTestCompleteVsIncompletePi(Binary64Field.getInstance());
    }

    @Test
    public void testNomeMediumParameter() {
        doTestNomeMediumParameter(Binary64Field.getInstance());
    }

    @Test
    public void testNomeSmallParameter() {
        doTestNomeSmallParameter(Binary64Field.getInstance());
    }

    @Test
    public void testPrecomputedDelta() {
        doTestPrecomputedDelta(Binary64Field.getInstance());
    }

    @Test
    public void testIntegralsSmallParameter() {
        doTestIntegralsSmallParameter(Binary64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestBigE(final Field<T> field, final double m,
                                                                final double expected, final double tol) {
        Assert.assertEquals(expected,
                            LegendreEllipticIntegral.bigE(field.getZero().newInstance(m)).getReal(),
                            tol);
    }

    private <T extends CalculusFieldElement<T>> void doTestBigE(final Field<T> field,
                                                                final double phi, final double m,
                                                                final double expected, final double tol) {
        Assert.assertEquals(expected,
                            LegendreEllipticIntegral.bigE(field.getZero().newInstance(phi),
                                                          field.getZero().newInstance(m)).getReal(),
                            tol);
    }

    private <T extends CalculusFieldElement<T>> void doTestBigF(final Field<T> field,
                                                                final double phi, final double m,
                                                                final double expected, final double tol) {
        Assert.assertEquals(expected,
                            LegendreEllipticIntegral.bigF(field.getZero().newInstance(phi),
                                                          field.getZero().newInstance(m)).getReal(),
                            tol);
    }

    private <T extends CalculusFieldElement<T>> void doTestBigPi(final Field<T> field,
                                                                 final double phi, final double alpha2, final double m,
                                                                 final double expected, final double tol) {
        Assert.assertEquals(expected,
                            LegendreEllipticIntegral.bigPi(field.getZero().newInstance(alpha2),
                                                           field.getZero().newInstance(phi),
                                                           field.getZero().newInstance(m)).getReal(),
                            tol);
    }

    private <T extends CalculusFieldElement<T>> void doTestCompleteVsIncompleteF(final Field<T> field) {
        for (double m = 0.01; m < 1; m += 0.01) {
            double complete   = LegendreEllipticIntegral.bigK(field.getZero().newInstance(m)).getReal();
            double incomplete = LegendreEllipticIntegral.bigF(field.getZero().newInstance(MathUtils.SEMI_PI),
                                                              field.getZero().newInstance(m)).getReal();
            Assert.assertEquals(complete, incomplete, FastMath.ulp(complete));
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestCompleteVsIncompleteE(final Field<T> field) {
        for (double m = 0.01; m < 1; m += 0.01) {
            double complete   = LegendreEllipticIntegral.bigE(field.getZero().newInstance(m)).getReal();
            double incomplete = LegendreEllipticIntegral.bigE(field.getZero().newInstance(MathUtils.SEMI_PI),
                                                              field.getZero().newInstance(m)).getReal();
            Assert.assertEquals(complete, incomplete, 4 * FastMath.ulp(complete));
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestCompleteVsIncompleteD(final Field<T> field) {
        for (double m = 0.01; m < 1; m += 0.01) {
            double complete   = LegendreEllipticIntegral.bigD(field.getZero().newInstance(m)).getReal();
            double incomplete = LegendreEllipticIntegral.bigD(field.getZero().newInstance(MathUtils.SEMI_PI),
                                                              field.getZero().newInstance(m)).getReal();
            Assert.assertEquals(complete, incomplete, FastMath.ulp(complete));
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestCompleteVsIncompletePi(final Field<T> field) {
        for (double alpha2 = 0.01; alpha2 < 1; alpha2 += 0.01) {
            for (double m = 0.01; m < 1; m += 0.01) {
                double complete   = LegendreEllipticIntegral.bigPi(field.getZero().newInstance(alpha2),
                                                                   field.getZero().newInstance(m)).getReal();
                double incomplete = LegendreEllipticIntegral.bigPi(field.getZero().newInstance(alpha2),
                                                                   field.getZero().newInstance(MathUtils.SEMI_PI),
                                                                   field.getZero().newInstance(m)).getReal();
                Assert.assertEquals(complete, incomplete, FastMath.ulp(complete));
            }
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestNomeMediumParameter(final Field<T> field) {
        Assert.assertEquals(0.0857957337021947665168, LegendreEllipticIntegral.nome(field.getZero().newInstance(0.75)).getReal(), 1.0e-15);
    }

    private <T extends CalculusFieldElement<T>> void doTestNomeSmallParameter(final Field<T> field) {
        Assert.assertEquals(5.9375e-18, LegendreEllipticIntegral.nome(field.getZero().newInstance(0.95e-16)).getReal(), 1.0e-22);
    }

    private <T extends CalculusFieldElement<T>> void doTestIntegralsSmallParameter(final Field<T> field) {
        Assert.assertEquals(7.8539816428e-10,
                            LegendreEllipticIntegral.bigK(field.getZero().newInstance(2.0e-9)).getReal() - MathUtils.SEMI_PI,
                            1.0e-15);
    }

    private <T extends CalculusFieldElement<T>> void doTestPrecomputedDelta(final Field<T> field) {

        T n   = field.getZero().newInstance(0.7);
        T m   = field.getZero().newInstance(0.2);
        T phi = field.getZero().newInstance(1.2);
        T ref = field.getZero().newInstance(1.8264362537906997);
        Assert.assertEquals(0.0, LegendreEllipticIntegral.bigPi(n, phi, m).subtract(ref).getReal(), 1.0e-15);

        // no argument reduction and no precomputed delta
        final T csc     = phi.sin().reciprocal();
        final T csc2    = csc.multiply(csc);
        final T cM1     = csc2.subtract(1);
        final T cMm     = csc2.subtract(m);
        final T cMn     = csc2.subtract(n);
        final T pinphim = CarlsonEllipticIntegral.rF(cM1, cMm, csc2).
                          add(CarlsonEllipticIntegral.rJ(cM1, cMm, csc2, cMn).multiply(n).divide(3));
        Assert.assertEquals(0.0, pinphim.subtract(ref).getReal(), 1.0e-15);

    }

}
