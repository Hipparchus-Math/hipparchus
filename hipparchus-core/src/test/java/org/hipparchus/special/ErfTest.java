/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 *  contributor license agreements.  See the NOTICE file distributed with
 *  this work for additional information regarding copyright ownership.
 *  The ASF licenses this file to You under the Apache License, Version 2.0
 *  (the "License"); you may not use this file except in compliance with
 *  the License.  You may obtain a copy of the License at
 *
 *     https://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

package org.hipparchus.special;

import org.hipparchus.Field;
import org.hipparchus.UnitTestUtils;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 */
class ErfTest {
    final Field<Binary64> field = Binary64Field.getInstance();
    final Binary64 zero = field.getZero();
    final Binary64 one = field.getOne();

    @Test
    void testErf0() {
        double actual = Erf.erf(0.0);
        double expected = 0.0;
        assertEquals(expected, actual, 1.0e-15);
        assertEquals(1 - expected, Erf.erfc(0.0), 1.0e-15);
    }

    @Test
    void testErf0Field() {
        Binary64 actual   = Erf.erf(zero);
        Binary64 expected = zero;
        assertEquals(zero.getReal(), actual.getReal(), 1.0e-15);
        assertEquals(one.subtract(expected).getReal(), Erf.erfc(zero).getReal(), 1.0e-15);
    }

    @Test
    void testErf1960() {
        double x = 1.960 / FastMath.sqrt(2.0);
        double actual = Erf.erf(x);
        double expected = 0.95;
        assertEquals(expected, actual, 1.0e-5);
        assertEquals(1 - actual, Erf.erfc(x), 1.0e-15);

        actual = Erf.erf(-x);
        expected = -expected;
        assertEquals(expected, actual, 1.0e-5);
        assertEquals(1 - actual, Erf.erfc(-x), 1.0e-15);
    }

    @Test
    void testErf1960Field() {
        Binary64 x = one.multiply(1.960).divide(FastMath.sqrt(2.0));
        Binary64 actual = Erf.erf(x);
        Binary64 expected = one.multiply(0.95);
        assertEquals(expected.getReal(), actual.getReal(), 1.0e-5);
        assertEquals(one.subtract(actual).getReal(), Erf.erfc(x).getReal(), 1.0e-15);

        actual = Erf.erf(x.negate());
        expected = expected.negate();
        assertEquals(expected.getReal(), actual.getReal(), 1.0e-5);
        assertEquals(one.subtract(actual).getReal(), Erf.erfc(x.negate()).getReal(), 1.0e-15);
    }

    @Test
    void testErf2576() {
        double x = 2.576 / FastMath.sqrt(2.0);
        double actual = Erf.erf(x);
        double expected = 0.99;
        assertEquals(expected, actual, 1.0e-5);
        assertEquals(1 - actual, Erf.erfc(x), 1e-15);

        actual = Erf.erf(-x);
        expected = -expected;
        assertEquals(expected, actual, 1.0e-5);
        assertEquals(1 - actual, Erf.erfc(-x), 1.0e-15);
    }

    @Test
    void testErf2576Field() {
        Binary64 x = one.multiply(2.576).divide(FastMath.sqrt(2.0));
        Binary64 actual = Erf.erf(x);
        Binary64 expected = one.multiply(0.99);
        assertEquals(expected.getReal(), actual.getReal(), 1.0e-5);
        assertEquals(one.subtract(actual).getReal(), Erf.erfc(x).getReal(), 1.0e-15);

        actual = Erf.erf(x.negate());
        expected = expected.negate();
        assertEquals(expected.getReal(), actual.getReal(), 1.0e-5);
        assertEquals(one.subtract(actual).getReal(), Erf.erfc(x.negate()).getReal(), 1.0e-15);
    }

    @Test
    void testErf2807() {
        double x = 2.807 / FastMath.sqrt(2.0);
        double actual = Erf.erf(x);
        double expected = 0.995;
        assertEquals(expected, actual, 1.0e-5);
        assertEquals(1 - actual, Erf.erfc(x), 1.0e-15);

        actual = Erf.erf(-x);
        expected = -expected;
        assertEquals(expected, actual, 1.0e-5);
        assertEquals(1 - actual, Erf.erfc(-x), 1.0e-15);
    }

    @Test
    void testErf2807Field() {
        Binary64 x = one.multiply(2.807).divide(FastMath.sqrt(2.0));
        Binary64 actual = Erf.erf(x);
        Binary64 expected = one.multiply(0.995);
        assertEquals(expected.getReal(), actual.getReal(), 1.0e-5);
        assertEquals(one.subtract(actual).getReal(), Erf.erfc(x).getReal(), 1.0e-15);

        actual = Erf.erf(x.negate());
        expected = expected.negate();
        assertEquals(expected.getReal(), actual.getReal(), 1.0e-5);
        assertEquals(one.subtract(actual).getReal(), Erf.erfc(x.negate()).getReal(), 1.0e-15);
    }

    @Test
    void testErf3291() {
        double x = 3.291 / FastMath.sqrt(2.0);
        double actual = Erf.erf(x);
        double expected = 0.999;
        assertEquals(expected, actual, 1.0e-5);
        assertEquals(1 - expected, Erf.erfc(x), 1.0e-5);

        actual = Erf.erf(-x);
        expected = -expected;
        assertEquals(expected, actual, 1.0e-5);
        assertEquals(1 - expected, Erf.erfc(-x), 1.0e-5);
    }

    @Test
    void testErf3291Field() {
        Binary64 x = one.multiply(3.291).divide(FastMath.sqrt(2.0));
        Binary64 actual = Erf.erf(x);
        Binary64 expected = one.multiply(0.999);
        assertEquals(expected.getReal(), actual.getReal(), 1.0e-5);
        assertEquals(one.subtract(actual).getReal(), Erf.erfc(x).getReal(), 1.0e-15);

        actual = Erf.erf(x.negate());
        expected = expected.negate();
        assertEquals(expected.getReal(), actual.getReal(), 1.0e-5);
        assertEquals(one.subtract(actual).getReal(), Erf.erfc(x.negate()).getReal(), 1.0e-15);
    }

    /**
     * MATH-301, MATH-456
     */
    @Test
    void testLargeValues() {
        for (int i = 1; i < 200; i*=10) {
            double result = Erf.erf(i);
            assertFalse(Double.isNaN(result));
            assertTrue(result > 0 && result <= 1);
            result = Erf.erf(-i);
            assertFalse(Double.isNaN(result));
            assertTrue(result >= -1 && result < 0);
            result = Erf.erfc(i);
            assertFalse(Double.isNaN(result));
            assertTrue(result >= 0 && result < 1);
            result = Erf.erfc(-i);
            assertFalse(Double.isNaN(result));
            assertTrue(result >= 1 && result <= 2);
        }
        assertEquals(-1, Erf.erf(Double.NEGATIVE_INFINITY), 0);
        assertEquals(1, Erf.erf(Double.POSITIVE_INFINITY), 0);
        assertEquals(2, Erf.erfc(Double.NEGATIVE_INFINITY), 0);
        assertEquals(0, Erf.erfc(Double.POSITIVE_INFINITY), 0);
    }

    /**
     * MATH-301, MATH-456
     */
    @Test
    void testLargeValuesField() {
        for (int i = 1; i < 200; i*=10) {
            final Binary64 iField = new Binary64(i);
            Binary64 result = Erf.erf(iField);
            assertFalse(result.isNaN());
            assertTrue(result.getReal() > 0 && result.getReal() <= 1);
            result = Erf.erf(iField.negate());
            assertFalse(result.isNaN());
            assertTrue(result.getReal() >= -1 && result.getReal() < 0);
            result = Erf.erfc(iField);
            assertFalse(result.isNaN());
            assertTrue(result.getReal() >= 0 && result.getReal() < 1);
            result = Erf.erfc(iField.negate());
            assertFalse(result.isNaN());
            assertTrue(result.getReal() >= 1 && result.getReal() <= 2);
        }
        assertEquals(one.negate().getReal(), Erf.erf(new Binary64(Double.NEGATIVE_INFINITY)).getReal(), 0);
        assertEquals(one.getReal(), Erf.erf(new Binary64(Double.POSITIVE_INFINITY)).getReal(), 0);
        assertEquals(one.multiply(2).getReal(), new Binary64(Erf.erfc(Double.NEGATIVE_INFINITY)).getReal(), 0);
        assertEquals(zero.getReal(), Erf.erfc(new Binary64(Double.POSITIVE_INFINITY)).getReal(), 0);
    }

    /**
     * Compare Erf.erf against reference values computed using GCC 4.2.1 (Apple OSX packaged version)
     * erfl (extended precision erf).
     */
    @Test
    void testErfGnu() {
        final double tol = 1E-15;
        final double[] gnuValues = new double[] {-1, -1, -1, -1, -1,
        -1, -1, -1, -0.99999999999999997848,
        -0.99999999999999264217, -0.99999999999846254017, -0.99999999980338395581, -0.99999998458274209971,
        -0.9999992569016276586, -0.99997790950300141459, -0.99959304798255504108, -0.99532226501895273415,
        -0.96610514647531072711, -0.84270079294971486948, -0.52049987781304653809,  0,
         0.52049987781304653809, 0.84270079294971486948, 0.96610514647531072711, 0.99532226501895273415,
         0.99959304798255504108, 0.99997790950300141459, 0.9999992569016276586, 0.99999998458274209971,
         0.99999999980338395581, 0.99999999999846254017, 0.99999999999999264217, 0.99999999999999997848,
         1,  1,  1,  1,
         1,  1,  1,  1};
        double x = -10d;
        for (int i = 0; i < 41; i++) {
            assertEquals(gnuValues[i], Erf.erf(x), tol);
            x += 0.5d;
        }
    }

    /**
     * Compare Erf.erf against reference values computed using GCC 4.2.1 (Apple OSX packaged version)
     * erfl (extended precision erf).
     */
    @Test
    void testErfGnuField() {
        final double tol = 1E-15;
        final Binary64[] gnuValues = new Binary64[] {one.negate(), one.negate(), one.negate(), one.negate(), one.negate(),
        one.negate(), one.negate(), one.negate(), new Binary64(-0.99999999999999997848),
        new Binary64(-0.99999999999999264217), new Binary64(-0.99999999999846254017), new Binary64(-0.99999999980338395581), new Binary64(-0.99999998458274209971),
        new Binary64(-0.9999992569016276586), new Binary64(-0.99997790950300141459), new Binary64(-0.99959304798255504108), new Binary64(-0.99532226501895273415),
        new Binary64(-0.96610514647531072711), new Binary64(-0.84270079294971486948), new Binary64(-0.52049987781304653809),  zero,
         new Binary64(0.52049987781304653809), new Binary64(0.84270079294971486948), new Binary64(0.96610514647531072711), new Binary64(0.99532226501895273415),
         new Binary64(0.99959304798255504108), new Binary64(0.99997790950300141459), new Binary64(0.9999992569016276586), new Binary64(0.99999998458274209971),
         new Binary64(0.99999999980338395581), new Binary64(0.99999999999846254017), new Binary64(0.99999999999999264217), new Binary64(0.99999999999999997848),
         one,  one,  one,  one,
         one,  one,  one,  one};
        Binary64 x = one.multiply(-10d);
        for (int i = 0; i < 41; i++) {
            assertEquals(gnuValues[i].getReal(), Erf.erf(x).getReal(), tol);
            x = x.add(0.5d);
        }
    }

    /**
     * Compare Erf.erfc against reference values computed using GCC 4.2.1 (Apple OSX packaged version)
     * erfcl (extended precision erfc).
     */
    @Test
    void testErfcGnu() {
        final double tol = 1E-15;
        final double[] gnuValues = new double[] { 2,  2,  2,  2,  2,
        2,  2,  2, 1.9999999999999999785,
        1.9999999999999926422, 1.9999999999984625402, 1.9999999998033839558, 1.9999999845827420998,
        1.9999992569016276586, 1.9999779095030014146, 1.9995930479825550411, 1.9953222650189527342,
        1.9661051464753107271, 1.8427007929497148695, 1.5204998778130465381,  1,
        0.47950012218695346194, 0.15729920705028513051, 0.033894853524689272893, 0.0046777349810472658333,
        0.00040695201744495893941, 2.2090496998585441366E-05, 7.4309837234141274516E-07, 1.5417257900280018858E-08,
        1.966160441542887477E-10, 1.5374597944280348501E-12, 7.3578479179743980661E-15, 2.1519736712498913103E-17,
        3.8421483271206474691E-20, 4.1838256077794144006E-23, 2.7766493860305691016E-26, 1.1224297172982927079E-29,
        2.7623240713337714448E-33, 4.1370317465138102353E-37, 3.7692144856548799402E-41, 2.0884875837625447567E-45};
        double x = -10d;
        for (int i = 0; i < 41; i++) {
            assertEquals(gnuValues[i], Erf.erfc(x), tol);
            x += 0.5d;
        }
    }

    /**
     * Compare Erf.erfc against reference values computed using GCC 4.2.1 (Apple OSX packaged version)
     * erfcl (extended precision erfc).
     */
    @Test
    void testErfcGnuField() {
        final double tol = 1E-15;
        final Binary64[] gnuValues = new Binary64[] { new Binary64(2),  new Binary64(2),  new Binary64(2),  new Binary64(2),  new Binary64(2),
        new Binary64(2),  new Binary64(2),  new Binary64(2), new Binary64(1.9999999999999999785),
        new Binary64(1.9999999999999926422), new Binary64(1.9999999999984625402), new Binary64(1.9999999998033839558), new Binary64(1.9999999845827420998),
        new Binary64(1.9999992569016276586), new Binary64(1.9999779095030014146), new Binary64(1.9995930479825550411), new Binary64(1.9953222650189527342),
        new Binary64(1.9661051464753107271), new Binary64(1.8427007929497148695), new Binary64(1.5204998778130465381),  one,
        new Binary64(0.47950012218695346194), new Binary64(0.15729920705028513051), new Binary64(0.033894853524689272893), new Binary64(0.0046777349810472658333),
        new Binary64(0.00040695201744495893941), new Binary64(2.2090496998585441366E-05), new Binary64(7.4309837234141274516E-07), new Binary64(1.5417257900280018858E-08),
        new Binary64(1.966160441542887477E-10), new Binary64(1.5374597944280348501E-12), new Binary64(7.3578479179743980661E-15), new Binary64(2.1519736712498913103E-17),
        new Binary64(3.8421483271206474691E-20), new Binary64(4.1838256077794144006E-23), new Binary64(2.7766493860305691016E-26), new Binary64(1.1224297172982927079E-29),
        new Binary64(2.7623240713337714448E-33), new Binary64(4.1370317465138102353E-37), new Binary64(3.7692144856548799402E-41),new Binary64( 2.0884875837625447567E-45)};
        Binary64 x = new Binary64(-10d);
        for (int i = 0; i < 41; i++) {
            assertEquals(gnuValues[i].getReal(), Erf.erfc(x).getReal(), tol);
            x = x.add(0.5d);
        }
    }

    /**
     * Tests erfc against reference data computed using Maple reported in Marsaglia, G,,
     * "Evaluating the Normal Distribution," Journal of Statistical Software, July, 2004.
     * http//www.jstatsoft.org/v11/a05/paper
     */
    @Test
    void testErfcMaple() {
        double[][] ref = new double[][]
                        {{0.1, 4.60172162722971e-01},
                         {1.2, 1.15069670221708e-01},
                         {2.3, 1.07241100216758e-02},
                         {3.4, 3.36929265676881e-04},
                         {4.5, 3.39767312473006e-06},
                         {5.6, 1.07175902583109e-08},
                         {6.7, 1.04209769879652e-11},
                         {7.8, 3.09535877195870e-15},
                         {8.9, 2.79233437493966e-19},
                         {10.0, 7.61985302416053e-24},
                         {11.1, 6.27219439321703e-29},
                         {12.2, 1.55411978638959e-34},
                         {13.3, 1.15734162836904e-40},
                         {14.4, 2.58717592540226e-47},
                         {15.5, 1.73446079179387e-54},
                         {16.6, 3.48454651995041e-62}
        };
        for (int i = 0; i < 15; i++) {
            final double result = 0.5*Erf.erfc(ref[i][0]/FastMath.sqrt(2));
            assertEquals(ref[i][1], result, 1E-15);
            UnitTestUtils.customAssertRelativelyEquals(ref[i][1], result, 1E-13);
        }
    }

    /**
     * Tests erfc against reference data computed using Maple reported in Marsaglia, G,,
     * "Evaluating the Normal Distribution," Journal of Statistical Software, July, 2004.
     * http//www.jstatsoft.org/v11/a05/paper
     */
    @Test
    void testErfcMapleField() {
        Binary64[][] ref = new Binary64[][]
                {{new Binary64(0.1), new Binary64(4.60172162722971e-01)},
                 {new Binary64(1.2), new Binary64(1.15069670221708e-01)},
                 {new Binary64(2.3), new Binary64(1.07241100216758e-02)},
                 {new Binary64(3.4), new Binary64(3.36929265676881e-04)},
                 {new Binary64(4.5), new Binary64(3.39767312473006e-06)},
                 {new Binary64(5.6), new Binary64(1.07175902583109e-08)},
                 {new Binary64(6.7), new Binary64(1.04209769879652e-11)},
                 {new Binary64(7.8), new Binary64(3.09535877195870e-15)},
                 {new Binary64(8.9), new Binary64(2.79233437493966e-19)},
                 {new Binary64(10.0), new Binary64(7.61985302416053e-24)},
                 {new Binary64(11.1), new Binary64(6.27219439321703e-29)},
                 {new Binary64(12.2), new Binary64(1.55411978638959e-34)},
                 {new Binary64(13.3), new Binary64(1.15734162836904e-40)},
                 {new Binary64(14.4), new Binary64(2.58717592540226e-47)},
                 {new Binary64(15.5), new Binary64(1.73446079179387e-54)},
                 {new Binary64(16.6), new Binary64(3.48454651995041e-62)}
                };
        for (int i = 0; i < 15; i++) {
            final Binary64 result = Erf.erfc(ref[i][0].divide(FastMath.sqrt(2))).multiply(0.5);
            assertEquals(ref[i][1].getReal(), result.getReal(), 1E-15);
            UnitTestUtils.customAssertRelativelyEquals(ref[i][1].getReal(), result.getReal(), 1E-13);
        }
    }

    /**
     * Test the implementation of Erf.erf(double, double) for consistency with results
     * obtained from Erf.erf(double) and Erf.erfc(double).
     */
    @Test
    void testTwoArgumentErf() {
        double[] xi = new double[]{-2.0, -1.0, -0.9, -0.1, 0.0, 0.1, 0.9, 1.0, 2.0};
        for(double x1 : xi) {
            for(double x2 : xi) {
                double a = Erf.erf(x1, x2);
                double b = Erf.erf(x2) - Erf.erf(x1);
                double c = Erf.erfc(x1) - Erf.erfc(x2);
                assertEquals(a, b, 1E-15);
                assertEquals(a, c, 1E-15);
            }
        }
    }

    /**
     * Test the implementation of Erf.erf(double, double) for consistency with results
     * obtained from Erf.erf(double) and Erf.erfc(double).
     */
    @Test
    void testTwoArgumentErfField() {
        Binary64[] xi = new Binary64[]{new Binary64(-2.0), new Binary64(-1.0), new Binary64(-0.9), new Binary64(-0.1), new Binary64(0.0), new Binary64(0.1), new Binary64(0.9), new Binary64(1.0), new Binary64(2.0)};
        for(Binary64 x1 : xi) {
            for(Binary64 x2 : xi) {
                Binary64 a = Erf.erf(x1, x2);
                Binary64 b = Erf.erf(x2).subtract(Erf.erf(x1));
                Binary64 c = Erf.erfc(x1).subtract(Erf.erfc(x2));
                assertEquals(a.getReal(), b.getReal(), 1E-15);
                assertEquals(a.getReal(), c.getReal(), 1E-15);
            }
        }
    }

    @Test
    void testErfInvNaN() {
        assertTrue(Double.isNaN(Erf.erfInv(-1.001)));
        assertTrue(Double.isNaN(Erf.erfInv(+1.001)));
    }

    @Test
    void testErfInvNaNField() {
        assertTrue((Erf.erfInv(new Binary64(-1.001))).isNaN());
        assertTrue(Erf.erfInv(new Binary64(+1.001)).isNaN());
    }

    @Test
    void testErfInvInfinite() {
        assertTrue(Double.isInfinite(Erf.erfInv(-1)));
        assertTrue(Erf.erfInv(-1) < 0);
        assertTrue(Double.isInfinite(Erf.erfInv(+1)));
        assertTrue(Erf.erfInv(+1) > 0);
    }

    @Test
    void testErfInvInfiniteField() {
        assertTrue(Double.isInfinite(Erf.erfInv(-1)));
        assertTrue(Erf.erfInv(-1) < 0);
        assertTrue(Double.isInfinite(Erf.erfInv(+1)));
        assertTrue(Erf.erfInv(+1) > 0);
    }

    @Test
    void testErfInv() {
        for (double x = -5.9; x < 5.9; x += 0.01) {
            final double y = Erf.erf(x);
            final double dydx = 2 * FastMath.exp(-x * x) / FastMath.sqrt(FastMath.PI);
            assertEquals(x, Erf.erfInv(y), 1.0e-15 / dydx);
        }
    }

    @Test
    void testErfInvField() {
        for (Binary64 x = new Binary64(-5.9); x.getReal() < 5.9; x = x.add(0.01)) {
            final Binary64 y = Erf.erf(x);
            final Binary64 dydx = x.square().negate().exp().multiply(2/FastMath.sqrt(FastMath.PI));
            assertEquals(x.getReal(), Erf.erfInv(y).getReal(), 1.0e-15 / dydx.getReal());
        }
    }

    @Test
    void testErfcInvNaN() {
        assertTrue(Double.isNaN(Erf.erfcInv(-0.001)));
        assertTrue(Double.isNaN(Erf.erfcInv(+2.001)));
    }

    @Test
    void testErfcInvNaNField() {
        assertTrue(Erf.erfcInv(new Binary64(-0.001)).isNaN());
        assertTrue(Erf.erfcInv(new Binary64(+2.001)).isNaN());
    }

    @Test
    void testErfcInvInfinite() {
        assertTrue(Double.isInfinite(Erf.erfcInv(-0)));
        assertTrue(Erf.erfcInv( 0) > 0);
        assertTrue(Double.isInfinite(Erf.erfcInv(+2)));
        assertTrue(Erf.erfcInv(+2) < 0);
    }

    @Test
    void testErfcInvInfiniteField() {
        assertTrue(Erf.erfcInv(new Binary64(-0)).isInfinite());
        assertTrue(Erf.erfcInv( zero).getReal() > 0);
        assertTrue(Erf.erfcInv(new Binary64(+2)).isInfinite());
        assertTrue(Erf.erfcInv(new Binary64(+2)).getReal() < 0);
    }

    @Test
    void testErfcInv() {
        for (double x = -5.85; x < 5.9; x += 0.01) {
            final double y = Erf.erfc(x);
            final double dydxAbs = 2 * FastMath.exp(-x * x) / FastMath.sqrt(FastMath.PI);
            assertEquals(x, Erf.erfcInv(y), 1.0e-15 / dydxAbs);
        }
    }

    @Test
    void testErfcInvField() {
        for (Binary64 x = new Binary64(-5.85); x.getReal() < 5.9; x = x.add(0.01)) {
            final Binary64 y = Erf.erfc(x);
            final Binary64 dydxAbs = x.square().negate().exp().multiply(2/FastMath.sqrt(FastMath.PI));
            assertEquals(x.getReal(), Erf.erfcInv(y).getReal(), 1.0e-15 / dydxAbs.getReal());
        }
    }
}
