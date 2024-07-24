/* Copyright 2018 Ulf Adams
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

/*
 * This is not the original file distributed by Ulf Adams in project
 * https://github.com/ulfjack/ryu
 * It has been modified by the Hipparchus project.
 */
package org.hipparchus.util;

import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well19937a;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class RyuDoubleTest {

    private void customAssertD2sEquals(String expected, double f) {
        assertEquals(expected, RyuDouble.doubleToString(f));
    }

    @Test
    void simpleCases() {
        customAssertD2sEquals("0.0", 0);
        customAssertD2sEquals("-0.0", Double.longBitsToDouble(0x8000000000000000L));
        customAssertD2sEquals("1.0", 1.0d);
        customAssertD2sEquals("-1.0", -1.0d);
        customAssertD2sEquals("NaN", Double.NaN);
        customAssertD2sEquals("Infinity", Double.POSITIVE_INFINITY);
        customAssertD2sEquals("-Infinity", Double.NEGATIVE_INFINITY);
    }

    @Test
    void switchToSubnormal() {
        customAssertD2sEquals("2.2250738585072014E-308", Double.longBitsToDouble(0x0010000000000000L));
    }

    /**
     * Floating point values in the range 1.0E-3 <= x < 1.0E7 have to be printed
     * without exponent. This test checks the values at those boundaries.
     */
    @Test
    void boundaryConditions() {
        // x = 1.0E7
        customAssertD2sEquals("1.0E7", 1.0E7d);
        // x < 1.0E7
        customAssertD2sEquals("9999999.999999998", 9999999.999999998d);
        // x = 1.0E-3
        customAssertD2sEquals("0.001", 0.001d);
        // x < 1.0E-3
        customAssertD2sEquals("9.999999999999998E-4", 0.0009999999999999998d);
    }

    @Test
    void test10Power() {
        for (int e = -20; e < -3; ++e) {
            customAssertD2sEquals("1.0E" + e, FastMath.pow(10.0, e));
        }
        customAssertD2sEquals("0.001", FastMath.pow(10.0, -3));
        customAssertD2sEquals("0.01", FastMath.pow(10.0, -2));
        customAssertD2sEquals("0.1", FastMath.pow(10.0, -1));
        customAssertD2sEquals("1.0", FastMath.pow(10.0, 0));
        customAssertD2sEquals("10.0", FastMath.pow(10.0, 1));
        customAssertD2sEquals("100.0", FastMath.pow(10.0, 2));
        customAssertD2sEquals("1000.0", FastMath.pow(10.0, 3));
        customAssertD2sEquals("10000.0", FastMath.pow(10.0, 4));
        customAssertD2sEquals("100000.0", FastMath.pow(10.0, 5));
        customAssertD2sEquals("1000000.0", FastMath.pow(10.0, 6));
        for (int e = 7; e < 20; ++e) {
            customAssertD2sEquals("1.0E" + e, FastMath.pow(10.0, e));
        }
    }

    @Test
    void minAndMax() {
        customAssertD2sEquals("1.7976931348623157E308", Double.longBitsToDouble(0x7fefffffffffffffL));
        customAssertD2sEquals("4.9E-324", Double.longBitsToDouble(1));
    }

    @Test
    void roundingModeEven() {
        customAssertD2sEquals("-2.109808898695963E16", -2.109808898695963E16);
    }

    @Test
    void regressionTest() {
        customAssertD2sEquals("4.940656E-318", 4.940656E-318d);
        customAssertD2sEquals("1.18575755E-316", 1.18575755E-316d);
        customAssertD2sEquals("2.989102097996E-312", 2.989102097996E-312d);
        customAssertD2sEquals("9.0608011534336E15", 9.0608011534336E15d);
        customAssertD2sEquals("4.708356024711512E18", 4.708356024711512E18);
        customAssertD2sEquals("9.409340012568248E18", 9.409340012568248E18);
        // This number naively requires 65 bit for the intermediate results if we reduce the lookup
        // table by half. This checks that we don't loose any information in that case.
        customAssertD2sEquals("1.8531501765868567E21", 1.8531501765868567E21);
        customAssertD2sEquals("-3.347727380279489E33", -3.347727380279489E33);
        // Discovered by Andriy Plokhotnyuk, see #29.
        customAssertD2sEquals("1.9430376160308388E16", 1.9430376160308388E16);
        customAssertD2sEquals("-6.9741824662760956E19", -6.9741824662760956E19);
        customAssertD2sEquals("4.3816050601147837E18", 4.3816050601147837E18);
    }

    @Test
    void testMantissaSize() {
        assertEquals("1.0",                    RyuDouble.doubleToString(1.0,                -20, 20));
        assertEquals("21.0",                   RyuDouble.doubleToString(21.0,               -20, 20));
        assertEquals("321.0",                  RyuDouble.doubleToString(321.0,              -20, 20));
        assertEquals("4321.0",                 RyuDouble.doubleToString(4321.0,             -20, 20));
        assertEquals("54321.0",                RyuDouble.doubleToString(54321.0,            -20, 20));
        assertEquals("654321.0",               RyuDouble.doubleToString(654321.0,           -20, 20));
        assertEquals("7654321.0",              RyuDouble.doubleToString(7654321.0,          -20, 20));
        assertEquals("87654321.0",             RyuDouble.doubleToString(87654321.0,         -20, 20));
        assertEquals("987654321.0",            RyuDouble.doubleToString(987654321.0,        -20, 20));
        assertEquals("1987654321.0",           RyuDouble.doubleToString(1987654321.0,       -20, 20));
        assertEquals("21987654321.0",          RyuDouble.doubleToString(21987654321.0,      -20, 20));
        assertEquals("321987654321.0",         RyuDouble.doubleToString(321987654321.0,     -20, 20));
        assertEquals("4321987654321.0",        RyuDouble.doubleToString(4321987654321.0,    -20, 20));
        assertEquals("54321987654321.0",       RyuDouble.doubleToString(54321987654321.0,   -20, 20));
        assertEquals("654321987654321.0",      RyuDouble.doubleToString(654321987654321.0,  -20, 20));
        assertEquals("7654321987654321.0",     RyuDouble.doubleToString(7654321987654321.0, -20, 20));
    }

    @Test
    void testStandard() {
        RandomGenerator random = new Well19937a(0xca939d6d82eff2d6l);
        for (int i = 0; i < 1000000; ++i) {
            final long   l = random.nextLong();
            final double d = Double.longBitsToDouble(l);
            final String s1 = Double.toString(d);
            final String s2 = RyuDouble.doubleToString(d);
            if (!s1.equals(s2)) {
                // in about 3% cases, standard Double.toString and Ryū give different results
                // in theses cases, Ryū finds either an equal or shorter representation
                // and both representations are information preserving (i.e. they parse to the same number)
                assertTrue(s2.length() <= s1.length());
                assertEquals(l, Double.doubleToRawLongBits(Double.parseDouble(s1)));
                assertEquals(l, Double.doubleToRawLongBits(Double.parseDouble(s2)));
            }
        }
    }

}
