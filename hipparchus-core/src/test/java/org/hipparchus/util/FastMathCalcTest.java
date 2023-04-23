/*
 * Licensed to the Hipparchus project under one or more
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
package org.hipparchus.util;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;

public class FastMathCalcTest {

    @Test
    public void testExpIntTables() {

        final double[] fmTableA   = getD1("ExpIntTable", "EXP_INT_TABLE_A");
        final double[] fmTableB   = getD1("ExpIntTable", "EXP_INT_TABLE_B");
        final int      len        = getInt("EXP_INT_TABLE_LEN");
        final int      max        = getInt("EXP_INT_TABLE_MAX_INDEX");
        Assert.assertEquals(len, fmTableA.length);
        Assert.assertEquals(len, fmTableB.length);

        final double[] tmp   = new double[2];
        final double[] recip = new double[2];
        for (int i = 0; i < max; i++) {
            FastMathCalc.expint(i, tmp);
            if (i == 0) {
                Assert.assertEquals(fmTableA[max], tmp[0], FastMath.ulp(fmTableA[i]));
                Assert.assertEquals(fmTableB[max], tmp[1], FastMath.ulp(fmTableB[i]));
            } else {
                FastMathCalc.splitReciprocal(tmp, recip);
                Assert.assertEquals(fmTableA[max - i], recip[0], FastMath.ulp(fmTableA[i]));
                Assert.assertEquals(fmTableB[max - i], recip[1], FastMath.ulp(fmTableB[i]));
            }
        }

    }

    @Test
    public void testExpFracTables() {

        final double[] fmTableA   = getD1("ExpFracTable", "EXP_FRAC_TABLE_A");
        final double[] fmTableB   = getD1("ExpFracTable", "EXP_FRAC_TABLE_B");
        final int      len        = getInt("EXP_FRAC_TABLE_LEN");
        Assert.assertEquals(len, fmTableA.length);
        Assert.assertEquals(len, fmTableB.length);

        final double factor = 1d / (len - 1);
        final double[] tmp = new double[2];
        for (int i = 0; i < len; i++) {
            FastMathCalc.slowexp(i * factor, tmp);
            Assert.assertEquals(fmTableA[i], tmp[0], FastMath.ulp(fmTableA[i]));
            Assert.assertEquals(fmTableB[i], tmp[1], FastMath.ulp(fmTableB[i]));
        }

    }

    @Test
    public void testLnMantTables() {
        final double[][] fmTable  = getD2("lnMant", "LN_MANT");
        final int      len        = getInt("LN_MANT_LEN");
        Assert.assertEquals(len, fmTable.length);

        for (int i = 0; i < len; i++) {
            final double d = Double.longBitsToDouble( (((long) i) << 42) | 0x3ff0000000000000L );
            final double[] tmp = FastMathCalc.slowLog(d);
            Assert.assertEquals(fmTable[i].length, tmp.length);
            for (int j = 0; j < fmTable[i].length; ++j) {
                Assert.assertEquals(fmTable[i][j], tmp[j], FastMath.ulp(fmTable[i][j]));
            }
        }

    }

    @Test
    public void testSplit() {
        checkSplit(0x3ffe0045dab7321fl, 0x3ffe0045c0000000l, 0x3e7ab7321f000000l);
        checkSplit(0x3ffe0045fab7321fl, 0x3ffe004600000000l, 0xbe55233784000000l);
        checkSplit(0x7dfedcba9876543fl, 0x7dfedcba80000000l, 0x7c7876543f000000l);
        checkSplit(0x7dfedcbaf876543fl, 0x7dfedcbb00000000l, 0xfc5e26af04000000l);
        checkSplit(0xfdfedcba9876543fl, 0xfdfedcba80000000l, 0xfc7876543f000000l);
        checkSplit(0xfdfedcbaf876543fl, 0xfdfedcbb00000000l, 0x7c5e26af04000000l);
    }

    private void checkSplit(final long bits, final long high, final long low) {
        try {
            Method split = FastMathCalc.class.getDeclaredMethod("split", Double.TYPE, double[].class);
            split.setAccessible(true);
            double   d      = Double.longBitsToDouble(bits);
            double[] result = new double[2];
            split.invoke(null, d, result);
            Assert.assertEquals(bits, Double.doubleToRawLongBits(result[0] + result[1]));
            Assert.assertEquals(high, Double.doubleToRawLongBits(result[0]));
            Assert.assertEquals(low,  Double.doubleToRawLongBits(result[1]));

        } catch (NoSuchMethodException | SecurityException | IllegalArgumentException |
                 IllegalAccessException | InvocationTargetException e) {
            Assert.fail(e.getLocalizedMessage());
        }
    }

    @Test
    public void testSinCosTanTables() {
        try {
            final double[] sinA = getFastMathTable("SINE_TABLE_A");
            final double[] sinB = getFastMathTable("SINE_TABLE_B");
            final double[] cosA = getFastMathTable("COSINE_TABLE_A");
            final double[] cosB = getFastMathTable("COSINE_TABLE_B");
            final double[] tanA = getFastMathTable("TANGENT_TABLE_A");
            final double[] tanB = getFastMathTable("TANGENT_TABLE_B");
            Method buildSinCosTables = FastMathCalc.class.getDeclaredMethod("buildSinCosTables",
                                                                            double[].class, double[].class, double[].class, double[].class,
                                                                            Integer.TYPE,
                                                                            double[].class, double[].class);
            buildSinCosTables.setAccessible(true);
            final double[] calcSinA = new double[sinA.length];
            final double[] calcSinB = new double[sinB.length];
            final double[] calcCosA = new double[cosA.length];
            final double[] calcCosB = new double[cosB.length];
            final double[] calcTanA = new double[tanA.length];
            final double[] calcTanB = new double[tanB.length];
            buildSinCosTables.invoke(null, calcSinA, calcSinB, calcCosA, calcCosB, sinA.length, calcTanA, calcTanB);
            checkTable(sinA, calcSinA, 0);
            checkTable(sinB, calcSinB, 0);
            checkTable(cosA, calcCosA, 0);
            checkTable(cosB, calcCosB, 0);
            checkTable(tanA, calcTanA, 0);
            checkTable(tanB, calcTanB, 0);

        } catch (NoSuchMethodException | SecurityException | IllegalArgumentException |
                 IllegalAccessException | InvocationTargetException e) {
            Assert.fail(e.getLocalizedMessage());
        }
    }

    private double[] getFastMathTable(final String name) {
        try {
            final Field field = FastMath.class.getDeclaredField(name);
            field.setAccessible(true);
            return (double[]) field.get(null);
        } catch (NoSuchFieldException | SecurityException |
                 IllegalArgumentException | IllegalAccessException e) {
            Assert.fail(e.getLocalizedMessage());
            return null;
        }
    }

    private void checkTable(final double[] reference, final double[] actual, int maxUlps) {
        Assert.assertEquals(reference.length, actual.length);
        for (int i = 0; i < reference.length; ++i) {
            Assert.assertTrue(Precision.equals(reference[i], actual[i], maxUlps));
        }
    }

    @Test
    public void testPrintArray1() {
        try (ByteArrayOutputStream bos = new ByteArrayOutputStream();
             PrintStream ps = new PrintStream(bos, true, StandardCharsets.UTF_8.name())) {
            Method printArray = FastMathCalc.class.getDeclaredMethod("printarray", PrintStream.class,
                                                                     String.class, Integer.TYPE, double[].class);
            printArray.setAccessible(true);
            printArray.invoke(null, ps, "name", 2, new double[] { 1.25, -0.5 });
            Assert.assertEquals(String.format("name=%n" +
                                              "    {%n" +
                                              "        +1.25d,%n" +
                                              "        -0.5d,%n" +
                                              "    };%n"),
                                bos.toString(StandardCharsets.UTF_8.name()));
        } catch (IOException | NoSuchMethodException | IllegalAccessException |
                 IllegalArgumentException | InvocationTargetException e) {
            e.printStackTrace();
            Assert.fail(e.getLocalizedMessage());
        }
    }

    @Test
    public void testPrintArray2() {
        try (ByteArrayOutputStream bos = new ByteArrayOutputStream();
             PrintStream ps = new PrintStream(bos, true, StandardCharsets.UTF_8.name())) {
            Method printArray = FastMathCalc.class.getDeclaredMethod("printarray", PrintStream.class,
                                                                     String.class, Integer.TYPE, double[][].class);
            printArray.setAccessible(true);
            printArray.invoke(null, ps, "name", 2, new double[][] { { 1.25, -0.5 }, { 0.0, 3.0 } });
            Assert.assertEquals(String.format("name%n" + 
                                              "    { %n" + 
                                              "        {+1.25d,                  -0.5d,                   }, // 0%n" +
                                              "        {+0.0d,                   +3.0d,                   }, // 1%n" +
                                              "    };%n"),
                                bos.toString(StandardCharsets.UTF_8.name()));
        } catch (IOException | NoSuchMethodException | IllegalAccessException |
                 IllegalArgumentException | InvocationTargetException e) {
            e.printStackTrace();
            Assert.fail(e.getLocalizedMessage());
        }
    }

    private double[] getD1(final String innerClassName, final String tableName) {
        try {

            final Class<?> inerClass = Arrays.stream(FastMath.class.getDeclaredClasses()).
                                                     filter(c -> c.getName().endsWith("$" + innerClassName)).
                                                     findFirst().
                                                     get();
            final Field fmTableField = inerClass.getDeclaredField(tableName);
            fmTableField.setAccessible(true);
            return (double[]) fmTableField.get(null);

        } catch (NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            Assert.fail(e.getLocalizedMessage());
            return null;
        }
    }

    private double[][] getD2(final String innerClassName, final String tableName) {
        try {

            final Class<?> inerClass = Arrays.stream(FastMath.class.getDeclaredClasses()).
                                                     filter(c -> c.getName().endsWith("$" + innerClassName)).
                                                     findFirst().
                                                     get();
            final Field fmTableField = inerClass.getDeclaredField(tableName);
            fmTableField.setAccessible(true);
            return (double[][]) fmTableField.get(null);

        } catch (NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            Assert.fail(e.getLocalizedMessage());
            return null;
        }
    }

    private int getInt(String lenName) {
        try {

            final Field fmLen = FastMath.class.getDeclaredField(lenName);
            fmLen.setAccessible(true);
            return ((Integer) fmLen.get(null)).intValue();

        } catch (NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            Assert.fail(e.getLocalizedMessage());
            return -1;
        }
    }

}
