/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.util;

import java.lang.reflect.Field;
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
