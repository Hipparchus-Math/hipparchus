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
package org.hipparchus.stat.correlation;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.linear.BlockRealMatrix;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Test cases for Kendall's Tau rank correlation.
 */
class KendallsCorrelationTest extends PearsonsCorrelationTest {

    private KendallsCorrelation correlation;

    @BeforeEach
    void setUp() {
        correlation = new KendallsCorrelation();
    }

    /**
     * Test Longley dataset against R.
     */
    @Override
    @Test
    public void testLongly() {
        RealMatrix matrix = createRealMatrix(longleyData, 16, 7);
        KendallsCorrelation corrInstance = new KendallsCorrelation(matrix);
        RealMatrix correlationMatrix = corrInstance.getCorrelationMatrix();
        double[] rData = new double[] {
                1, 0.9166666666666666, 0.9333333333333332, 0.3666666666666666, 0.05, 0.8999999999999999,
                0.8999999999999999, 0.9166666666666666, 1, 0.9833333333333333, 0.45, 0.03333333333333333,
                0.9833333333333333, 0.9833333333333333, 0.9333333333333332, 0.9833333333333333, 1,
                0.4333333333333333, 0.05, 0.9666666666666666, 0.9666666666666666, 0.3666666666666666,
                0.45, 0.4333333333333333, 1, -0.2166666666666666, 0.4666666666666666, 0.4666666666666666, 0.05,
                0.03333333333333333, 0.05, -0.2166666666666666, 1, 0.05, 0.05, 0.8999999999999999, 0.9833333333333333,
                0.9666666666666666, 0.4666666666666666, 0.05, 1, 0.9999999999999999, 0.8999999999999999,
                0.9833333333333333, 0.9666666666666666, 0.4666666666666666, 0.05, 0.9999999999999999, 1
        };
        UnitTestUtils.customAssertEquals("Kendall's correlation matrix", createRealMatrix(rData, 7, 7), correlationMatrix, 10E-15);
    }

    /**
     * Test R swiss fertility dataset.
     */
    @Test
    void testSwiss() {
        RealMatrix matrix = createRealMatrix(swissData, 47, 5);
        KendallsCorrelation corrInstance = new KendallsCorrelation(matrix);
        RealMatrix correlationMatrix = corrInstance.getCorrelationMatrix();
        double[] rData = new double[] {
                1, 0.1795465254708308, -0.4762437404200669, -0.3306111613580587, 0.2453703703703704,
                0.1795465254708308, 1, -0.4505221560842292, -0.4761645631778491, 0.2054604569820847,
                -0.4762437404200669, -0.4505221560842292, 1, 0.528943683925829, -0.3212755391722673,
                -0.3306111613580587, -0.4761645631778491, 0.528943683925829, 1, -0.08479652265379604,
                0.2453703703703704, 0.2054604569820847, -0.3212755391722673, -0.08479652265379604, 1
        };
        UnitTestUtils.customAssertEquals("Kendall's correlation matrix", createRealMatrix(rData, 5, 5), correlationMatrix, 10E-15);
    }

    @Test
    void testSimpleOrdered() {
        final int length = 10;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        for (int i = 0; i < length; i++) {
            xArray[i] = i;
            yArray[i] = i;
        }
        assertEquals(1.0, correlation.correlation(xArray, yArray), Double.MIN_VALUE);
    }

    @Test
    void testSimpleReversed() {
        final int length = 10;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        for (int i = 0; i < length; i++) {
            xArray[length - i - 1] = i;
            yArray[i] = i;
        }
        assertEquals(-1.0, correlation.correlation(xArray, yArray), Double.MIN_VALUE);
    }

    @Test
    void testSimpleOrderedPowerOf2() {
        final int length = 16;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        for (int i = 0; i < length; i++) {
            xArray[i] = i;
            yArray[i] = i;
        }
        assertEquals(1.0, correlation.correlation(xArray, yArray), Double.MIN_VALUE);
    }

    @Test
    void testSimpleReversedPowerOf2() {
        final int length = 16;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        for (int i = 0; i < length; i++) {
            xArray[length - i - 1] = i;
            yArray[i] = i;
        }
        assertEquals(-1.0, correlation.correlation(xArray, yArray), Double.MIN_VALUE);
    }

    @Test
    void testSimpleJumble() {
        //                                     A    B    C    D
        final double[] xArray = new double[] {1.0, 2.0, 3.0, 4.0};
        final double[] yArray = new double[] {1.0, 3.0, 2.0, 4.0};

        // 6 pairs: (A,B) (A,C) (A,D) (B,C) (B,D) (C,D)
        // (B,C) is discordant, the other 5 are concordant

        assertEquals((5 - 1) / (double) 6,
                correlation.correlation(xArray, yArray),
                Double.MIN_VALUE);
    }

    @Test
    void testBalancedJumble() {
        //                                     A    B    C    D
        final double[] xArray = new double[] {1.0, 2.0, 3.0, 4.0};
        final double[] yArray = new double[] {1.0, 4.0, 3.0, 2.0};

        // 6 pairs: (A,B) (A,C) (A,D) (B,C) (B,D) (C,D)
        // (A,B) (A,C), (A,D) are concordant, the other 3 are discordant

        assertEquals(0.0,
                correlation.correlation(xArray, yArray),
                Double.MIN_VALUE);
    }

    @Test
    void testOrderedTies() {
        final int length = 10;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        for (int i = 0; i < length; i++) {
            xArray[i] = i / 2;
            yArray[i] = i / 2;
        }
        // 5 pairs of points that are tied in both values.
        // 16 + 12 + 8 + 4 = 40 concordant
        // (40 - 0) / Math.sqrt((45 - 5) * (45 - 5)) = 1
        assertEquals(1.0, correlation.correlation(xArray, yArray), Double.MIN_VALUE);
    }


    @Test
    void testAllTiesInBoth() {
        final int length = 10;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        assertEquals(Double.NaN, correlation.correlation(xArray, yArray), 0);
    }

    @Test
    void testAllTiesInX() {
        final int length = 10;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        for (int i = 0; i < length; i++) {
            xArray[i] = i;
        }
        assertEquals(Double.NaN, correlation.correlation(xArray, yArray), 0);
    }

    @Test
    void testAllTiesInY() {
        final int length = 10;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        for (int i = 0; i < length; i++) {
            yArray[i] = i;
        }
        assertEquals(Double.NaN, correlation.correlation(xArray, yArray), 0);
    }

    @Test
    void testSingleElement() {
        final int length = 1;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        assertEquals(Double.NaN, correlation.correlation(xArray, yArray), 0);
    }

    @Test
    void testTwoElements() {
        final double[] xArray = new double[] {2.0, 1.0};
        final double[] yArray = new double[] {1.0, 2.0};
        assertEquals(-1.0, correlation.correlation(xArray, yArray), Double.MIN_VALUE);
    }

    @Test
    void test2dDoubleArray() {
        final double[][] input = new double[][] {
                new double[] {2.0, 1.0, 2.0},
                new double[] {1.0, 2.0, 1.0},
                new double[] {0.0, 0.0, 0.0}
        };

        final double[][] expected = new double[][] {
                new double[] {1.0, 1.0 / 3.0, 1.0},
                new double[] {1.0 / 3.0, 1.0, 1.0 / 3.0},
                new double[] {1.0, 1.0 / 3.0, 1.0}};

        assertEquals(correlation.computeCorrelationMatrix(input),
                new BlockRealMatrix(expected));

    }

    @Test
    void testBlockMatrix() {
        final double[][] input = new double[][] {
                new double[] {2.0, 1.0, 2.0},
                new double[] {1.0, 2.0, 1.0},
                new double[] {0.0, 0.0, 0.0}
        };

        final double[][] expected = new double[][] {
                new double[] {1.0, 1.0 / 3.0, 1.0},
                new double[] {1.0 / 3.0, 1.0, 1.0 / 3.0},
                new double[] {1.0, 1.0 / 3.0, 1.0}};

        assertEquals(
                correlation.computeCorrelationMatrix(new BlockRealMatrix(input)),
                new BlockRealMatrix(expected));
    }

    @Test
    void testLargeArray() {
        // test integer overflow detected in MATH-1068
        double[] xArray = new double[100000];
        Arrays.fill(xArray, 0, 2500, 1.0);

        assertEquals(1.0, correlation.correlation(xArray, xArray), 1e-6);
    }

    @Test
    void testMath1277() {
        // example that led to a correlation coefficient outside of [-1, 1]
        // due to a bug reported in MATH-1277
        RandomGenerator rng = new Well1024a(0);
        double[] xArray = new double[120000];
        double[] yArray = new double[120000];
        for (int i = 0; i < xArray.length; ++i) {
            xArray[i] =  rng.nextDouble();
        }
        for (int i = 0; i < yArray.length; ++i) {
            yArray[i] =  rng.nextDouble();
        }
        double coefficient = correlation.correlation(xArray, yArray);
        assertTrue(1.0 >= coefficient && -1.0 <= coefficient);
    }
}
