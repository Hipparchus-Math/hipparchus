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
package org.hipparchus.stat.inference;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the MannWhitneyUTest class.
 */

class MannWhitneyUTestTest {

    protected MannWhitneyUTest testStatistic = new MannWhitneyUTest();

    /**
     * Target values for most tests below were computed using R, version 3.4.4.
     * The wilcox.test function was used throughout, passing (along with x and y)
     *   alternative = "two.sized"
     *   mu = 0
     *   paired = FALSE
     * For the exact tests,
     *   exact = TRUE
     * For normal approximation (exact = false in MannWhitneyTest API)
     *   exact = FALSE
     *   correct = TRUE (tells R to do the continuity correction).
     * The target for the U statistic value was obtained by running
     * the R test twice, the second time switching the order of x and y and
     * taking the minimum value. This is correct because the statistic that R
     * reports as W is the Wilcoxon Rank Sum with x before y and U is the
     * minimum of this value and the rank sum with y before x.
     */

    @Test
    void testMannWhitneyUSimple() {
        /*
         * Target values computed using R version 3.4.4.
         */
        final double[] x = {
            19, 22, 16, 29, 24
        };
        final double[] y = {
            20, 11, 17, 12
        };

        assertEquals(3, testStatistic.mannWhitneyU(x, y), 1e-10);
        assertEquals(0.11134688653,
                            testStatistic.mannWhitneyUTest(x, y, false), 1e-5);
        assertEquals(0.11111111111,
                            testStatistic.mannWhitneyUTest(x, y, true), 1e-5);
    }

    @Test
    void testDisjoint() {
        /*
         * Target values computed using R version 3.4.4.
         */
        final double[] x = {
            1, 2, 3, 4, 5
        };
        final double[] y = {
            6, 7, 8, 9, 10, 11
        };
        assertEquals(0, testStatistic.mannWhitneyU(x, y), 1e-10);
        assertEquals(0.0081131172656,
                            testStatistic.mannWhitneyUTest(x, y, false), 1e-5);
        assertEquals(0.004329004329,
                            testStatistic.mannWhitneyUTest(x, y, true), 1e-5);

    }

    @Test
    void testMannWhitneyUInputValidation() {
        /*
         * Samples must be present, i.e. length > 0
         */
        try {
            testStatistic.mannWhitneyUTest(new double[] {}, new double[] {
                1.0
            });
            fail("x does not contain samples (exact), MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            testStatistic.mannWhitneyUTest(new double[] {
                1.0
            }, new double[] {});
            fail("y does not contain samples (exact), MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        /*
         * x and y is null
         */
        try {
            testStatistic.mannWhitneyUTest(null, null);
            fail("x and y is null (exact), NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }

        try {
            testStatistic.mannWhitneyUTest(null, null);
            fail("x and y is null (asymptotic), NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }

        /*
         * x or y is null
         */
        try {
            testStatistic.mannWhitneyUTest(null, new double[] {
                1.0
            });
            fail("x is null (exact), NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }

        try {
            testStatistic.mannWhitneyUTest(new double[] {
                1.0
            }, null);
            fail("y is null (exact), NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }
    }

    @Test
    void testLargeDatasetExact() {
        /**
         * Expected values computed using R 3.4.4
         */
        final double[] x = {
            11, 22, 19, 22.3, 16, 29, 24, 5.2, 7, 3, 44, 72, 43, 18, 65
        };
        final double[] y = {
            15, 32, 38, 5, 6, 29.1, 31, 73, 88, 70, 50, 60, 93, 112, 190
        };

        assertEquals(59, testStatistic.mannWhitneyU(x, y), 1e-10);
        assertEquals(0.027925057353,
                            testStatistic.mannWhitneyUTest(x, y, false), 1e-9);
        assertEquals(0.02635404434,
                            testStatistic.mannWhitneyUTest(x, y, true), 1e-9);
        // Should default to exact test
        assertEquals(0.02635404434, testStatistic.mannWhitneyUTest(x, y),
                            1e-9);
    }

    @Test
    void testDatasetTooLargeForExact() {
        /**
         * Expected values computed using R 3.4.4
         */
        final double[] x = {
            11, 22, 19, 22.3, 16, 29, 24, 5.2, 7, 3, 44, 72, 43, 18, 65, 69, 71,
            115, 117, 119, 121, 123, 124, 125, 126, 127
        };
        final double[] y = {
            15, 32, 38, 5, 6, 29.1, 31, 73, 88, 70, 50, 60, 93, 112, 190, 200,
            201, 202, 203, 204, 205, 207, 209, 210, 212
        };

        assertEquals(204, testStatistic.mannWhitneyU(x, y), 1e-10);
        assertEquals(0.023177956065,
                            testStatistic.mannWhitneyUTest(x, y, false), 1e-9);
        // Should default to normal approximation
        assertEquals(0.023177956065,
                            testStatistic.mannWhitneyUTest(x, y), 1e-9);
        assertEquals(0.022259264963,
                            testStatistic.mannWhitneyUTest(x, y, true), 1e-9);
    }

    @Test
    void testExactHiP() {
        /**
         * Expected values computed using R 3.4.4
         */
        final double[] x = {
            0, 2, 4, 6, 8, 10, 12, 14, 16, 18
        };
        final double[] y = {
            1, 3, 5, 7, 9, 11, 13, 15, 17, 19
        };
        assertEquals(45, testStatistic.mannWhitneyU(x, y), 1e-10);
        assertEquals(0.7337299957,
                            testStatistic.mannWhitneyUTest(x, y, false), 1e-9);
        assertEquals(0.73936435082,
                            testStatistic.mannWhitneyUTest(x, y, true), 1e-9);

    }

    @Test
    void testExactLarge() {
        /**
         * Expected values computed using R 3.4.4
         */
        final double[] x = {
            0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34,
            36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60
        };
        final double[] y = {
            1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35,
            37, 39, 43, 45, 49, 51, 55, 61, 63, 65, 67, 69
        };
        assertEquals(441, testStatistic.mannWhitneyU(x, y), 1e-10);
        assertEquals(0.73459710599,
                            testStatistic.mannWhitneyUTest(x, y, false), 1e-9);
        assertEquals(0.73642668965,
                            testStatistic.mannWhitneyUTest(x, y, true), 1e-9);
    }

    @Test
    void testExactVerySmall() {
        /**
         * Expected values computed using R 3.4.4
         */
        final double[] x = {
            1,2
        };
        final double[] y = {
            1.5, 2.5
        };
        assertEquals(1, testStatistic.mannWhitneyU(x, y), 1e-10);
        assertEquals(0.66666666667,
                            testStatistic.mannWhitneyUTest(x, y, true), 1e-9);
        assertEquals(0.6985353583,
                            testStatistic.mannWhitneyUTest(x, y, false), 1e-9);

    }

    @Test
    void testExactDegenerate() {
        final double[] x = {
            1
        };
        final double[] y = {
            1.5
        };
        assertEquals(0, testStatistic.mannWhitneyU(x, y), 1e-10);
        assertEquals(1.0,
                            testStatistic.mannWhitneyUTest(x, y, true), 1e-9);
        assertEquals(1.0,
                            testStatistic.mannWhitneyUTest(x, y, false), 1e-9);
    }

    @Test
    void testApproximateWithTies() {
        /**
         * Expected values computed using R 3.4.4
         */
        final double[] x = {
            0, 2, 4, 6, 8, 10, 12, 14, 16, 18
        };
        final double[] y = {
            1, 3, 5, 7, 8, 10, 10, 13, 15, 17, 19
        };
        try {
            testStatistic.mannWhitneyUTest(x, y, true);
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        assertEquals(50.5, testStatistic.mannWhitneyU(x, y), 1e-10);
        assertEquals(0.77784391371,
                            testStatistic.mannWhitneyUTest(x, y, false), 1e-9);

    }

    @Test
    void testApproximateWithTies2() {
        /**
         * Expected values computed using R 3.4.4
         */
        final double[] x = {
            2, 3, 7, 11, 23, 45, 48, 55, 70, 81, 92, 95, 97, 100, 110, 123, 125
        };
        final double[] y = {
            3.5, 4, 8, 12, 25, 46, 49, 56, 70, 81, 92, 95, 97, 97, 100, 112,
            125, 127
        };
        try {
            testStatistic.mannWhitneyUTest(x, y, true);
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        assertEquals(142, testStatistic.mannWhitneyU(x, y), 1e-10);
        assertEquals(0.72874565042,
                            testStatistic.mannWhitneyUTest(x, y, false), 1e-9);

    }

    @Test
    void testIdenticalArrays() {
        /**
         * Expected values computed using R 3.4.4
         */
        final double[] x = {
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10
        };
        final double[] y = {
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10
        };
        try {
            testStatistic.mannWhitneyUTest(x, y, true);
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        assertEquals(50, testStatistic.mannWhitneyU(x, y), 1e-10);
        assertEquals(1.0, testStatistic.mannWhitneyUTest(x, y, false),
                            1e-10);

    }

    @Test
    void testExactThrowsOnTies() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final double[] x = {
                1, 5, 7
            };
            final double[] y = {
                2, 3, 1
            };
            testStatistic.mannWhitneyUTest(x, y, true);
        });
    }

    @Test
    void testBigDataSet() {
        double[] d1 = new double[1500];
        double[] d2 = new double[1500];
        for (int i = 0; i < 1500; i++) {
            d1[i] = 2 * i;
            d2[i] = 2 * i + 1;
        }
        double result = testStatistic.mannWhitneyUTest(d1, d2);
        assertTrue(result > 0.1);
    }

    @Test
    void testBigDataSetOverflow() {
        // MATH-1145
        double[] d1 = new double[110000];
        double[] d2 = new double[110000];
        for (int i = 0; i < 110000; i++) {
            d1[i] = i;
            d2[i] = i;
        }
        assertEquals(1.0, testStatistic.mannWhitneyUTest(d1, d2, false),
                            1E-7);
    }
}
