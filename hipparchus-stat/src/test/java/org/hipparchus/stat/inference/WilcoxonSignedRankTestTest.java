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
import org.junit.Assert;
import org.junit.Test;

/**
 * Test cases for the WilcoxonSignedRangTest class.
 */

public class WilcoxonSignedRankTestTest {

    protected WilcoxonSignedRankTest testStatistic = new WilcoxonSignedRankTest();

    @Test
    public void testWilcoxonSignedRankSimple() {
        /*
         * Hollandar and Wolfe data from R docs.
         * Target values computed using R version 3.4.4.
         * x <- c(1.83, 0.50, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
         * y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
         */
        final double x[] = {
            1.83, 0.50, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.30
        };
        final double y[] = {
            0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29
        };

        /*
         * EXACT:
         * wilcox.test(x, y, alternative = "two.sided", mu = 0, paired = TRUE
         * exact = TRUE, correct = FALSE)
         * V = 40, p-value = 0.03906
         * Expected values are from R, version 3.4.4.
         */
        Assert.assertEquals(40, testStatistic.wilcoxonSignedRank(x, y), 1e-10);
        Assert.assertEquals(0.03906,
                            testStatistic.wilcoxonSignedRankTest(x, y, true),
                            1e-5);

        /*
         * ASYMPTOTIC:
         * wilcox.test(x, y, alternative = "two.sided", mu = 0,
         * paired = TRUE, exact = FALSE, correct = TRUE)
         * V = 40, p-value = 0.044010984013
         */
        Assert.assertEquals(40, testStatistic.wilcoxonSignedRank(x, y), 1e-10);
        Assert.assertEquals(0.044010984013,
                            testStatistic.wilcoxonSignedRankTest(x, y, false),
                            1e-10);
    }

    @Test
    public void testWilcoxonSignedRankSimple2() {
        /*
         * x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91)
         * y <- c(1.15, 0.88, 0.90, 0.74, 1.21, 2.0, 1.72)
         * Expected values are from R, version 3.4.4.
         */
        final double[] x = {0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91};
        final double[] y = {1.15, 0.88, 0.90, 0.74, 1.21, 2.0, 1.72};
        Assert.assertEquals(16,  testStatistic.wilcoxonSignedRank(x, y), 0);
        // Exact
        Assert.assertEquals(0.8125,
                            testStatistic.wilcoxonSignedRankTest(x, y, true),
                            1e-10);
        // Asymptotic
        Assert.assertEquals(0.79984610566,
                            testStatistic.wilcoxonSignedRankTest(x, y, false),
                            1e-10);
    }


    @Test
    public void testWilcoxonSignedRankTiesDiscarded() {
        /*
         * Verify that tied pairs are discarded.
         */
        final double x[] = {
            1.83, 0.50, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.30
        };
        final double y[] = {
            1.83, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29
        };
        final double xp[] = {
            0.50, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.30
        };
        final double yp[] = {
            0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29
        };
        Assert.assertEquals(testStatistic.wilcoxonSignedRank(xp, yp),
                            testStatistic.wilcoxonSignedRank(x, y),
                            0);
        Assert.assertEquals(testStatistic.wilcoxonSignedRankTest(xp, yp, true),
                            testStatistic.wilcoxonSignedRankTest(x, y, true),
                            0);
        Assert.assertEquals(testStatistic.wilcoxonSignedRankTest(xp, yp, false),
                            testStatistic.wilcoxonSignedRankTest(x, y, false),
                            0);
    }

    @Test
    public void testWilcoxonSignedRankInputValidation() {
        /*
         * Exact only for sample size <= 30
         */
        final double[] x1 = new double[30];
        final double[] x2 = new double[31];
        final double[] y1 = new double[30];
        final double[] y2 = new double[31];
        for (int i = 0; i < 30; ++i) {
            x1[i] = x2[i] = y1[i] = y2[i] = i;
        }

        // Exactly 30 is okay
        // testStatistic.wilcoxonSignedRankTest(x1, y1, true);

        try {
            testStatistic.wilcoxonSignedRankTest(x2, y2, true);
            Assert
                .fail("More than 30 samples and exact chosen, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        /*
         * Samples must be present, i.e. length > 0
         */
        try {
            testStatistic.wilcoxonSignedRankTest(new double[] {}, new double[] {
                1.0
            }, true);
            Assert
                .fail("x does not contain samples (exact), MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            testStatistic.wilcoxonSignedRankTest(new double[] {}, new double[] {
                1.0
            }, false);
            Assert
                .fail("x does not contain samples (asymptotic), MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            testStatistic.wilcoxonSignedRankTest(new double[] {
                1.0
            }, new double[] {}, true);
            Assert
                .fail("y does not contain samples (exact), MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            testStatistic.wilcoxonSignedRankTest(new double[] {
                1.0
            }, new double[] {}, false);
            Assert
                .fail("y does not contain samples (asymptotic), MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        /*
         * Samples not same size, i.e. cannot be paired
         */
        try {
            testStatistic.wilcoxonSignedRankTest(new double[] {
                1.0, 2.0
            }, new double[] {
                3.0
            }, true);
            Assert
                .fail("x and y not same size (exact), MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            testStatistic.wilcoxonSignedRankTest(new double[] {
                1.0, 2.0
            }, new double[] {
                3.0
            }, false);
            Assert
                .fail("x and y not same size (asymptotic), MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        /*
         * x and y is null
         */
        try {
            testStatistic.wilcoxonSignedRankTest(null, null, true);
            Assert
                .fail("x and y is null (exact), NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }

        try {
            testStatistic.wilcoxonSignedRankTest(null, null, false);
            Assert
                .fail("x and y is null (asymptotic), NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }

        /*
         * x or y is null
         */
        try {
            testStatistic.wilcoxonSignedRankTest(null, new double[] {
                1.0
            }, true);
            Assert.fail("x is null (exact), NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }

        try {
            testStatistic.wilcoxonSignedRankTest(null, new double[] {
                1.0
            }, false);
            Assert
                .fail("x is null (asymptotic), NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }

        try {
            testStatistic.wilcoxonSignedRankTest(new double[] {
                1.0
            }, null, true);
            Assert.fail("y is null (exact), NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }

        try {
            testStatistic.wilcoxonSignedRankTest(new double[] {
                1.0
            }, null, false);
            Assert
                .fail("y is null (asymptotic), NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }
    }

    @Test(expected = MathIllegalArgumentException.class)
    public void testBadInputAllTies() {
        testStatistic.wilcoxonSignedRankTest(new double[] {
            1.0, 2.0, 3.0
        }, new double[] {
            1.0, 2.0, 3.0
        }, true);
    }

    @Test
    public void testDegenerateOnePair() {
        final double[] x = {1};
        final double[] y = {2};
        Assert.assertEquals(1.0, testStatistic.wilcoxonSignedRank(x,y), 0);
        Assert.assertEquals(1.0, testStatistic.wilcoxonSignedRankTest(x,y, true), 0);
    }

}
