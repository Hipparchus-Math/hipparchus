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

import java.util.ArrayList;
import java.util.List;

import org.hipparchus.distribution.continuous.NormalDistribution;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.stat.ranking.NaNStrategy;
import org.hipparchus.stat.ranking.NaturalRanking;
import org.hipparchus.stat.ranking.TiesStrategy;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;

/**
 * An implementation of the Wilcoxon signed-rank test.
 *
 * This implementation currently handles only paired (equal length) samples
 * and discards tied pairs from the analysis. The latter behavior differs from
 * the R implementation of wilcox.test and corresponds to the "wilcox"
 * zero_method configurable in scipy.stats.wilcoxon.
 */
public class WilcoxonSignedRankTest { // NOPMD - this is not a Junit test class, PMD false positive here

    /** Ranking algorithm. */
    private final NaturalRanking naturalRanking;

    /**
     * Create a test instance where NaN's are left in place and ties get the
     * average of applicable ranks.
     */
    public WilcoxonSignedRankTest() {
        naturalRanking = new NaturalRanking(NaNStrategy.FIXED,
                                            TiesStrategy.AVERAGE);
    }

    /**
     * Create a test instance using the given strategies for NaN's and ties.
     *
     * @param nanStrategy specifies the strategy that should be used for
     *        Double.NaN's
     * @param tiesStrategy specifies the strategy that should be used for ties
     */
    public WilcoxonSignedRankTest(final NaNStrategy nanStrategy,
                                  final TiesStrategy tiesStrategy) {
        naturalRanking = new NaturalRanking(nanStrategy, tiesStrategy);
    }

    /**
     * Ensures that the provided arrays fulfills the assumptions. Also computes
     * and returns the number of tied pairs (i.e., zero differences).
     *
     * @param x first sample
     * @param y second sample
     * @return the number of indices where x[i] == y[i]
     * @throws NullArgumentException if {@code x} or {@code y} are {@code null}.
     * @throws MathIllegalArgumentException if {@code x} or {@code y} are
     *         zero-length
     * @throws MathIllegalArgumentException if {@code x} and {@code y} do not
     *         have the same length.
     * @throws MathIllegalArgumentException if all pairs are tied (i.e., if no
     *         data remains when tied pairs have been removed.
     */
    private int ensureDataConformance(final double[] x, final double[] y)
        throws MathIllegalArgumentException, NullArgumentException {

        if (x == null || y == null) {
            throw new NullArgumentException();
        }
        if (x.length == 0 || y.length == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NO_DATA);
        }
        MathArrays.checkEqualLength(y, x);
        int nTies = 0;
        for (int i = 0; i < x.length; i++) {
            if (x[i] == y[i]) {
                nTies++;
            }
        }
        if (x.length - nTies == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.INSUFFICIENT_DATA);
        }
        return nTies;
    }

    /**
     * Calculates y[i] - x[i] for all i, discarding ties.
     *
     * @param x first sample
     * @param y second sample
     * @return z = y - x (minus tied values)
     */
    private double[] calculateDifferences(final double[] x, final double[] y) {
        final List<Double> differences = new ArrayList<>();
        for (int i = 0; i < x.length; ++i) {
            if (y[i] != x[i]) {
                differences.add(y[i] - x[i]);
            }
        }
        final int nDiff = differences.size();
        final double[] z = new double[nDiff];
        for (int i = 0; i < nDiff; i++) {
            z[i] = differences.get(i);
        }
        return z;
    }

    /**
     * Calculates |z[i]| for all i
     *
     * @param z sample
     * @return |z|
     * @throws NullArgumentException if {@code z} is {@code null}
     * @throws MathIllegalArgumentException if {@code z} is zero-length.
     */
    private double[] calculateAbsoluteDifferences(final double[] z)
        throws MathIllegalArgumentException, NullArgumentException {

        if (z == null) {
            throw new NullArgumentException();
        }

        if (z.length == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NO_DATA);
        }

        final double[] zAbs = new double[z.length];

        for (int i = 0; i < z.length; ++i) {
            zAbs[i] = FastMath.abs(z[i]);
        }

        return zAbs;
    }

    /**
     * Computes the
     * <a href="http://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test">
     * Wilcoxon signed ranked statistic</a> comparing means for two related
     * samples or repeated measurements on a single sample.
     * <p>
     * This statistic can be used to perform a Wilcoxon signed ranked test
     * evaluating the null hypothesis that the two related samples or repeated
     * measurements on a single sample have equal mean.
     * </p>
     * <p>
     * Let X<sub>i</sub> denote the i'th individual of the first sample and
     * Y<sub>i</sub> the related i'th individual in the second sample. Let
     * Z<sub>i</sub> = Y<sub>i</sub> - X<sub>i</sub>.
     * </p>
     * <p>* <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The differences Z<sub>i</sub> must be independent.</li>
     * <li>Each Z<sub>i</sub> comes from a continuous population (they must be
     * identical) and is symmetric about a common median.</li>
     * <li>The values that X<sub>i</sub> and Y<sub>i</sub> represent are
     * ordered, so the comparisons greater than, less than, and equal to are
     * meaningful.</li>
     * </ul>
     *
     * @param x the first sample
     * @param y the second sample
     * @return wilcoxonSignedRank statistic (the larger of W+ and W-)
     * @throws NullArgumentException if {@code x} or {@code y} are {@code null}.
     * @throws MathIllegalArgumentException if {@code x} or {@code y} are
     *         zero-length.
     * @throws MathIllegalArgumentException if {@code x} and {@code y} do not
     *         have the same length.
     */
    public double wilcoxonSignedRank(final double[] x, final double[] y)
        throws MathIllegalArgumentException, NullArgumentException {

        ensureDataConformance(x, y);

        final double[] z = calculateDifferences(x, y);
        final double[] zAbs = calculateAbsoluteDifferences(z);

        final double[] ranks = naturalRanking.rank(zAbs);

        double Wplus = 0;

        for (int i = 0; i < z.length; ++i) {
            if (z[i] > 0) {
                Wplus += ranks[i];
            }
        }

        final int n = z.length;
        final double Wminus = ((n * (n + 1)) / 2.0) - Wplus;

        return FastMath.max(Wplus, Wminus);
    }

    /**
     * Calculates the p-value associated with a Wilcoxon signed rank statistic
     * by enumerating all possible rank sums and counting the number that exceed
     * the given value.
     *
     * @param stat Wilcoxon signed rank statistic value
     * @param n number of subjects (corresponding to x.length)
     * @return two-sided exact p-value
     */
    private double calculateExactPValue(final double stat, final int n) {
        final int m = 1 << n;
        int largerRankSums = 0;
        for (int i = 0; i < m; ++i) {
            int rankSum = 0;

            // Generate all possible rank sums
            for (int j = 0; j < n; ++j) {

                // (i >> j) & 1 extract i's j-th bit from the right
                if (((i >> j) & 1) == 1) {
                    rankSum += j + 1;
                }
            }

            if (rankSum >= stat) {
                ++largerRankSums;
            }
        }

        /*
         * largerRankSums / m gives the one-sided p-value, so it's multiplied
         * with 2 to get the two-sided p-value
         */
        return 2 * ((double) largerRankSums) / m;
    }

    /**
     * Computes an estimate of the (2-sided) p-value using the normal
     * approximation. Includes a continuity correction in computing the
     * correction factor.
     *
     * @param stat Wilcoxon rank sum statistic
     * @param n number of subjects (corresponding to x.length minus any tied ranks)
     * @return two-sided asymptotic p-value
     */
    private double calculateAsymptoticPValue(final double stat, final int n) {

        final double ES = n * (n + 1) / 4.0;

        /*
         * Same as (but saves computations): final double VarW = ((double) (N *
         * (N + 1) * (2*N + 1))) / 24;
         */
        final double VarS = ES * ((2 * n + 1) / 6.0);

        double z = stat - ES;
        final double t = FastMath.signum(z);
        z = (z - t * 0.5) / FastMath.sqrt(VarS);

        // want 2-sided tail probability, so make sure z < 0
        if (z > 0) {
            z = -z;
        }
        final NormalDistribution standardNormal = new NormalDistribution(0, 1);
        return 2 * standardNormal.cumulativeProbability(z);
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <a href= "http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue">
     * p-value</a>, associated with a
     * <a href="http://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test">
     * Wilcoxon signed ranked statistic</a> comparing mean for two related
     * samples or repeated measurements on a single sample.
     * <p>
     * Let X<sub>i</sub> denote the i'th individual of the first sample and
     * Y<sub>i</sub> the related i'th individual in the second sample. Let
     * Z<sub>i</sub> = Y<sub>i</sub> - X<sub>i</sub>.
     * </p>
     * <p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The differences Z<sub>i</sub> must be independent.</li>
     * <li>Each Z<sub>i</sub> comes from a continuous population (they must be
     * identical) and is symmetric about a common median.</li>
     * <li>The values that X<sub>i</sub> and Y<sub>i</sub> represent are
     * ordered, so the comparisons greater than, less than, and equal to are
     * meaningful.</li>
     * </ul>
     * <p><strong>Implementation notes</strong>:</p>
     * <ul>
     * <li>Tied pairs are discarded from the data.</li>
     * <li>When {@code exactPValue} is false, the normal approximation is used
     * to estimate the p-value including a continuity correction factor.
     * {@code wilcoxonSignedRankTest(x, y, true)} should give the same results
     * as {@code  wilcox.test(x, y, alternative = "two.sided", mu = 0,
     *     paired = TRUE, exact = FALSE, correct = TRUE)} in R (as long as
     * there are no tied pairs in the data).</li>
     * </ul>
     *
     * @param x the first sample
     * @param y the second sample
     * @param exactPValue if the exact p-value is wanted (only works for
     *        x.length &lt;= 30, if true and x.length &gt; 30, MathIllegalArgumentException is thrown)
     * @return p-value
     * @throws NullArgumentException if {@code x} or {@code y} are {@code null}.
     * @throws MathIllegalArgumentException if {@code x} or {@code y} are
     *         zero-length or for all i, x[i] == y[i]
     * @throws MathIllegalArgumentException if {@code x} and {@code y} do not
     *         have the same length.
     * @throws MathIllegalArgumentException if {@code exactPValue} is
     *         {@code true} and {@code x.length} &gt; 30
     * @throws MathIllegalStateException if the p-value can not be computed due
     *         to a convergence error
     * @throws MathIllegalStateException if the maximum number of iterations is
     *         exceeded
     */
    public double wilcoxonSignedRankTest(final double[] x, final double[] y,
                                         final boolean exactPValue)
        throws MathIllegalArgumentException, NullArgumentException,
        MathIllegalStateException {

        final int nTies = ensureDataConformance(x, y);

        final int n = x.length - nTies;
        final double stat = wilcoxonSignedRank(x, y);

        if (exactPValue && n > 30) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_LARGE,
                                                   n, 30);
        }

        if (exactPValue) {
            return calculateExactPValue(stat, n);
        } else {
            return calculateAsymptoticPValue(stat, n);
        }
    }
}
