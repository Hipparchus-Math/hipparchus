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

import java.util.Map;
import java.util.TreeMap;
import java.util.stream.LongStream;

import org.hipparchus.distribution.continuous.NormalDistribution;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.stat.LocalizedStatFormats;
import org.hipparchus.stat.ranking.NaNStrategy;
import org.hipparchus.stat.ranking.NaturalRanking;
import org.hipparchus.stat.ranking.TiesStrategy;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/**
 * An implementation of the Mann-Whitney U test.
 * <p>
 * The definitions and computing formulas used in this implementation follow
 * those in the article,
 * <a href="http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U"> Mann-Whitney U
 * Test</a>
 * <p>
 * In general, results correspond to (and have been tested against) the R
 * wilcox.test function, with {@code exact} meaning the same thing in both APIs
 * and {@code CORRECT} uniformly true in this implementation. For example,
 * wilcox.test(x, y, alternative = "two.sided", mu = 0, paired = FALSE, exact = FALSE
 * correct = TRUE) will return the same p-value as mannWhitneyUTest(x, y,
 * false). The minimum of the W value returned by R for wilcox.test(x, y...) and
 * wilcox.test(y, x...) should equal mannWhitneyU(x, y...).
 */
public class MannWhitneyUTest { // NOPMD - this is not a Junit test class, PMD false positive here

    /**
     * If the combined dataset contains no more values than this, test defaults to
     * exact test.
     */
    private static final int SMALL_SAMPLE_SIZE = 50;

    /** Ranking algorithm. */
    private final NaturalRanking naturalRanking;

    /** Normal distribution */
    private final NormalDistribution standardNormal;

    /**
     * Create a test instance using where NaN's are left in place and ties get
     * the average of applicable ranks.
     */
    public MannWhitneyUTest() {
        naturalRanking = new NaturalRanking(NaNStrategy.FIXED,
                                            TiesStrategy.AVERAGE);
        standardNormal = new NormalDistribution(0, 1);
    }

    /**
     * Create a test instance using the given strategies for NaN's and ties.
     *
     * @param nanStrategy specifies the strategy that should be used for
     *        Double.NaN's
     * @param tiesStrategy specifies the strategy that should be used for ties
     */
    public MannWhitneyUTest(final NaNStrategy nanStrategy,
                            final TiesStrategy tiesStrategy) {
        naturalRanking = new NaturalRanking(nanStrategy, tiesStrategy);
        standardNormal = new NormalDistribution(0, 1);
    }

    /**
     * Computes the
     * <a href="http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U">
     * Mann-Whitney U statistic</a> comparing means for two independent samples
     * possibly of different lengths.
     * <p>
     * This statistic can be used to perform a Mann-Whitney U test evaluating
     * the null hypothesis that the two independent samples have equal mean.
     * <p>
     * Let X<sub>i</sub> denote the i'th individual of the first sample and
     * Y<sub>j</sub> the j'th individual in the second sample. Note that the
     * samples can have different lengths.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>All observations in the two samples are independent.</li>
     * <li>The observations are at least ordinal (continuous are also
     * ordinal).</li>
     * </ul>
     *
     * @param x the first sample
     * @param y the second sample
     * @return Mann-Whitney U statistic (minimum of U<sup>x</sup> and
     *         U<sup>y</sup>)
     * @throws NullArgumentException if {@code x} or {@code y} are {@code null}.
     * @throws MathIllegalArgumentException if {@code x} or {@code y} are
     *         zero-length.
     */
    public double mannWhitneyU(final double[] x, final double[] y)
        throws MathIllegalArgumentException, NullArgumentException {

        ensureDataConformance(x, y);

        final double[] z = concatenateSamples(x, y);
        final double[] ranks = naturalRanking.rank(z);

        double sumRankX = 0;

        /*
         * The ranks for x is in the first x.length entries in ranks because x
         * is in the first x.length entries in z
         */
        for (int i = 0; i < x.length; ++i) {
            sumRankX += ranks[i];
        }

        /*
         * U1 = R1 - (n1 * (n1 + 1)) / 2 where R1 is sum of ranks for sample 1,
         * e.g. x, n1 is the number of observations in sample 1.
         */
        final double U1 = sumRankX - ((long) x.length * (x.length + 1)) / 2;

        /*
         * U1 + U2 = n1 * n2
         */
        final double U2 = (long) x.length * y.length - U1;

        return FastMath.min(U1, U2);
    }

    /**
     * Concatenate the samples into one array.
     *
     * @param x first sample
     * @param y second sample
     * @return concatenated array
     */
    private double[] concatenateSamples(final double[] x, final double[] y) {
        final double[] z = new double[x.length + y.length];

        System.arraycopy(x, 0, z, 0, x.length);
        System.arraycopy(y, 0, z, x.length, y.length);

        return z;
    }

    /**
     * Returns the asymptotic <i>observed significance level</i>, or
     * <a href="http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue">
     * p-value</a>, associated with a <a href=
     * "http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U">Mann-Whitney U
     * Test</a> comparing means for two independent samples.
     * <p>
     * Let X<sub>i</sub> denote the i'th individual of the first sample and
     * Y<sub>j</sub> the j'th individual in the second sample.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>All observations in the two samples are independent.</li>
     * <li>The observations are at least ordinal.</li>
     * </ul>
     * <p>
     * If there are no ties in the data and both samples are small (less than or
     * equal to 50 values in the combined dataset), an exact test is performed;
     * otherwise the test uses the normal approximation (with continuity
     * correction).
     * <p>
     * If the combined dataset contains ties, the variance used in the normal
     * approximation is bias-adjusted using the formula in the reference above.
     *
     * @param x the first sample
     * @param y the second sample
     * @return approximate 2-sized p-value
     * @throws NullArgumentException if {@code x} or {@code y} are {@code null}.
     * @throws MathIllegalArgumentException if {@code x} or {@code y} are
     *         zero-length
     */
    public double mannWhitneyUTest(final double[] x, final double[] y)
        throws MathIllegalArgumentException, NullArgumentException {
        ensureDataConformance(x, y);

        // If samples are both small and there are no ties, perform exact test
        if (x.length + y.length <= SMALL_SAMPLE_SIZE &&
            tiesMap(x, y).isEmpty()) {
            return mannWhitneyUTest(x, y, true);
        } else { // Normal approximation
            return mannWhitneyUTest(x, y, false);
        }
    }

    /**
     * Returns the asymptotic <i>observed significance level</i>, or
     * <a href="http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue">
     * p-value</a>, associated with a <a href=
     * "http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U">Mann-Whitney U
     * Test</a> comparing means for two independent samples.
     * <p>
     * Let X<sub>i</sub> denote the i'th individual of the first sample and
     * Y<sub>j</sub> the j'th individual in the second sample.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>All observations in the two samples are independent.</li>
     * <li>The observations are at least ordinal.</li>
     * </ul>
     * <p>
     * If {@code exact} is {@code true}, the p-value reported is exact, computed
     * using the exact distribution of the U statistic. The computation in this
     * case requires storage on the order of the product of the two sample
     * sizes, so this should not be used for large samples.
     * <p>
     * If {@code exact} is {@code false}, the normal approximation is used to
     * estimate the p-value.
     * <p>
     * If the combined dataset contains ties and {@code exact} is {@code true},
     * MathIllegalArgumentException is thrown. If {@code exact} is {@code false}
     * and the ties are present, the variance used to compute the approximate
     * p-value in the normal approximation is bias-adjusted using the formula in
     * the reference above.
     *
     * @param x the first sample
     * @param y the second sample
     * @param exact true means compute the p-value exactly, false means use the
     *        normal approximation
     * @return approximate 2-sided p-value
     * @throws NullArgumentException if {@code x} or {@code y} are {@code null}.
     * @throws MathIllegalArgumentException if {@code x} or {@code y} are
     *         zero-length or if {@code exact} is {@code true} and ties are
     *         present in the data
     */
    public double mannWhitneyUTest(final double[] x, final double[] y,
                                   final boolean exact)
        throws MathIllegalArgumentException, NullArgumentException {
        ensureDataConformance(x, y);
        final Map<Double, Integer> tiesMap = tiesMap(x, y);
        final double u = mannWhitneyU(x, y);
        if (exact) {
            if (!tiesMap.isEmpty()) {
                throw new MathIllegalArgumentException(LocalizedStatFormats.TIES_ARE_NOT_ALLOWED);
            }
            return exactP(x.length, y.length, u);
        }

        return approximateP(u, x.length, y.length,
                            varU(x.length, y.length, tiesMap));
    }

    /**
     * Ensures that the provided arrays fulfills the assumptions.
     *
     * @param x first sample
     * @param y second sample
     * @throws NullArgumentException if {@code x} or {@code y} are {@code null}.
     * @throws MathIllegalArgumentException if {@code x} or {@code y} are
     *         zero-length.
     */
    private void ensureDataConformance(final double[] x, final double[] y)
        throws MathIllegalArgumentException, NullArgumentException {

        if (x == null || y == null) {
            throw new NullArgumentException();
        }
        if (x.length == 0 || y.length == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NO_DATA);
        }
    }

    /**
     * Estimates the 2-sided p-value associated with a Mann-Whitney U statistic
     * value using the normal approximation.
     * <p>
     * The variance passed in is assumed to be corrected for ties. Continuity
     * correction is applied to the normal approximation.
     *
     * @param u Mann-Whitney U statistic
     * @param n1 number of subjects in first sample
     * @param n2 number of subjects in second sample
     * @param varU variance of U (corrected for ties if these exist)
     * @return two-sided asymptotic p-value
     * @throws MathIllegalStateException if the p-value can not be computed due
     *         to a convergence error
     * @throws MathIllegalStateException if the maximum number of iterations is
     *         exceeded
     */
    private double approximateP(final double u, final int n1, final int n2,
                                final double varU)
        throws MathIllegalStateException {

        final double mu = (long) n1 * n2 / 2.0;

        // If u == mu, return 1
        if (Precision.equals(mu, u)) {
            return 1;
        }

        // Force z <= 0 so we get tail probability. Also apply continuity
        // correction
        final double z = -Math.abs((u - mu) + 0.5) / FastMath.sqrt(varU);

        return 2 * standardNormal.cumulativeProbability(z);
    }

    /**
     * Calculates the (2-sided) p-value associated with a Mann-Whitney U
     * statistic.
     * <p>
     * To compute the p-value, the probability densities for each value of U up
     * to and including u are summed and the resulting tail probability is
     * multiplied by 2.
     * <p>
     * The result of this computation is only valid when the combined n + m
     * sample has no tied values.
     * <p>
     * This method should not be used for large values of n or m as it maintains
     * work arrays of size n*m.
     *
     * @param u Mann-Whitney U statistic value
     * @param n first sample size
     * @param m second sample size
     * @return two-sided exact p-value
     */
    private double exactP(final int n, final int m, final double u) {
        final double nm = m * n;
        if (u > nm) { // Quick exit if u is out of range
            return 1;
        }
        // Need to convert u to a mean deviation, so cumulative probability is
        // tail probability
        final double crit = u < nm / 2 ? u : nm / 2 - u;

        double cum = 0d;
        for (int ct = 0; ct <= crit; ct++) {
            cum += uDensity(n, m, ct);
        }
        return 2 * cum;
    }

    /**
     * Computes the probability density function for the Mann-Whitney U
     * statistic.
     * <p>
     * This method should not be used for large values of n or m as it maintains
     * work arrays of size n*m.
     *
     * @param n first sample size
     * @param m second sample size
     * @param u U-statistic value
     * @return the probability that a U statistic derived from random samples of
     *         size n and m (containing no ties) equals u
     */
    private double uDensity(final int n, final int m, double u) {
        if (u < 0 || u > m * n) {
            return 0;
        }
        final long[] freq = uFrequencies(n, m);
        return freq[(int) FastMath.round(u + 1)] /
               (double) LongStream.of(freq).sum();
    }

    /**
     * Computes frequency counts for values of the Mann-Whitney U statistc. If
     * freq[] is the returned array, freq[u + 1] counts the frequency of U = u
     * among all possible n-m orderings. Therefore, P(u = U) = freq[u + 1] / sum
     * where sum is the sum of the values in the returned array.
     * <p>
     * Implements the algorithm presented in "Algorithm AS 62: A Generator for
     * the Sampling Distribution of the Mann-Whitney U Statistic", L. C. Dinneen
     * and B. C. Blakesley Journal of the Royal Statistical Society. Series C
     * (Applied Statistics) Vol. 22, No. 2 (1973), pp. 269-273.
     *
     * @param n first sample size
     * @param m second sample size
     * @return array of U statistic value frequencies
     */
    private long[] uFrequencies(final int n, final int m) {
        final int max = FastMath.max(m, n);
        if (max > 100) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_LARGE,
                                                   max, 100);
        }
        final int min = FastMath.min(m, n);
        final long[] out = new long[n * m + 2];
        final long[] work = new long[n * m + 2];
        for (int i = 1; i < out.length; i++) {
            out[i] = (i <= (max + 1)) ? 1 : 0;
        }
        work[1] = 0;
        int in = max;
        for (int i = 2; i <= min; i++) {
            work[i] = 0;
            in = in + max;
            int n1 = in + 2;
            long l = 1 + in / 2;
            int k = i;
            for (int j = 1; j <= l; j++) {
                k++;
                n1 = n1 - 1;
                final long sum = out[j] + work[j];
                out[j] = sum;
                work[k] = sum - out[n1];
                out[n1] = sum;
            }
        }
        return out;
    }

    /**
     * Computes the variance for a U-statistic associated with samples of
     * sizes{@code n} and {@code m} and ties described by {@code tiesMap}. If
     * {@code tiesMap} is non-empty, the multiplicity counts in its values set
     * are used to adjust the variance.
     *
     * @param n first sample size
     * @param m second sample size
     * @param tiesMap map of <value, multiplicity>
     * @return ties-adjusted variance
     */
    private double varU(final int n, final int m,
                        Map<Double, Integer> tiesMap) {
        final double nm = (long) n * m;
        if (tiesMap.isEmpty()) {
            return nm * (n + m + 1) / 12.0;
        }
        final long tSum = tiesMap.entrySet().stream()
            .mapToLong(e -> e.getValue() * e.getValue() * e.getValue() -
                            e.getValue())
            .sum();
        final double totalN = n + m;
        return (nm / 12) * (totalN + 1 - tSum / (totalN * (totalN - 1)));

    }

    /**
     * Creates a map whose keys are values occurring more than once in the
     * combined dataset formed from x and y. Map entry values are the number of
     * occurrences. The returned map is empty iff there are no ties in the data.
     *
     * @param x first dataset
     * @param y second dataset
     * @return map of <value, number of times it occurs> for values occurring
     *         more than once or an empty map if there are no ties (the returned
     *         map is <em>not</em> thread-safe, which is OK in the context of the callers)
     */
    private Map<Double, Integer> tiesMap(final double[] x, final double[] y) {
        final Map<Double, Integer> tiesMap = new TreeMap<>(); // NOPMD - no concurrent access in the callers context
        for (int i = 0; i < x.length; i++) {
            tiesMap.merge(x[i], 1, Integer::sum);
        }
        for (int i = 0; i < y.length; i++) {
            tiesMap.merge(y[i], 1, Integer::sum);
        }
        tiesMap.entrySet().removeIf(e -> e.getValue() == 1);
        return tiesMap;
    }
}
