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

import java.util.Collection;

import org.hipparchus.distribution.RealDistribution;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.stat.descriptive.StatisticalSummary;

/**
 * A collection of static methods to create inference test instances or to
 * perform inference tests.
 */
public class InferenceTestUtils  {

    /** Singleton TTest instance. */
    private static final TTest T_TEST = new TTest();

    /** Singleton ChiSquareTest instance. */
    private static final ChiSquareTest CHI_SQUARE_TEST = new ChiSquareTest();

    /** Singleton OneWayAnova instance. */
    private static final OneWayAnova ONE_WAY_ANANOVA = new OneWayAnova();

    /** Singleton G-Test instance. */
    private static final GTest G_TEST = new GTest();

    /** Singleton K-S test instance */
    private static final KolmogorovSmirnovTest KS_TEST = new KolmogorovSmirnovTest();

    /**
     * Prevent instantiation.
     */
    private InferenceTestUtils() {
        super();
    }

    /**
     * Computes a 2-sample t statistic,  under the hypothesis of equal
     * subpopulation variances.  To compute a t-statistic without the
     * equal variances hypothesis, use {@link #t(double[], double[])}.
     * <p>
     * This statistic can be used to perform a (homoscedastic) two-sample
     * t-test to compare sample means.</p>
     * <p>
     * The t-statistic is</p>
     * <p>
     * &nbsp;&nbsp;<code>  t = (m1 - m2) / (sqrt(1/n1 +1/n2) sqrt(var))</code>
     * </p><p>
     * where <strong><code>n1</code></strong> is the size of first sample;
     * <strong><code> n2</code></strong> is the size of second sample;
     * <strong><code> m1</code></strong> is the mean of first sample;
     * <strong><code> m2</code></strong> is the mean of second sample
     * and <strong><code>var</code></strong> is the pooled variance estimate:
     * </p><p>
     * <code>var = sqrt(((n1 - 1)var1 + (n2 - 1)var2) / ((n1-1) + (n2-1)))</code>
     * </p><p>
     * with <strong><code>var1</code></strong> the variance of the first sample and
     * <strong><code>var2</code></strong> the variance of the second sample.
     * </p><p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The observed array lengths must both be at least 2.
     * </li></ul>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @return t statistic
     * @throws NullArgumentException if the arrays are <code>null</code>
     * @throws MathIllegalArgumentException if the length of the arrays is &lt; 2
     */
    public static double homoscedasticT(final double[] sample1, final double[] sample2)
        throws MathIllegalArgumentException, NullArgumentException {
        return T_TEST.homoscedasticT(sample1, sample2);
    }

    /**
     * Computes a 2-sample t statistic, comparing the means of the datasets
     * described by two {@link StatisticalSummary} instances, under the
     * assumption of equal subpopulation variances.  To compute a t-statistic
     * without the equal variances assumption, use
     * {@link #t(StatisticalSummary, StatisticalSummary)}.
     * <p>
     * This statistic can be used to perform a (homoscedastic) two-sample
     * t-test to compare sample means.</p>
     * <p>
     * The t-statistic returned is</p>
     * <p>
     * &nbsp;&nbsp;<code>  t = (m1 - m2) / (sqrt(1/n1 +1/n2) sqrt(var))</code>
     * </p><p>
     * where <strong><code>n1</code></strong> is the size of first sample;
     * <strong><code> n2</code></strong> is the size of second sample;
     * <strong><code> m1</code></strong> is the mean of first sample;
     * <strong><code> m2</code></strong> is the mean of second sample
     * and <strong><code>var</code></strong> is the pooled variance estimate:
     * </p><p>
     * <code>var = sqrt(((n1 - 1)var1 + (n2 - 1)var2) / ((n1-1) + (n2-1)))</code>
     * </p><p>
     * with <strong><code>var1</code></strong> the variance of the first sample and
     * <strong><code>var2</code></strong> the variance of the second sample.
     * </p><p>
     * <strong>Preconditions</strong>:</p><ul>
     * <li>The datasets described by the two Univariates must each contain
     * at least 2 observations.
     * </li></ul>
     *
     * @param sampleStats1 StatisticalSummary describing data from the first sample
     * @param sampleStats2 StatisticalSummary describing data from the second sample
     * @return t statistic
     * @throws NullArgumentException if the sample statistics are <code>null</code>
     * @throws MathIllegalArgumentException if the number of samples is &lt; 2
     */
    public static double homoscedasticT(final StatisticalSummary sampleStats1,
                                        final StatisticalSummary sampleStats2)
        throws MathIllegalArgumentException, NullArgumentException {
        return T_TEST.homoscedasticT(sampleStats1, sampleStats2);
    }

    /**
     * Performs a
     * <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm">
     * two-sided t-test</a> evaluating the null hypothesis that <code>sample1</code>
     * and <code>sample2</code> are drawn from populations with the same mean,
     * with significance level <code>alpha</code>,  assuming that the
     * subpopulation variances are equal.  Use
     * {@link #tTest(double[], double[], double)} to perform the test without
     * the assumption of equal variances.
     * <p>
     * Returns <code>true</code> iff the null hypothesis that the means are
     * equal can be rejected with confidence <code>1 - alpha</code>.  To
     * perform a 1-sided test, use <code>alpha * 2.</code>  To perform the test
     * without the assumption of equal subpopulation variances, use
     * {@link #tTest(double[], double[], double)}.</p>
     * <p>
     * A pooled variance estimate is used to compute the t-statistic. See
     * {@link #t(double[], double[])} for the formula. The sum of the sample
     * sizes minus 2 is used as the degrees of freedom.</p>
     * <p>
     * <strong>Examples:</strong></p><ol>
     * <li>To test the (2-sided) hypothesis <code>mean 1 = mean 2 </code> at
     * the 95% level, use <br><code>tTest(sample1, sample2, 0.05). </code>
     * </li>
     * <li>To test the (one-sided) hypothesis <code> mean 1 &lt; mean 2, </code>
     * at the 99% level, first verify that the measured mean of
     * <code>sample 1</code> is less than the mean of <code>sample 2</code>
     * and then use
     * <br><code>tTest(sample1, sample2, 0.02) </code>
     * </li></ol>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The observed array lengths must both be at least 2.
     * </li>
     * <li> <code> 0 &lt; alpha &lt; 0.5 </code>
     * </li></ul>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @param alpha significance level of the test
     * @return true if the null hypothesis can be rejected with
     * confidence 1 - alpha
     * @throws NullArgumentException if the arrays are <code>null</code>
     * @throws MathIllegalArgumentException if the length of the arrays is &lt; 2
     * @throws MathIllegalArgumentException if <code>alpha</code> is not in the range (0, 0.5]
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static boolean homoscedasticTTest(final double[] sample1, final double[] sample2,
                                             final double alpha)
        throws MathIllegalArgumentException, NullArgumentException, MathIllegalStateException {
        return T_TEST.homoscedasticTTest(sample1, sample2, alpha);
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <i>p-value</i>, associated with a two-sample, two-tailed t-test
     * comparing the means of the input arrays, under the assumption that
     * the two samples are drawn from subpopulations with equal variances.
     * To perform the test without the equal variances assumption, use
     * {@link #tTest(double[], double[])}.
     * <p>
     * The number returned is the smallest significance level
     * at which one can reject the null hypothesis that the two means are
     * equal in favor of the two-sided alternative that they are different.
     * For a one-sided test, divide the returned value by 2.</p>
     * <p>
     * A pooled variance estimate is used to compute the t-statistic.  See
     * {@link #homoscedasticT(double[], double[])}. The sum of the sample sizes
     * minus 2 is used as the degrees of freedom.</p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the p-value depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The observed array lengths must both be at least 2.
     * </li></ul>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @return p-value for t-test
     * @throws NullArgumentException if the arrays are <code>null</code>
     * @throws MathIllegalArgumentException if the length of the arrays is &lt; 2
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static double homoscedasticTTest(final double[] sample1, final double[] sample2)
        throws MathIllegalArgumentException, NullArgumentException, MathIllegalStateException {
        return T_TEST.homoscedasticTTest(sample1, sample2);
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <i>p-value</i>, associated with a two-sample, two-tailed t-test
     * comparing the means of the datasets described by two StatisticalSummary
     * instances, under the hypothesis of equal subpopulation variances. To
     * perform a test without the equal variances assumption, use
     * {@link #tTest(StatisticalSummary, StatisticalSummary)}.
     * <p>
     * The number returned is the smallest significance level
     * at which one can reject the null hypothesis that the two means are
     * equal in favor of the two-sided alternative that they are different.
     * For a one-sided test, divide the returned value by 2.</p>
     * <p>
     * See {@link #homoscedasticT(double[], double[])} for the formula used to
     * compute the t-statistic. The sum of the  sample sizes minus 2 is used as
     * the degrees of freedom.</p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the p-value depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">here</a>
     * </p><p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The datasets described by the two Univariates must each contain
     * at least 2 observations.
     * </li></ul>
     *
     * @param sampleStats1  StatisticalSummary describing data from the first sample
     * @param sampleStats2  StatisticalSummary describing data from the second sample
     * @return p-value for t-test
     * @throws NullArgumentException if the sample statistics are <code>null</code>
     * @throws MathIllegalArgumentException if the number of samples is &lt; 2
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static double homoscedasticTTest(final StatisticalSummary sampleStats1,
                                            final StatisticalSummary sampleStats2)
        throws MathIllegalArgumentException, NullArgumentException, MathIllegalStateException {
        return T_TEST.homoscedasticTTest(sampleStats1, sampleStats2);
    }

    /**
     * Computes a paired, 2-sample t-statistic based on the data in the input
     * arrays.  The t-statistic returned is equivalent to what would be returned by
     * computing the one-sample t-statistic {@link #t(double, double[])}, with
     * <code>mu = 0</code> and the sample array consisting of the (signed)
     * differences between corresponding entries in <code>sample1</code> and
     * <code>sample2.</code>
     * <p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The input arrays must have the same length and their common length
     * must be at least 2.
     * </li></ul>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @return t statistic
     * @throws NullArgumentException if the arrays are <code>null</code>
     * @throws MathIllegalArgumentException if the arrays are empty
     * @throws MathIllegalArgumentException if the length of the arrays is not equal
     * @throws MathIllegalArgumentException if the length of the arrays is &lt; 2
     */
    public static double pairedT(final double[] sample1, final double[] sample2)
        throws MathIllegalArgumentException, NullArgumentException {
        return T_TEST.pairedT(sample1, sample2);
    }

    /**
     * Performs a paired t-test evaluating the null hypothesis that the
     * mean of the paired differences between <code>sample1</code> and
     * <code>sample2</code> is 0 in favor of the two-sided alternative that the
     * mean paired difference is not equal to 0, with significance level
     * <code>alpha</code>.
     * <p>
     * Returns <code>true</code> iff the null hypothesis can be rejected with
     * confidence <code>1 - alpha</code>.  To perform a 1-sided test, use
     * <code>alpha * 2</code></p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The input array lengths must be the same and their common length
     * must be at least 2.
     * </li>
     * <li> <code> 0 &lt; alpha &lt; 0.5 </code>
     * </li></ul>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @param alpha significance level of the test
     * @return true if the null hypothesis can be rejected with
     * confidence 1 - alpha
     * @throws NullArgumentException if the arrays are <code>null</code>
     * @throws MathIllegalArgumentException if the arrays are empty
     * @throws MathIllegalArgumentException if the length of the arrays is not equal
     * @throws MathIllegalArgumentException if the length of the arrays is &lt; 2
     * @throws MathIllegalArgumentException if <code>alpha</code> is not in the range (0, 0.5]
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static boolean pairedTTest(final double[] sample1, final double[] sample2,
                                      final double alpha)
        throws MathIllegalArgumentException, NullArgumentException, MathIllegalStateException {
        return T_TEST.pairedTTest(sample1, sample2, alpha);
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <i> p-value</i>, associated with a paired, two-sample, two-tailed t-test
     * based on the data in the input arrays.
     * <p>
     * The number returned is the smallest significance level
     * at which one can reject the null hypothesis that the mean of the paired
     * differences is 0 in favor of the two-sided alternative that the mean paired
     * difference is not equal to 0. For a one-sided test, divide the returned
     * value by 2.</p>
     * <p>
     * This test is equivalent to a one-sample t-test computed using
     * {@link #tTest(double, double[])} with <code>mu = 0</code> and the sample
     * array consisting of the signed differences between corresponding elements of
     * <code>sample1</code> and <code>sample2.</code></p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the p-value depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The input array lengths must be the same and their common length must
     * be at least 2.
     * </li></ul>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @return p-value for t-test
     * @throws NullArgumentException if the arrays are <code>null</code>
     * @throws MathIllegalArgumentException if the arrays are empty
     * @throws MathIllegalArgumentException if the length of the arrays is not equal
     * @throws MathIllegalArgumentException if the length of the arrays is &lt; 2
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static double pairedTTest(final double[] sample1, final double[] sample2)
        throws MathIllegalArgumentException, NullArgumentException, MathIllegalStateException {
        return T_TEST.pairedTTest(sample1, sample2);
    }

    /**
     * Computes a <a href="http://www.itl.nist.gov/div898/handbook/prc/section2/prc22.htm#formula">
     * t statistic </a> given observed values and a comparison constant.
     * <p>
     * This statistic can be used to perform a one sample t-test for the mean.
     * </p><p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The observed array length must be at least 2.
     * </li></ul>
     *
     * @param mu comparison constant
     * @param observed array of values
     * @return t statistic
     * @throws NullArgumentException if <code>observed</code> is <code>null</code>
     * @throws MathIllegalArgumentException if the length of <code>observed</code> is &lt; 2
     */
    public static double t(final double mu, final double[] observed)
        throws MathIllegalArgumentException, NullArgumentException {
        return T_TEST.t(mu, observed);
    }

    /**
     * Computes a <a href="http://www.itl.nist.gov/div898/handbook/prc/section2/prc22.htm#formula">
     * t statistic </a> to use in comparing the mean of the dataset described by
     * <code>sampleStats</code> to <code>mu</code>.
     * <p>
     * This statistic can be used to perform a one sample t-test for the mean.
     * </p><p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li><code>observed.getN() &ge; 2</code>.
     * </li></ul>
     *
     * @param mu comparison constant
     * @param sampleStats DescriptiveStatistics holding sample summary statitstics
     * @return t statistic
     * @throws NullArgumentException if <code>sampleStats</code> is <code>null</code>
     * @throws MathIllegalArgumentException if the number of samples is &lt; 2
     */
    public static double t(final double mu, final StatisticalSummary sampleStats)
        throws MathIllegalArgumentException, NullArgumentException {
        return T_TEST.t(mu, sampleStats);
    }

    /**
     * Computes a 2-sample t statistic, without the hypothesis of equal
     * subpopulation variances.  To compute a t-statistic assuming equal
     * variances, use {@link #homoscedasticT(double[], double[])}.
     * <p>
     * This statistic can be used to perform a two-sample t-test to compare
     * sample means.</p>
     * <p>
     * The t-statistic is</p>
     * <p>
     * &nbsp;&nbsp; <code>  t = (m1 - m2) / sqrt(var1/n1 + var2/n2)</code>
     * </p><p>
     *  where <strong><code>n1</code></strong> is the size of the first sample
     * <strong><code> n2</code></strong> is the size of the second sample;
     * <strong><code> m1</code></strong> is the mean of the first sample;
     * <strong><code> m2</code></strong> is the mean of the second sample;
     * <strong><code> var1</code></strong> is the variance of the first sample;
     * <strong><code> var2</code></strong> is the variance of the second sample;
     * </p><p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The observed array lengths must both be at least 2.
     * </li></ul>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @return t statistic
     * @throws NullArgumentException if the arrays are <code>null</code>
     * @throws MathIllegalArgumentException if the length of the arrays is &lt; 2
     */
    public static double t(final double[] sample1, final double[] sample2)
        throws MathIllegalArgumentException, NullArgumentException {
        return T_TEST.t(sample1, sample2);
    }

    /**
     * Computes a 2-sample t statistic, comparing the means of the datasets
     * described by two {@link StatisticalSummary} instances, without the
     * assumption of equal subpopulation variances.  Use
     * {@link #homoscedasticT(StatisticalSummary, StatisticalSummary)} to
     * compute a t-statistic under the equal variances assumption.
     * <p>
     * This statistic can be used to perform a two-sample t-test to compare
     * sample means.</p>
     * <p>
      * The returned  t-statistic is</p>
     * <p>
     * &nbsp;&nbsp; <code>  t = (m1 - m2) / sqrt(var1/n1 + var2/n2)</code>
     * </p><p>
     * where <strong><code>n1</code></strong> is the size of the first sample;
     * <strong><code> n2</code></strong> is the size of the second sample;
     * <strong><code> m1</code></strong> is the mean of the first sample;
     * <strong><code> m2</code></strong> is the mean of the second sample
     * <strong><code> var1</code></strong> is the variance of the first sample;
     * <strong><code> var2</code></strong> is the variance of the second sample
     * </p><p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The datasets described by the two Univariates must each contain
     * at least 2 observations.
     * </li></ul>
     *
     * @param sampleStats1 StatisticalSummary describing data from the first sample
     * @param sampleStats2 StatisticalSummary describing data from the second sample
     * @return t statistic
     * @throws NullArgumentException if the sample statistics are <code>null</code>
     * @throws MathIllegalArgumentException if the number of samples is &lt; 2
     */
    public static double t(final StatisticalSummary sampleStats1,
                           final StatisticalSummary sampleStats2)
        throws MathIllegalArgumentException, NullArgumentException {
        return T_TEST.t(sampleStats1, sampleStats2);
    }

    /**
     * Performs a <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm">
     * two-sided t-test</a> evaluating the null hypothesis that the mean of the population from
     * which <code>sample</code> is drawn equals <code>mu</code>.
     * <p>
     * Returns <code>true</code> iff the null hypothesis can be
     * rejected with confidence <code>1 - alpha</code>.  To
     * perform a 1-sided test, use <code>alpha * 2</code></p>
     * <p>
     * <strong>Examples:</strong></p><ol>
     * <li>To test the (2-sided) hypothesis <code>sample mean = mu </code> at
     * the 95% level, use <br><code>tTest(mu, sample, 0.05) </code>
     * </li>
     * <li>To test the (one-sided) hypothesis <code> sample mean &lt; mu </code>
     * at the 99% level, first verify that the measured sample mean is less
     * than <code>mu</code> and then use
     * <br><code>tTest(mu, sample, 0.02) </code>
     * </li></ol>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the one-sample
     * parametric t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/sg_glos.html#one-sample">here</a>
     * </p><p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The observed array length must be at least 2.
     * </li></ul>
     *
     * @param mu constant value to compare sample mean against
     * @param sample array of sample data values
     * @param alpha significance level of the test
     * @return p-value
     * @throws NullArgumentException if the sample array is <code>null</code>
     * @throws MathIllegalArgumentException if the length of the array is &lt; 2
     * @throws MathIllegalArgumentException if <code>alpha</code> is not in the range (0, 0.5]
     * @throws MathIllegalStateException if an error computing the p-value
     */
    public static boolean tTest(final double mu, final double[] sample, final double alpha)
        throws MathIllegalArgumentException, NullArgumentException, MathIllegalStateException {
        return T_TEST.tTest(mu, sample, alpha);
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <i>p-value</i>, associated with a one-sample, two-tailed t-test
     * comparing the mean of the input array with the constant <code>mu</code>.
     * <p>
     * The number returned is the smallest significance level
     * at which one can reject the null hypothesis that the mean equals
     * <code>mu</code> in favor of the two-sided alternative that the mean
     * is different from <code>mu</code>. For a one-sided test, divide the
     * returned value by 2.</p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">here</a>
     * </p><p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The observed array length must be at least 2.
     * </li></ul>
     *
     * @param mu constant value to compare sample mean against
     * @param sample array of sample data values
     * @return p-value
     * @throws NullArgumentException if the sample array is <code>null</code>
     * @throws MathIllegalArgumentException if the length of the array is &lt; 2
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static double tTest(final double mu, final double[] sample)
        throws MathIllegalArgumentException, NullArgumentException,
        MathIllegalStateException {
        return T_TEST.tTest(mu, sample);
    }

    /**
     * Performs a <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm">
     * two-sided t-test</a> evaluating the null hypothesis that the mean of the
     * population from which the dataset described by <code>stats</code> is
     * drawn equals <code>mu</code>.
     * <p>
     * Returns <code>true</code> iff the null hypothesis can be rejected with
     * confidence <code>1 - alpha</code>.  To  perform a 1-sided test, use
     * <code>alpha * 2.</code></p>
     * <p>
     * <strong>Examples:</strong></p><ol>
     * <li>To test the (2-sided) hypothesis <code>sample mean = mu </code> at
     * the 95% level, use <br><code>tTest(mu, sampleStats, 0.05) </code>
     * </li>
     * <li>To test the (one-sided) hypothesis <code> sample mean &lt; mu </code>
     * at the 99% level, first verify that the measured sample mean is less
     * than <code>mu</code> and then use
     * <br><code>tTest(mu, sampleStats, 0.02) </code>
     * </li></ol>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the one-sample
     * parametric t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/sg_glos.html#one-sample">here</a>
     * </p><p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The sample must include at least 2 observations.
     * </li></ul>
     *
     * @param mu constant value to compare sample mean against
     * @param sampleStats StatisticalSummary describing sample data values
     * @param alpha significance level of the test
     * @return p-value
     * @throws NullArgumentException if <code>sampleStats</code> is <code>null</code>
     * @throws MathIllegalArgumentException if the number of samples is &lt; 2
     * @throws MathIllegalArgumentException if <code>alpha</code> is not in the range (0, 0.5]
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static boolean tTest(final double mu, final StatisticalSummary sampleStats,
                                final double alpha)
        throws MathIllegalArgumentException, NullArgumentException, MathIllegalStateException {
        return T_TEST.tTest(mu, sampleStats, alpha);
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <i>p-value</i>, associated with a one-sample, two-tailed t-test
     * comparing the mean of the dataset described by <code>sampleStats</code>
     * with the constant <code>mu</code>.
     * <p>
     * The number returned is the smallest significance level
     * at which one can reject the null hypothesis that the mean equals
     * <code>mu</code> in favor of the two-sided alternative that the mean
     * is different from <code>mu</code>. For a one-sided test, divide the
     * returned value by 2.</p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The sample must contain at least 2 observations.
     * </li></ul>
     *
     * @param mu constant value to compare sample mean against
     * @param sampleStats StatisticalSummary describing sample data
     * @return p-value
     * @throws NullArgumentException if <code>sampleStats</code> is <code>null</code>
     * @throws MathIllegalArgumentException if the number of samples is &lt; 2
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static double tTest(final double mu, final StatisticalSummary sampleStats)
        throws MathIllegalArgumentException, NullArgumentException,
        MathIllegalStateException {
        return T_TEST.tTest(mu, sampleStats);
    }

    /**
     * Performs a
     * <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm">
     * two-sided t-test</a> evaluating the null hypothesis that <code>sample1</code>
     * and <code>sample2</code> are drawn from populations with the same mean,
     * with significance level <code>alpha</code>.  This test does not assume
     * that the subpopulation variances are equal.  To perform the test assuming
     * equal variances, use
     * {@link #homoscedasticTTest(double[], double[], double)}.
     * <p>
     * Returns <code>true</code> iff the null hypothesis that the means are
     * equal can be rejected with confidence <code>1 - alpha</code>.  To
     * perform a 1-sided test, use <code>alpha * 2</code></p>
     * <p>
     * See {@link #t(double[], double[])} for the formula used to compute the
     * t-statistic.  Degrees of freedom are approximated using the
     * <a href="http://www.itl.nist.gov/div898/handbook/prc/section3/prc31.htm">
     * Welch-Satterthwaite approximation.</a></p>
     * <p>
     * <strong>Examples:</strong></p><ol>
     * <li>To test the (2-sided) hypothesis <code>mean 1 = mean 2 </code> at
     * the 95% level,  use
     * <br><code>tTest(sample1, sample2, 0.05). </code>
     * </li>
     * <li>To test the (one-sided) hypothesis <code> mean 1 &lt; mean 2 </code>,
     * at the 99% level, first verify that the measured  mean of <code>sample 1</code>
     * is less than the mean of <code>sample 2</code> and then use
     * <br><code>tTest(sample1, sample2, 0.02) </code>
     * </li></ol>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The observed array lengths must both be at least 2.
     * </li>
     * <li> <code> 0 &lt; alpha &lt; 0.5 </code>
     * </li></ul>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @param alpha significance level of the test
     * @return true if the null hypothesis can be rejected with
     * confidence 1 - alpha
     * @throws NullArgumentException if the arrays are <code>null</code>
     * @throws MathIllegalArgumentException if the length of the arrays is &lt; 2
     * @throws MathIllegalArgumentException if <code>alpha</code> is not in the range (0, 0.5]
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static boolean tTest(final double[] sample1, final double[] sample2,
                                final double alpha)
        throws MathIllegalArgumentException, NullArgumentException, MathIllegalStateException {
        return T_TEST.tTest(sample1, sample2, alpha);
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <i>p-value</i>, associated with a two-sample, two-tailed t-test
     * comparing the means of the input arrays.
     * <p>
     * The number returned is the smallest significance level
     * at which one can reject the null hypothesis that the two means are
     * equal in favor of the two-sided alternative that they are different.
     * For a one-sided test, divide the returned value by 2.</p>
     * <p>
     * The test does not assume that the underlying popuation variances are
     * equal  and it uses approximated degrees of freedom computed from the
     * sample data to compute the p-value.  The t-statistic used is as defined in
     * {@link #t(double[], double[])} and the Welch-Satterthwaite approximation
     * to the degrees of freedom is used,
     * as described
     * <a href="http://www.itl.nist.gov/div898/handbook/prc/section3/prc31.htm">
     * here.</a>  To perform the test under the assumption of equal subpopulation
     * variances, use {@link #homoscedasticTTest(double[], double[])}.</p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the p-value depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The observed array lengths must both be at least 2.
     * </li></ul>
     *
     * @param sample1 array of sample data values
     * @param sample2 array of sample data values
     * @return p-value for t-test
     * @throws NullArgumentException if the arrays are <code>null</code>
     * @throws MathIllegalArgumentException if the length of the arrays is &lt; 2
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static double tTest(final double[] sample1, final double[] sample2)
        throws MathIllegalArgumentException, NullArgumentException,
        MathIllegalStateException {
        return T_TEST.tTest(sample1, sample2);
    }

    /**
     * Performs a
     * <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm">
     * two-sided t-test</a> evaluating the null hypothesis that
     * <code>sampleStats1</code> and <code>sampleStats2</code> describe
     * datasets drawn from populations with the same mean, with significance
     * level <code>alpha</code>.   This test does not assume that the
     * subpopulation variances are equal.  To perform the test under the equal
     * variances assumption, use
     * {@link #homoscedasticTTest(StatisticalSummary, StatisticalSummary)}.
     * <p>
     * Returns <code>true</code> iff the null hypothesis that the means are
     * equal can be rejected with confidence <code>1 - alpha</code>.  To
     * perform a 1-sided test, use <code>alpha * 2</code></p>
     * <p>
     * See {@link #t(double[], double[])} for the formula used to compute the
     * t-statistic.  Degrees of freedom are approximated using the
     * <a href="http://www.itl.nist.gov/div898/handbook/prc/section3/prc31.htm">
     * Welch-Satterthwaite approximation.</a></p>
     * <p>
     * <strong>Examples:</strong></p><ol>
     * <li>To test the (2-sided) hypothesis <code>mean 1 = mean 2 </code> at
     * the 95%, use
     * <br><code>tTest(sampleStats1, sampleStats2, 0.05) </code>
     * </li>
     * <li>To test the (one-sided) hypothesis <code> mean 1 &lt; mean 2 </code>
     * at the 99% level,  first verify that the measured mean of
     * <code>sample 1</code> is less than  the mean of <code>sample 2</code>
     * and then use
     * <br><code>tTest(sampleStats1, sampleStats2, 0.02) </code>
     * </li></ol>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the test depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The datasets described by the two Univariates must each contain
     * at least 2 observations.
     * </li>
     * <li> <code> 0 &lt; alpha &lt; 0.5 </code>
     * </li></ul>
     *
     * @param sampleStats1 StatisticalSummary describing sample data values
     * @param sampleStats2 StatisticalSummary describing sample data values
     * @param alpha significance level of the test
     * @return true if the null hypothesis can be rejected with
     * confidence 1 - alpha
     * @throws NullArgumentException if the sample statistics are <code>null</code>
     * @throws MathIllegalArgumentException if the number of samples is &lt; 2
     * @throws MathIllegalArgumentException if <code>alpha</code> is not in the range (0, 0.5]
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static boolean tTest(final StatisticalSummary sampleStats1,
                                final StatisticalSummary sampleStats2,
                                final double alpha)
        throws MathIllegalArgumentException, NullArgumentException, MathIllegalStateException {
        return T_TEST.tTest(sampleStats1, sampleStats2, alpha);
    }

    /**
     * Returns the <i>observed significance level</i>, or
     * <i>p-value</i>, associated with a two-sample, two-tailed t-test
     * comparing the means of the datasets described by two StatisticalSummary
     * instances.
     * <p>
     * The number returned is the smallest significance level
     * at which one can reject the null hypothesis that the two means are
     * equal in favor of the two-sided alternative that they are different.
     * For a one-sided test, divide the returned value by 2.</p>
     * <p>
     * The test does not assume that the underlying population variances are
     * equal  and it uses approximated degrees of freedom computed from the
     * sample data to compute the p-value.   To perform the test assuming
     * equal variances, use
     * {@link #homoscedasticTTest(StatisticalSummary, StatisticalSummary)}.</p>
     * <p>
     * <strong>Usage Note:</strong><br>
     * The validity of the p-value depends on the assumptions of the parametric
     * t-test procedure, as discussed
     * <a href="http://www.basic.nwu.edu/statguidefiles/ttest_unpaired_ass_viol.html">
     * here</a></p>
     * <p>
     * <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The datasets described by the two Univariates must each contain
     * at least 2 observations.
     * </li></ul>
     *
     * @param sampleStats1  StatisticalSummary describing data from the first sample
     * @param sampleStats2  StatisticalSummary describing data from the second sample
     * @return p-value for t-test
     * @throws NullArgumentException if the sample statistics are <code>null</code>
     * @throws MathIllegalArgumentException if the number of samples is &lt; 2
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static double tTest(final StatisticalSummary sampleStats1,
                               final StatisticalSummary sampleStats2)
        throws MathIllegalArgumentException, NullArgumentException,
        MathIllegalStateException {
        return T_TEST.tTest(sampleStats1, sampleStats2);
    }

    /**
     * Computes the <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda35f.htm">
     * Chi-Square statistic</a> comparing <code>observed</code> and <code>expected</code>
     * frequency counts.
     * <p>
     * This statistic can be used to perform a Chi-Square test evaluating the null
     * hypothesis that the observed counts follow the expected distribution.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>Expected counts must all be positive.</li>
     * <li>Observed counts must all be &ge; 0.</li>
     * <li>The observed and expected arrays must have the same length and
     * their common length must be at least 2.</li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>IllegalArgumentException</code> is thrown.
     * <p>
     * <strong>Note: </strong>This implementation rescales the
     * <code>expected</code> array if necessary to ensure that the sum of the
     * expected and observed counts are equal.
     *
     * @param observed array of observed frequency counts
     * @param expected array of expected frequency counts
     * @return chiSquare test statistic
     * @throws MathIllegalArgumentException if <code>observed</code> has negative entries
     * @throws MathIllegalArgumentException if <code>expected</code> has entries that are
     * not strictly positive
     * @throws MathIllegalArgumentException if the arrays length is less than 2
     */
    public static double chiSquare(final double[] expected, final long[] observed)
        throws MathIllegalArgumentException {
        return CHI_SQUARE_TEST.chiSquare(expected, observed);
    }

    /**
     * Computes the Chi-Square statistic associated with a
     * <a href="http://www.itl.nist.gov/div898/handbook/prc/section4/prc45.htm">
     * chi-square test of independence</a> based on the input <code>counts</code>
     * array, viewed as a two-way table.
     * <p>
     * The rows of the 2-way table are
     * <code>count[0], ... , count[count.length - 1] </code>
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>All counts must be &ge; 0.</li>
     * <li>The count array must be rectangular (i.e. all count[i] subarrays
     * must have the same length).</li>
     * <li>The 2-way table represented by <code>counts</code> must have at
     * least 2 columns and at least 2 rows.</li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>IllegalArgumentException</code> is thrown.
     *
     * @param counts array representation of 2-way table
     * @return chiSquare test statistic
     * @throws NullArgumentException if the array is null
     * @throws MathIllegalArgumentException if the array is not rectangular
     * @throws MathIllegalArgumentException if {@code counts} has negative entries
     */
    public static double chiSquare(final long[][] counts)
        throws MathIllegalArgumentException, NullArgumentException {
        return CHI_SQUARE_TEST.chiSquare(counts);
    }

    /**
     * Performs a <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda35f.htm">
     * Chi-square goodness of fit test</a> evaluating the null hypothesis that the
     * observed counts conform to the frequency distribution described by the expected
     * counts, with significance level <code>alpha</code>.  Returns true iff the null
     * hypothesis can be rejected with 100 * (1 - alpha) percent confidence.
     * <p>
     * <strong>Example:</strong><br>
     * To test the hypothesis that <code>observed</code> follows
     * <code>expected</code> at the 99% level, use
     * <code>chiSquareTest(expected, observed, 0.01)</code>
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>Expected counts must all be positive.</li>
     * <li>Observed counts must all be &ge; 0.</li>
     * <li>The observed and expected arrays must have the same length and
     * their common length must be at least 2.</li>
     * <li><code> 0 &lt; alpha &lt; 0.5</code></li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>IllegalArgumentException</code> is thrown.
     * <p>
     * <strong>Note: </strong>This implementation rescales the
     * <code>expected</code> array if necessary to ensure that the sum of the
     * expected and observed counts are equal.
     *
     * @param observed array of observed frequency counts
     * @param expected array of expected frequency counts
     * @param alpha significance level of the test
     * @return true iff null hypothesis can be rejected with confidence
     * 1 - alpha
     * @throws MathIllegalArgumentException if <code>observed</code> has negative entries
     * @throws MathIllegalArgumentException if <code>expected</code> has entries that are
     * not strictly positive
     * @throws MathIllegalArgumentException if the arrays length is less than 2
     * @throws MathIllegalArgumentException if <code>alpha</code> is not in the range (0, 0.5]
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static boolean chiSquareTest(final double[] expected, final long[] observed,
                                        final double alpha)
        throws MathIllegalArgumentException, MathIllegalStateException {
        return CHI_SQUARE_TEST.chiSquareTest(expected, observed, alpha);
    }

    /**
     * Returns the <i>observed significance level</i>, or <a href=
     * "http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue">
     * p-value</a>, associated with a
     * <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda35f.htm">
     * Chi-square goodness of fit test</a> comparing the <code>observed</code>
     * frequency counts to those in the <code>expected</code> array.
     * <p>
     * The number returned is the smallest significance level at which one can reject
     * the null hypothesis that the observed counts conform to the frequency distribution
     * described by the expected counts.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>Expected counts must all be positive.</li>
     * <li>Observed counts must all be &ge; 0.</li>
     * <li>The observed and expected arrays must have the same length and
     * their common length must be at least 2.</li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>IllegalArgumentException</code> is thrown.
     * <p>
     * <strong>Note: </strong>This implementation rescales the
     * <code>expected</code> array if necessary to ensure that the sum of the
     * expected and observed counts are equal.
     *
     * @param observed array of observed frequency counts
     * @param expected array of expected frequency counts
     * @return p-value
     * @throws MathIllegalArgumentException if <code>observed</code> has negative entries
     * @throws MathIllegalArgumentException if <code>expected</code> has entries that are
     * not strictly positive
     * @throws MathIllegalArgumentException if the arrays length is less than 2
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static double chiSquareTest(final double[] expected, final long[] observed)
        throws MathIllegalArgumentException, MathIllegalStateException {
        return CHI_SQUARE_TEST.chiSquareTest(expected, observed);
    }

    /**
     * Performs a <a href="http://www.itl.nist.gov/div898/handbook/prc/section4/prc45.htm">
     * chi-square test of independence</a> evaluating the null hypothesis that the
     * classifications represented by the counts in the columns of the input 2-way table
     * are independent of the rows, with significance level <code>alpha</code>.
     * Returns true iff the null hypothesis can be rejected with 100 * (1 - alpha) percent
     * confidence.
     * <p>
     * The rows of the 2-way table are
     * <code>count[0], ... , count[count.length - 1] </code>
     * <p>
     * <strong>Example:</strong><br>
     * To test the null hypothesis that the counts in
     * <code>count[0], ... , count[count.length - 1] </code>
     * all correspond to the same underlying probability distribution at the 99% level,
     * use <code>chiSquareTest(counts, 0.01)</code>.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>All counts must be &ge; 0.</li>
     * <li>The count array must be rectangular (i.e. all count[i] subarrays must have the
     * same length).</li>
     * <li>The 2-way table represented by <code>counts</code> must have at least 2 columns and
     * at least 2 rows.</li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>IllegalArgumentException</code> is thrown.
     *
     * @param counts array representation of 2-way table
     * @param alpha significance level of the test
     * @return true iff null hypothesis can be rejected with confidence
     * 1 - alpha
     * @throws NullArgumentException if the array is null
     * @throws MathIllegalArgumentException if the array is not rectangular
     * @throws MathIllegalArgumentException if {@code counts} has any negative entries
     * @throws MathIllegalArgumentException if <code>alpha</code> is not in the range (0, 0.5]
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static boolean chiSquareTest(final long[][] counts, final double alpha)
        throws MathIllegalArgumentException, NullArgumentException, MathIllegalStateException {
        return CHI_SQUARE_TEST.chiSquareTest(counts, alpha);
    }

    /**
     * Returns the <i>observed significance level</i>, or <a href=
     * "http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue">
     * p-value</a>, associated with a
     * <a href="http://www.itl.nist.gov/div898/handbook/prc/section4/prc45.htm">
     * chi-square test of independence</a> based on the input <code>counts</code>
     * array, viewed as a two-way table.
     * <p>
     * The rows of the 2-way table are
     * <code>count[0], ... , count[count.length - 1] </code>
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>All counts must be &ge; 0.</li>
     * <li>The count array must be rectangular (i.e. all count[i] subarrays must have
     * the same length).</li>
     * <li>The 2-way table represented by <code>counts</code> must have at least 2
     * columns and at least 2 rows.</li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>IllegalArgumentException</code> is thrown.
     *
     * @param counts array representation of 2-way table
     * @return p-value
     * @throws NullArgumentException if the array is null
     * @throws MathIllegalArgumentException if the array is not rectangular
     * @throws MathIllegalArgumentException if {@code counts} has negative entries
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static double chiSquareTest(final long[][] counts)
        throws MathIllegalArgumentException, NullArgumentException, MathIllegalStateException {
        return CHI_SQUARE_TEST.chiSquareTest(counts);
    }

    /**
     * Computes a
     * <a href="http://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/chi2samp.htm">
     * Chi-Square two sample test statistic</a> comparing bin frequency counts
     * in <code>observed1</code> and <code>observed2</code>.
     * <p>
     * The sums of frequency counts in the two samples are not required to be the
     * same. The formula used to compute the test statistic is
     * </p>
     * <code>
     * &sum;[(K * observed1[i] - observed2[i]/K)<sup>2</sup> / (observed1[i] + observed2[i])]
     * </code>
     * <p>
     * where
     * </p>
     * <code>K = √[&sum;(observed2 / &sum;(observed1)]</code>
     * <p>
     * This statistic can be used to perform a Chi-Square test evaluating the
     * null hypothesis that both observed counts follow the same distribution.
     * </p>
     * <p><strong>Preconditions</strong>:</p>
     * <ul>
     * <li>Observed counts must be non-negative.</li>
     * <li>Observed counts for a specific bin must not both be zero.</li>
     * <li>Observed counts for a specific sample must not all be 0.</li>
     * <li>The arrays <code>observed1</code> and <code>observed2</code> must have
     * the same length and their common length must be at least 2.</li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>IllegalArgumentException</code> is thrown.
     * </p>
     *
     * @param observed1 array of observed frequency counts of the first data set
     * @param observed2 array of observed frequency counts of the second data set
     * @return chiSquare test statistic
     * @throws MathIllegalArgumentException the the length of the arrays does not match
     * @throws MathIllegalArgumentException if any entries in <code>observed1</code> or
     * <code>observed2</code> are negative
     * @throws MathIllegalArgumentException if either all counts of <code>observed1</code> or
     * <code>observed2</code> are zero, or if the count at some index is zero
     * for both arrays
     */
    public static double chiSquareDataSetsComparison(final long[] observed1,
                                                     final long[] observed2)
        throws MathIllegalArgumentException {
        return CHI_SQUARE_TEST.chiSquareDataSetsComparison(observed1, observed2);
    }

    /**
     * Returns the <i>observed significance level</i>, or <a href=
     * "http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue">
     * p-value</a>, associated with a Chi-Square two sample test comparing
     * bin frequency counts in <code>observed1</code> and
     * <code>observed2</code>.
     * <p>
     * The number returned is the smallest significance level at which one
     * can reject the null hypothesis that the observed counts conform to the
     * same distribution.
     * <p>
     * See {@link #chiSquareDataSetsComparison(long[], long[])} for details
     * on the formula used to compute the test statistic. The degrees of
     * of freedom used to perform the test is one less than the common length
     * of the input observed count arrays.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>Observed counts must be non-negative.</li>
     * <li>Observed counts for a specific bin must not both be zero.</li>
     * <li>Observed counts for a specific sample must not all be 0.</li>
     * <li>The arrays <code>observed1</code> and <code>observed2</code> must
     * have the same length and their common length must be at least 2.</li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>IllegalArgumentException</code> is thrown.
     *
     * @param observed1 array of observed frequency counts of the first data set
     * @param observed2 array of observed frequency counts of the second data set
     * @return p-value
     * @throws MathIllegalArgumentException the the length of the arrays does not match
     * @throws MathIllegalArgumentException if any entries in <code>observed1</code> or
     * <code>observed2</code> are negative
     * @throws MathIllegalArgumentException if either all counts of <code>observed1</code> or
     * <code>observed2</code> are zero, or if the count at the same index is zero
     * for both arrays
     * @throws MathIllegalStateException if an error occurs computing the p-value
     */
    public static double chiSquareTestDataSetsComparison(final long[] observed1,
                                                         final long[] observed2)
        throws MathIllegalArgumentException,
        MathIllegalStateException {
        return CHI_SQUARE_TEST.chiSquareTestDataSetsComparison(observed1, observed2);
    }

    /**
     * Performs a Chi-Square two sample test comparing two binned data
     * sets. The test evaluates the null hypothesis that the two lists of
     * observed counts conform to the same frequency distribution, with
     * significance level <code>alpha</code>.  Returns true iff the null
     * hypothesis can be rejected with 100 * (1 - alpha) percent confidence.
     * <p>
     * See {@link #chiSquareDataSetsComparison(long[], long[])} for
     * details on the formula used to compute the Chisquare statistic used
     * in the test. The degrees of of freedom used to perform the test is
     * one less than the common length of the input observed count arrays.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>Observed counts must be non-negative.</li>
     * <li>Observed counts for a specific bin must not both be zero.</li>
     * <li>Observed counts for a specific sample must not all be 0.</li>
     * <li>The arrays <code>observed1</code> and <code>observed2</code> must
     * have the same length and their common length must be at least 2.</li>
     * <li><code> 0 &lt; alpha &lt; 0.5</code></li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, an
     * <code>IllegalArgumentException</code> is thrown.
     *
     * @param observed1 array of observed frequency counts of the first data set
     * @param observed2 array of observed frequency counts of the second data set
     * @param alpha significance level of the test
     * @return true iff null hypothesis can be rejected with confidence
     * 1 - alpha
     * @throws MathIllegalArgumentException the the length of the arrays does not match
     * @throws MathIllegalArgumentException if any entries in <code>observed1</code> or
     * <code>observed2</code> are negative
     * @throws MathIllegalArgumentException if either all counts of <code>observed1</code> or
     * <code>observed2</code> are zero, or if the count at the same index is zero
     * for both arrays
     * @throws MathIllegalArgumentException if <code>alpha</code> is not in the range (0, 0.5]
     * @throws MathIllegalStateException if an error occurs performing the test
     */
    public static boolean chiSquareTestDataSetsComparison(final long[] observed1,
                                                          final long[] observed2,
                                                          final double alpha)
        throws MathIllegalArgumentException, MathIllegalStateException {
        return CHI_SQUARE_TEST.chiSquareTestDataSetsComparison(observed1, observed2, alpha);
    }

    /**
     * Computes the ANOVA F-value for a collection of <code>double[]</code>
     * arrays.
     *
     * <p><strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The categoryData <code>Collection</code> must contain
     * <code>double[]</code> arrays.</li>
     * <li> There must be at least two <code>double[]</code> arrays in the
     * <code>categoryData</code> collection and each of these arrays must
     * contain at least two values.</li></ul>
     * <p>
     * This implementation computes the F statistic using the definitional
     * formula</p>
     * <pre>
     *   F = msbg/mswg</pre>
     * <p>where</p>
     * <pre>
     *  msbg = between group mean square
     *  mswg = within group mean square</pre>
     * <p>
     * are as defined <a href="http://faculty.vassar.edu/lowry/ch13pt1.html">
     * here</a></p>
     *
     * @param categoryData <code>Collection</code> of <code>double[]</code>
     * arrays each containing data for one category
     * @return Fvalue
     * @throws NullArgumentException if <code>categoryData</code> is <code>null</code>
     * @throws MathIllegalArgumentException if the length of the <code>categoryData</code>
     * array is less than 2 or a contained <code>double[]</code> array does not have
     * at least two values
     */
    public static double oneWayAnovaFValue(final Collection<double[]> categoryData)
        throws MathIllegalArgumentException, NullArgumentException {
        return ONE_WAY_ANANOVA.anovaFValue(categoryData);
    }

    /**
     * Computes the ANOVA P-value for a collection of <code>double[]</code>
     * arrays.
     *
     * <p><strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The categoryData <code>Collection</code> must contain
     * <code>double[]</code> arrays.</li>
     * <li> There must be at least two <code>double[]</code> arrays in the
     * <code>categoryData</code> collection and each of these arrays must
     * contain at least two values.</li></ul>
     * <p>
     * This implementation uses the
     * {@link org.hipparchus.distribution.continuous.FDistribution
     * Hipparchus F Distribution implementation} to estimate the exact
     * p-value, using the formula</p>
     * <pre>
     *   p = 1 - cumulativeProbability(F)</pre>
     * <p>
     * where <code>F</code> is the F value and <code>cumulativeProbability</code>
     * is the Hipparchus implementation of the F distribution.</p>
     *
     * @param categoryData <code>Collection</code> of <code>double[]</code>
     * arrays each containing data for one category
     * @return Pvalue
     * @throws NullArgumentException if <code>categoryData</code> is <code>null</code>
     * @throws MathIllegalArgumentException if the length of the <code>categoryData</code>
     * array is less than 2 or a contained <code>double[]</code> array does not have
     * at least two values
     * @throws MathIllegalStateException if the p-value can not be computed due to a convergence error
     * @throws MathIllegalStateException if the maximum number of iterations is exceeded
     */
    public static double oneWayAnovaPValue(final Collection<double[]> categoryData)
        throws MathIllegalArgumentException, NullArgumentException,
        MathIllegalStateException {
        return ONE_WAY_ANANOVA.anovaPValue(categoryData);
    }

    /**
     * Performs an ANOVA test, evaluating the null hypothesis that there
     * is no difference among the means of the data categories.
     *
     * <p><strong>Preconditions</strong>:</p>
     * <ul>
     * <li>The categoryData <code>Collection</code> must contain
     * <code>double[]</code> arrays.</li>
     * <li> There must be at least two <code>double[]</code> arrays in the
     * <code>categoryData</code> collection and each of these arrays must
     * contain at least two values.</li>
     * <li>alpha must be strictly greater than 0 and less than or equal to 0.5.
     * </li></ul>
     * <p>
     * This implementation uses the
     * {@link org.hipparchus.distribution.continuous.FDistribution
     * Hipparchus F Distribution implementation} to estimate the exact
     * p-value, using the formula</p><pre>
     *   p = 1 - cumulativeProbability(F)</pre>
     * <p>where <code>F</code> is the F value and <code>cumulativeProbability</code>
     * is the Hipparchus implementation of the F distribution.</p>
     * <p>True is returned iff the estimated p-value is less than alpha.</p>
     *
     * @param categoryData <code>Collection</code> of <code>double[]</code>
     * arrays each containing data for one category
     * @param alpha significance level of the test
     * @return true if the null hypothesis can be rejected with
     * confidence 1 - alpha
     * @throws NullArgumentException if <code>categoryData</code> is <code>null</code>
     * @throws MathIllegalArgumentException if the length of the <code>categoryData</code>
     * array is less than 2 or a contained <code>double[]</code> array does not have
     * at least two values
     * @throws MathIllegalArgumentException if <code>alpha</code> is not in the range (0, 0.5]
     * @throws MathIllegalStateException if the p-value can not be computed due to a convergence error
     * @throws MathIllegalStateException if the maximum number of iterations is exceeded
     */
    public static boolean oneWayAnovaTest(final Collection<double[]> categoryData,
                                          final double alpha)
        throws MathIllegalArgumentException, NullArgumentException, MathIllegalStateException {
        return ONE_WAY_ANANOVA.anovaTest(categoryData, alpha);
    }

    /**
     * Computes the <a href="http://en.wikipedia.org/wiki/G-test">G statistic
     * for Goodness of Fit</a> comparing {@code observed} and {@code expected}
     * frequency counts.
     * <p>
     * This statistic can be used to perform a G test (Log-Likelihood Ratio
     * Test) evaluating the null hypothesis that the observed counts follow the
     * expected distribution.
     * <p>
     * <strong>Preconditions</strong>:
     * <ul>
     * <li>Expected counts must all be positive.</li>
     * <li>Observed counts must all be &ge; 0.</li>
     * <li>The observed and expected arrays must have the same length and their
     * common length must be at least 2. </li>
     * </ul>
     * <p>
     * If any of the preconditions are not met, a
     * {@code MathIllegalArgumentException} is thrown.
     * <p>
     * <strong>Note:</strong>This implementation rescales the
     * {@code expected} array if necessary to ensure that the sum of the
     * expected and observed counts are equal.
     *
     * @param observed array of observed frequency counts
     * @param expected array of expected frequency counts
     * @return G-Test statistic
     * @throws MathIllegalArgumentException if {@code observed} has negative entries
     * @throws MathIllegalArgumentException if {@code expected} has entries that
     * are not strictly positive
     * @throws MathIllegalArgumentException if the array lengths do not match or
     * are less than 2.
     */
    public static double g(final double[] expected, final long[] observed)
        throws MathIllegalArgumentException {
        return G_TEST.g(expected, observed);
    }

    /**
     * Returns the <i>observed significance level</i>, or <a href=
     * "http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue"> p-value</a>,
     * associated with a G-Test for goodness of fit comparing the
     * {@code observed} frequency counts to those in the {@code expected} array.
     *
     * <p>The number returned is the smallest significance level at which one
     * can reject the null hypothesis that the observed counts conform to the
     * frequency distribution described by the expected counts.</p>
     *
     * <p>The probability returned is the tail probability beyond
     * {@link #g(double[], long[]) g(expected, observed)}
     * in the ChiSquare distribution with degrees of freedom one less than the
     * common length of {@code expected} and {@code observed}.</p>
     *
     * <p> <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>Expected counts must all be positive. </li>
     * <li>Observed counts must all be &ge; 0. </li>
     * <li>The observed and expected arrays must have the
     * same length and their common length must be at least 2.</li>
     * </ul>
     *
     * <p>If any of the preconditions are not met, a
     * {@code MathIllegalArgumentException} is thrown.</p>
     *
     * <p><strong>Note:</strong>This implementation rescales the
     * {@code expected} array if necessary to ensure that the sum of the
     *  expected and observed counts are equal.</p>
     *
     * @param observed array of observed frequency counts
     * @param expected array of expected frequency counts
     * @return p-value
     * @throws MathIllegalArgumentException if {@code observed} has negative entries
     * @throws MathIllegalArgumentException if {@code expected} has entries that
     * are not strictly positive
     * @throws MathIllegalArgumentException if the array lengths do not match or
     * are less than 2.
     * @throws MathIllegalStateException if an error occurs computing the
     * p-value.
     */
    public static double gTest(final double[] expected, final long[] observed)
        throws MathIllegalArgumentException, MathIllegalStateException {
        return G_TEST.gTest(expected, observed);
    }

    /**
     * Returns the intrinsic (Hardy-Weinberg proportions) p-Value, as described
     * in p64-69 of McDonald, J.H. 2009. Handbook of Biological Statistics
     * (2nd ed.). Sparky House Publishing, Baltimore, Maryland.
     *
     * <p> The probability returned is the tail probability beyond
     * {@link #g(double[], long[]) g(expected, observed)}
     * in the ChiSquare distribution with degrees of freedom two less than the
     * common length of {@code expected} and {@code observed}.</p>
     *
     * @param observed array of observed frequency counts
     * @param expected array of expected frequency counts
     * @return p-value
     * @throws MathIllegalArgumentException if {@code observed} has negative entries
     * @throws MathIllegalArgumentException {@code expected} has entries that are
     * not strictly positive
     * @throws MathIllegalArgumentException if the array lengths do not match or
     * are less than 2.
     * @throws MathIllegalStateException if an error occurs computing the
     * p-value.
     */
    public static double gTestIntrinsic(final double[] expected, final long[] observed)
        throws MathIllegalArgumentException, MathIllegalStateException {
        return G_TEST.gTestIntrinsic(expected, observed);
    }

    /**
     * Performs a G-Test (Log-Likelihood Ratio Test) for goodness of fit
     * evaluating the null hypothesis that the observed counts conform to the
     * frequency distribution described by the expected counts, with
     * significance level {@code alpha}. Returns true iff the null
     * hypothesis can be rejected with {@code 100 * (1 - alpha)} percent confidence.
     *
     * <p><strong>Example:</strong><br> To test the hypothesis that
     * {@code observed} follows {@code expected} at the 99% level,
     * use </p><p>
     * {@code gTest(expected, observed, 0.01)}</p>
     *
     * <p>Returns true iff {@link #gTest(double[], long[])
     *  gTestGoodnessOfFitPValue(expected, observed)} &gt; alpha</p>
     *
     * <p><strong>Preconditions</strong>:</p>
     * <ul>
     * <li>Expected counts must all be positive. </li>
     * <li>Observed counts must all be &ge; 0. </li>
     * <li>The observed and expected arrays must have the same length and their
     * common length must be at least 2.
     * <li> {@code 0 < alpha < 0.5} </li></ul>
     *
     * <p>If any of the preconditions are not met, a
     * {@code MathIllegalArgumentException} is thrown.</p>
     *
     * <p><strong>Note:</strong>This implementation rescales the
     * {@code expected} array if necessary to ensure that the sum of the
     * expected and observed counts are equal.</p>
     *
     * @param observed array of observed frequency counts
     * @param expected array of expected frequency counts
     * @param alpha significance level of the test
     * @return true iff null hypothesis can be rejected with confidence 1 -
     * alpha
     * @throws MathIllegalArgumentException if {@code observed} has negative entries
     * @throws MathIllegalArgumentException if {@code expected} has entries that
     * are not strictly positive
     * @throws MathIllegalArgumentException if the array lengths do not match or
     * are less than 2.
     * @throws MathIllegalStateException if an error occurs computing the
     * p-value.
     * @throws MathIllegalArgumentException if alpha is not strictly greater than zero
     * and less than or equal to 0.5
     */
    public static boolean gTest(final double[] expected, final long[] observed,
                                final double alpha)
        throws MathIllegalArgumentException, MathIllegalStateException {
        return G_TEST.gTest(expected, observed, alpha);
    }

    /**
     * <p>Computes a G (Log-Likelihood Ratio) two sample test statistic for
     * independence comparing frequency counts in
     * {@code observed1} and {@code observed2}. The sums of frequency
     * counts in the two samples are not required to be the same. The formula
     * used to compute the test statistic is </p>
     *
     * <p>{@code 2 * totalSum * [H(rowSums) + H(colSums) - H(k)]}</p>
     *
     * <p> where {@code H} is the
     * <a href="http://en.wikipedia.org/wiki/Entropy_%28information_theory%29">
     * Shannon Entropy</a> of the random variable formed by viewing the elements
     * of the argument array as incidence counts; <br>
     * {@code k} is a matrix with rows {@code [observed1, observed2]}; <br>
     * {@code rowSums, colSums} are the row/col sums of {@code k}; <br>
     * and {@code totalSum} is the overall sum of all entries in {@code k}.</p>
     *
     * <p>This statistic can be used to perform a G test evaluating the null
     * hypothesis that both observed counts are independent </p>
     *
     * <p> <strong>Preconditions</strong>:</p>
     * <ul>
     * <li>Observed counts must be non-negative. </li>
     * <li>Observed counts for a specific bin must not both be zero. </li>
     * <li>Observed counts for a specific sample must not all be  0. </li>
     * <li>The arrays {@code observed1} and {@code observed2} must have
     * the same length and their common length must be at least 2. </li></ul>
     *
     * <p>If any of the preconditions are not met, a
     * {@code MathIllegalArgumentException} is thrown.</p>
     *
     * @param observed1 array of observed frequency counts of the first data set
     * @param observed2 array of observed frequency counts of the second data
     * set
     * @return G-Test statistic
     * @throws MathIllegalArgumentException the the lengths of the arrays do not
     * match or their common length is less than 2
     * @throws MathIllegalArgumentException if any entry in {@code observed1} or
     * {@code observed2} is negative
     * @throws MathIllegalArgumentException if either all counts of
     * {@code observed1} or {@code observed2} are zero, or if the count
     * at the same index is zero for both arrays.
     */
    public static double gDataSetsComparison(final long[] observed1,
                                                  final long[] observed2)
        throws MathIllegalArgumentException {
        return G_TEST.gDataSetsComparison(observed1, observed2);
    }

    /**
     * Calculates the root log-likelihood ratio for 2 state Datasets. See
     * {@link #gDataSetsComparison(long[], long[] )}.
     *
     * <p>Given two events A and B, let k11 be the number of times both events
     * occur, k12 the incidence of B without A, k21 the count of A without B,
     * and k22 the number of times neither A nor B occurs.  What is returned
     * by this method is </p>
     *
     * <p>{@code (sgn) sqrt(gValueDataSetsComparison({k11, k12}, {k21, k22})}</p>
     *
     * <p>where {@code sgn} is -1 if {@code k11 / (k11 + k12) < k21 / (k21 + k22))};<br>
     * 1 otherwise.</p>
     *
     * <p>Signed root LLR has two advantages over the basic LLR: a) it is positive
     * where k11 is bigger than expected, negative where it is lower b) if there is
     * no difference it is asymptotically normally distributed. This allows one
     * to talk about "number of standard deviations" which is a more common frame
     * of reference than the chi^2 distribution.</p>
     *
     * @param k11 number of times the two events occurred together (AB)
     * @param k12 number of times the second event occurred WITHOUT the
     * first event (notA,B)
     * @param k21 number of times the first event occurred WITHOUT the
     * second event (A, notB)
     * @param k22 number of times something else occurred (i.e. was neither
     * of these events (notA, notB)
     * @return root log-likelihood ratio
     *
     */
    public static double rootLogLikelihoodRatio(final long k11, final long k12, final long k21, final long k22)
        throws MathIllegalArgumentException {
        return G_TEST.rootLogLikelihoodRatio(k11, k12, k21, k22);
    }


    /**
     * <p>Returns the <i>observed significance level</i>, or <a href=
     * "http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue">
     * p-value</a>, associated with a G-Value (Log-Likelihood Ratio) for two
     * sample test comparing bin frequency counts in {@code observed1} and
     * {@code observed2}.</p>
     *
     * <p>The number returned is the smallest significance level at which one
     * can reject the null hypothesis that the observed counts conform to the
     * same distribution. </p>
     *
     * <p>See {@link #gTest(double[], long[])} for details
     * on how the p-value is computed.  The degrees of of freedom used to
     * perform the test is one less than the common length of the input observed
     * count arrays.</p>
     *
     * <p><strong>Preconditions</strong>:</p>
     * <ul> <li>Observed counts must be non-negative. </li>
     * <li>Observed counts for a specific bin must not both be zero. </li>
     * <li>Observed counts for a specific sample must not all be 0. </li>
     * <li>The arrays {@code observed1} and {@code observed2} must
     * have the same length and their common length must be at least 2. </li>
     * </ul>
     * <p> If any of the preconditions are not met, a
     * {@code MathIllegalArgumentException} is thrown.</p>
     *
     * @param observed1 array of observed frequency counts of the first data set
     * @param observed2 array of observed frequency counts of the second data
     * set
     * @return p-value
     * @throws MathIllegalArgumentException the the length of the arrays does not
     * match or their common length is less than 2
     * @throws MathIllegalArgumentException if any of the entries in {@code observed1} or
     * {@code observed2} are negative
     * @throws MathIllegalArgumentException if either all counts of {@code observed1} or
     * {@code observed2} are zero, or if the count at some index is
     * zero for both arrays
     * @throws MathIllegalStateException if an error occurs computing the
     * p-value.
     */
    public static double gTestDataSetsComparison(final long[] observed1,
                                                        final long[] observed2)
        throws MathIllegalArgumentException,
        MathIllegalStateException {
        return G_TEST.gTestDataSetsComparison(observed1, observed2);
    }

    /**
     * <p>Performs a G-Test (Log-Likelihood Ratio Test) comparing two binned
     * data sets. The test evaluates the null hypothesis that the two lists
     * of observed counts conform to the same frequency distribution, with
     * significance level {@code alpha}. Returns true iff the null
     * hypothesis can be rejected  with 100 * (1 - alpha) percent confidence.
     * </p>
     * <p>See {@link #gDataSetsComparison(long[], long[])} for details
     * on the formula used to compute the G (LLR) statistic used in the test and
     * {@link #gTest(double[], long[])} for information on how
     * the observed significance level is computed. The degrees of of freedom used
     * to perform the test is one less than the common length of the input observed
     * count arrays. </p>
     *
     * <p><strong>Preconditions</strong>:</p>
     * <ul>
     * <li>Observed counts must be non-negative. </li>
     * <li>Observed counts for a specific bin must not both be zero. </li>
     * <li>Observed counts for a specific sample must not all be 0. </li>
     * <li>The arrays {@code observed1} and {@code observed2} must
     * have the same length and their common length must be at least 2. </li>
     * <li>{@code 0 < alpha < 0.5} </li></ul>
     *
     * <p>If any of the preconditions are not met, a
     * {@code MathIllegalArgumentException} is thrown.</p>
     *
     * @param observed1 array of observed frequency counts of the first data set
     * @param observed2 array of observed frequency counts of the second data
     * set
     * @param alpha significance level of the test
     * @return true iff null hypothesis can be rejected with confidence 1 -
     * alpha
     * @throws MathIllegalArgumentException the the length of the arrays does not
     * match
     * @throws MathIllegalArgumentException if any of the entries in {@code observed1} or
     * {@code observed2} are negative
     * @throws MathIllegalArgumentException if either all counts of {@code observed1} or
     * {@code observed2} are zero, or if the count at some index is
     * zero for both arrays
     * @throws MathIllegalArgumentException if {@code alpha} is not in the range
     * (0, 0.5]
     * @throws MathIllegalStateException if an error occurs performing the test
     */
    public static boolean gTestDataSetsComparison(final long[] observed1,
                                                  final long[] observed2,
                                                  final double alpha)
        throws MathIllegalArgumentException, MathIllegalStateException {
        return G_TEST.gTestDataSetsComparison(observed1, observed2, alpha);
    }

    /**
     * Computes the one-sample Kolmogorov-Smirnov test statistic, \(D_n=\sup_x |F_n(x)-F(x)|\) where
     * \(F\) is the distribution (cdf) function associated with {@code distribution}, \(n\) is the
     * length of {@code data} and \(F_n\) is the empirical distribution that puts mass \(1/n\) at
     * each of the values in {@code data}.
     *
     * @param dist reference distribution
     * @param data sample being evaluated
     * @return Kolmogorov-Smirnov statistic \(D_n\)
     * @throws MathIllegalArgumentException if {@code data} does not have length at least 2
     * @throws org.hipparchus.exception.NullArgumentException if {@code data} is null
     */
    public static double kolmogorovSmirnovStatistic(RealDistribution dist, double[] data)
            throws MathIllegalArgumentException, NullArgumentException {
        return KS_TEST.kolmogorovSmirnovStatistic(dist, data);
    }

    /**
     * Computes the <i>p-value</i>, or <i>observed significance level</i>, of a one-sample <a
     * href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"> Kolmogorov-Smirnov test</a>
     * evaluating the null hypothesis that {@code data} conforms to {@code distribution}.
     *
     * @param dist reference distribution
     * @param data sample being being evaluated
     * @return the p-value associated with the null hypothesis that {@code data} is a sample from
     *         {@code distribution}
     * @throws MathIllegalArgumentException if {@code data} does not have length at least 2
     * @throws org.hipparchus.exception.NullArgumentException if {@code data} is null
     */
    public static double kolmogorovSmirnovTest(RealDistribution dist, double[] data)
            throws MathIllegalArgumentException, NullArgumentException {
        return KS_TEST.kolmogorovSmirnovTest(dist, data);
    }

    /**
     * Computes the <i>p-value</i>, or <i>observed significance level</i>, of a one-sample <a
     * href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"> Kolmogorov-Smirnov test</a>
     * evaluating the null hypothesis that {@code data} conforms to {@code distribution}. If
     * {@code exact} is true, the distribution used to compute the p-value is computed using
     * extended precision. See {@link KolmogorovSmirnovTest#cdfExact(double, int)}.
     *
     * @param dist reference distribution
     * @param data sample being being evaluated
     * @param strict whether or not to force exact computation of the p-value
     * @return the p-value associated with the null hypothesis that {@code data} is a sample from
     *         {@code distribution}
     * @throws MathIllegalArgumentException if {@code data} does not have length at least 2
     * @throws org.hipparchus.exception.NullArgumentException if {@code data} is null
     */
    public static double kolmogorovSmirnovTest(RealDistribution dist, double[] data, boolean strict)
            throws MathIllegalArgumentException, NullArgumentException {
        return KS_TEST.kolmogorovSmirnovTest(dist, data, strict);
    }

    /**
     * Performs a <a href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"> Kolmogorov-Smirnov
     * test</a> evaluating the null hypothesis that {@code data} conforms to {@code distribution}.
     *
     * @param dist reference distribution
     * @param data sample being being evaluated
     * @param alpha significance level of the test
     * @return true iff the null hypothesis that {@code data} is a sample from {@code distribution}
     *         can be rejected with confidence 1 - {@code alpha}
     * @throws MathIllegalArgumentException if {@code data} does not have length at least 2
     * @throws org.hipparchus.exception.NullArgumentException if {@code data} is null
     */
    public static boolean kolmogorovSmirnovTest(RealDistribution dist, double[] data, double alpha)
            throws MathIllegalArgumentException, NullArgumentException {
        return KS_TEST.kolmogorovSmirnovTest(dist, data, alpha);
    }

    /**
     * Computes the two-sample Kolmogorov-Smirnov test statistic, \(D_{n,m}=\sup_x |F_n(x)-F_m(x)|\)
     * where \(n\) is the length of {@code x}, \(m\) is the length of {@code y}, \(F_n\) is the
     * empirical distribution that puts mass \(1/n\) at each of the values in {@code x} and \(F_m\)
     * is the empirical distribution of the {@code y} values.
     *
     * @param x first sample
     * @param y second sample
     * @return test statistic \(D_{n,m}\) used to evaluate the null hypothesis that {@code x} and
     *         {@code y} represent samples from the same underlying distribution
     * @throws MathIllegalArgumentException if either {@code x} or {@code y} does not have length at
     *         least 2
     * @throws org.hipparchus.exception.NullArgumentException if either {@code x} or {@code y} is null
     */
    public static double kolmogorovSmirnovStatistic(double[] x, double[] y)
            throws MathIllegalArgumentException, NullArgumentException {
        return KS_TEST.kolmogorovSmirnovStatistic(x, y);
    }

    /**
     * Computes the <i>p-value</i>, or <i>observed significance level</i>, of a two-sample <a
     * href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"> Kolmogorov-Smirnov test</a>
     * evaluating the null hypothesis that {@code x} and {@code y} are samples drawn from the same
     * probability distribution. Assumes the strict form of the inequality used to compute the
     * p-value. See {@link KolmogorovSmirnovTest#kolmogorovSmirnovTest(RealDistribution, double[], boolean)}.
     *
     * @param x first sample dataset
     * @param y second sample dataset
     * @return p-value associated with the null hypothesis that {@code x} and {@code y} represent
     *         samples from the same distribution
     * @throws MathIllegalArgumentException if either {@code x} or {@code y} does not have length at
     *         least 2
     * @throws org.hipparchus.exception.NullArgumentException if either {@code x} or {@code y} is null
     */
    public static double kolmogorovSmirnovTest(double[] x, double[] y)
            throws MathIllegalArgumentException, NullArgumentException {
        return KS_TEST.kolmogorovSmirnovTest(x, y);
    }

    /**
     * Computes the <i>p-value</i>, or <i>observed significance level</i>, of a two-sample <a
     * href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"> Kolmogorov-Smirnov test</a>
     * evaluating the null hypothesis that {@code x} and {@code y} are samples drawn from the same
     * probability distribution. Specifically, what is returned is an estimate of the probability
     * that the {@link KolmogorovSmirnovTest#kolmogorovSmirnovStatistic(double[], double[])} associated with a randomly
     * selected partition of the combined sample into subsamples of sizes {@code x.length} and
     * {@code y.length} will strictly exceed (if {@code strict} is {@code true}) or be at least as
     * large as {@code strict = false}) as {@code kolmogorovSmirnovStatistic(x, y)}.
     * <ul>
     * <li>For small samples (where the product of the sample sizes is less than
     * {@link KolmogorovSmirnovTest#LARGE_SAMPLE_PRODUCT}), the exact p-value is computed using the method presented
     * in [4], implemented in {@link #exactP(double, int, int, boolean)}. </li>
     * <li>When the product of the sample sizes exceeds {@link KolmogorovSmirnovTest#LARGE_SAMPLE_PRODUCT}, the
     * asymptotic distribution of \(D_{n,m}\) is used. See {@link #approximateP(double, int, int)}
     * for details on the approximation.</li>
     * </ul><p>
     * If {@code x.length * y.length} &lt; {@link KolmogorovSmirnovTest#LARGE_SAMPLE_PRODUCT} and the combined set of values in
     * {@code x} and {@code y} contains ties, random jitter is added to {@code x} and {@code y} to
     * break ties before computing \(D_{n,m}\) and the p-value. The jitter is uniformly distributed
     * on (-minDelta / 2, minDelta / 2) where minDelta is the smallest pairwise difference between
     * values in the combined sample.</p>
     * <p>
     * If ties are known to be present in the data, {@link KolmogorovSmirnovTest#bootstrap(double[], double[], int, boolean)}
     * may be used as an alternative method for estimating the p-value.</p>
     *
     * @param x first sample dataset
     * @param y second sample dataset
     * @param strict whether or not the probability to compute is expressed as a strict inequality
     *        (ignored for large samples)
     * @return p-value associated with the null hypothesis that {@code x} and {@code y} represent
     *         samples from the same distribution
     * @throws MathIllegalArgumentException if either {@code x} or {@code y} does not have length at
     *         least 2
     * @throws org.hipparchus.exception.NullArgumentException if either {@code x} or {@code y} is null
     * @see KolmogorovSmirnovTest#bootstrap(double[], double[], int, boolean)
     */
    public static double kolmogorovSmirnovTest(double[] x, double[] y, boolean strict)
            throws MathIllegalArgumentException, NullArgumentException  {
        return KS_TEST.kolmogorovSmirnovTest(x, y, strict);
    }

    /**
     * Computes \(P(D_{n,m} &gt; d)\) if {@code strict} is {@code true}; otherwise \(P(D_{n,m} \ge
     * d)\), where \(D_{n,m}\) is the 2-sample Kolmogorov-Smirnov statistic. See
     * {@link KolmogorovSmirnovTest#kolmogorovSmirnovStatistic(double[], double[])} for the definition of \(D_{n,m}\).
     * <p>
     * The returned probability is exact, implemented by unwinding the recursive function
     * definitions presented in [4] from the class javadoc.
     * </p>
     *
     * @param d D-statistic value
     * @param n first sample size
     * @param m second sample size
     * @param strict whether or not the probability to compute is expressed as a strict inequality
     * @return probability that a randomly selected m-n partition of m + n generates \(D_{n,m}\)
     *         greater than (resp. greater than or equal to) {@code d}
     */
    public static double exactP(double d, int m, int n, boolean strict) {
        return KS_TEST.exactP(d, n, m, strict);
    }

    /**
     * Uses the Kolmogorov-Smirnov distribution to approximate \(P(D_{n,m} &gt; d)\) where \(D_{n,m}\)
     * is the 2-sample Kolmogorov-Smirnov statistic. See
     * {@link KolmogorovSmirnovTest#kolmogorovSmirnovStatistic(double[], double[])} for the definition of \(D_{n,m}\).
     * <p>
     * Specifically, what is returned is \(1 - k(d \sqrt{mn / (m + n)})\) where \(k(t) = 1 + 2
     * \sum_{i=1}^\infty (-1)^i e^{-2 i^2 t^2}\). See {@link KolmogorovSmirnovTest#ksSum(double, double, int)} for
     * details on how convergence of the sum is determined. This implementation passes {@code ksSum}
     * {@link KolmogorovSmirnovTest#KS_SUM_CAUCHY_CRITERION} as {@code tolerance} and
     * {@link KolmogorovSmirnovTest#MAXIMUM_PARTIAL_SUM_COUNT} as {@code maxIterations}.
     * </p>
     *
     * @param d D-statistic value
     * @param n first sample size
     * @param m second sample size
     * @return approximate probability that a randomly selected m-n partition of m + n generates
     *         \(D_{n,m}\) greater than {@code d}
     */
    public static double approximateP(double d, int n, int m) {
        return KS_TEST.approximateP(d, n, m);
    }

}
