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

package org.hipparchus.random;

import org.hipparchus.distribution.RealDistribution;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.fraction.BigFraction;
import org.hipparchus.fraction.BigFractionField;
import org.hipparchus.linear.Array2DRowFieldMatrix;
import org.hipparchus.linear.FieldMatrix;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

import java.math.RoundingMode;
import java.util.Arrays;

/**
 * This is a stripped-down copy of the class from hipparchus-stat, fir test purposes only.
 */
class KolmogorovSmirnovTest { // NOPMD - this is not a Junit test class, PMD false positive here

    /**
     * Bound on the number of partial sums
     */
    protected static final int MAXIMUM_PARTIAL_SUM_COUNT = 100000;

    /** Convergence criterion for the sums in #pelzGood(double, double, int)} */
    protected static final double PG_SUM_RELATIVE_ERROR = 1.0e-10;

    /**
     * Computes the one-sample Kolmogorov-Smirnov test statistic, \(D_n=\sup_x |F_n(x)-F(x)|\) where
     * \(F\) is the distribution (cdf) function associated with {@code distribution}, \(n\) is the
     * length of {@code data} and \(F_n\) is the empirical distribution that puts mass \(1/n\) at
     * each of the values in {@code data}.
     *
     * @param distribution reference distribution
     * @param data sample being evaluated
     * @return Kolmogorov-Smirnov statistic \(D_n\)
     * @throws MathIllegalArgumentException if {@code data} does not have length at least 2
     * @throws org.hipparchus.exception.NullArgumentException if {@code data} is null
     */
    public double kolmogorovSmirnovStatistic(RealDistribution distribution, double[] data) {
        checkArray(data);
        final int n = data.length;
        final double nd = n;
        final double[] dataCopy = new double[n];
        System.arraycopy(data, 0, dataCopy, 0, n);
        Arrays.sort(dataCopy);
        double d = 0d;
        for (int i = 1; i <= n; i++) {
            final double yi = distribution.cumulativeProbability(dataCopy[i - 1]);
            final double currD = FastMath.max(yi - (i - 1) / nd, i / nd - yi);
            if (currD > d) {
                d = currD;
            }
        }
        return d;
    }

    /**
     * Calculates {@code P(D_n < d)} using the method described in [1] with quick decisions for extreme
     * values given in [2] (see above).
     *
     * @param d statistic
     * @param n sample size
     * @return \(P(D_n &lt; d)\)
     * @throws MathRuntimeException if algorithm fails to convert {@code h} to a
     *         {@link BigFraction} in expressing {@code d} as \((k
     *         - h) / m\) for integer {@code k, m} and \(0 &lt;= h &lt; 1\)
     */
    public double cdf(double d, int n)
        throws MathRuntimeException {
        return cdf(d, n, false);
    }

    /**
     * Calculates {@code P(D_n < d)} using method described in [1] with quick decisions for extreme
     * values given in [2] (see above).
     *
     * @param d statistic
     * @param n sample size
     * @param exact whether the probability should be calculated exact using
     *        {@link BigFraction} everywhere at the expense of
     *        very slow execution time, or if {@code double} should be used convenient places to
     *        gain speed. Almost never choose {@code true} in real applications unless you are very
     *        sure; {@code true} is almost solely for verification purposes.
     * @return \(P(D_n &lt; d)\)
     * @throws MathRuntimeException if algorithm fails to convert {@code h} to a
     *         {@link BigFraction} in expressing {@code d} as \((k
     *         - h) / m\) for integer {@code k, m} and \(0 \lt;= h &lt; 1\).
     */
    public double cdf(double d, int n, boolean exact)
        throws MathRuntimeException {

        final double ninv = 1 / ((double) n);
        final double ninvhalf = 0.5 * ninv;

        if (d <= ninvhalf) {
            return 0;
        } else if (ninvhalf < d && d <= ninv) {
            double res = 1;
            final double f = 2 * d - ninv;
            // n! f^n = n*f * (n-1)*f * ... * 1*x
            for (int i = 1; i <= n; ++i) {
                res *= i * f;
            }
            return res;
        } else if (1 - ninv <= d && d < 1) {
            return 1 - 2 * Math.pow(1 - d, n);
        } else if (1 <= d) {
            return 1;
        }
        if (exact) {
            return exactK(d, n);
        }
        if (n <= 140) {
            return roundedK(d, n);
        }
        return pelzGood(d, n);
    }

    /**
     * Calculates the exact value of {@code P(D_n < d)} using the method described in [1] (reference
     * in class javadoc above) and {@link BigFraction} (see
     * above).
     *
     * @param d statistic
     * @param n sample size
     * @return the two-sided probability of \(P(D_n < d)\)
     * @throws MathRuntimeException if algorithm fails to convert {@code h} to a
     *         {@link BigFraction} in expressing {@code d} as \((k
     *         - h) / m\) for integer {@code k, m} and \(0 \le h < 1\).
     */
    private double exactK(double d, int n)
        throws MathRuntimeException {

        final int k = (int) Math.ceil(n * d);

        final FieldMatrix<BigFraction> H = this.createExactH(d, n);
        final FieldMatrix<BigFraction> Hpower = H.power(n);

        BigFraction pFrac = Hpower.getEntry(k - 1, k - 1);

        for (int i = 1; i <= n; ++i) {
            pFrac = pFrac.multiply(i).divide(n);
        }

        /*
         * BigFraction.doubleValue converts numerator to double and the denominator to double and
         * divides afterwards. That gives NaN quite easy. This does not (scale is the number of
         * digits):
         */
        return pFrac.bigDecimalValue(20, RoundingMode.HALF_UP).doubleValue();
    }

    /**
     * Calculates {@code P(D_n < d)} using method described in [1] and doubles (see above).
     *
     * @param d statistic
     * @param n sample size
     * @return \(P(D_n < d)\)
     */
    private double roundedK(double d, int n) {

        final int k = (int) Math.ceil(n * d);
        final RealMatrix H = this.createRoundedH(d, n);
        final RealMatrix Hpower = H.power(n);

        double pFrac = Hpower.getEntry(k - 1, k - 1);
        for (int i = 1; i <= n; ++i) {
            pFrac *= ((double) i) / n;
        }

        return pFrac;
    }

    /**
     * Computes the Pelz-Good approximation for \(P(D_n &lt; d)\) as described in [2] in the class javadoc.
     *
     * @param d value of d-statistic (x in [2])
     * @param n sample size
     * @return \(P(D_n &lt; d)\)
     */
    public double pelzGood(double d, int n) {
        // Change the variable since approximation is for the distribution evaluated at d / sqrt(n)
        final double sqrtN = FastMath.sqrt(n);
        final double z = d * sqrtN;
        final double z2 = d * d * n;
        final double z4 = z2 * z2;
        final double z6 = z4 * z2;
        final double z8 = z4 * z4;

        // Compute K_0(z)
        double sum = 0;
        double z2Term = MathUtils.PI_SQUARED / (8 * z2);
        int k = 1;
        for (; k < MAXIMUM_PARTIAL_SUM_COUNT; k++) {
            final double kTerm = 2 * k - 1;
            final double increment = FastMath.exp(-z2Term * kTerm * kTerm);
            sum += increment;
            if (increment <= PG_SUM_RELATIVE_ERROR * sum) {
                break;
            }
        }
        if (k == MAXIMUM_PARTIAL_SUM_COUNT) {
            throw new MathIllegalStateException(LocalizedCoreFormats.MAX_COUNT_EXCEEDED, MAXIMUM_PARTIAL_SUM_COUNT);
        }
        double ret = sum * FastMath.sqrt(2 * FastMath.PI) / z;

        // K_1(z)
        // Sum is -inf to inf, but k term is always (k + 1/2) ^ 2, so really have
        // twice the sum from k = 0 to inf (k = -1 is same as 0, -2 same as 1, ...)
        final double twoZ2 = 2 * z2;
        sum = 0;
        for (k = 0; k < MAXIMUM_PARTIAL_SUM_COUNT; k++) {
            final double kTerm = k + 0.5;
            final double kTerm2 = kTerm * kTerm;
            final double increment = (MathUtils.PI_SQUARED * kTerm2 - z2) * FastMath.exp(-MathUtils.PI_SQUARED * kTerm2 / twoZ2);
            sum += increment;
            if (FastMath.abs(increment) < PG_SUM_RELATIVE_ERROR * FastMath.abs(sum)) {
                break;
            }
        }
        if (k == MAXIMUM_PARTIAL_SUM_COUNT) {
            throw new MathIllegalStateException(LocalizedCoreFormats.MAX_COUNT_EXCEEDED, MAXIMUM_PARTIAL_SUM_COUNT);
        }
        final double sqrtHalfPi = FastMath.sqrt(FastMath.PI / 2);
        // Instead of doubling sum, divide by 3 instead of 6
        ret += sum * sqrtHalfPi / (3 * z4 * sqrtN);

        // K_2(z)
        // Same drill as K_1, but with two doubly infinite sums, all k terms are even powers.
        final double z4Term = 2 * z4;
        final double z6Term = 6 * z6;
        z2Term = 5 * z2;
        final double pi4 = MathUtils.PI_SQUARED * MathUtils.PI_SQUARED;
        sum = 0;
        for (k = 0; k < MAXIMUM_PARTIAL_SUM_COUNT; k++) {
            final double kTerm = k + 0.5;
            final double kTerm2 = kTerm * kTerm;
            final double increment =  (z6Term + z4Term + MathUtils.PI_SQUARED * (z4Term - z2Term) * kTerm2 +
                    pi4 * (1 - twoZ2) * kTerm2 * kTerm2) * FastMath.exp(-MathUtils.PI_SQUARED * kTerm2 / twoZ2);
            sum += increment;
            if (FastMath.abs(increment) < PG_SUM_RELATIVE_ERROR * FastMath.abs(sum)) {
                break;
            }
        }
        if (k == MAXIMUM_PARTIAL_SUM_COUNT) {
            throw new MathIllegalStateException(LocalizedCoreFormats.MAX_COUNT_EXCEEDED, MAXIMUM_PARTIAL_SUM_COUNT);
        }
        double sum2 = 0;
        for (k = 1; k < MAXIMUM_PARTIAL_SUM_COUNT; k++) {
            final double kTerm2 = k * k;
            final double increment = MathUtils.PI_SQUARED * kTerm2 * FastMath.exp(-MathUtils.PI_SQUARED * kTerm2 / twoZ2);
            sum2 += increment;
            if (FastMath.abs(increment) < PG_SUM_RELATIVE_ERROR * FastMath.abs(sum2)) {
                break;
            }
        }
        if (k == MAXIMUM_PARTIAL_SUM_COUNT) {
            throw new MathIllegalStateException(LocalizedCoreFormats.MAX_COUNT_EXCEEDED, MAXIMUM_PARTIAL_SUM_COUNT);
        }
        // Again, adjust coefficients instead of doubling sum, sum2
        ret += (sqrtHalfPi / n) * (sum / (36 * z2 * z2 * z2 * z) - sum2 / (18 * z2 * z));

        // K_3(z) One more time with feeling - two doubly infinite sums, all k powers even.
        // Multiply coefficient denominators by 2, so omit doubling sums.
        final double pi6 = pi4 * MathUtils.PI_SQUARED;
        sum = 0;
        for (k = 0; k < MAXIMUM_PARTIAL_SUM_COUNT; k++) {
            final double kTerm = k + 0.5;
            final double kTerm2 = kTerm * kTerm;
            final double kTerm4 = kTerm2 * kTerm2;
            final double kTerm6 = kTerm4 * kTerm2;
            final double increment = (pi6 * kTerm6 * (5 - 30 * z2) + pi4 * kTerm4 * (-60 * z2 + 212 * z4) +
                            MathUtils.PI_SQUARED * kTerm2 * (135 * z4 - 96 * z6) - 30 * z6 - 90 * z8) *
                    FastMath.exp(-MathUtils.PI_SQUARED * kTerm2 / twoZ2);
            sum += increment;
            if (FastMath.abs(increment) < PG_SUM_RELATIVE_ERROR * FastMath.abs(sum)) {
                break;
            }
        }
        if (k == MAXIMUM_PARTIAL_SUM_COUNT) {
            throw new MathIllegalStateException(LocalizedCoreFormats.MAX_COUNT_EXCEEDED, MAXIMUM_PARTIAL_SUM_COUNT);
        }
        sum2 = 0;
        for (k = 1; k < MAXIMUM_PARTIAL_SUM_COUNT; k++) {
            final double kTerm2 = k * k;
            final double kTerm4 = kTerm2 * kTerm2;
            final double increment = (-pi4 * kTerm4 + 3 * MathUtils.PI_SQUARED * kTerm2 * z2) *
                    FastMath.exp(-MathUtils.PI_SQUARED * kTerm2 / twoZ2);
            sum2 += increment;
            if (FastMath.abs(increment) < PG_SUM_RELATIVE_ERROR * FastMath.abs(sum2)) {
                break;
            }
        }
        if (k == MAXIMUM_PARTIAL_SUM_COUNT) {
            throw new MathIllegalStateException(LocalizedCoreFormats.MAX_COUNT_EXCEEDED, MAXIMUM_PARTIAL_SUM_COUNT);
        }
        return ret + (sqrtHalfPi / (sqrtN * n)) * (sum / (3240 * z6 * z4) + sum2 / (108 * z6));

    }

    /***
     * Creates {@code H} of size {@code m x m} as described in [1] (see above).
     *
     * @param d statistic
     * @param n sample size
     * @return H matrix
     * @throws MathIllegalArgumentException if fractional part is greater than 1
     * @throws MathIllegalStateException if algorithm fails to convert {@code h} to a
     *         {@link BigFraction} in expressing {@code d} as \((k
     *         - h) / m\) for integer {@code k, m} and \(0 <= h < 1\).
     */
    private FieldMatrix<BigFraction> createExactH(double d, int n)
        throws MathIllegalArgumentException, MathIllegalStateException {

        final int k = (int) Math.ceil(n * d);
        final int m = 2 * k - 1;
        final double hDouble = k - n * d;
        if (hDouble >= 1) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_LARGE_BOUND_EXCLUDED,
                                                   hDouble, 1.0);
        }
        BigFraction h;
        try {
            h = new BigFraction(hDouble, 1.0e-20, 10000);
        } catch (final MathIllegalStateException e1) {
            try {
                h = new BigFraction(hDouble, 1.0e-10, 10000);
            } catch (final MathIllegalStateException e2) {
                h = new BigFraction(hDouble, 1.0e-5, 10000);
            }
        }
        final BigFraction[][] Hdata = new BigFraction[m][m];

        /*
         * Start by filling everything with either 0 or 1.
         */
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                if (i - j + 1 < 0) {
                    Hdata[i][j] = BigFraction.ZERO;
                } else {
                    Hdata[i][j] = BigFraction.ONE;
                }
            }
        }

        /*
         * Setting up power-array to avoid calculating the same value twice: hPowers[0] = h^1 ...
         * hPowers[m-1] = h^m
         */
        final BigFraction[] hPowers = new BigFraction[m];
        hPowers[0] = h;
        for (int i = 1; i < m; ++i) {
            hPowers[i] = h.multiply(hPowers[i - 1]);
        }

        /*
         * First column and last row has special values (each other reversed).
         */
        for (int i = 0; i < m; ++i) {
            Hdata[i][0] = Hdata[i][0].subtract(hPowers[i]);
            Hdata[m - 1][i] = Hdata[m - 1][i].subtract(hPowers[m - i - 1]);
        }

        /*
         * [1] states: "For 1/2 < h < 1 the bottom left element of the matrix should be (1 - 2*h^m +
         * (2h - 1)^m )/m!" Since 0 <= h < 1, then if h > 1/2 is sufficient to check:
         */
        if (h.compareTo(BigFraction.ONE_HALF) == 1) {
            Hdata[m - 1][0] = Hdata[m - 1][0].add(h.multiply(2).subtract(1).pow(m));
        }

        /*
         * Aside from the first column and last row, the (i, j)-th element is 1/(i - j + 1)! if i -
         * j + 1 >= 0, else 0. 1's and 0's are already put, so only division with (i - j + 1)! is
         * needed in the elements that have 1's. There is no need to calculate (i - j + 1)! and then
         * divide - small steps avoid overflows. Note that i - j + 1 > 0 <=> i + 1 > j instead of
         * j'ing all the way to m. Also note that it is started at g = 2 because dividing by 1 isn't
         * really necessary.
         */
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < i + 1; ++j) {
                if (i - j + 1 > 0) {
                    for (int g = 2; g <= i - j + 1; ++g) {
                        Hdata[i][j] = Hdata[i][j].divide(g);
                    }
                }
            }
        }
        return new Array2DRowFieldMatrix<>(BigFractionField.getInstance(), Hdata);
    }

    /***
     * Creates {@code H} of size {@code m x m} as described in [1] (see above)
     * using double-precision.
     *
     * @param d statistic
     * @param n sample size
     * @return H matrix
     * @throws MathIllegalArgumentException if fractional part is greater than 1
     */
    private RealMatrix createRoundedH(double d, int n)
        throws MathIllegalArgumentException {

        final int k = (int) Math.ceil(n * d);
        final int m = 2 * k - 1;
        final double h = k - n * d;
        if (h >= 1) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_LARGE_BOUND_EXCLUDED,
                                                   h, 1.0);
        }
        final double[][] Hdata = new double[m][m];

        /*
         * Start by filling everything with either 0 or 1.
         */
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                if (i - j + 1 < 0) {
                    Hdata[i][j] = 0;
                } else {
                    Hdata[i][j] = 1;
                }
            }
        }

        /*
         * Setting up power-array to avoid calculating the same value twice: hPowers[0] = h^1 ...
         * hPowers[m-1] = h^m
         */
        final double[] hPowers = new double[m];
        hPowers[0] = h;
        for (int i = 1; i < m; ++i) {
            hPowers[i] = h * hPowers[i - 1];
        }

        /*
         * First column and last row has special values (each other reversed).
         */
        for (int i = 0; i < m; ++i) {
            Hdata[i][0] = Hdata[i][0] - hPowers[i];
            Hdata[m - 1][i] -= hPowers[m - i - 1];
        }

        /*
         * [1] states: "For 1/2 < h < 1 the bottom left element of the matrix should be (1 - 2*h^m +
         * (2h - 1)^m )/m!" Since 0 <= h < 1, then if h > 1/2 is sufficient to check:
         */
        if (Double.compare(h, 0.5) > 0) {
            Hdata[m - 1][0] += FastMath.pow(2 * h - 1, m);
        }

        /*
         * Aside from the first column and last row, the (i, j)-th element is 1/(i - j + 1)! if i -
         * j + 1 >= 0, else 0. 1's and 0's are already put, so only division with (i - j + 1)! is
         * needed in the elements that have 1's. There is no need to calculate (i - j + 1)! and then
         * divide - small steps avoid overflows. Note that i - j + 1 > 0 <=> i + 1 > j instead of
         * j'ing all the way to m. Also note that it is started at g = 2 because dividing by 1 isn't
         * really necessary.
         */
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < i + 1; ++j) {
                if (i - j + 1 > 0) {
                    for (int g = 2; g <= i - j + 1; ++g) {
                        Hdata[i][j] /= g;
                    }
                }
            }
        }
        return MatrixUtils.createRealMatrix(Hdata);
    }

    /**
     * Verifies that {@code array} has length at least 2.
     *
     * @param array array to test
     * @throws org.hipparchus.exception.NullArgumentException if array is null
     * @throws MathIllegalArgumentException if array is too short
     */
    private void checkArray(double[] array) {
        MathUtils.checkNotNull(array);
        if (array.length < 2) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.INSUFFICIENT_OBSERVED_POINTS_IN_SAMPLE,
                                                   array.length, 2);
        }
    }

}
