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
package org.hipparchus.util;

import java.lang.reflect.Array;
import java.util.Iterator;
import java.util.List;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.special.Gamma;

/**
 * Combinatorial utilities.
 */
public final class CombinatoricsUtils {

    /** Maximum index of Bell number that fits into a long.
     * @since 2.2
     */
    public static final int MAX_BELL = 25;

    /** All long-representable factorials */
    static final long[] FACTORIALS = {
                       1L,                  1L,                   2L,
                       6L,                 24L,                 120L,
                     720L,               5040L,               40320L,
                  362880L,            3628800L,            39916800L,
               479001600L,         6227020800L,         87178291200L,
           1307674368000L,     20922789888000L,     355687428096000L,
        6402373705728000L, 121645100408832000L, 2432902008176640000L };

    /** Stirling numbers of the second kind. */
    static final AtomicReference<long[][]> STIRLING_S2 = new AtomicReference<> (null);

    /**
     * Default implementation of {@link #factorialLog(int)} method:
     * <ul>
     *  <li>No pre-computation</li>
     *  <li>No cache allocation</li>
     * </ul>
     */
    private static final FactorialLog FACTORIAL_LOG_NO_CACHE = FactorialLog.create();

    /** Bell numbers.
     * @since 2.2
     */
    private static final AtomicReference<long[]> BELL = new AtomicReference<> (null);

    /** Private constructor (class contains only static methods). */
    private CombinatoricsUtils() {}


    /**
     * Returns an exact representation of the <a
     * href="http://mathworld.wolfram.com/BinomialCoefficient.html"> Binomial
     * Coefficient</a>, "{@code n choose k}", the number of
     * {@code k}-element subsets that can be selected from an
     * {@code n}-element set.
     * <p><Strong>Preconditions</strong>:</p>
     * <ul>
     * <li> {@code 0 <= k <= n } (otherwise
     * {@code MathIllegalArgumentException} is thrown)</li>
     * <li> The result is small enough to fit into a {@code long}. The
     * largest value of {@code n} for which all coefficients are
     * {@code  < Long.MAX_VALUE} is 66. If the computed value exceeds
     * {@code Long.MAX_VALUE} a {@code MathRuntimeException} is
     * thrown.</li>
     * </ul>
     *
     * @param n the size of the set
     * @param k the size of the subsets to be counted
     * @return {@code n choose k}
     * @throws MathIllegalArgumentException if {@code n < 0}.
     * @throws MathIllegalArgumentException if {@code k > n}.
     * @throws MathRuntimeException if the result is too large to be
     * represented by a long integer.
     */
    public static long binomialCoefficient(final int n, final int k)
        throws MathIllegalArgumentException, MathRuntimeException {
        checkBinomial(n, k);
        if ((n == k) || (k == 0)) {
            return 1;
        }
        if ((k == 1) || (k == n - 1)) {
            return n;
        }
        // Use symmetry for large k
        if (k > n / 2) {
            return binomialCoefficient(n, n - k);
        }

        // We use the formula
        // (n choose k) = n! / (n-k)! / k!
        // (n choose k) == ((n-k+1)*...*n) / (1*...*k)
        // which could be written
        // (n choose k) == (n-1 choose k-1) * n / k
        long result = 1;
        if (n <= 61) {
            // For n <= 61, the naive implementation cannot overflow.
            int i = n - k + 1;
            for (int j = 1; j <= k; j++) {
                result = result * i / j;
                i++;
            }
        } else if (n <= 66) {
            // For n > 61 but n <= 66, the result cannot overflow,
            // but we must take care not to overflow intermediate values.
            int i = n - k + 1;
            for (int j = 1; j <= k; j++) {
                // We know that (result * i) is divisible by j,
                // but (result * i) may overflow, so we split j:
                // Filter out the gcd, d, so j/d and i/d are integer.
                // result is divisible by (j/d) because (j/d)
                // is relative prime to (i/d) and is a divisor of
                // result * (i/d).
                final long d = ArithmeticUtils.gcd(i, j);
                result = (result / (j / d)) * (i / d);
                i++;
            }
        } else {
            // For n > 66, a result overflow might occur, so we check
            // the multiplication, taking care to not overflow
            // unnecessary.
            int i = n - k + 1;
            for (int j = 1; j <= k; j++) {
                final long d = ArithmeticUtils.gcd(i, j);
                result = ArithmeticUtils.mulAndCheck(result / (j / d), i / d);
                i++;
            }
        }
        return result;
    }

    /**
     * Returns a {@code double} representation of the <a
     * href="http://mathworld.wolfram.com/BinomialCoefficient.html"> Binomial
     * Coefficient</a>, "{@code n choose k}", the number of
     * {@code k}-element subsets that can be selected from an
     * {@code n}-element set.
     * <p>* <Strong>Preconditions</strong>:</p>
     * <ul>
     * <li> {@code 0 <= k <= n } (otherwise
     * {@code IllegalArgumentException} is thrown)</li>
     * <li> The result is small enough to fit into a {@code double}. The
     * largest value of {@code n} for which all coefficients are &lt;
     * Double.MAX_VALUE is 1029. If the computed value exceeds Double.MAX_VALUE,
     * Double.POSITIVE_INFINITY is returned</li>
     * </ul>
     *
     * @param n the size of the set
     * @param k the size of the subsets to be counted
     * @return {@code n choose k}
     * @throws MathIllegalArgumentException if {@code n < 0}.
     * @throws MathIllegalArgumentException if {@code k > n}.
     * @throws MathRuntimeException if the result is too large to be
     * represented by a long integer.
     */
    public static double binomialCoefficientDouble(final int n, final int k)
        throws MathIllegalArgumentException, MathRuntimeException {
        checkBinomial(n, k);
        if ((n == k) || (k == 0)) {
            return 1d;
        }
        if ((k == 1) || (k == n - 1)) {
            return n;
        }
        if (k > n/2) {
            return binomialCoefficientDouble(n, n - k);
        }
        if (n < 67) {
            return binomialCoefficient(n,k);
        }

        double result = 1d;
        for (int i = 1; i <= k; i++) {
             result *= ((double) (n - k + i)) / i;
        }

        return FastMath.floor(result + 0.5);
    }

    /**
     * Returns the natural {@code log} of the <a
     * href="http://mathworld.wolfram.com/BinomialCoefficient.html"> Binomial
     * Coefficient</a>, "{@code n choose k}", the number of
     * {@code k}-element subsets that can be selected from an
     * {@code n}-element set.
     * <p>* <Strong>Preconditions</strong>:</p>
     * <ul>
     * <li> {@code 0 <= k <= n } (otherwise
     * {@code MathIllegalArgumentException} is thrown)</li>
     * </ul>
     *
     * @param n the size of the set
     * @param k the size of the subsets to be counted
     * @return {@code n choose k}
     * @throws MathIllegalArgumentException if {@code n < 0}.
     * @throws MathIllegalArgumentException if {@code k > n}.
     * @throws MathRuntimeException if the result is too large to be
     * represented by a long integer.
     */
    public static double binomialCoefficientLog(final int n, final int k)
        throws MathIllegalArgumentException, MathRuntimeException {
        checkBinomial(n, k);
        if ((n == k) || (k == 0)) {
            return 0;
        }
        if ((k == 1) || (k == n - 1)) {
            return FastMath.log(n);
        }

        /*
         * For values small enough to do exact integer computation,
         * return the log of the exact value
         */
        if (n < 67) {
            return FastMath.log(binomialCoefficient(n,k));
        }

        /*
         * Return the log of binomialCoefficientDouble for values that will not
         * overflow binomialCoefficientDouble
         */
        if (n < 1030) {
            return FastMath.log(binomialCoefficientDouble(n, k));
        }

        if (k > n / 2) {
            return binomialCoefficientLog(n, n - k);
        }

        /*
         * Sum logs for values that could overflow
         */
        double logSum = 0;

        // n!/(n-k)!
        for (int i = n - k + 1; i <= n; i++) {
            logSum += FastMath.log(i);
        }

        // divide by k!
        for (int i = 2; i <= k; i++) {
            logSum -= FastMath.log(i);
        }

        return logSum;
    }

    /**
     * Returns n!. Shorthand for {@code n} <a
     * href="http://mathworld.wolfram.com/Factorial.html"> Factorial</a>, the
     * product of the numbers {@code 1,...,n}.
     * <p>* <Strong>Preconditions</strong>:</p>
     * <ul>
     * <li> {@code n >= 0} (otherwise
     * {@code MathIllegalArgumentException} is thrown)</li>
     * <li> The result is small enough to fit into a {@code long}. The
     * largest value of {@code n} for which {@code n!} does not exceed
     * Long.MAX_VALUE} is 20. If the computed value exceeds {@code Long.MAX_VALUE}
     * an {@code MathRuntimeException } is thrown.</li>
     * </ul>
     *
     * @param n argument
     * @return {@code n!}
     * @throws MathRuntimeException if the result is too large to be represented
     * by a {@code long}.
     * @throws MathIllegalArgumentException if {@code n < 0}.
     * @throws MathIllegalArgumentException if {@code n > 20}: The factorial value is too
     * large to fit in a {@code long}.
     */
    public static long factorial(final int n) throws MathIllegalArgumentException {
        if (n < 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.FACTORIAL_NEGATIVE_PARAMETER, n);
        }
        if (n > 20) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_LARGE, n, 20);
        }
        return FACTORIALS[n];
    }

    /**
     * Compute n!, the<a href="http://mathworld.wolfram.com/Factorial.html">
     * factorial</a> of {@code n} (the product of the numbers 1 to n), as a
     * {@code double}.
     * The result should be small enough to fit into a {@code double}: The
     * largest {@code n} for which {@code n!} does not exceed
     * {@code Double.MAX_VALUE} is 170. If the computed value exceeds
     * {@code Double.MAX_VALUE}, {@code Double.POSITIVE_INFINITY} is returned.
     *
     * @param n Argument.
     * @return {@code n!}
     * @throws MathIllegalArgumentException if {@code n < 0}.
     */
    public static double factorialDouble(final int n) throws MathIllegalArgumentException {
        if (n < 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.FACTORIAL_NEGATIVE_PARAMETER,
                                           n);
        }
        if (n < 21) {
            return FACTORIALS[n];
        }
        return FastMath.floor(FastMath.exp(factorialLog(n)) + 0.5);
    }

    /**
     * Compute the natural logarithm of the factorial of {@code n}.
     *
     * @param n Argument.
     * @return {@code log(n!)}
     * @throws MathIllegalArgumentException if {@code n < 0}.
     */
    public static double factorialLog(final int n) throws MathIllegalArgumentException {
        return FACTORIAL_LOG_NO_CACHE.value(n);
    }

    /**
     * Returns the <a
     * href="http://mathworld.wolfram.com/StirlingNumberoftheSecondKind.html">
     * Stirling number of the second kind</a>, "{@code S(n,k)}", the number of
     * ways of partitioning an {@code n}-element set into {@code k} non-empty
     * subsets.
     * <p>
     * The preconditions are {@code 0 <= k <= n } (otherwise
     * {@code MathIllegalArgumentException} is thrown)
     * </p>
     * @param n the size of the set
     * @param k the number of non-empty subsets
     * @return {@code S(n,k)}
     * @throws MathIllegalArgumentException if {@code k < 0}.
     * @throws MathIllegalArgumentException if {@code k > n}.
     * @throws MathRuntimeException if some overflow happens, typically for n exceeding 25 and
     * k between 20 and n-2 (S(n,n-1) is handled specifically and does not overflow)
     */
    public static long stirlingS2(final int n, final int k)
        throws MathIllegalArgumentException, MathRuntimeException {
        if (k < 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL, k, 0);
        }
        if (k > n) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_LARGE,
                                                   k, n);
        }

        long[][] stirlingS2 = STIRLING_S2.get();

        if (stirlingS2 == null) {
            // the cache has never been initialized, compute the first numbers
            // by direct recurrence relation

            // as S(26,9) = 11201516780955125625 is larger than Long.MAX_VALUE
            // we must stop computation at row 26
            final int maxIndex = 26;
            stirlingS2 = new long[maxIndex][];
            stirlingS2[0] = new long[] { 1l };
            for (int i = 1; i < stirlingS2.length; ++i) {
                stirlingS2[i] = new long[i + 1];
                stirlingS2[i][0] = 0;
                stirlingS2[i][1] = 1;
                stirlingS2[i][i] = 1;
                for (int j = 2; j < i; ++j) {
                    stirlingS2[i][j] = j * stirlingS2[i - 1][j] + stirlingS2[i - 1][j - 1];
                }
            }

            // atomically save the cache
            STIRLING_S2.compareAndSet(null, stirlingS2);

        }

        if (n < stirlingS2.length) {
            // the number is in the small cache
            return stirlingS2[n][k];
        } else {
            // use explicit formula to compute the number without caching it
            if (k == 0) {
                return 0;
            } else if (k == 1 || k == n) {
                return 1;
            } else if (k == 2) {
                return (1l << (n - 1)) - 1l;
            } else if (k == n - 1) {
                return binomialCoefficient(n, 2);
            } else {
                // definition formula: note that this may trigger some overflow
                long sum = 0;
                long sign = ((k & 0x1) == 0) ? 1 : -1;
                for (int j = 1; j <= k; ++j) {
                    sign = -sign;
                    sum += sign * binomialCoefficient(k, j) * ArithmeticUtils.pow(j, n);
                    if (sum < 0) {
                        // there was an overflow somewhere
                        throw new MathRuntimeException(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE,
                                                          n, 0, stirlingS2.length - 1);
                    }
                }
                return sum / factorial(k);
            }
        }

    }

    /**
     * Returns an iterator whose range is the k-element subsets of {0, ..., n - 1}
     * represented as {@code int[]} arrays.
     * <p>
     * The arrays returned by the iterator are sorted in descending order and
     * they are visited in lexicographic order with significance from right to
     * left. For example, combinationsIterator(4, 2) returns an Iterator that
     * will generate the following sequence of arrays on successive calls to
     * {@code next()}:</p><p>
     * {@code [0, 1], [0, 2], [1, 2], [0, 3], [1, 3], [2, 3]}
     * </p><p>
     * If {@code k == 0} an Iterator containing an empty array is returned and
     * if {@code k == n} an Iterator containing [0, ..., n -1] is returned.</p>
     *
     * @param n Size of the set from which subsets are selected.
     * @param k Size of the subsets to be enumerated.
     * @return an {@link Iterator iterator} over the k-sets in n.
     * @throws MathIllegalArgumentException if {@code n < 0}.
     * @throws MathIllegalArgumentException if {@code k > n}.
     */
    public static Iterator<int[]> combinationsIterator(int n, int k) {
        return new Combinations(n, k).iterator();
    }

    /**
     * Check binomial preconditions.
     *
     * @param n Size of the set.
     * @param k Size of the subsets to be counted.
     * @throws MathIllegalArgumentException if {@code n < 0}.
     * @throws MathIllegalArgumentException if {@code k > n}.
     */
    public static void checkBinomial(final int n,
                                     final int k)
        throws MathIllegalArgumentException {
        if (n < k) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.BINOMIAL_INVALID_PARAMETERS_ORDER,
                                                k, n, true);
        }
        if (n < 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.BINOMIAL_NEGATIVE_PARAMETER, n);
        }
    }

    /** Compute the Bell number (number of partitions of a set).
     * @param n number of elements of the set
     * @return Bell number Bₙ
     * @since 2.2
     */
    public static long bellNumber(final int n) {

        if (n < 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL, n, 0);
        }
        if (n > MAX_BELL) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_LARGE, n, MAX_BELL);
        }

        long[] bell = BELL.get();

        if (bell == null) {

            // the cache has never been initialized, compute the numbers using the Bell triangle
            // storage for one line of the Bell triangle
            bell = new long[MAX_BELL];
            bell[0] = 1l;

            final long[] row = new long[bell.length];
            for (int k = 1; k < row.length; ++k) {

                // first row, with one element
                row[0] = 1l;

                // iterative computation of rows
                for (int i = 1; i < k; ++i) {
                    long previous = row[0];
                    row[0] = row[i - 1];
                    for (int j = 1; j <= i; ++j) {
                        long rj = row[j - 1] + previous;
                        previous = row[j];
                        row[j] = rj;
                    }
                }

                bell[k] = row[k - 1];

            }

            BELL.compareAndSet(null, bell);

        }

        return bell[n];

    }

    /** Generate a stream of partitions of a list.
     * <p>
     * This method implements the iterative algorithm described in
     * <a href="https://academic.oup.com/comjnl/article/32/3/281/331557">Short Note:
     * A Fast Iterative Algorithm for Generating Set Partitions</a>
     * by B. Djokić, M. Miyakawa, S. Sekiguchi, I. Semba, and I. Stojmenović
     * (The Computer Journal, Volume 32, Issue 3, 1989, Pages 281–282,
     * <a href="https://doi.org/10.1093/comjnl/32.3.281">https://doi.org/10.1093/comjnl/32.3.281</a>
     * </p>
     * @param <T> type of the list elements
     * @param list list to partition
     * @return stream of partitions of the list, each partition is an array or parts
     * and each part is a list of elements
     * @since 2.2
     */
    public static <T> Stream<List<T>[]> partitions(final List<T> list) {

        // handle special cases of empty and singleton lists
        if (list.size() < 2) {
            @SuppressWarnings("unchecked")
            final List<T>[] partition = (List<T>[]) Array.newInstance(List.class, 1);
            partition[0] = list;
            final Stream.Builder<List<T>[]> builder = Stream.builder();
            return builder.add(partition).build();
        }

        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(new PartitionsIterator<>(list),
                                                                        Spliterator.DISTINCT | Spliterator.NONNULL |
                                                                        Spliterator.IMMUTABLE | Spliterator.ORDERED),
                                    false);

    }

    /** Generate a stream of permutations of a list.
     * <p>
     * This method implements the Steinhaus–Johnson–Trotter algorithm
     * with Even's speedup
     * <a href="https://en.wikipedia.org/wiki/Steinhaus%E2%80%93Johnson%E2%80%93Trotter_algorithm">Steinhaus–Johnson–Trotter algorithm</a>
     * @param <T> type of the list elements
     * @param list list to permute
     * @return stream of permutations of the list
     * @since 2.2
     */
    public static <T> Stream<List<T>> permutations(final List<T> list) {

        // handle special cases of empty and singleton lists
        if (list.size() < 2) {
            return Stream.of(list);
        }

        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(new PermutationsIterator<>(list),
                                                                        Spliterator.DISTINCT | Spliterator.NONNULL |
                                                                        Spliterator.IMMUTABLE | Spliterator.ORDERED),
                                    false);

    }

    /**
     * Class for computing the natural logarithm of the factorial of {@code n}.
     * It allows to allocate a cache of precomputed values.
     * In case of cache miss, computation is preformed by a call to
     * {@link Gamma#logGamma(double)}.
     */
    public static final class FactorialLog {
        /**
         * Precomputed values of the function:
         * {@code LOG_FACTORIALS[i] = log(i!)}.
         */
        private final double[] LOG_FACTORIALS;

        /**
         * Creates an instance, reusing the already computed values if available.
         *
         * @param numValues Number of values of the function to compute.
         * @param cache Existing cache.
         * @throw MathIllegalArgumentException if {@code n < 0}.
         */
        private FactorialLog(int numValues,
                             double[] cache) {
            if (numValues < 0) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL, numValues, 0);
            }

            LOG_FACTORIALS = new double[numValues];

            final int beginCopy = 2;
            final int endCopy = cache == null || cache.length <= beginCopy ?
                beginCopy : cache.length <= numValues ?
                cache.length : numValues;

            // Copy available values.
            for (int i = beginCopy; i < endCopy; i++) {
                LOG_FACTORIALS[i] = cache[i];
            }

            // Precompute.
            for (int i = endCopy; i < numValues; i++) {
                LOG_FACTORIALS[i] = LOG_FACTORIALS[i - 1] + FastMath.log(i);
            }
        }

        /**
         * Creates an instance with no precomputed values.
         * @return instance with no precomputed values
         */
        public static FactorialLog create() {
            return new FactorialLog(0, null);
        }

        /**
         * Creates an instance with the specified cache size.
         *
         * @param cacheSize Number of precomputed values of the function.
         * @return a new instance where {@code cacheSize} values have been
         * precomputed.
         * @throws MathIllegalArgumentException if {@code n < 0}.
         */
        public FactorialLog withCache(final int cacheSize) {
            return new FactorialLog(cacheSize, LOG_FACTORIALS);
        }

        /**
         * Computes {@code log(n!)}.
         *
         * @param n Argument.
         * @return {@code log(n!)}.
         * @throws MathIllegalArgumentException if {@code n < 0}.
         */
        public double value(final int n) {
            if (n < 0) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.FACTORIAL_NEGATIVE_PARAMETER,
                                               n);
            }

            // Use cache of precomputed values.
            if (n < LOG_FACTORIALS.length) {
                return LOG_FACTORIALS[n];
            }

            // Use cache of precomputed factorial values.
            if (n < FACTORIALS.length) {
                return FastMath.log(FACTORIALS[n]);
            }

            // Delegate.
            return Gamma.logGamma(n + 1);
        }
    }
}
