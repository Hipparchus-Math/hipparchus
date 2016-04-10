/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
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
package org.hipparchus.random;

import java.io.Serializable;
import java.util.Arrays;
import java.util.HashSet;

import org.hipparchus.distribution.IntegerDistribution;
import org.hipparchus.distribution.RealDistribution;
import org.hipparchus.util.CombinatoricsUtils;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.hipparchus.util.ResizableDoubleArray;

/**
 *
 * A {@link RandomGenerator} that provides additional methods to generate deviates
 * following probability distributions.
 *
 */
public class RandomData implements RandomGenerator, Serializable {

    /**
     * Cached random normal value.  The default implementation for
     * {@link #nextGaussian} generates pairs of values and this field caches the
     * second value so that the full algorithm is not executed for every
     * activation.  The value {@code Double.NaN} signals that there is
     * no cached value.  Use {@link #clear} to clear the cached value.
     */
    private double cachedNormalDeviate = Double.NaN;

    /**
     * Used when generating Exponential samples.
     * Table containing the constants
     * q_i = sum_{j=1}^i (ln 2)^j/j! = ln 2 + (ln 2)^2/2 + ... + (ln 2)^i/i!
     * until the largest representable fraction below 1 is exceeded.
     *
     * Note that
     * 1 = 2 - 1 = exp(ln 2) - 1 = sum_{n=1}^infty (ln 2)^n / n!
     * thus q_i -> 1 as i -> +inf,
     * so the higher i, the closer to one we get (the series is not alternating).
     *
     * By trying, n = 16 in Java is enough to reach 1.0.
     */
    private static final double[] EXPONENTIAL_SA_QI;

    /** Serializable version identifier. */
    private static final long serialVersionUID = 2306581345647615033L;

    /** Source of random data */
    private final RandomGenerator randomGenerator;

    private static final HashSet<String> CONTINUOUS_NAMES = new HashSet<String>();
    private static final HashSet<String> DISCRETE_NAMES = new HashSet<String>();

    /**
     * Initialize tables.
     */
    static {
        /**
         * Filling EXPONENTIAL_SA_QI table.
         * Note that we don't want qi = 0 in the table.
         */
        final double LN2 = FastMath.log(2);
        double qi = 0;
        int i = 1;

        /**
         * ArithmeticUtils provides factorials up to 20, so let's use that
         * limit together with Precision.EPSILON to generate the following
         * code (a priori, we know that there will be 16 elements, but it is
         * better to not hardcode it).
         */
        final ResizableDoubleArray ra = new ResizableDoubleArray(20);

        while (qi < 1) {
            qi += FastMath.pow(LN2, i) / CombinatoricsUtils.factorial(i);
            ra.addElement(qi);
            ++i;
        }

        EXPONENTIAL_SA_QI = ra.getElements();

        final String[] discreteCns = {"list", "them", "all"};
        DISCRETE_NAMES.addAll(Arrays.asList(discreteCns));
        final String[] continuousCns = {"list", "them", "all"};
        CONTINUOUS_NAMES.addAll(Arrays.asList(continuousCns));
    }

    /**
     * Prevent instantiation without a generator argument
     */
    @SuppressWarnings("unused")
    private RandomData() { randomGenerator = null; }

    /**
     * Construct a RandomAdaptor wrapping the supplied RandomGenerator.
     *
     * @param randomGenerator  the wrapped generator
     */
    public RandomData(RandomGenerator randomGenerator) {
        this.randomGenerator = randomGenerator;
    }

    /**
     * Factory method to create a {@codeRandomData} instance using the supplied
     * {@code RandomGenerator}.
     *
     * @param randomGenerator source of random bits
     * @return a RandomData using the given RandomGenerator to source bits
     */
    public static RandomData of(RandomGenerator randomGenerator) {
        return new RandomData(randomGenerator);
    }

    /**
     * Returns the next pseudo-random, uniformly distributed
     * <code>boolean</code> value from this random number generator's
     * sequence.
     *
     * @return  the next pseudo-random, uniformly distributed
     * <code>boolean</code> value from this random number generator's
     * sequence
     */
    @Override
    public boolean nextBoolean() {
        return randomGenerator.nextBoolean();
    }

     /**
     * Generates random bytes and places them into a user-supplied
     * byte array.  The number of random bytes produced is equal to
     * the length of the byte array.
     *
     * @param bytes the non-null byte array in which to put the
     * random bytes
     */
    @Override
    public void nextBytes(byte[] bytes) {
        randomGenerator.nextBytes(bytes);
    }

     /**
     * Returns the next pseudo-random, uniformly distributed
     * <code>double</code> value between <code>0.0</code> and
     * <code>1.0</code> from this random number generator's sequence.
     *
     * @return  the next pseudo-random, uniformly distributed
     *  <code>double</code> value between <code>0.0</code> and
     *  <code>1.0</code> from this random number generator's sequence
     */
    @Override
    public double nextDouble() {
        return randomGenerator.nextDouble();
    }

    /**
     * Returns the next pseudo-random, uniformly distributed <code>float</code>
     * value between <code>0.0</code> and <code>1.0</code> from this random
     * number generator's sequence.
     *
     * @return  the next pseudo-random, uniformly distributed <code>float</code>
     * value between <code>0.0</code> and <code>1.0</code> from this
     * random number generator's sequence
     */
    @Override
    public float nextFloat() {
        return randomGenerator.nextFloat();
    }

     /**
     * Returns the next pseudo-random, uniformly distributed <code>int</code>
     * value from this random number generator's sequence.
     * All 2<font size="-1"><sup>32</sup></font> possible {@code int} values
     * should be produced with  (approximately) equal probability.
     *
     * @return the next pseudo-random, uniformly distributed <code>int</code>
     *  value from this random number generator's sequence
     */
    @Override
    public int nextInt() {
        return randomGenerator.nextInt();
    }

    /**
     * Returns a pseudo-random, uniformly distributed {@code int} value
     * between 0 (inclusive) and the specified value (exclusive), drawn from
     * this random number generator's sequence.
     *
     * @param n the bound on the random number to be returned.  Must be
     * positive.
     * @return  a pseudo-random, uniformly distributed {@code int}
     * value between 0 (inclusive) and n (exclusive).
     * @throws IllegalArgumentException  if n is not positive.
     */
    @Override
    public int nextInt(int n) {
        return randomGenerator.nextInt(n);
    }

    /**
     * Returns the next pseudo-random, uniformly distributed <code>long</code>
     * value from this random number generator's sequence.  All
     * 2<font size="-1"><sup>64</sup></font> possible {@code long} values
     * should be produced with (approximately) equal probability.
     *
     * @return  the next pseudo-random, uniformly distributed <code>long</code>
     *value from this random number generator's sequence
     */
    @Override
    public long nextLong() {
        return randomGenerator.nextLong();
    }

    public void setSeed(int seed) {
        if (randomGenerator != null) {  // required to avoid NPE in constructor
            randomGenerator.setSeed(seed);
        }
    }

    public void setSeed(int[] seed) {
        if (randomGenerator != null) {  // required to avoid NPE in constructor
            randomGenerator.setSeed(seed);
        }
    }

    /** {@inheritDoc} */
    @Override
    public void setSeed(long seed) {
        if (randomGenerator != null) {  // required to avoid NPE in constructor
            randomGenerator.setSeed(seed);
        }
    }

    /**
     * Returns the next pseudo-random, Gaussian ("normally") distributed
     * {@code double} value with mean {@code 0.0} and standard
     * deviation {@code 1.0} from this random number generator's sequence.
     * <p>
     * The default implementation uses the <em>Polar Method</em>
     * due to G.E.P. Box, M.E. Muller and G. Marsaglia, as described in
     * D. Knuth, <u>The Art of Computer Programming</u>, 3.4.1C.</p>
     * <p>
     * The algorithm generates a pair of independent random values.  One of
     * these is cached for reuse, so the full algorithm is not executed on each
     * activation.  Implementations that do not override this method should
     * make sure to call {@link #clear} to clear the cached value in the
     * implementation of {@link #setSeed(long)}.</p>
     *
     * @return  the next pseudorandom, Gaussian ("normally") distributed
     * {@code double} value with mean {@code 0.0} and
     * standard deviation {@code 1.0} from this random number
     *  generator's sequence
     */
    @Override
    public double nextGaussian() {
        if (!Double.isNaN(cachedNormalDeviate)) {
            double dev = cachedNormalDeviate;
            cachedNormalDeviate = Double.NaN;
            return dev;
        }
        double v1 = 0;
        double v2 = 0;
        double s = 1;
        while (s >=1 ) {
            v1 = 2 * randomGenerator.nextDouble() - 1;
            v2 = 2 * randomGenerator.nextDouble() - 1;
            s = v1 * v1 + v2 * v2;
        }
        if (s != 0) {
            s = FastMath.sqrt(-2 * FastMath.log(s) / s);
        }
        cachedNormalDeviate = v2 * s;
        return v1 * s;
    }

    /**
     * Returns the next pseudo-random beta-distributed value with the given shape and scale parameters.
     *
     * @param alpha First shape parameter (must be positive).
     * @param beta Second shape parameter (must be positive).
     * @return beta-distributed random deviate
     */
    public double nextBeta(double alpha, double beta) {
        // TODO: validate parameters
        final double a = FastMath.min(alpha, beta);
        final double b = FastMath.max(alpha, beta);

        if (a > 1) {
            return algorithmBB(alpha, a, b);
        } else {
            return algorithmBC(alpha, b, a);
        }
    }

    /**
     * Returns the next pseudo-random, exponentially distributed deviate.
     *
     * @param mean mean of the exponential distribution
     * @return exponentially distributed deviate about the given mean
     */
    public double nextExponential(double mean) {
        // TODO: validate parameter
        // Step 1:
        double a = 0;
        double u = randomGenerator.nextDouble();

        // Step 2 and 3:
        while (u < 0.5) {
            a += EXPONENTIAL_SA_QI[0];
            u *= 2;
        }

        // Step 4 (now u >= 0.5):
        u += u - 1;

        // Step 5:
        if (u <= EXPONENTIAL_SA_QI[0]) {
            return mean * (a + u);
        }

        // Step 6:
        int i = 0; // Should be 1, be we iterate before it in while using 0
        double u2 = randomGenerator.nextDouble();
        double umin = u2;

        // Step 7 and 8:
        do {
            ++i;
            u2 = randomGenerator.nextDouble();

            if (u2 < umin) {
                umin = u2;
            }

            // Step 8:
        } while (u > EXPONENTIAL_SA_QI[i]); // Ensured to exit since EXPONENTIAL_SA_QI[MAX] = 1

        return mean * (a + umin * EXPONENTIAL_SA_QI[0]);
    }

    /**
     * Returns the next pseudo-random gamma-distributed value with the given shape and scale parameters.
     *
     * @param shape shape parameter of the distribution
     * @param scale scale parameter of the distribution
     * @return gamma-distributed random deviate
     */
    public double nextGamma(double shape, double scale) {
        if (shape < 1) {
            // [1]: p. 228, Algorithm GS

            while (true) {
                // Step 1:
                final double u = randomGenerator.nextDouble();
                final double bGS = 1 + shape / FastMath.E;
                final double p = bGS * u;

                if (p <= 1) {
                    // Step 2:

                    final double x = FastMath.pow(p, 1 / shape);
                    final double u2 = randomGenerator.nextDouble();

                    if (u2 > FastMath.exp(-x)) {
                        // Reject
                        continue;
                    } else {
                        return scale * x;
                    }
                } else {
                    // Step 3:

                    final double x = -1 * FastMath.log((bGS - p) / shape);
                    final double u2 = randomGenerator.nextDouble();

                    if (u2 > FastMath.pow(x, shape - 1)) {
                        // Reject
                        continue;
                    } else {
                        return scale * x;
                    }
                }
            }
        }

        // Now shape >= 1

        final double d = shape - 0.333333333333333333;
        final double c = 1 / (3 * FastMath.sqrt(d));

        while (true) {
            final double x = randomGenerator.nextGaussian();
            final double v = (1 + c * x) * (1 + c * x) * (1 + c * x);

            if (v <= 0) {
                continue;
            }

            final double x2 = x * x;
            final double u = randomGenerator.nextDouble();

            // Squeeze
            if (u < 1 - 0.0331 * x2 * x2) {
                return scale * d * v;
            }

            if (FastMath.log(u) < 0.5 * x2 + d * (1 - v + FastMath.log(v))) {
                return scale * d * v;
            }
        }
    }

    /**
     * Returns the next normally-distributed pseudo-random deviate.
     *
     * @param mean mean of the normal distribution
     * @param standardDeviation standard deviation of the normal distribution
     * @return a random value, normally distributed with the given mean and standard deviation
     */
    public double nextNormal(double mean, double standardDeviation) {
        //TODO: Check parameters
        return standardDeviation * nextGaussian() + mean;
    }

    /**
     * Returns a random deviate from the given distribution.
     *
     * @param dist the distribution to sample from
     * @return a random value following the given distribution
     */
    public double nextDeviate(RealDistribution dist) {
        if (CONTINUOUS_NAMES.contains(dist.getClass().getName())) {
            return callContinuousMethod(dist);
        }
        return dist.inverseCumulativeProbability(randomGenerator.nextDouble());
    }

    /**
     * Returns a random deviate from the given distribution.
     *
     * @param dist the distribution to sample from
     * @return a random value following the given distribution
     */
    public int nextDeviate(IntegerDistribution dist) {
        if (DISCRETE_NAMES.contains(dist.getClass().getName())) {
            return callDiscreteMethod(dist);
        }
        return dist.inverseCumulativeProbability(randomGenerator.nextDouble());
    }

    /**
     * Returns one Beta sample using Cheng's BB algorithm, when both &alpha; and &beta; are greater than 1.
     * @param random random generator to use
     * @param a0 distribution first shape parameter (&alpha;)
     * @param a min(&alpha;, &beta;) where &alpha;, &beta; are the two distribution shape parameters
     * @param b max(&alpha;, &beta;) where &alpha;, &beta; are the two distribution shape parameters
     * @return sampled value
     */
    private double algorithmBB(
                                      final double a0,
                                      final double a,
                                      final double b) {
        final double alpha = a + b;
        final double beta = FastMath.sqrt((alpha - 2.) / (2. * a * b - alpha));
        final double gamma = a + 1. / beta;

        double r;
        double w;
        double t;
        do {
            final double u1 = randomGenerator.nextDouble();
            final double u2 = randomGenerator.nextDouble();
            final double v = beta * (FastMath.log(u1) - FastMath.log1p(-u1));
            w = a * FastMath.exp(v);
            final double z = u1 * u1 * u2;
            r = gamma * v - 1.3862944;
            final double s = a + r - w;
            if (s + 2.609438 >= 5 * z) {
                break;
            }

            t = FastMath.log(z);
            if (s >= t) {
                break;
            }
        } while (r + alpha * (FastMath.log(alpha) - FastMath.log(b + w)) < t);

        w = FastMath.min(w, Double.MAX_VALUE);
        return Precision.equals(a, a0) ? w / (b + w) : b / (b + w);
    }

    /**
     * Returns one Beta sample using Cheng's BC algorithm, when at least one of &alpha; and &beta;
     * is smaller than 1.
     *
     * @param random random generator to use
     * @param a0 distribution first shape parameter (&alpha;)
     * @param a max(&alpha;, &beta;) where &alpha;, &beta; are the two distribution shape parameters
     * @param b min(&alpha;, &beta;) where &alpha;, &beta; are the two distribution shape parameters
     * @return sampled value
     */
    private double algorithmBC(
                                      final double a0,
                                      final double a,
                                      final double b) {
        final double alpha = a + b;
        final double beta = 1. / b;
        final double delta = 1. + a - b;
        final double k1 = delta * (0.0138889 + 0.0416667 * b) / (a * beta - 0.777778);
        final double k2 = 0.25 + (0.5 + 0.25 / delta) * b;

        double w;
        for (;;) {
            final double u1 = randomGenerator.nextDouble();
            final double u2 = randomGenerator.nextDouble();
            final double y = u1 * u2;
            final double z = u1 * y;
            if (u1 < 0.5) {
                if (0.25 * u2 + z - y >= k1) {
                    continue;
                }
            } else {
                if (z <= 0.25) {
                    final double v = beta * (FastMath.log(u1) - FastMath.log1p(-u1));
                    w = a * FastMath.exp(v);
                    break;
                }

                if (z >= k2) {
                    continue;
                }
            }

            final double v = beta * (FastMath.log(u1) - FastMath.log1p(-u1));
            w = a * FastMath.exp(v);
            if (alpha * (FastMath.log(alpha) - FastMath.log(b + w) + v) - 1.3862944 >= FastMath.log(z)) {
                break;
            }
        }

        w = FastMath.min(w, Double.MAX_VALUE);
        return Precision.equals(a, a0) ? w / (b + w) : b / (b + w);
    }

    private double callContinuousMethod(RealDistribution dist) {
        // Ugly if with InstanceOf checks and delegation
        return 0;
    }

    private int callDiscreteMethod(IntegerDistribution dist) {
     // Ugly if with InstanceOf checks and delegation
        return 0;
    }

}
