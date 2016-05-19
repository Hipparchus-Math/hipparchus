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
import java.util.Collection;
import java.util.HashMap;

import org.hipparchus.distribution.IntegerDistribution;
import org.hipparchus.distribution.RealDistribution;
import org.hipparchus.distribution.continuous.BetaDistribution;
import org.hipparchus.distribution.continuous.GammaDistribution;
import org.hipparchus.distribution.continuous.LogNormalDistribution;
import org.hipparchus.distribution.continuous.NormalDistribution;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.CombinatoricsUtils;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.Precision;
import org.hipparchus.util.ResizableDoubleArray;

/**
 *
 * A class for generating random data.
 *
 */
public class RandomDataGenerator implements RandomGenerator, Serializable {

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

    /** Map of <classname, switch constant> for continuous distributions */
    private static final HashMap<String, Integer> CONTINUOUS_NAMES = new HashMap<String, Integer>();
    /** Map of <classname, switch constant> for discrete distributions */
    private static final HashMap<String, Integer> DISCRETE_NAMES = new HashMap<String, Integer>();
    /** beta distribution */
    private static final int BETA = 0;
    /** gamma distribution */
    private static final int GAMMA = 1;
    /** exponential distribution */
    private static final int EXPONENTIAL = 2;
    /** normal distribution */
    private static final int NORMAL = 3;
    /** poisson distribution */
    private static final int POISSON = 4;
    /** uniform integer distribution */
    private static final int UNIFORM_INT = 5;
    /** uniform real distribution */
    private static final int UNIFORM_REAL = 6;
    /** log normal distribution */
    private static final int LOG_NORMAL = 7;

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

        CONTINUOUS_NAMES.put("BetaDistribution", BETA);
        CONTINUOUS_NAMES.put("GammaDistribution", GAMMA);
        CONTINUOUS_NAMES.put("ExponentialDistribution", EXPONENTIAL);
        CONTINUOUS_NAMES.put("NormalDistribution", NORMAL);
        CONTINUOUS_NAMES.put("LogNormalDistribution", LOG_NORMAL);
        CONTINUOUS_NAMES.put("UniformDistribution", UNIFORM_REAL);
        DISCRETE_NAMES.put("PoissonDistribution", POISSON);
        DISCRETE_NAMES.put("UniformIntegerDistribution", UNIFORM_INT);
    }

    /**
     * Cached random normal value.  The default implementation for
     * {@link #nextGaussian} generates pairs of values and this field caches the
     * second value so that the full algorithm is not executed for every
     * activation.  The value {@code Double.NaN} signals that there is
     * no cached value.
     */
    private double cachedNormalDeviate = Double.NaN;

    /** Source of random data */
    private final RandomGenerator randomGenerator;

    /**
     * Construct a RandomDataGenerator with a default RandomGenerator as its source of random data.
     */
    public RandomDataGenerator() {
        this(new Well19937c());
    }

    /**
     * Construct a RandomDataGenerator using the given RandomGenerator as its source of random data.
     *
     * @param randomGenerator the underlying PRNG
     */
    public RandomDataGenerator(RandomGenerator randomGenerator) {
        this.randomGenerator = randomGenerator;
    }

    /**
     * Construct a RandomDataGenerator with a default RandomGenerator as its source of random data, initialized
     * with the given seed value.
     *
     * @param seed seed value
     */
    public RandomDataGenerator(long seed) {
        this(new Well19937c(seed));
    }

    /**
     * Construct a RandomDataGenerator using the given RandomGenerator as its source of random data, initialized
     * with the given seed.
     *
     * @param randomGenerator the underlying PRNG
     * @param seed seed for the PRNG
     */
    public RandomDataGenerator(long seed, RandomGenerator randomGenerator) {
        this.randomGenerator = randomGenerator;
        randomGenerator.setSeed(seed);
    }

    /**
     * Factory method to create a {@code RandomData} instance using the supplied
     * {@code RandomGenerator}.
     *
     * @param randomGenerator source of random bits
     * @return a RandomData using the given RandomGenerator to source bits
     */
    public static RandomDataGenerator of(RandomGenerator randomGenerator) {
        return new RandomDataGenerator(randomGenerator);
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

    /** {@inheritDoc} */
    @Override
    public void setSeed(int seed) {
        if (randomGenerator != null) {
            randomGenerator.setSeed(seed);
        }
    }

    /** {@inheritDoc} */
    @Override
    public void setSeed(int[] seed) {
        if (randomGenerator != null) {
            randomGenerator.setSeed(seed);
        }
    }

    /** {@inheritDoc} */
    @Override
    public void setSeed(long seed) {
        if (randomGenerator != null) {
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
     * activation.</p>
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
        if (mean <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.MEAN, mean);
        }
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
        if (standardDeviation <= 0) {
            throw new MathIllegalArgumentException (LocalizedCoreFormats.NUMBER_TOO_SMALL, mean, 0);
        }
        return standardDeviation * nextGaussian() + mean;
    }

    /**
     * Returns the next log-normally-distributed pseudo-random deviate.
     *
     * @param shape shape parameter of the log-normal distribution
     * @param scale scale parameter of the log-normal distribution
     * @return a random value, normally distributed with the given mean and standard deviation
     */
    public double nextLogNormal(double shape, double scale) {
        if (shape <= 0) {
            throw new MathIllegalArgumentException (LocalizedCoreFormats.NUMBER_TOO_SMALL, shape, 0);
        }
        return FastMath.exp(scale + shape * nextGaussian());
    }

    /**
     * Returns a poisson-distributed deviate with the given mean.
     *
     * @param mean expected value
     * @return poisson deviate
     * @throws MathIllegalArgumentException if mean is not strictly positive
     */
    public int nextPoisson(double mean) {
        if (mean <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL, mean, 0);
        }
        final double pivot = 40.0d;
        if (mean < pivot) {
            double p = FastMath.exp(-mean);
            long n = 0;
            double r = 1.0d;
            double rnd = 1.0d;

            while (n < 1000 * mean) {
                rnd = randomGenerator.nextDouble();
                r *= rnd;
                if (r >= p) {
                    n++;
                } else {
                    return (int) FastMath.min(n, Integer.MAX_VALUE);
                }
            }
            return (int) FastMath.min(n, Integer.MAX_VALUE);
        } else {
            final double lambda = FastMath.floor(mean);
            final double lambdaFractional = mean - lambda;
            final double logLambda = FastMath.log(lambda);
            final double logLambdaFactorial = CombinatoricsUtils.factorialLog((int) lambda);
            final long y2 = lambdaFractional < Double.MIN_VALUE ? 0 : nextPoisson(lambdaFractional);
            final double delta = FastMath.sqrt(lambda * FastMath.log(32 * lambda / FastMath.PI + 1));
            final double halfDelta = delta / 2;
            final double twolpd = 2 * lambda + delta;
            final double a1 = FastMath.sqrt(FastMath.PI * twolpd) * FastMath.exp(1 / (8 * lambda));
            final double a2 = (twolpd / delta) * FastMath.exp(-delta * (1 + delta) / twolpd);
            final double aSum = a1 + a2 + 1;
            final double p1 = a1 / aSum;
            final double p2 = a2 / aSum;
            final double c1 = 1 / (8 * lambda);

            double x = 0;
            double y = 0;
            double v = 0;
            int a = 0;
            double t = 0;
            double qr = 0;
            double qa = 0;
            for (;;) {
                final double u = randomGenerator.nextDouble();
                if (u <= p1) {
                    final double n = randomGenerator.nextGaussian();
                    x = n * FastMath.sqrt(lambda + halfDelta) - 0.5d;
                    if (x > delta || x < -lambda) {
                        continue;
                    }
                    y = x < 0 ? FastMath.floor(x) : FastMath.ceil(x);
                    final double e = nextExponential(1);
                    v = -e - (n * n / 2) + c1;
                } else {
                    if (u > p1 + p2) {
                        y = lambda;
                        break;
                    } else {
                        x = delta + (twolpd / delta) * nextExponential(1);
                        y = FastMath.ceil(x);
                        v = -nextExponential(1) - delta * (x + 1) / twolpd;
                    }
                }
                a = x < 0 ? 1 : 0;
                t = y * (y + 1) / (2 * lambda);
                if (v < -t && a == 0) {
                    y = lambda + y;
                    break;
                }
                qr = t * ((2 * y + 1) / (6 * lambda) - 1);
                qa = qr - (t * t) / (3 * (lambda + a * (y + 1)));
                if (v < qa) {
                    y = lambda + y;
                    break;
                }
                if (v > qr) {
                    continue;
                }
                if (v < y * logLambda - CombinatoricsUtils.factorialLog((int) (y + lambda)) + logLambdaFactorial) {
                    y = lambda + y;
                    break;
                }
            }
            return (int) FastMath.min(y2 + (long) y, Integer.MAX_VALUE);
        }
    }

    /**
     * Returns a random deviate from the given distribution.
     *
     * @param dist the distribution to sample from
     * @return a random value following the given distribution
     */
    public double nextDeviate(RealDistribution dist) {
        String className = dist.getClass().getSimpleName();
        Integer val = CONTINUOUS_NAMES.get(className);
        if (val != null) {
            switch(val) {
                case EXPONENTIAL:
                    return nextExponential(dist.getNumericalMean());
                case GAMMA:
                    GammaDistribution gammaDist = (GammaDistribution) dist;
                    return nextGamma(gammaDist.getShape(), gammaDist.getScale());
                case NORMAL:
                    NormalDistribution normalDist = (NormalDistribution) dist;
                    return nextNormal(normalDist.getMean(), normalDist.getStandardDeviation());
                case LOG_NORMAL:
                    LogNormalDistribution logNormalDist = (LogNormalDistribution) dist;
                    return nextLogNormal(logNormalDist.getShape(), logNormalDist.getScale());
                case BETA:
                    BetaDistribution betaDist = (BetaDistribution) dist;
                    return nextBeta(betaDist.getAlpha(), betaDist.getBeta());
                case UNIFORM_REAL:
                    return nextUniform(dist.getSupportLowerBound(), dist.getSupportUpperBound());
                default:
                    throw new MathRuntimeException(LocalizedCoreFormats.INTERNAL_ERROR);
            }
        }
        return dist.inverseCumulativeProbability(randomGenerator.nextDouble());
    }

    /**
     * Returns an array of random deviates from the given distribution.
     *
     * @param dist the distribution to sample from
     * @param size the number of values to return
     *
     * @return an array of {@code size }values following the given distribution
     */
    public double[] nextDeviates(RealDistribution dist, int size) {
        //TODO: check parameters
        double[] out = new double[size];
        for (int i = 0; i < size; i++) {
            out[i] = nextDeviate(dist);
        }
        return out;
    }

    /**
     * Returns a random deviate from the given distribution.
     *
     * @param dist the distribution to sample from
     * @return a random value following the given distribution
     */
    public int nextDeviate(IntegerDistribution dist) {
        String className = dist.getClass().getSimpleName();
        Integer val = DISCRETE_NAMES.get(className);
        if (val != null) {
            switch(val) {
                case POISSON:
                    return nextPoisson(dist.getNumericalMean());
                case UNIFORM_INT:
                    return nextInt(dist.getSupportLowerBound(), dist.getSupportUpperBound());

                default:
                    throw new MathRuntimeException(LocalizedCoreFormats.INTERNAL_ERROR);
            }
        }

        return dist.inverseCumulativeProbability(randomGenerator.nextDouble());
    }

    /**
     * Returns an array of random deviates from the given distribution.
     *
     * @param dist the distribution to sample from
     * @param size the number of values to return
     *
     * @return an array of {@code size }values following the given distribution
     */
    public int[] nextDeviates(IntegerDistribution dist, int size) {
        //TODO: check parameters
        int[] out = new int[size];
        for (int i = 0; i < size; i++) {
            out[i] = nextDeviate(dist);
        }
        return out;
    }

    /**
     * Returns one Beta sample using Cheng's BB algorithm, when both &alpha; and &beta; are greater than 1.
     *
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
     * Returns a Beta-distribute value using Cheng's BC algorithm, when at least one of &alpha; and &beta;
     * is smaller than 1.
     *
     * @param a0 distribution first shape parameter (&alpha;)
     * @param a max(&alpha;, &beta;) where &alpha;, &beta; are the two distribution shape parameters
     * @param b min(&alpha;, &beta;) where &alpha;, &beta; are the two distribution shape parameters
     * @return sampled value
     */
    private double algorithmBC(final double a0,
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

    /**
     * Returns a uniformly distributed random integer between lower and upper (inclusive).
     *
     * @param lower lower bound for the generated value
     * @param upper upper bound for the generated value
     * @return a random integer value within the given bounds
     * @throws MathIllegalArgumentException if lower is not strictly less than or equal to upper
     */
    public int nextInt(int lower, int upper) {
        if (lower >= upper) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.LOWER_BOUND_NOT_BELOW_UPPER_BOUND,
                                                lower, upper);
        }
        final int max = (upper - lower) + 1;
        if (max <= 0) {
            // The range is too wide to fit in a positive int (larger
            // than 2^31); as it covers more than half the integer range,
            // we use a simple rejection method.
            while (true) {
                final int r = nextInt();
                if (r >= lower &&
                    r <= upper) {
                    return r;
                }
            }
        } else {
            // We can shift the range and directly generate a positive int.
            return lower + nextInt(max);
        }
    }

    /**
     * Returns a uniformly distributed random long integer between lower and upper (inclusive).
     *
     * @param lower lower bound for the generated value
     * @param upper upper bound for the generated value
     * @return a random long integer value within the given bounds
     * @throws MathIllegalArgumentException if lower is not strictly less than or equal to upper
     */
    public long nextLong(final long lower, final long upper) throws MathIllegalArgumentException {
        if (lower >= upper) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.LOWER_BOUND_NOT_BELOW_UPPER_BOUND,
                                                lower, upper);
        }
        final long max = (upper - lower) + 1;
        if (max <= 0) {
            // the range is too wide to fit in a positive long (larger than 2^63); as it covers
            // more than half the long range, we use directly a simple rejection method
            while (true) {
                final long r = randomGenerator.nextLong();
                if (r >= lower && r <= upper) {
                    return r;
                }
            }
        } else if (max < Integer.MAX_VALUE){
            // we can shift the range and generate directly a positive int
            return lower + randomGenerator.nextInt((int) max);
        } else {
            // we can shift the range and generate directly a positive long
            return lower + nextLong(max);
        }
    }

    /**
     * Returns a double value uniformly distributed over [lower, upper]
     * @param lower lower bound
     * @param upper upper bound
     * @return uniform deviate
     * @throws MathIllegalArgumentException if upper is less than or equal to upper
     */
    public double nextUniform(double lower, double upper) {
        if (upper <= lower) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.LOWER_BOUND_NOT_BELOW_UPPER_BOUND, lower, upper);
        }
        if (Double.isInfinite(lower) || Double.isInfinite(upper)) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.INFINITE_BOUND);
        }
        if (Double.isNaN(lower) || Double.isNaN(upper)) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NAN_NOT_ALLOWED);
        }
        final double u = randomGenerator.nextDouble();
        return u * upper + (1 - u) * lower;
    }

    /**
     * Generates a random string of hex characters of length {@code len}.
     * <p>
     * The generated string will be random, but not cryptographically secure. To generate
     * cryptographically secure strings, use {@link #nextSecureHexString(int)}.
     * <p>
     * <strong>Algorithm Description:</strong> hex strings are generated using a
     * 2-step process.
     * <ol>
     * <li>{@code len / 2 + 1} binary bytes are generated using the underlying
     * Random</li>
     * <li>Each binary byte is translated into 2 hex digits</li>
     * </ol>
     * </p>
     *
     * @param len the desired string length.
     * @return the random string.
     * @throws MathIllegalArgumentException if {@code len <= 0}.
     */
    public String nextHexString(int len) throws MathIllegalArgumentException {
        if (len <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.LENGTH, len);
        }

        // Initialize output buffer
        StringBuilder outBuffer = new StringBuilder();

        // Get int(len/2)+1 random bytes
        byte[] randomBytes = new byte[(len / 2) + 1];
        randomGenerator.nextBytes(randomBytes);

        // Convert each byte to 2 hex digits
        for (int i = 0; i < randomBytes.length; i++) {
            Integer c = Integer.valueOf(randomBytes[i]);

            /*
             * Add 128 to byte value to make interval 0-255 before doing hex
             * conversion. This guarantees <= 2 hex digits from toHexString()
             * toHexString would otherwise add 2^32 to negative arguments.
             */
            String hex = Integer.toHexString(c.intValue() + 128);

            // Make sure we add 2 hex digits for each byte
            if (hex.length() == 1) {
                hex = "0" + hex;
            }
            outBuffer.append(hex);
        }
        return outBuffer.toString().substring(0, len);
    }

    /**
     * Generates an integer array of length {@code k} whose entries are selected
     * randomly, without repetition, from the integers {@code 0, ..., n - 1}
     * (inclusive).
     * <p>
     * Generated arrays represent permutations of {@code n} taken {@code k} at a
     * time.</p>
     * This method calls {@link MathArrays#shuffle(int[],RandomGenerator)
     * MathArrays.shuffle} in order to create a random shuffle of the set
     * of natural numbers {@code { 0, 1, ..., n - 1 }}.
     *
     *
     * @param n the domain of the permutation
     * @param k the size of the permutation
     * @return a random {@code k}-permutation of {@code n}, as an array of
     * integers
     * @throws MathIllegalArgumentException if {@code k > n}.
     * @throws MathIllegalArgumentException if {@code k <= 0}.
     */
    public int[] nextPermutation(int n, int k)
        throws MathIllegalArgumentException {
        if (k > n) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.PERMUTATION_EXCEEDS_N,
                                                k, n, true);
        }
        if (k <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.PERMUTATION_SIZE,
                                                   k);
        }

        final int[] index = MathArrays.natural(n);
        MathArrays.shuffle(index, randomGenerator);

        // Return a new array containing the first "k" entries of "index".
        return Arrays.copyOf(index, k);
    }

    /**
     * Returns an array of {@code k} objects selected randomly from the
     * Collection {@code c}.
     * <p>
     * Sampling from {@code c} is without replacement; but if {@code c} contains
     * identical objects, the sample may include repeats.  If all elements of
     * {@code c} are distinct, the resulting object array represents a
     * <a href="http://rkb.home.cern.ch/rkb/AN16pp/node250.html#SECTION0002500000000000000000">
     * Simple Random Sample</a> of size {@code k} from the elements of
     * {@code c}.</p>
     * <p>This method calls {@link #nextPermutation(int,int) nextPermutation(c.size(), k)}
     * in order to sample the collection.
     * </p>
     *
     * @param c the collection to be sampled
     * @param k the size of the sample
     * @return a random sample of {@code k} elements from {@code c}
     * @throws MathIllegalArgumentException if {@code k > c.size()}.
     * @throws MathIllegalArgumentException if {@code k <= 0}.
     */
    public Object[] nextSample(Collection<?> c, int k) throws MathIllegalArgumentException {

        int len = c.size();
        if (k > len) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.SAMPLE_SIZE_EXCEEDS_COLLECTION_SIZE,
                                                k, len, true);
        }
        if (k <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_OF_SAMPLES, k);
        }

        Object[] objects = c.toArray();
        int[] index = nextPermutation(len, k);
        Object[] result = new Object[k];
        for (int i = 0; i < k; i++) {
            result[i] = objects[index[i]];
        }
        return result;
    }

    /**
     * Returns an array of {@code k} double values selected randomly from the
     * double array {@code a}.
     * <p>
     * Sampling from {@code a} is without replacement; but if {@code a} contains
     * identical objects, the sample may include repeats.  If all elements of
     * {@code a} are distinct, the resulting object array represents a
     * <a href="http://rkb.home.cern.ch/rkb/AN16pp/node250.html#SECTION0002500000000000000000">
     * Simple Random Sample</a> of size {@code k} from the elements of
     * {@code a}.</p>
     *
     * @param a the array to be sampled
     * @param k the size of the sample
     * @return a random sample of {@code k} elements from {@code a}
     * @throws MathIllegalArgumentException if {@code k > c.size()}.
     * @throws MathIllegalArgumentException if {@code k <= 0}.
     */
    public double[] nextSample(double[] a, int k) throws MathIllegalArgumentException {
        int len = a.length;
        if (k > len) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.SAMPLE_SIZE_EXCEEDS_COLLECTION_SIZE,
                                                k, len, true);
        }
        if (k <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_OF_SAMPLES, k);
        }
        int[] index = nextPermutation(len, k);
        double[] result = new double[k];
        for (int i = 0; i < k; i++) {
            result[i] = a[index[i]];
        }
        return result;
    }

    /**
     * Returns a pseudorandom, uniformly distributed {@code long} value
     * between 0 (inclusive) and the specified value (exclusive), drawn from
     * this random number generator's sequence.
     *
     * @param n the bound on the random number to be returned.  Must be
     * positive.
     * @return  a pseudorandom, uniformly distributed {@code long}
     * value between 0 (inclusive) and n (exclusive).
     * @throws MathIllegalArgumentException  if n is not positive.
     */
    private long nextLong(final long n) throws MathIllegalArgumentException {
        if (n > 0) {
            final byte[] byteArray = new byte[8];
            long bits;
            long val;
            do {
                randomGenerator.nextBytes(byteArray);
                bits = 0;
                for (final byte b : byteArray) {
                    bits = (bits << 8) | (((long) b) & 0xffL);
                }
                bits &= 0x7fffffffffffffffL;
                val  = bits % n;
            } while (bits - val + (n - 1) < 0);
            return val;
        }
        throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL, n);
    }

}
