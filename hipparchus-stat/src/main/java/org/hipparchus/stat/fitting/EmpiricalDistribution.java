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

package org.hipparchus.stat.fitting;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.hipparchus.distribution.RealDistribution;
import org.hipparchus.distribution.continuous.AbstractRealDistribution;
import org.hipparchus.distribution.continuous.ConstantRealDistribution;
import org.hipparchus.distribution.continuous.NormalDistribution;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.stat.descriptive.StatisticalSummary;
import org.hipparchus.stat.descriptive.StreamingStatistics;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/**
 * <p>Represents an <a href="http://http://en.wikipedia.org/wiki/Empirical_distribution_function">
 * empirical probability distribution</a> -- a probability distribution derived
 * from observed data without making any assumptions about the functional form
 * of the population distribution that the data come from.</p>
 *
 * <p>An <code>EmpiricalDistribution</code> maintains data structures, called
 * <i>distribution digests</i>, that describe empirical distributions and
 * support the following operations:
 * </p>
 * <ul>
 * <li>loading the distribution from a file of observed data values</li>
 * <li>dividing the input data into "bin ranges" and reporting bin frequency
 *     counts (data for histogram)</li>
 * <li>reporting univariate statistics describing the full set of data values
 *     as well as the observations within each bin</li>
 * <li>generating random values from the distribution</li>
 * </ul>
 * <p>
 * Applications can use <code>EmpiricalDistribution</code> to build grouped
 * frequency histograms representing the input data or to generate random values
 * "like" those in the input file -- i.e., the values generated will follow the
 * distribution of the values in the file.
 * </p>
 *
 * <p>The implementation uses what amounts to the
 * <a href="http://nedwww.ipac.caltech.edu/level5/March02/Silverman/Silver2_6.html">
 * Variable Kernel Method</a> with Gaussian smoothing:</p>
 * <p><strong>Digesting the input file</strong></p>
 * <ol><li>Pass the file once to compute min and max.</li>
 * <li>Divide the range from min-max into <code>binCount</code> "bins."</li>
 * <li>Pass the data file again, computing bin counts and univariate
 *     statistics (mean, std dev.) for each of the bins </li>
 * <li>Divide the interval (0,1) into subintervals associated with the bins,
 *     with the length of a bin's subinterval proportional to its count.</li></ol>
 * <strong>Generating random values from the distribution</strong><ol>
 * <li>Generate a uniformly distributed value in (0,1) </li>
 * <li>Select the subinterval to which the value belongs.
 * <li>Generate a random Gaussian value with mean = mean of the associated
 *     bin and std dev = std dev of associated bin.</li></ol>
 *
 * <p>EmpiricalDistribution implements the {@link RealDistribution} interface
 * as follows.  Given x within the range of values in the dataset, let B
 * be the bin containing x and let K be the within-bin kernel for B.  Let P(B-)
 * be the sum of the probabilities of the bins below B and let K(B) be the
 * mass of B under K (i.e., the integral of the kernel density over B).  Then
 * set P(X &lt; x) = P(B-) + P(B) * K(x) / K(B) where K(x) is the kernel distribution
 * evaluated at x. This results in a cdf that matches the grouped frequency
 * distribution at the bin endpoints and interpolates within bins using
 * within-bin kernels.</p>
 *
 *<p><strong>USAGE NOTES:</strong></p>
 *<ul>
 *<li>The <code>binCount</code> is set by default to 1000.  A good rule of thumb
 *    is to set the bin count to approximately the length of the input file divided
 *    by 10. </li>
 *<li>The input file <i>must</i> be a plain text file containing one valid numeric
 *    entry per line.</li>
 * </ul>
 *
 */
public class EmpiricalDistribution extends AbstractRealDistribution {

    /** Default bin count */
    public static final int DEFAULT_BIN_COUNT = 1000;

    /** Character set for file input */
    private static final String FILE_CHARSET = "US-ASCII";

    /** Serializable version identifier */
    private static final long serialVersionUID = 5729073523949762654L;

    /** RandomDataGenerator instance to use in repeated calls to getNext() */
    protected final RandomDataGenerator randomData;

    /** List of SummaryStatistics objects characterizing the bins */
    private final List<StreamingStatistics> binStats;

    /** Sample statistics */
    private StreamingStatistics sampleStats;

    /** Max loaded value */
    private double max = Double.NEGATIVE_INFINITY;

    /** Min loaded value */
    private double min = Double.POSITIVE_INFINITY;

    /** Grid size */
    private double delta;

    /** number of bins */
    private final int binCount;

    /** is the distribution loaded? */
    private boolean loaded;

    /** upper bounds of subintervals in (0,1) "belonging" to the bins */
    private double[] upperBounds;

    /**
     * Creates a new EmpiricalDistribution with the default bin count.
     */
    public EmpiricalDistribution() {
        this(DEFAULT_BIN_COUNT);
    }

    /**
     * Creates a new EmpiricalDistribution with the specified bin count.
     *
     * @param binCount number of bins. Must be strictly positive.
     * @throws MathIllegalArgumentException if {@code binCount <= 0}.
     */
    public EmpiricalDistribution(int binCount) {
        this(binCount, new RandomDataGenerator());
    }

    /**
     * Creates a new EmpiricalDistribution with the specified bin count using the
     * provided {@link RandomGenerator} as the source of random data.
     *
     * @param binCount number of bins. Must be strictly positive.
     * @param generator random data generator (may be null, resulting in default JDK generator)
     * @throws MathIllegalArgumentException if {@code binCount <= 0}.
     */
    public EmpiricalDistribution(int binCount, RandomGenerator generator) {
        this(binCount, RandomDataGenerator.of(generator));
    }

    /**
     * Creates a new EmpiricalDistribution with default bin count using the
     * provided {@link RandomGenerator} as the source of random data.
     *
     * @param generator random data generator (may be null, resulting in default JDK generator)
     */
    public EmpiricalDistribution(RandomGenerator generator) {
        this(DEFAULT_BIN_COUNT, generator);
    }

    /**
     * Private constructor to allow lazy initialisation of the RNG contained
     * in the {@link #randomData} instance variable.
     *
     * @param binCount number of bins. Must be strictly positive.
     * @param randomData Random data generator.
     * @throws MathIllegalArgumentException if {@code binCount <= 0}.
     */
    private EmpiricalDistribution(int binCount,
                                  RandomDataGenerator randomData) {
        if (binCount <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL_BOUND_EXCLUDED,
                                                   binCount, 0);
        }
        this.binCount = binCount;
        this.randomData = randomData;
        binStats = new ArrayList<>();
    }

    /**
     * Computes the empirical distribution from the provided
     * array of numbers.
     *
     * @param in the input data array
     * @exception NullArgumentException if in is null
     */
    public void load(double[] in) throws NullArgumentException {
        try (DataAdapter da = new ArrayDataAdapter(in)) {
            da.computeStats();
            // new adapter for the second pass
            fillBinStats(new ArrayDataAdapter(in));
        } catch (IOException ex) {
            // Can't happen
            throw MathRuntimeException.createInternalError(); // NOPMD - cannot happen
        }
        loaded = true;

    }

    /**
     * Computes the empirical distribution using data read from a URL.
     *
     * <p>The input file <i>must</i> be an ASCII text file containing one
     * valid numeric entry per line.</p>
     *
     * @param url url of the input file
     *
     * @throws IOException if an IO error occurs
     * @throws NullArgumentException if url is null
     * @throws MathIllegalArgumentException if URL contains no data
     */
    public void load(URL url) throws IOException, MathIllegalArgumentException, NullArgumentException {
        MathUtils.checkNotNull(url);
        Charset charset = Charset.forName(FILE_CHARSET);
        try (InputStream       is1  = url.openStream();
             InputStreamReader isr1 = new InputStreamReader(is1, charset);
             BufferedReader    br1  = new BufferedReader(isr1);
             DataAdapter       da1  = new StreamDataAdapter(br1)) {
            da1.computeStats();
            if (sampleStats.getN() == 0) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.URL_CONTAINS_NO_DATA, url);
            }
            // new adapter for the second pass
            try (InputStream       is2  = url.openStream();
                 InputStreamReader isr2 = new InputStreamReader(is2, charset);
                 BufferedReader    br2  = new BufferedReader(isr2);
                 DataAdapter       da2  = new StreamDataAdapter(br2)) {
                fillBinStats(da2);
                loaded = true;
            }
        }
    }

    /**
     * Computes the empirical distribution from the input file.
     *
     * <p>The input file <i>must</i> be an ASCII text file containing one
     * valid numeric entry per line.</p>
     *
     * @param file the input file
     * @throws IOException if an IO error occurs
     * @throws NullArgumentException if file is null
     */
    public void load(File file) throws IOException, NullArgumentException {
        MathUtils.checkNotNull(file);
        Charset charset = Charset.forName(FILE_CHARSET);
        try (InputStream    is1 = Files.newInputStream(file.toPath());
             BufferedReader br1 = new BufferedReader(new InputStreamReader(is1, charset));
             DataAdapter    da1 = new StreamDataAdapter(br1)) {
            da1.computeStats();
            // new adapter for second pass
            try (InputStream    is2 = Files.newInputStream(file.toPath());
                 BufferedReader in2 = new BufferedReader(new InputStreamReader(is2, charset));
                 DataAdapter    da2 = new StreamDataAdapter(in2)) {
                fillBinStats(da2);
            }
            loaded = true;
        }
    }

    /**
     * Provides methods for computing <code>sampleStats</code> and
     * <code>beanStats</code> abstracting the source of data.
     */
    private abstract static class DataAdapter implements Closeable {

        /**
         * Compute bin stats.
         *
         * @throws IOException  if an error occurs computing bin stats
         */
        public abstract void computeBinStats() throws IOException;

        /**
         * Compute sample statistics.
         *
         * @throws IOException if an error occurs computing sample stats
         */
        public abstract void computeStats() throws IOException;

    }

    /**
     * <code>DataAdapter</code> for data provided through some input stream
     */
    private class StreamDataAdapter extends DataAdapter {

        /** Input stream providing access to the data */
        private BufferedReader inputStream;

        /**
         * Create a StreamDataAdapter from a BufferedReader
         *
         * @param in BufferedReader input stream
         */
        StreamDataAdapter(BufferedReader in) {
            inputStream = in;
        }

        /** {@inheritDoc} */
        @Override
        public void computeBinStats() throws IOException {
            for (String str = inputStream.readLine(); str != null; str = inputStream.readLine()) {
                final double val = Double.parseDouble(str);
                StreamingStatistics stats = binStats.get(findBin(val));
                stats.addValue(val);
            }
        }

        /** {@inheritDoc} */
        @Override
        public void computeStats() throws IOException {
            sampleStats = new StreamingStatistics();
            for (String str = inputStream.readLine(); str != null; str = inputStream.readLine()) {
                final double val = Double.parseDouble(str);
                sampleStats.addValue(val);
            }
        }

        /** {@inheritDoc} */
        @Override
        public void close() throws IOException {
            if (inputStream != null) {
                inputStream.close();
                inputStream = null;
            }
        }

    }

    /**
     * <code>DataAdapter</code> for data provided as array of doubles.
     */
    private class ArrayDataAdapter extends DataAdapter {

        /** Array of input  data values */
        private final double[] inputArray;

        /**
         * Construct an ArrayDataAdapter from a double[] array
         *
         * @param in double[] array holding the data, a reference to the array will be stored
         * @throws NullArgumentException if in is null
         */
        ArrayDataAdapter(double[] in) throws NullArgumentException {
            super();
            MathUtils.checkNotNull(in);
            inputArray = in; // NOPMD - storing a reference to the array is intentional and documented here
        }

        /** {@inheritDoc} */
        @Override
        public void computeStats() {
            sampleStats = new StreamingStatistics();
            for (double v : inputArray) {
                sampleStats.addValue(v);
            }
        }

        /** {@inheritDoc} */
        @Override
        public void computeBinStats() {
            for (double v : inputArray) {
                StreamingStatistics stats =
                        binStats.get(findBin(v));
                stats.addValue(v);
            }
        }

        /** {@inheritDoc} */
        @Override
        public void close() {
            // nothing to do
        }

    }

    /**
     * Fills binStats array (second pass through data file).
     *
     * @param da object providing access to the data
     * @throws IOException  if an IO error occurs
     */
    private void fillBinStats(final DataAdapter da)
        throws IOException {
        // Set up grid
        min = sampleStats.getMin();
        max = sampleStats.getMax();
        delta = (max - min)/binCount;

        // Initialize binStats ArrayList
        if (!binStats.isEmpty()) {
            binStats.clear();
        }
        for (int i = 0; i < binCount; i++) {
            StreamingStatistics stats = new StreamingStatistics();
            binStats.add(i,stats);
        }

        // Filling data in binStats Array
        da.computeBinStats();

        // Assign upperBounds based on bin counts
        upperBounds = new double[binCount];
        upperBounds[0] =
        ((double) binStats.get(0).getN()) / sampleStats.getN();
        for (int i = 1; i < binCount-1; i++) {
            upperBounds[i] = upperBounds[i-1] +
            ((double) binStats.get(i).getN()) / sampleStats.getN();
        }
        upperBounds[binCount-1] = 1.0d;
    }

    /**
     * Returns the index of the bin to which the given value belongs
     *
     * @param value  the value whose bin we are trying to find
     * @return the index of the bin containing the value
     */
    private int findBin(double value) {
        return FastMath.min(
                FastMath.max((int) FastMath.ceil((value - min) / delta) - 1, 0),
                binCount - 1);
    }

    /**
     * Generates a random value from this distribution.
     * <strong>Preconditions:</strong><ul>
     * <li>the distribution must be loaded before invoking this method</li></ul>
     * @return the random value.
     * @throws MathIllegalStateException if the distribution has not been loaded
     */
    public double getNextValue() throws MathIllegalStateException {

        if (!loaded) {
            throw new MathIllegalStateException(LocalizedCoreFormats.DISTRIBUTION_NOT_LOADED);
        }

        return inverseCumulativeProbability(randomData.nextDouble());
    }

    /**
     * Returns a {@link StatisticalSummary} describing this distribution.
     * <strong>Preconditions:</strong><ul>
     * <li>the distribution must be loaded before invoking this method</li></ul>
     *
     * @return the sample statistics
     * @throws IllegalStateException if the distribution has not been loaded
     */
    public StatisticalSummary getSampleStats() {
        return sampleStats;
    }

    /**
     * Returns the number of bins.
     *
     * @return the number of bins.
     */
    public int getBinCount() {
        return binCount;
    }

    /**
     * Returns a List of {@link StreamingStatistics} instances containing
     * statistics describing the values in each of the bins.  The list is
     * indexed on the bin number.
     *
     * @return List of bin statistics.
     */
    public List<StreamingStatistics> getBinStats() {
        return binStats;
    }

    /**
     * <p>Returns a fresh copy of the array of upper bounds for the bins.
     * Bins are: <br>
     * [min,upperBounds[0]],(upperBounds[0],upperBounds[1]],...,
     *  (upperBounds[binCount-2], upperBounds[binCount-1] = max].</p>
     *
     * @return array of bin upper bounds
     */
    public double[] getUpperBounds() {
        double[] binUpperBounds = new double[binCount];
        for (int i = 0; i < binCount - 1; i++) {
            binUpperBounds[i] = min + delta * (i + 1);
        }
        binUpperBounds[binCount - 1] = max;
        return binUpperBounds;
    }

    /**
     * <p>Returns a fresh copy of the array of upper bounds of the subintervals
     * of [0,1] used in generating data from the empirical distribution.
     * Subintervals correspond to bins with lengths proportional to bin counts.</p>
     *
     * <strong>Preconditions:</strong><ul>
     * <li>the distribution must be loaded before invoking this method</li></ul>
     *
     * @return array of upper bounds of subintervals used in data generation
     * @throws NullPointerException unless a {@code load} method has been
     * called beforehand.
     */
    public double[] getGeneratorUpperBounds() {
        int len = upperBounds.length;
        double[] out = new double[len];
        System.arraycopy(upperBounds, 0, out, 0, len);
        return out;
    }

    /**
     * Property indicating whether or not the distribution has been loaded.
     *
     * @return true if the distribution has been loaded
     */
    public boolean isLoaded() {
        return loaded;
    }

    /**
     * Reseeds the random number generator used by {@link #getNextValue()}.
     *
     * @param seed random generator seed
     */
    public void reSeed(long seed) {
        randomData.setSeed(seed);
    }

    // Distribution methods ---------------------------

    /**
     * {@inheritDoc}
     *
     * <p>Returns the kernel density normalized so that its integral over each bin
     * equals the bin mass.</p>
     *
     * <p>Algorithm description:</p>
     * <ol>
     * <li>Find the bin B that x belongs to.</li>
     * <li>Compute K(B) = the mass of B with respect to the within-bin kernel (i.e., the
     * integral of the kernel density over B).</li>
     * <li>Return k(x) * P(B) / K(B), where k is the within-bin kernel density
     * and P(B) is the mass of B.</li></ol>
     */
    @Override
    public double density(double x) {
        if (x < min || x > max) {
            return 0d;
        }
        final int binIndex = findBin(x);
        final RealDistribution kernel = getKernel(binStats.get(binIndex));
        double kernelDensity = kernel.density(x);
        if (kernelDensity == 0d) {
            return 0d;
        }
        return kernelDensity * pB(binIndex) / kB(binIndex);
    }

    /**
     * {@inheritDoc}
     *
     * <p>Algorithm description:</p>
     * <ol>
     * <li>Find the bin B that x belongs to.</li>
     * <li>Compute P(B) = the mass of B and P(B-) = the combined mass of the bins below B.</li>
     * <li>Compute K(B) = the probability mass of B with respect to the within-bin kernel
     * and K(B-) = the kernel distribution evaluated at the lower endpoint of B</li>
     * <li>Return P(B-) + P(B) * [K(x) - K(B-)] / K(B) where
     * K(x) is the within-bin kernel distribution function evaluated at x.</li></ol>
     * <p>If K is a constant distribution, we return P(B-) + P(B) (counting the full
     * mass of B).</p>
     *
     */
    @Override
    public double cumulativeProbability(double x) {
        if (x < min) {
            return 0d;
        } else if (x >= max) {
            return 1d;
        }
        final int binIndex = findBin(x);
        final double pBminus = pBminus(binIndex);
        final double pB = pB(binIndex);
        final RealDistribution kernel = k(x);
        if (kernel instanceof ConstantRealDistribution) {
            if (x < kernel.getNumericalMean()) {
                return pBminus;
            } else {
                return pBminus + pB;
            }
        }
        final double[] binBounds = getUpperBounds();
        final double kB = kB(binIndex);
        final double lower = binIndex == 0 ? min : binBounds[binIndex - 1];
        final double withinBinCum =
            (kernel.cumulativeProbability(x) -  kernel.cumulativeProbability(lower)) / kB;
        return pBminus + pB * withinBinCum;
    }

    /**
     * {@inheritDoc}
     *
     * <p>Algorithm description:</p>
     * <ol>
     * <li>Find the smallest i such that the sum of the masses of the bins
     *  through i is at least p.</li>
     * <li>
     *   Let K be the within-bin kernel distribution for bin i.<br>
     *   Let K(B) be the mass of B under K. <br>
     *   Let K(B-) be K evaluated at the lower endpoint of B (the combined
     *   mass of the bins below B under K).<br>
     *   Let P(B) be the probability of bin i.<br>
     *   Let P(B-) be the sum of the bin masses below bin i. <br>
     *   Let pCrit = p - P(B-)<br>
     * <li>Return the inverse of K evaluated at <br>
     *    K(B-) + pCrit * K(B) / P(B) </li>
     *  </ol>
     *
     */
    @Override
    public double inverseCumulativeProbability(final double p) throws MathIllegalArgumentException {
        MathUtils.checkRangeInclusive(p, 0, 1);

        if (p == 0.0) {
            return getSupportLowerBound();
        }

        if (p == 1.0) {
            return getSupportUpperBound();
        }

        int i = 0;
        while (cumBinP(i) < p) {
            i++;
        }

        final RealDistribution kernel = getKernel(binStats.get(i));
        final double kB = kB(i);
        final double[] binBounds = getUpperBounds();
        final double lower = i == 0 ? min : binBounds[i - 1];
        final double kBminus = kernel.cumulativeProbability(lower);
        final double pB = pB(i);
        final double pBminus = pBminus(i);
        final double pCrit = p - pBminus;
        if (pCrit <= 0) {
            return lower;
        }
        return kernel.inverseCumulativeProbability(kBminus + pCrit * kB / pB);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getNumericalMean() {
       return sampleStats.getMean();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getNumericalVariance() {
        return sampleStats.getVariance();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getSupportLowerBound() {
       return min;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getSupportUpperBound() {
        return max;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean isSupportConnected() {
        return true;
    }

    /**
     * Reseed the underlying PRNG.
     *
     * @param seed new seed value
     */
    public void reseedRandomGenerator(long seed) {
        randomData.setSeed(seed);
    }

    /**
     * The probability of bin i.
     *
     * @param i the index of the bin
     * @return the probability that selection begins in bin i
     */
    private double pB(int i) {
        return i == 0 ? upperBounds[0] :
            upperBounds[i] - upperBounds[i - 1];
    }

    /**
     * The combined probability of the bins up to but not including bin i.
     *
     * @param i the index of the bin
     * @return the probability that selection begins in a bin below bin i.
     */
    private double pBminus(int i) {
        return i == 0 ? 0 : upperBounds[i - 1];
    }

    /**
     * Mass of bin i under the within-bin kernel of the bin.
     *
     * @param i index of the bin
     * @return the difference in the within-bin kernel cdf between the
     * upper and lower endpoints of bin i
     */
    private double kB(int i) {
        final double[] binBounds = getUpperBounds();
        final RealDistribution kernel = getKernel(binStats.get(i));
        return i == 0 ? kernel.probability(min, binBounds[0]) :
            kernel.probability(binBounds[i - 1], binBounds[i]);
    }

    /**
     * The within-bin kernel of the bin that x belongs to.
     *
     * @param x the value to locate within a bin
     * @return the within-bin kernel of the bin containing x
     */
    private RealDistribution k(double x) {
        final int binIndex = findBin(x);
        return getKernel(binStats.get(binIndex));
    }

    /**
     * The combined probability of the bins up to and including binIndex.
     *
     * @param binIndex maximum bin index
     * @return sum of the probabilities of bins through binIndex
     */
    private double cumBinP(int binIndex) {
        return upperBounds[binIndex];
    }

    /**
     * The within-bin smoothing kernel. Returns a Gaussian distribution
     * parameterized by {@code bStats}, unless the bin contains less than 2
     * observations, in which case a constant distribution is returned.
     *
     * @param bStats summary statistics for the bin
     * @return within-bin kernel parameterized by bStats
     */
    protected RealDistribution getKernel(StreamingStatistics bStats) {
        if (bStats.getN() < 2 || bStats.getVariance() == 0) {
            return new ConstantRealDistribution(bStats.getMean());
        } else {
            return new NormalDistribution(bStats.getMean(),
                                          bStats.getStandardDeviation());
        }
    }
}
