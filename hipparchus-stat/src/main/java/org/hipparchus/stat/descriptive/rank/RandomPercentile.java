/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.stat.descriptive.rank;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well19937c;
import org.hipparchus.stat.StatUtils;
import org.hipparchus.stat.descriptive.AbstractStorelessUnivariateStatistic;
import org.hipparchus.stat.descriptive.StorelessUnivariateStatistic;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;

/**
 * A {@link StorelessUnivariateStatistic} estimating percentiles using the
 * <a href=http://dimacs.rutgers.edu/~graham/pubs/papers/nquantiles.pdf>RANDOM</a>
 * Algorithm.
 * <p>
 * Storage requirements for the RANDOM algorithm depend on the desired
 * accuracy of quantile estimates. Quantile estimate accuracy is defined as follows.
 * <p>
 * Let \(X\) be the set of all data values consumed from the stream.
 * If \(\epsilon\) is the configured accuracy, \(\hat{q}\) is an estimated
 * quantile (what is returned by {@link #getResult()} or {@link #getResult(double)}),
 * \(rank(\hat{q}) = |{x \in X : x < \hat{q}}|\) is the rank of \(\hat{q}\) in the data and
 * \(n = |X|\) is the number of observations, then we can expect
 * \(rank(\hat{q}) / n < \epsilon\).
 * <p>
 * If the {@code epsilon} configuration parameter is set to \(\epsilon\), then
 * the algorithm maintains \(\left\lceil {log_{2}(1/\epsilon)}\right\rceil + 1\) buffers
 * of size \(\left\lceil {1/\epsilon \sqrt{log_2(1/\epsilon)}}\right\rceil\).  When
 * {@code epsilon} is set to the default value of \(10^{-4}\), this makes 15 buffers
 * of size 36,453.
 * <p>
 * The algorithm uses the buffers to maintain samples of data from the stream.  Up until
 * the point where all buffers are full, the entire sample is stored in the buffers.
 * If one of the {@code getResult} methods is called when all data are available in memory
 * and there is room to make a copy of the data (meaning the combined set of buffers is
 * less than half full), the {@code getResult} method delegates to a {@link Percentile}
 * instance to compute and return the exact value for the desired quantile.
 * For default {@code epsilon}, this means exact values will be returned whenever fewer than
 * \(\left\lceil {15 \times 36453 / 2} \right\rceil = 273,398\) values have been consumed
 * from the data stream.
 * <p>
 * When buffers become full, the algorithm merges buffers so that they effectively represent
 * a larger set of values than they can hold. Subsequently, data values are sampled from the
 * stream to fill buffers freed by merge operations.  Both the merging and the sampling
 * require random selection, which is done using a {@code RandomGenerator}.  To get
 * repeatable results for large data streams, users should provide {@code RandomGenerator}
 * instances with fixed seeds. {@code RandomPercentile} itself does not reseed or othewise
 * initialize the {@code RandomGenerator} provided to it.  By default, it uses a
 * {@link Well19937c} generator with the default seed.
 * <p>
 * Note: This implementation is not thread-safe.
 */
public class RandomPercentile
    extends AbstractStorelessUnivariateStatistic implements StorelessUnivariateStatistic, Serializable {

    private static final long serialVersionUID = 1L;
    /** Storage size of each buffer */
    private final int s;
    /** Number of buffers minus 1 */
    private final int h;
    /** Data structure used to manage buffers */
    private final BufferMap bufferPool;
    /** Default quantile (what {@link #getResult()} returns) - scaled 0-1 */
    private final double quantile;
    /** Bound on the quantile estimation error */
    private final double epsilon;
    /** Source of random data */
    private final RandomGenerator randomGenerator;
    /** Number of elements consumed from the input data stream */
    private long n = 0;
    /** Buffer currently being filled */
    private Buffer currentBuffer = null;
    /** Default quantile estimation error setting */
    public static final double DEFAULT_EPSILON = 1e-4;

    /**
     * Constructs a {@code RandomPercentile} with quantile estimation error
     * {@code epsilon} and default percentile {@code percentile}, using
     * {@code randomGenerator} as its source of random data.
     *
     * @param epsilon bound on quantile estimation error (see class javadoc)
     * @param randomGenerator PRNG used in sampling and merge operations
     * @param percentile default percentile (what is returned by {@link #getResult()}
     * @throws MathIllegalArgumentException if percentile is not in the range [0, 100]
     */
    public RandomPercentile(double epsilon, RandomGenerator randomGenerator, double percentile) {
        if (percentile > 100 || percentile < 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.OUT_OF_RANGE,
                                                   percentile, 0, 100);
        }
        if (epsilon <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL,
                                                   epsilon, 0);
        }
        this.h = (int) FastMath.ceil(log2(1/epsilon));
        this.s = (int) FastMath.ceil(FastMath.sqrt(log2(1/epsilon)) / epsilon);
        this.randomGenerator = randomGenerator;
        bufferPool = new BufferMap(h + 1, s, randomGenerator);
        currentBuffer = bufferPool.create(0);
        this.quantile = percentile / 100;
        this.epsilon = epsilon;
    }

    /**
     * Constructs a {@code RandomPercentile} with quantile estimation error
     * {@code epsilon} and default percentile {@code percentile}, using
     * the default PRNG as source of random data.
     *
     * @param epsilon bound on quantile estimation error (see class javadoc)
     * @param percentile default percentile (what is returned by {@link #getResult()}
     * @throws MathIllegalArgumentException if percentile is not in the range [0, 100]
     */
    public RandomPercentile(double epsilon, double percentile) {
        this(epsilon, new Well19937c(), percentile);
    }

    /**
     * Constructs a {@code RandomPercentile} with quantile estimation error
     * set to the default ({@link #DEFAULT_EPSILON}) and default percentile
     * {@code percentile}, using the default PRNG as source of random data.
     *
     * @param percentile default percentile (what is returned by {@link #getResult()}
     * @throws MathIllegalArgumentException if percentile is not in the range [0, 100]
     */
    public RandomPercentile(double percentile) {
        this(DEFAULT_EPSILON, new Well19937c(), percentile);
    }

    /**
     * Constructs a {@code RandomPercentile} with quantile estimation error
     * set to the default ({@link #DEFAULT_EPSILON}) and default percentile
     * {@code percentile}, using the {@code randomGenerator} as source of
     * random data.
     *
     * @param randomGenerator PRNG used in sampling and merge operations
     * @param percentile default percentile (what is returned by {@link #getResult()}
     * @throws MathIllegalArgumentException if percentile is not in the range [0, 100]
     */
    public RandomPercentile(double percentile, RandomGenerator randomGenerator) {
        this(DEFAULT_EPSILON, randomGenerator, percentile);
    }

    /**
     * Constructs a {@code RandomPercentile} with quantile estimation error
     * set to the default ({@link #DEFAULT_EPSILON}) and default percentile
     * equal to the median, using the default PRNG as source of random data.
     */
    public RandomPercentile() {
        this(DEFAULT_EPSILON, new Well19937c(), 50d);
    }

    /**
     * Copy constructor, creates a new {@code RandomPercentile} identical
     * to the {@code original}.  Note: the RandomGenerator used by the new
     * instance is referenced, not copied - i.e., the new instance shares
     * a generator with the original.
     *
     * @param original the {@code PSquarePercentile} instance to copy
     */
    public RandomPercentile(RandomPercentile original) {
        super();
        this.h = original.h;
        this.n = original.n;
        this.quantile = original.quantile;
        this.s = original.s;
        this.epsilon = original.epsilon;
        this.bufferPool = new BufferMap(original.bufferPool);
        this.randomGenerator = original.randomGenerator;
        Iterator<Buffer> iterator = bufferPool.iterator();
        Buffer current = null;
        Buffer curr = null;
        // See if there is a partially filled buffer - that will be currentBuffer
        while (current == null && iterator.hasNext()) {
            curr = iterator.next();
            if (curr.hasCapacity()) {
                current = curr;
            }
        }
        // If there is no partially filled buffer, just assign the last one.
        // Next increment() will find no capacity and create a new one or trigger
        // a merge.
        this.currentBuffer = current == null ? curr : current;
    }

    @Override
    public long getN() {
        return n;
    }

    public double getQuantile() {
        return quantile;
    }

    @Override
    public double evaluate(final double[] values, final int begin, final int length)
        throws MathIllegalArgumentException {
        if (MathArrays.verifyValues(values, begin, length)) {
            RandomPercentile randomPercentile = new RandomPercentile(this.epsilon,
                                                                     this.randomGenerator,
                                                                     this.quantile * 100);
            randomPercentile.incrementAll(values, begin, length);
            return randomPercentile.getResult();
        }
        return Double.NaN;
    }

    @Override
    public StorelessUnivariateStatistic copy() {
        return new RandomPercentile(this);
    }

    @Override
    public void clear() {
        n = 0;
        bufferPool.clear();
        currentBuffer = bufferPool.create(0);
    }

    @Override
    public double getResult() {
        return getResult(this.quantile * 100);

    }

    /**
     * Returns an estimate of the given percentile.
     *
     * @param percentile desired percentile (scaled 0 - 100)
     * @return estimated percentile
     * @throws MathIllegalArgumentException if percentile is out of the range [0, 100]
     */
    public double getResult(double percentile) {
        if (percentile > 100 || percentile < 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.OUT_OF_RANGE,
                                                   percentile, 0, 100);
        }
        // Convert to internal quantile scale
        final double q = percentile / 100;
        // First get global min and max to bound search.
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        double bMin;
        double bMax;
        Iterator<Buffer> bufferIterator = bufferPool.iterator();
        while (bufferIterator.hasNext()) {
            Buffer buffer = bufferIterator.next();
            bMin = StatUtils.min(buffer.getData());
            if (bMin < min) {
                min = bMin;
            }
            bMax = StatUtils.max(buffer.getData());
            if (bMax > max) {
                max = bMax;
            }
        }

        // Handle degenerate cases
        if (Double.compare(q, 0d) == 0 || n == 1) {
            return min;
        }
        if (Double.compare(q, 1) == 0) {
            return max;
        }
        if (n == 0) {
            return Double.NaN;
        }

        // See if we have all data in memory and enough free memory to copy.
        // If so, use Percentile to perform exact computation.
        if (bufferPool.halfEmpty()) {
            return new Percentile(percentile).evaluate(bufferPool.levelZeroData());
        }

        // Compute target rank
        final double targetRank = q * n;

        // Start with initial guess min + quantile * (max - min).
        double estimate = min + q * (max - min);
        double estimateRank = getRank(estimate);
        double lower;
        double upper;
        if (estimateRank == targetRank) {
            return estimate;
        }
        if (estimateRank > targetRank) {
            upper = estimate;
            lower = min;
        } else {
            lower = estimate;
            upper = max;
        }
        final double eps = epsilon / 2;
        double intervalWidth = FastMath.abs(upper - lower);
        while (FastMath.abs(estimateRank / n - q) > eps && intervalWidth > eps / n) {
            if (estimateRank == targetRank) {
                return estimate;
            }
            if (estimateRank > targetRank) {
                upper = estimate;
            } else {
                lower = estimate;
            }
            intervalWidth = FastMath.abs(upper - lower);
            estimate = lower + intervalWidth / 2;
            estimateRank = getRank(estimate);
        }
        return estimate;
    }

    /**
     * Gets the estimated rank of {@code value}, i.e.  \(|\{x \in X : x < value\}|\)
     * where \(X\) is the set of values that have been consumed from the stream.
     *
     * @param value value whose overall rank is sought
     * @return estimated number of sample values that are strictly less than {@code value}
     */
    public double getRank(double value) {
        double rankSum = 0;
        Iterator<Buffer> bufferIterator = bufferPool.iterator();
        while (bufferIterator.hasNext()) {
            Buffer buffer = bufferIterator.next();
            rankSum += buffer.rankOf(value) * FastMath.pow(2, buffer.level);
        }
        return rankSum;
    }

    /**
     * Returns the estimated quantile position of value in the dataset.
     * Specifically, what is returned is an estimate of \(|\{x \in X : x < value\}| / |X|\)
     * where \(X\) is the set of values that have been consumed from the stream.
     *
     * @param value value whose quantile rank is sought.
     * @return estimated proportion of sample values that are strictly less than {@code value}
     */
    public double getQuantileRank(double value) {
        return getRank(value) / getN();
    }

    @Override
    public void increment(double d) {
        n++;
        if (!currentBuffer.hasCapacity()) { // Need to get a new buffer to fill
            // First see if we have not yet created all the buffers
            if (bufferPool.canCreate()) {
                final int level = (int) Math.ceil(Math.max(0, log2(n/(s * FastMath.pow(2, h - 1)))));
                currentBuffer = bufferPool.create(level);
            } else { // All buffers have been created - need to merge to free one
                currentBuffer = bufferPool.merge();
            }
        }
        currentBuffer.consume(d);
    }

    private static class Buffer implements Serializable {
        private static final long serialVersionUID = 1L;
        /** Number of values actually stored in the buffer */
        private final int size;
        /** Data sampled from the stream */
        private final double[] data;
        /** PRNG used for merges and stream sampling */
        private final RandomGenerator randomGenerator;
        /** Level of the buffer */
        private int level = 0;
        /** Block size  = 2^level */
        private long blockSize;
        /** Next location in backing array for stored (taken) value */
        private int next = 0;
        /** Number of values consumed in current 2^level block of values from the stream */
        private long consumed = 0;
        /** Index of next value to take in current 2^level block */
        private long nextToTake = 0;

        public Buffer(int size, int level, RandomGenerator randomGenerator) {
            this.size = size;
            data = new double[size];
            this.level = level;
            this.randomGenerator = randomGenerator;
            computeBlockSize();
        }

        /**
         * Sets blockSize and nextToTake based on level.
         */
        private void computeBlockSize() {
            if (level == 0) {
                blockSize = 1;
            } else {
                long product = 1;
                for (int i = 0; i < level; i++) {
                    product *= 2;
                }
                blockSize = product;
            }
            if (blockSize > 1) {
                nextToTake = randomGenerator.nextLong(blockSize);
            }
        }

        /**
         * Consumes a value from the input stream.
         * <p>
         * For each 2^level values consumed, one is added to the buffer.
         * The buffer is not considered full until 2^level * size values
         * have been consumed.
         * <p>
         * Sorts the data array if the consumption renders the buffer full.
         *
         * @param value value to consume from the stream
         */
        public void consume(double value) {
            if (consumed == nextToTake) {
                data[next] = value;
                next++;
            }
            consumed++;
            if (consumed == blockSize) {
                if (next == size) {   // Buffer is full
                    Arrays.sort(data);
                } else {              // Reset in-block counter and nextToTake
                    consumed = 0;
                    if (blockSize > 1) {
                        nextToTake = randomGenerator.nextLong(blockSize);
                    }
                }
            }
        }

        /**
         * Merge this with other.  After the merge, this will be the merged
         * buffer and other will be free.  Both will have level+1. Post-merge,
         * other can be used to accept new data.
         *
         * @param other initially full other buffer at the same level as this.
         */
        public void mergeWith(Buffer other) {
            // Make sure both this and other are full and have the same level
            if (this.hasCapacity() || other.hasCapacity() || other.level != this.level) {
                throw new MathIllegalStateException(LocalizedCoreFormats.INTERNAL_ERROR);
            }
            // Randomly select one of the two entries for each slot
            for (int i = 0; i < size; i++) {
                if (randomGenerator.nextBoolean()) {
                    data[i] = other.data[i];
                }
            }
            // Re-sort data
            Arrays.sort(data);
            // Bump level of both buffers
            other.setLevel(level + 1);
            this.setLevel(level + 1);
            // Clear the free one (and compute new blocksize)
            other.clear();
        }

        /**
         * @return true if the buffer has capacity; false if it is full
         */
        public boolean hasCapacity() {
            // Buffer has capacity if it has not yet set all of its data
            // values or if it has but still has not finished its last block
            return next < size || consumed < blockSize;
        }

        /**
         * Sets the level of the buffer.
         *
         * @param level new level value
         */
        public void setLevel(int level) {
            this.level = level;
        }

        /**
         * Clears data, recomputes blockSize and resets consumed and nextToTake.
         */
        public void clear() {
            consumed = 0;
            next = 0;
            computeBlockSize();
        }

        /**
         * Returns a copy of the data that has been added to the buffer
         *
         * @return possibly unsorted copy of the portion of the buffer that has been filled
         */
        public double[] getData() {
            final double[] out = new double[next];
            System.arraycopy(data, 0, out, 0, next);
            return out;
        }

        /**
         * Returns the rank of value among the sampled values in this buffer.
         *
         * @param value value whose rank is sought
         * @return |{v in data : v < value}|
         */
        public int rankOf(double value) {
            int ret = 0;
            if (!hasCapacity()) { // Full sorted buffer, can do binary search
                ret = Arrays.binarySearch(data, value);
                if (ret < 0) {
                    return -ret - 1;
                } else {
                    return ret;
                }
            } else { // have to count - not sorted yet and can't sort yet
                for (int i = 0; i < next; i++) {
                    if (data[i] < value) {
                        ret++;
                    }
                }
                return ret;
            }
        }

        /**
         * @return the level of this buffer
         */
        public int getLevel() {
            return level;
        }
    }

    private static double log2(double x) {
        return Math.log(x) / Math.log(2);
    }

    /**
     * A map structure to hold the buffers.  Keys are levels
     * and values are lists of buffers at the given level.
     * Overall capacity is limited by the total number of buffers.
     */
    private static class BufferMap implements Iterable<Buffer>, Serializable {
        private static final long serialVersionUID = 1L;
        /** Total number of buffers that can be created - cap for count */
        private final int capacity;
        /** PRNG used in merges */
        private final RandomGenerator randomGenerator;
        /** Total count of all buffers */
        private int count = 0;
        /** Uniform buffer size */
        private final int bufferSize;
        /** Backing store for the buffer map */
        final HashMap<Integer,List<Buffer>> bufferMap = new HashMap<>();

        /**
         * Creates a BufferMap that can manage up to capacity buffers.
         * Buffers created by the pool with have size = buffersize.
         *
         * @param capacity cap on the number of buffers
         * @param bufferSize size of each buffer
         * @param randomGenerator RandomGenerator to use in merges
         */
        public BufferMap(int capacity, int bufferSize, RandomGenerator randomGenerator) {
            this.bufferSize = bufferSize;
            this.capacity = capacity;
            this.randomGenerator = randomGenerator;
        }

        /**
         * Copy constructor.
         *
         * @param BufferMap to copy
         */
        public BufferMap(BufferMap original) {
            super();
            this.bufferSize = original.bufferSize;
            this.capacity = original.capacity;
            this.count = 0;
            this.randomGenerator = original.randomGenerator;
            Iterator<Buffer> iterator = original.iterator();
            Buffer current = null;
            Buffer newCopy = null;
            while (iterator.hasNext()) {
                current = iterator.next();
                // Create and register a new buffer at the same level
                newCopy = create(current.getLevel());
                // Consume the data
                final double[] data = current.getData();
                for (double value : data) {
                    newCopy.consume(value);
                }
            }
        }

        /**
         * Tries to create a buffer with the given level.
         * <p>
         * If there is capacity to create a new buffer (i.e., fewer than
         * count have been created), a new buffer is created with the given
         * level, registered and returned.  If capacity has been reached,
         * null is returned.
         *
         * @return an empty buffer or null if a buffer can't be provided
         */
        public Buffer create(int level) {
            if (!canCreate()) {
                return null;
            }
            Buffer buffer = new Buffer(bufferSize, level, randomGenerator);
            List<Buffer> bufferList = bufferMap.get(level);
            if (bufferList == null) {
                bufferList = new ArrayList<Buffer>();
                bufferMap.put(level, bufferList);
            }
            bufferList.add(buffer);
            count++;
            return buffer;
        }

        /**
         * Returns true if there is capacity to create a new buffer.
         *
         * @return true if fewer than capacity buffers have been created.
         */
        public boolean canCreate() {
            return count < capacity;
        }

        /**
         * Returns true if we have used less than half of the allocated storage.
         * <p>
         * Includes a check to make sure all buffers have level 0;
         * but this should always be the case.
         * <p>
         * When this method returns true, we have all consumed data in storage
         * and enough space to make a copy of the combined dataset.
         *
         * @return true if all buffers have level 0 and less than half of the
         * available storage has been used
         */
        public boolean halfEmpty() {
            return count * 2 < capacity &&
                    bufferMap.size() == 1 &&
                    bufferMap.containsKey(0);
        }

        /**
         * Returns a fresh copy of all data from level 0 buffers.
         *
         * @return combined data stored in all level 0 buffers
         */
        public double[] levelZeroData() {
            List<Buffer> levelZeroBuffers = bufferMap.get(0);
            // First determine the combined size of the data
            int length = 0;
            for (Buffer buffer : levelZeroBuffers) {
                if (!buffer.hasCapacity()) { // full buffer
                    length += buffer.size;
                } else {
                    length += buffer.next;  // filled amount
                }
            }
            // Copy the data
            int pos = 0;
            int currLen;
            final double[] out = new double[length];
            for (Buffer buffer : levelZeroBuffers) {
                if (!buffer.hasCapacity()) {
                    currLen = buffer.size;
                } else {
                    currLen =  buffer.next;
                }
                System.arraycopy(buffer.data, 0, out, pos, currLen);
                pos += currLen;
            }
            return out;
        }

        /**
         * Finds the lowest level l where there exist at least two buffers,
         * merges them to create a new buffer with level l+1 and returns
         * a free buffer with level l+1.
         *
         * @return free buffer that can accept data
         */
        public Buffer merge() {
            int l = 0;
            List<Buffer> mergeCandidates = null;
            // Find the lowest level containing at least two buffers
            while (mergeCandidates == null && l < bufferMap.size()) {
                if (bufferMap.get(l).size() > 1) {
                    mergeCandidates = bufferMap.get(l);
                } else {
                    l++;
                }
            }
            if (mergeCandidates == null) {
                // Should never happen
                throw new MathIllegalStateException(LocalizedCoreFormats.INTERNAL_ERROR);
            }
            Buffer buffer1 = mergeCandidates.get(0);
            Buffer buffer2 = mergeCandidates.get(1);
            // Remove buffers to be merged
            mergeCandidates.remove(0);
            mergeCandidates.remove(0);
            // Merge the buffers
            buffer1.mergeWith(buffer2);
            // Now both buffers have level l+1; buffer1 is full and buffer2 is free.
            // Register both buffers
            List<Buffer> bufferList = bufferMap.get(l + 1);
            if (bufferList == null) {
                bufferList = new ArrayList<Buffer>();
                bufferMap.put(l + 1, bufferList);
            }
            bufferList.add(buffer1);
            bufferList.add(buffer2);
            // Return the free one
            return buffer2;
        }

        /**
         * Clears the buffer map.
         */
        public void clear() {
            for (List<Buffer> bufferList : bufferMap.values()) {
                bufferList.clear();
            }
            bufferMap.clear();
            count = 0;
        }

        /**
         * Returns an iterator over all of the buffers.
         */
        @Override
        public Iterator<Buffer> iterator() {
            Iterator<Buffer> it = new Iterator<Buffer>() {
                final Iterator<List<Buffer>> levelIterator = bufferMap.values().iterator();
                Iterator<Buffer> bufferIterator = levelIterator.next().iterator();

                @Override
                public boolean hasNext() {
                    return levelIterator.hasNext() || bufferIterator.hasNext();
                }

                @Override
                public Buffer next() {
                    if (!bufferIterator.hasNext()) {
                        bufferIterator = levelIterator.next().iterator();
                    }
                    return bufferIterator.next();
                }

                @Override
                public void remove() {
                    throw new UnsupportedOperationException();
                }
            };
            return it;
        }
    }
}
