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
package org.hipparchus.stat.descriptive;

import java.io.Serializable;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.stat.descriptive.moment.Variance;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.ResizableDoubleArray;


/**
 * Maintains a dataset of values of a single variable and computes descriptive
 * statistics based on stored data.
 * <p>
 * The {@link #getWindowSize() windowSize} property sets a limit on the number
 * of values that can be stored in the dataset. The default value, INFINITE_WINDOW,
 * puts no limit on the size of the dataset. This value should be used with
 * caution, as the backing store will grow without bound in this case.  For very
 * large datasets, {@link SummaryStatistics}, which does not store the dataset,
 * should be used instead of this class. If <code>windowSize</code> is not
 * INFINITE_WINDOW and more values are added than can be stored in the dataset,
 * new values are added in a "rolling" manner, with new values replacing the
 * "oldest" values in the dataset.
 * <p>
 * Note: this class is not threadsafe.
 */
class DescriptiveStatisticsImpl implements DescriptiveStatistics, Serializable {

    /** Serialization UID */
    private static final long serialVersionUID = 20160411L;

    /** The statistic used to calculate the population variance - fixed. */
    private static final UnivariateStatistic populationVariance = new Variance(false);

    /** Stored data values. */
    private final ResizableDoubleArray eDA;

    /** Maximum statistic implementation. */
    private final UnivariateStatistic maxImpl;
    /** Minimum statistic implementation. */
    private final UnivariateStatistic minImpl;
    /** Sum statistic implementation. */
    private final UnivariateStatistic sumImpl;
    /** Sum of squares statistic implementation. */
    private final UnivariateStatistic sumOfSquaresImpl;
    /** Mean statistic implementation. */
    private final UnivariateStatistic meanImpl;
    /** Variance statistic implementation. */
    private final UnivariateStatistic varianceImpl;
    /** Geometric mean statistic implementation. */
    private final UnivariateStatistic geometricMeanImpl;
    /** Kurtosis statistic implementation. */
    private final UnivariateStatistic kurtosisImpl;
    /** Skewness statistic implementation. */
    private final UnivariateStatistic skewnessImpl;
    /** Percentile statistic implementation. */
    private final QuantiledUnivariateStatistic percentileImpl;

    /** holds the window size. **/
    private int windowSize;

    /**
     * Copy constructor.
     * <p>
     * Construct a new DescriptiveStatistics instance that
     * is a copy of original.
     *
     * @param original DescriptiveStatistics instance to copy
     * @throws NullArgumentException if original is null
     */
    protected DescriptiveStatisticsImpl(DescriptiveStatisticsImpl original) {
        MathUtils.checkNotNull(original);

        // Copy data and window size
        this.eDA        = original.eDA.copy();
        this.windowSize = original.windowSize;

        // Copy implementations
        this.maxImpl           = original.maxImpl.copy();
        this.minImpl           = original.minImpl.copy();
        this.meanImpl          = original.meanImpl.copy();
        this.sumImpl           = original.sumImpl.copy();
        this.sumOfSquaresImpl  = original.sumOfSquaresImpl.copy();
        this.varianceImpl      = original.varianceImpl.copy();
        this.geometricMeanImpl = original.geometricMeanImpl.copy();
        this.kurtosisImpl      = original.kurtosisImpl.copy();
        this.skewnessImpl      = original.skewnessImpl.copy();
        this.percentileImpl    = original.percentileImpl.copy();
    }

    /**
     * Construct a new DescriptiveStatistics instance based
     * on the data provided by a builder.
     *
     * @param builder the builder to use.
     */
    protected DescriptiveStatisticsImpl(Builder builder) {
        this.maxImpl           = builder.maxImpl;
        this.minImpl           = builder.minImpl;
        this.meanImpl          = builder.meanImpl;
        this.sumImpl           = builder.sumImpl;
        this.sumOfSquaresImpl  = builder.sumsqImpl;
        this.varianceImpl      = builder.varianceImpl;
        this.geometricMeanImpl = builder.geometricMeanImpl;
        this.kurtosisImpl      = builder.kurtosisImpl;
        this.skewnessImpl      = builder.skewnessImpl;
        this.percentileImpl    = builder.percentileImpl;

        this.windowSize        = builder.windowSize;
        int initialCapacity    = this.windowSize < 0 ? 100 : this.windowSize;
        this.eDA               = builder.initialValues != null ?
            new ResizableDoubleArray(builder.initialValues) :
            new ResizableDoubleArray(initialCapacity);
    }

    /** {@inheritDoc} */
    @Override
    public DescriptiveStatistics copy() {
        return new DescriptiveStatisticsImpl(this);
    }

    /** {@inheritDoc} */
    @Override
    public void addValue(double v) {
        if (windowSize != INFINITE_WINDOW) {
            if (getN() == windowSize) {
                eDA.addElementRolling(v);
            } else if (getN() < windowSize) {
                eDA.addElement(v);
            }
        } else {
            eDA.addElement(v);
        }
    }

    /** {@inheritDoc} */
    @Override
    public void removeMostRecentValue() throws MathIllegalStateException {
        try {
            eDA.discardMostRecentElements(1);
        } catch (MathIllegalArgumentException ex) {
            throw new MathIllegalStateException(LocalizedCoreFormats.NO_DATA);
        }
    }

    /** {@inheritDoc} */
    @Override
    public double replaceMostRecentValue(double v) throws MathIllegalStateException {
        return eDA.substituteMostRecentElement(v);
    }

    /** {@inheritDoc} */
    @Override
    public double apply(UnivariateStatistic stat) {
        // No try-catch or advertised exception here because arguments are guaranteed valid
        return eDA.compute(stat);
    }

    /** {@inheritDoc} */
    @Override
    public double getMean() {
        return apply(meanImpl);
    }

    /** {@inheritDoc} */
    @Override
    public double getGeometricMean() {
        return apply(geometricMeanImpl);
    }

    /** {@inheritDoc} */
    @Override
    public double getVariance() {
        return apply(varianceImpl);
    }

    /** {@inheritDoc} */
    @Override
    public double getPopulationVariance() {
        return apply(populationVariance);
    }

    /** {@inheritDoc} */
    @Override
    public double getSkewness() {
        return apply(skewnessImpl);
    }

    /** {@inheritDoc} */
    @Override
    public double getKurtosis() {
        return apply(kurtosisImpl);
    }

    /** {@inheritDoc} */
    @Override
    public double getMax() {
        return apply(maxImpl);
    }

    /** {@inheritDoc} */
    @Override
    public double getMin() {
        return apply(minImpl);
    }

    /** {@inheritDoc} */
    @Override
    public double getSum() {
        return apply(sumImpl);
    }

    /** {@inheritDoc} */
    @Override
    public double getSumOfSquares() {
        return apply(sumOfSquaresImpl);
    }

    /** {@inheritDoc} */
    @Override
    public double getPercentile(final double p)
        throws MathIllegalArgumentException {

        percentileImpl.setQuantile(p);
        return apply(percentileImpl);
    }

    /** {@inheritDoc} */
    @Override
    public long getN() {
        return eDA.getNumElements();
    }

    /** {@inheritDoc} */
    @Override
    public void clear() {
        eDA.clear();
    }

    /** {@inheritDoc} */
    @Override
    public int getWindowSize() {
        return windowSize;
    }

    /** {@inheritDoc} */
    @Override
    public void setWindowSize(int windowSize)
        throws MathIllegalArgumentException {

        if (windowSize < 1 && windowSize != INFINITE_WINDOW) {
            throw new MathIllegalArgumentException(
                    LocalizedCoreFormats.NOT_POSITIVE_WINDOW_SIZE, windowSize);
        }

        this.windowSize = windowSize;

        // We need to check to see if we need to discard elements
        // from the front of the array.  If the windowSize is less than
        // the current number of elements.
        if (windowSize != INFINITE_WINDOW && windowSize < eDA.getNumElements()) {
            eDA.discardFrontElements(eDA.getNumElements() - windowSize);
        }
    }

    /** {@inheritDoc} */
    @Override
    public double[] getValues() {
        return eDA.getElements();
    }

    /** {@inheritDoc} */
    @Override
    public double getElement(int index) {
        return eDA.getElement(index);
    }

    /**
     * Generates a text report displaying univariate statistics from values
     * that have been added.  Each statistic is displayed on a separate
     * line.
     *
     * @return String with line feeds displaying statistics
     */
    @Override
    public String toString() {
        final StringBuilder outBuffer = new StringBuilder();
        final String endl = "\n";
        outBuffer.append("DescriptiveStatistics:").append(endl);
        outBuffer.append("n: ").append(getN()).append(endl);
        outBuffer.append("min: ").append(getMin()).append(endl);
        outBuffer.append("max: ").append(getMax()).append(endl);
        outBuffer.append("mean: ").append(getMean()).append(endl);
        outBuffer.append("std dev: ").append(getStandardDeviation()).append(endl);
        try {
            // No catch for MIAE because actual parameter is valid below
            outBuffer.append("median: ").append(getPercentile(50)).append(endl);
        } catch (MathIllegalStateException ex) {
            outBuffer.append("median: unavailable").append(endl);
        }
        outBuffer.append("skewness: ").append(getSkewness()).append(endl);
        outBuffer.append("kurtosis: ").append(getKurtosis()).append(endl);
        return outBuffer.toString();
    }

}
