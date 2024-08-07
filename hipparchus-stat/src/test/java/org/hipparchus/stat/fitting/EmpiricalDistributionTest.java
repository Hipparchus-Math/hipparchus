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

import org.hipparchus.UnitTestUtils;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.integration.BaseAbstractUnivariateIntegrator;
import org.hipparchus.analysis.integration.IterativeLegendreGaussIntegrator;
import org.hipparchus.distribution.RealDistribution;
import org.hipparchus.distribution.continuous.ConstantRealDistribution;
import org.hipparchus.distribution.continuous.NormalDistribution;
import org.hipparchus.distribution.continuous.RealDistributionAbstractTest;
import org.hipparchus.distribution.continuous.UniformRealDistribution;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.stat.descriptive.StreamingStatistics;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the EmpiricalDistribution class
 */
public final class EmpiricalDistributionTest extends RealDistributionAbstractTest {

    protected EmpiricalDistribution empiricalDistribution = null;
    protected EmpiricalDistribution empiricalDistribution2 = null;
    protected File file = null;
    protected URL url = null;
    protected double[] dataArray = null;
    protected final int n = 10000;

    @Override
    @BeforeEach
    public void setUp() {
        super.setUp();
        empiricalDistribution = new EmpiricalDistribution(100);
        url = getClass().getResource("testData.txt");
        final ArrayList<Double> list = new ArrayList<>();
        assertDoesNotThrow(() -> {
            empiricalDistribution2 = new EmpiricalDistribution(100);
            BufferedReader in =
                new BufferedReader(new InputStreamReader(
                    url.openStream()));
            String str = null;
            while ((str = in.readLine()) != null) {
                list.add(Double.valueOf(str));
            }
            in.close();
            in = null;
        }, "IOException ");

        dataArray = new double[list.size()];
        int i = 0;
        for (Double data : list) {
            dataArray[i] = data.doubleValue();
            i++;
        }
    }

    // MATH-1279
    @Test
    void testPrecondition1() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new EmpiricalDistribution(0);
        });
    }

    /**
     * Test EmpiricalDistrbution.load() using sample data file.<br>
     * Check that the sampleCount, mu and sigma match data in
     * the sample data file. Also verify that load is idempotent.
     */
    @Test
    void testLoad() throws Exception {
        // Load from a URL
        empiricalDistribution.load(url);
        checkDistribution();

        // Load again from a file (also verifies idempotency of load)
        File file = new File(url.toURI());
        empiricalDistribution.load(file);
        checkDistribution();
    }

    private void checkDistribution() {
        // testData File has 10000 values, with mean ~ 5.0, std dev ~ 1
        // Make sure that loaded distribution matches this
        assertEquals(1000, empiricalDistribution.getSampleStats().getN(), 10E-7);
        //TODO: replace with statistical tests
        assertEquals(
            5.069831575018909, empiricalDistribution.getSampleStats().getMean(), 10E-7);
        assertEquals(
            1.0173699343977738, empiricalDistribution.getSampleStats().getStandardDeviation(), 10E-7);
    }

    @Test
    void testLoadURLError() throws IOException {
        try {
            URL existing = getClass().getResource("testData.txt");
            URL nonexistent = new URL(existing.toString() + "-nonexistent");
            empiricalDistribution.load(nonexistent);
            fail("an exception should have been thrown");
        } catch (IOException ioe) {
            // expected
        }
    }

    @Test
    void testLoadFileError() throws IOException, URISyntaxException {
        try {
            File existing = new File(getClass().getResource("testData.txt").toURI());
            File nonexistent = new File(existing.getAbsolutePath() + "-nonexistent");
            empiricalDistribution.load(nonexistent);
            fail("an exception should have been thrown");
        } catch (IOException ioe) {
            // expected
        }
    }

    /**
     * Test EmpiricalDistrbution.load(double[]) using data taken from
     * sample data file.<br>
     * Check that the sampleCount, mu and sigma match data in
     * the sample data file.
     */
    @Test
    void testDoubleLoad() throws Exception {
        empiricalDistribution2.load(dataArray);
        // testData File has 10000 values, with mean ~ 5.0, std dev ~ 1
        // Make sure that loaded distribution matches this
        assertEquals(1000, empiricalDistribution2.getSampleStats().getN(), 10E-7);
        //TODO: replace with statistical tests
        assertEquals(
            5.069831575018909, empiricalDistribution2.getSampleStats().getMean(), 10E-7);
        assertEquals(
            1.0173699343977738, empiricalDistribution2.getSampleStats().getStandardDeviation(), 10E-7);

        double[] bounds = empiricalDistribution2.getGeneratorUpperBounds();
        assertEquals(100, bounds.length);
        assertEquals(1.0, bounds[99], 10e-12);

    }

    /**
      * Generate 1000 random values and make sure they look OK.<br>
      * Note that there is a non-zero (but very small) probability that
      * these tests will fail even if the code is working as designed.
      */
    @Test
    void testNext() throws Exception {
        tstGen(0.1);
        tstDoubleGen(0.1);
    }

    /**
      * Make sure exception thrown if digest getNext is attempted
      * before loading empiricalDistribution.
     */
    @Test
    void testNexFail() {
        try {
            empiricalDistribution.getNextValue();
            empiricalDistribution2.getNextValue();
            fail("Expecting MathIllegalStateException");
        } catch (MathIllegalStateException ex) {
            // expected
        }
    }

    /**
     * Make sure we can handle a grid size that is too fine
     */
    @Test
    void testGridTooFine() throws Exception {
        empiricalDistribution = new EmpiricalDistribution(1001);
        tstGen(0.1);
        empiricalDistribution2 = new EmpiricalDistribution(1001);
        tstDoubleGen(0.1);
    }

    /**
     * How about too fat?
     */
    @Test
    void testGridTooFat() throws Exception {
        empiricalDistribution = new EmpiricalDistribution(1);
        tstGen(5); // ridiculous tolerance; but ridiculous grid size
                   // really just checking to make sure we do not bomb
        empiricalDistribution2 = new EmpiricalDistribution(1);
        tstDoubleGen(5);
    }

    /**
     * Test bin index overflow problem (BZ 36450)
     */
    @Test
    void testBinIndexOverflow() throws Exception {
        double[] x = new double[] {9474.94326071674, 2080107.8865462579};
        new EmpiricalDistribution().load(x);
    }

    @Test
    void testSerialization() {
        // Empty
        EmpiricalDistribution dist = new EmpiricalDistribution();
        EmpiricalDistribution dist2 = (EmpiricalDistribution) UnitTestUtils.serializeAndRecover(dist);
        verifySame(dist, dist2);

        // Loaded
        empiricalDistribution2.load(dataArray);
        dist2 = (EmpiricalDistribution) UnitTestUtils.serializeAndRecover(empiricalDistribution2);
        verifySame(empiricalDistribution2, dist2);
    }

    @Test
    void testLoadNullDoubleArray() {
        assertThrows(NullArgumentException.class, () -> {
            new EmpiricalDistribution().load((double[]) null);
        });
    }

    @Test
    void testLoadNullURL() throws Exception {
        assertThrows(NullArgumentException.class, () -> {
            new EmpiricalDistribution().load((URL) null);
        });
    }

    @Test
    void testLoadNullFile() throws Exception {
        assertThrows(NullArgumentException.class, () -> {
            new EmpiricalDistribution().load((File) null);
        });
    }

    /**
     * MATH-298
     */
    @Test
    void testGetBinUpperBounds() {
        double[] testData = {0, 1, 1, 2, 3, 4, 4, 5, 6, 7, 8, 9, 10};
        EmpiricalDistribution dist = new EmpiricalDistribution(5);
        dist.load(testData);
        double[] expectedBinUpperBounds = {2, 4, 6, 8, 10};
        double[] expectedGeneratorUpperBounds = {4d/13d, 7d/13d, 9d/13d, 11d/13d, 1};
        double tol = 10E-12;
        UnitTestUtils.customAssertEquals(expectedBinUpperBounds, dist.getUpperBounds(), tol);
        UnitTestUtils.customAssertEquals(expectedGeneratorUpperBounds, dist.getGeneratorUpperBounds(), tol);
    }

    @Test
    void testReSeed() throws Exception {
        empiricalDistribution.load(url);
        empiricalDistribution.reSeed(100);
        final double [] values = new double[10];
        for (int i = 0; i < 10; i++) {
            values[i] = empiricalDistribution.getNextValue();
        }
        empiricalDistribution.reSeed(100);
        for (int i = 0; i < 10; i++) {
            assertEquals(values[i],empiricalDistribution.getNextValue(), 0d);
        }
    }

    private void verifySame(EmpiricalDistribution d1, EmpiricalDistribution d2) {
        assertEquals(d1.isLoaded(), d2.isLoaded());
        assertEquals(d1.getBinCount(), d2.getBinCount());
        assertEquals(d1.getSampleStats(), d2.getSampleStats());
        if (d1.isLoaded()) {
            for (int i = 0;  i < d1.getUpperBounds().length; i++) {
                assertEquals(d1.getUpperBounds()[i], d2.getUpperBounds()[i], 0);
            }
            assertEquals(d1.getBinStats(), d2.getBinStats());
        }
    }

    private void tstGen(double tolerance)throws Exception {
        empiricalDistribution.load(url);
        empiricalDistribution.reSeed(1000);
        StreamingStatistics stats = new StreamingStatistics();
        for (int i = 1; i < 1000; i++) {
            stats.addValue(empiricalDistribution.getNextValue());
        }
        assertEquals(5.069831575018909, stats.getMean(),tolerance,"mean");
        assertEquals(1.0173699343977738, stats.getStandardDeviation(),tolerance,"std dev");
    }

    private void tstDoubleGen(double tolerance)throws Exception {
        empiricalDistribution2.load(dataArray);
        empiricalDistribution2.reSeed(1000);
        StreamingStatistics stats = new StreamingStatistics();
        for (int i = 1; i < 1000; i++) {
            stats.addValue(empiricalDistribution2.getNextValue());
        }
        assertEquals(5.069831575018909, stats.getMean(), tolerance, "mean");
        assertEquals(1.0173699343977738, stats.getStandardDeviation(), tolerance, "std dev");
    }

    //  Setup for distribution tests

    @Override
    public RealDistribution makeDistribution() {
        // Create a uniform distribution on [0, 10,000]
        final double[] sourceData = new double[n + 1];
        for (int i = 0; i < n + 1; i++) {
            sourceData[i] = i;
        }
        EmpiricalDistribution dist = new EmpiricalDistribution();
        dist.load(sourceData);
        return dist;
    }

    /** Uniform bin mass = 10/10001 == mass of all but the first bin */
    private final double binMass = 10d / (n + 1);

    /** Mass of first bin = 11/10001 */
    private final double firstBinMass = 11d / (n + 1);

    @Override
    public double[] makeCumulativeTestPoints() {
       final double[] testPoints = new double[] {9, 10, 15, 1000, 5004, 9999};
       return testPoints;
    }


    @Override
    public double[] makeCumulativeTestValues() {
        /*
         * Bins should be [0, 10], (10, 20], ..., (9990, 10000]
         * Kernels should be N(4.5, 3.02765), N(14.5, 3.02765)...
         * Each bin should have mass 10/10000 = .001
         */
        final double[] testPoints = getCumulativeTestPoints();
        final double[] cumValues = new double[testPoints.length];
        final EmpiricalDistribution empiricalDistribution = (EmpiricalDistribution) makeDistribution();
        final double[] binBounds = empiricalDistribution.getUpperBounds();
        for (int i = 0; i < testPoints.length; i++) {
            final int bin = findBin(testPoints[i]);
            final double lower = bin == 0 ? empiricalDistribution.getSupportLowerBound() :
                binBounds[bin - 1];
            final double upper = binBounds[bin];
            // Compute bMinus = sum or mass of bins below the bin containing the point
            // First bin has mass 11 / 10000, the rest have mass 10 / 10000.
            final double bMinus = bin == 0 ? 0 : (bin - 1) * binMass + firstBinMass;
            final RealDistribution kernel = findKernel(lower, upper);
            final double withinBinKernelMass = kernel.probability(lower, upper);
            final double kernelCum = kernel.probability(lower, testPoints[i]);
            cumValues[i] = bMinus + (bin == 0 ? firstBinMass : binMass) * kernelCum/withinBinKernelMass;
        }
        return cumValues;
    }

    @Override
    public double[] makeDensityTestValues() {
        final double[] testPoints = getCumulativeTestPoints();
        final double[] densityValues = new double[testPoints.length];
        final EmpiricalDistribution empiricalDistribution = (EmpiricalDistribution) makeDistribution();
        final double[] binBounds = empiricalDistribution.getUpperBounds();
        for (int i = 0; i < testPoints.length; i++) {
            final int bin = findBin(testPoints[i]);
            final double lower = bin == 0 ? empiricalDistribution.getSupportLowerBound() :
                binBounds[bin - 1];
            final double upper = binBounds[bin];
            final RealDistribution kernel = findKernel(lower, upper);
            final double withinBinKernelMass = kernel.probability(lower, upper);
            final double density = kernel.density(testPoints[i]);
            densityValues[i] = density * (bin == 0 ? firstBinMass : binMass) / withinBinKernelMass;
        }
        return densityValues;
    }

    /**
     * Modify test integration bounds from the default. Because the distribution
     * has discontinuities at bin boundaries, integrals spanning multiple bins
     * will face convergence problems.  Only test within-bin integrals and spans
     * across no more than 3 bin boundaries.
     */
    @Override
    @Test
    public void testDensityIntegrals() {
        final RealDistribution distribution = makeDistribution();
        final double tol = 1.0e-9;
        final BaseAbstractUnivariateIntegrator integrator =
            new IterativeLegendreGaussIntegrator(5, 1.0e-12, 1.0e-10);
        final UnivariateFunction d = new UnivariateFunction() {
            @Override
            public double value(double x) {
                return distribution.density(x);
            }
        };
        final double[] lower = {0, 5, 1000, 5001, 9995};
        final double[] upper = {5, 12, 1030, 5010, 10000};
        for (int i = 1; i < 5; i++) {
            assertEquals(
                    distribution.probability(
                            lower[i], upper[i]),
                            integrator.integrate(
                                    1000000, // Triangle integrals are very slow to converge
                                    d, lower[i], upper[i]), tol);
        }
    }

    /**
     * MATH-984
     * Verify that sampled values do not go outside of the range of the data.
     */
    @Test
    void testSampleValuesRange() {
        // Concentrate values near the endpoints of (0, 1).
        // Unconstrained Gaussian kernel would generate values outside the interval.
        final double[] data = new double[100];
        for (int i = 0; i < 50; i++) {
            data[i] = 1 / ((double) i + 1);
        }
        for (int i = 51; i < 100; i++) {
            data[i] = 1 - 1 / (100 - (double) i + 2);
        }
        EmpiricalDistribution dist = new EmpiricalDistribution(10);
        dist.load(data);
        dist.reseedRandomGenerator(1000);
        for (int i = 0; i < 1000; i++) {
            final double dev = dist.getNextValue();
            assertTrue(dev < 1);
            assertTrue(dev > 0);
        }
    }

    /**
     * MATH-1203, MATH-1208
     */
    @Test
    void testNoBinVariance() {
        final double[] data = {0, 0, 1, 1};
        EmpiricalDistribution dist = new EmpiricalDistribution(2);
        dist.load(data);
        dist.reseedRandomGenerator(1000);
        for (int i = 0; i < 1000; i++) {
            final double dev = dist.getNextValue();
            assertTrue(dev == 0 || dev == 1);
        }
        assertEquals(0.5, dist.cumulativeProbability(0), Double.MIN_VALUE);
        assertEquals(1.0, dist.cumulativeProbability(1), Double.MIN_VALUE);
        assertEquals(0.5, dist.cumulativeProbability(0.5), Double.MIN_VALUE);
        assertEquals(0.5, dist.cumulativeProbability(0.7), Double.MIN_VALUE);
    }

    /**
     * Find the bin that x belongs (relative to {@link #makeDistribution()}).
     */
    private int findBin(double x) {
        // Number of bins below x should be trunc(x/10)
        final double nMinus = FastMath.floor(x / 10);
        final int bin =  (int) FastMath.round(nMinus);
        // If x falls on a bin boundary, it is in the lower bin
        return FastMath.floor(x / 10) == x / 10 ? bin - 1 : bin;
    }

    /**
     * Find the within-bin kernel for the bin with lower bound lower
     * and upper bound upper. All bins other than the first contain 10 points
     * exclusive of the lower bound and are centered at (lower + upper + 1) / 2.
     * The first bin includes its lower bound, 0, so has different mean and
     * standard deviation.
     */
    private RealDistribution findKernel(double lower, double upper) {
        if (lower < 1) {
            return new NormalDistribution(5d, 3.3166247903554);
        } else {
            return new NormalDistribution((upper + lower + 1) / 2d, 3.0276503540974917);
        }
    }

    @Test
    void testKernelOverrideConstant() {
        final EmpiricalDistribution dist = new ConstantKernelEmpiricalDistribution(5);
        final double[] data = {1d,2d,3d, 4d,5d,6d, 7d,8d,9d, 10d,11d,12d, 13d,14d,15d};
        dist.load(data);
        // Bin masses concentrated on 2, 5, 8, 11, 14 <- effectively discrete uniform distribution over these
        double[] values = {2d, 5d, 8d, 11d, 14d};
        for (int i = 0; i < 20; i++) {
            assertTrue(Arrays.binarySearch(values, dist.getNextValue()) >= 0);
        }
        final double tol = 10E-12;
        assertEquals(0.0, dist.cumulativeProbability(1), tol);
        assertEquals(0.2, dist.cumulativeProbability(2), tol);
        assertEquals(0.6, dist.cumulativeProbability(10), tol);
        assertEquals(0.8, dist.cumulativeProbability(12), tol);
        assertEquals(0.8, dist.cumulativeProbability(13), tol);
        assertEquals(1.0, dist.cumulativeProbability(15), tol);

        assertEquals(2.0, dist.inverseCumulativeProbability(0.1), tol);
        assertEquals(2.0, dist.inverseCumulativeProbability(0.2), tol);
        assertEquals(5.0, dist.inverseCumulativeProbability(0.3), tol);
        assertEquals(5.0, dist.inverseCumulativeProbability(0.4), tol);
        assertEquals(8.0, dist.inverseCumulativeProbability(0.5), tol);
        assertEquals(8.0, dist.inverseCumulativeProbability(0.6), tol);
    }

    @Test
    void testKernelOverrideUniform() {
        final EmpiricalDistribution dist = new UniformKernelEmpiricalDistribution(5);
        final double[] data = {1d,2d,3d, 4d,5d,6d, 7d,8d,9d, 10d,11d,12d, 13d,14d,15d};
        dist.load(data);
        // Kernels are uniform distributions on [1,3], [4,6], [7,9], [10,12], [13,15]
        final double[] bounds = {3d, 6d, 9d, 12d};
        final double tol = 10E-12;
        for (int i = 0; i < 20; i++) {
            final double v = dist.getNextValue();
            // Make sure v is not in the excluded range between bins - that is (bounds[i], bounds[i] + 1)
            for (int j = 0; j < bounds.length; j++) {
                assertFalse(v > bounds[j] + tol && v < bounds[j] + 1 - tol);
            }
        }
        assertEquals(0.0, dist.cumulativeProbability(1), tol);
        assertEquals(0.1, dist.cumulativeProbability(2), tol);
        assertEquals(0.6, dist.cumulativeProbability(10), tol);
        assertEquals(0.8, dist.cumulativeProbability(12), tol);
        assertEquals(0.8, dist.cumulativeProbability(13), tol);
        assertEquals(1.0, dist.cumulativeProbability(15), tol);

        assertEquals(2.0, dist.inverseCumulativeProbability(0.1), tol);
        assertEquals(3.0, dist.inverseCumulativeProbability(0.2), tol);
        assertEquals(5.0, dist.inverseCumulativeProbability(0.3), tol);
        assertEquals(6.0, dist.inverseCumulativeProbability(0.4), tol);
        assertEquals(8.0, dist.inverseCumulativeProbability(0.5), tol);
        assertEquals(9.0, dist.inverseCumulativeProbability(0.6), tol);
    }

    @Test
    void testEmptyBins() {
        double[] data = new double[10];
        for (int i = 0; i < 10; ++i) {
            data[i] = i < 5 ? 0 : 1;
        }
        EmpiricalDistribution edist = new EmpiricalDistribution(100);
        edist.load(data);
        assertEquals(0.5, edist.cumulativeProbability(0), Double.MIN_VALUE);
        assertEquals(0.5, edist.cumulativeProbability(0.3), Double.MIN_VALUE);
        assertEquals(0.5, edist.cumulativeProbability(0.9), Double.MIN_VALUE);
        assertEquals(1, edist.cumulativeProbability(1), Double.MIN_VALUE);
        assertEquals(1, edist.cumulativeProbability(1.5), Double.MIN_VALUE);
    }


    /**
     * Empirical distribution using a constant smoothing kernel.
     */
    private class ConstantKernelEmpiricalDistribution extends EmpiricalDistribution {
        private static final long serialVersionUID = 1L;
        public ConstantKernelEmpiricalDistribution(int i) {
            super(i);
        }
        // Use constant distribution equal to bin mean within bin
        @Override
        protected RealDistribution getKernel(StreamingStatistics bStats) {
            return new ConstantRealDistribution(bStats.getMean());
        }
    }

    /**
     * Empirical distribution using a uniform smoothing kernel.
     */
    private class UniformKernelEmpiricalDistribution extends EmpiricalDistribution {
        private static final long serialVersionUID = 2963149194515159653L;
        public UniformKernelEmpiricalDistribution(int i) {
            super(i);
        }
        @Override
        protected RealDistribution getKernel(StreamingStatistics bStats) {
            return new UniformRealDistribution(bStats.getMin(), bStats.getMax());
        }
    }
}
