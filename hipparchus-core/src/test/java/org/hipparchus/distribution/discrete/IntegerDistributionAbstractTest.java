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
package org.hipparchus.distribution.discrete;

import org.hipparchus.distribution.IntegerDistribution;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Abstract base class for {@link IntegerDistribution} tests.
 * <p>
 * To create a concrete test class for an integer distribution implementation,
 * implement makeDistribution() to return a distribution instance to use in
 * tests and each of the test data generation methods below.  In each case, the
 * test points and test values arrays returned represent parallel arrays of
 * inputs and expected values for the distribution returned by makeDistribution().
 * <p>
 * makeDensityTestPoints() -- arguments used to test probability density calculation
 * makeDensityTestValues() -- expected probability densities
 * makeCumulativeTestPoints() -- arguments used to test cumulative probabilities
 * makeCumulativeTestValues() -- expected cumulative probabilities
 * makeInverseCumulativeTestPoints() -- arguments used to test inverse cdf evaluation
 * makeInverseCumulativeTestValues() -- expected inverse cdf values
 * <p>
 * To implement additional test cases with different distribution instances and test data,
 * use the setXxx methods for the instance data in test cases and call the verifyXxx methods
 * to verify results.
 */
public abstract class IntegerDistributionAbstractTest {

//-------------------- Private test instance data -------------------------
    /** Discrete distribution instance used to perform tests */
    private IntegerDistribution distribution;

    /** Tolerance used in comparing expected and returned values */
    private double tolerance = 1E-12;

    /** Arguments used to test probability density calculations */
    private int[] densityTestPoints;

    /** Values used to test probability density calculations */
    private double[] densityTestValues;

    /** Values used to test logarithmic probability density calculations */
    private double[] logDensityTestValues;

    /** Arguments used to test cumulative probability density calculations */
    private int[] cumulativeTestPoints;

    /** Values used to test cumulative probability density calculations */
    private double[] cumulativeTestValues;

    /** Arguments used to test inverse cumulative probability density calculations */
    private double[] inverseCumulativeTestPoints;

    /** Values used to test inverse cumulative probability density calculations */
    private int[] inverseCumulativeTestValues;

    //-------------------- Abstract methods -----------------------------------

    /** Creates the default discrete distribution instance to use in tests. */
    public abstract IntegerDistribution makeDistribution();

    /** Creates the default probability density test input values */
    public abstract int[] makeDensityTestPoints();

    /** Creates the default probability density test expected values */
    public abstract double[] makeDensityTestValues();

    /** Creates the default logarithmic probability density test expected values.
     *
     * The default implementation simply computes the logarithm of all the values in
     * {@link #makeDensityTestValues()}.
     *
     * @return double[] the default logarithmic probability density test expected values.
     */
    public double[] makeLogDensityTestValues() {
        final double[] densityTestValues = makeDensityTestValues();
        final double[] logDensityTestValues = new double[densityTestValues.length];
        for (int i = 0; i < densityTestValues.length; i++) {
            logDensityTestValues[i] = FastMath.log(densityTestValues[i]);
        }
        return logDensityTestValues;
    }

    /** Creates the default cumulative probability density test input values */
    public abstract int[] makeCumulativeTestPoints();

    /** Creates the default cumulative probability density test expected values */
    public abstract double[] makeCumulativeTestValues();

    /** Creates the default inverse cumulative probability test input values */
    public abstract double[] makeInverseCumulativeTestPoints();

    /** Creates the default inverse cumulative probability density test expected values */
    public abstract int[] makeInverseCumulativeTestValues();

    //-------------------- Setup / tear down ----------------------------------

    /**
     * Setup sets all test instance data to default values
     */
    @BeforeEach
    public void setUp() {
        distribution = makeDistribution();
        densityTestPoints = makeDensityTestPoints();
        densityTestValues = makeDensityTestValues();
        logDensityTestValues = makeLogDensityTestValues();
        cumulativeTestPoints = makeCumulativeTestPoints();
        cumulativeTestValues = makeCumulativeTestValues();
        inverseCumulativeTestPoints = makeInverseCumulativeTestPoints();
        inverseCumulativeTestValues = makeInverseCumulativeTestValues();
    }

    /**
     * Cleans up test instance data
     */
    @AfterEach
    public void tearDown() {
        distribution = null;
        densityTestPoints = null;
        densityTestValues = null;
        logDensityTestValues = null;
        cumulativeTestPoints = null;
        cumulativeTestValues = null;
        inverseCumulativeTestPoints = null;
        inverseCumulativeTestValues = null;
    }

    //-------------------- Verification methods -------------------------------

    /**
     * Verifies that probability density calculations match expected values
     * using current test instance data
     */
    protected void verifyDensities() {
        for (int i = 0; i < densityTestPoints.length; i++) {
            assertEquals(densityTestValues[i],
                    distribution.probability(densityTestPoints[i]), getTolerance(), "Incorrect density value returned for " + densityTestPoints[i]);
        }
    }

    /**
     * Verifies that logarithmic probability density calculations match expected values
     * using current test instance data.
     */
    protected void verifyLogDensities() {
        for (int i = 0; i < densityTestPoints.length; i++) {
            // FIXME: when logProbability methods are added to IntegerDistribution in 4.0, remove cast below
            assertEquals(logDensityTestValues[i],
                    ((AbstractIntegerDistribution) distribution).logProbability(densityTestPoints[i]), tolerance, "Incorrect log density value returned for " + densityTestPoints[i]);
        }
    }

    /**
     * Verifies that cumulative probability density calculations match expected values
     * using current test instance data
     */
    protected void verifyCumulativeProbabilities() {
        for (int i = 0; i < cumulativeTestPoints.length; i++) {
            assertEquals(cumulativeTestValues[i],
                    distribution.cumulativeProbability(cumulativeTestPoints[i]), getTolerance(), "Incorrect cumulative probability value returned for " + cumulativeTestPoints[i]);
        }
    }


    /**
     * Verifies that inverse cumulative probability density calculations match expected values
     * using current test instance data
     */
    protected void verifyInverseCumulativeProbabilities() {
        for (int i = 0; i < inverseCumulativeTestPoints.length; i++) {
            assertEquals(inverseCumulativeTestValues[i],
                    distribution.inverseCumulativeProbability(inverseCumulativeTestPoints[i]),
                    "Incorrect inverse cumulative probability value returned for "
                    + inverseCumulativeTestPoints[i]);
        }
    }

    //------------------------ Default test cases -----------------------------

    /**
     * Verifies that probability density calculations match expected values
     * using default test instance data
     */
    @Test
    public void testDensities() {
        verifyDensities();
    }

    /**
     * Verifies that logarithmic probability density calculations match expected values
     * using default test instance data
     */
    @Test
    public void testLogDensities() {
        verifyLogDensities();
    }

    /**
     * Verifies that cumulative probability density calculations match expected values
     * using default test instance data
     */
    @Test
    public void testCumulativeProbabilities() {
        verifyCumulativeProbabilities();
    }

    /**
     * Verifies that inverse cumulative probability density calculations match expected values
     * using default test instance data
     */
    @Test
    public void testInverseCumulativeProbabilities() {
        verifyInverseCumulativeProbabilities();
    }

    @Test
    public void testConsistencyAtSupportBounds() {
        final int lower = distribution.getSupportLowerBound();
        assertEquals(0.0, distribution.cumulativeProbability(lower - 1), 0.0, "Cumulative probability mmust be 0 below support lower bound.");
        assertEquals(distribution.probability(lower), distribution.cumulativeProbability(lower), getTolerance(), "Cumulative probability of support lower bound must be equal to probability mass at this point.");
        assertEquals(lower, distribution.inverseCumulativeProbability(0.0), "Inverse cumulative probability of 0 must be equal to support lower bound.");

        final int upper = distribution.getSupportUpperBound();
        if (upper != Integer.MAX_VALUE)
            assertEquals(1.0, distribution.cumulativeProbability(upper), 0.0, "Cumulative probability of support upper bound must be equal to 1.");
        assertEquals(upper, distribution.inverseCumulativeProbability(1.0), "Inverse cumulative probability of 1 must be equal to support upper bound.");
    }

    /**
     * Verifies that illegal arguments are correctly handled
     */
    @Test
    public void testIllegalArguments() {
        try {
            distribution.probability(1, 0);
            fail("Expecting MathIllegalArgumentException for bad cumulativeProbability interval");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            distribution.inverseCumulativeProbability(-1);
            fail("Expecting MathIllegalArgumentException for p = -1");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            distribution.inverseCumulativeProbability(2);
            fail("Expecting MathIllegalArgumentException for p = 2");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    //------------------ Getters / Setters for test instance data -----------
    /**
     * @return Returns the cumulativeTestPoints.
     */
    protected int[] getCumulativeTestPoints() {
        return cumulativeTestPoints;
    }

    /**
     * @param cumulativeTestPoints The cumulativeTestPoints to set.
     */
    protected void setCumulativeTestPoints(int[] cumulativeTestPoints) {
        this.cumulativeTestPoints = cumulativeTestPoints;
    }

    /**
     * @return Returns the cumulativeTestValues.
     */
    protected double[] getCumulativeTestValues() {
        return cumulativeTestValues;
    }

    /**
     * @param cumulativeTestValues The cumulativeTestValues to set.
     */
    protected void setCumulativeTestValues(double[] cumulativeTestValues) {
        this.cumulativeTestValues = cumulativeTestValues;
    }

    /**
     * @return Returns the densityTestPoints.
     */
    protected int[] getDensityTestPoints() {
        return densityTestPoints;
    }

    /**
     * @param densityTestPoints The densityTestPoints to set.
     */
    protected void setDensityTestPoints(int[] densityTestPoints) {
        this.densityTestPoints = densityTestPoints;
    }

    /**
     * @return Returns the densityTestValues.
     */
    protected double[] getDensityTestValues() {
        return densityTestValues;
    }

    /**
     * @param densityTestValues The densityTestValues to set.
     */
    protected void setDensityTestValues(double[] densityTestValues) {
        this.densityTestValues = densityTestValues;
    }

    /**
     * @return Returns the distribution.
     */
    protected IntegerDistribution getDistribution() {
        return distribution;
    }

    /**
     * @param distribution The distribution to set.
     */
    protected void setDistribution(IntegerDistribution distribution) {
        this.distribution = distribution;
    }

    /**
     * @return Returns the inverseCumulativeTestPoints.
     */
    protected double[] getInverseCumulativeTestPoints() {
        return inverseCumulativeTestPoints;
    }

    /**
     * @param inverseCumulativeTestPoints The inverseCumulativeTestPoints to set.
     */
    protected void setInverseCumulativeTestPoints(double[] inverseCumulativeTestPoints) {
        this.inverseCumulativeTestPoints = inverseCumulativeTestPoints;
    }

    /**
     * @return Returns the inverseCumulativeTestValues.
     */
    protected int[] getInverseCumulativeTestValues() {
        return inverseCumulativeTestValues;
    }

    /**
     * @param inverseCumulativeTestValues The inverseCumulativeTestValues to set.
     */
    protected void setInverseCumulativeTestValues(int[] inverseCumulativeTestValues) {
        this.inverseCumulativeTestValues = inverseCumulativeTestValues;
    }

    /**
     * @return Returns the tolerance.
     */
    protected double getTolerance() {
        return tolerance;
    }

    /**
     * @param tolerance The tolerance to set.
     */
    protected void setTolerance(double tolerance) {
        this.tolerance = tolerance;
    }

}
