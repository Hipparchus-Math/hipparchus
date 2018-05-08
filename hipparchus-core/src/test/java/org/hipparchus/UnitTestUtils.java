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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexFormat;
import org.hipparchus.distribution.RealDistribution;
import org.hipparchus.distribution.continuous.ChiSquaredDistribution;
import org.hipparchus.linear.BlockRealMatrix;
import org.hipparchus.linear.FieldMatrix;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.junit.Assert;

/**
 */
public class UnitTestUtils {
    /**
     * Collection of static methods used in math unit tests.
     */
    private UnitTestUtils() {
        super();
    }

    /**
     * Verifies that expected and actual are within delta, or are both NaN or
     * infinities of the same sign.
     */
    public static void assertEquals(double expected, double actual, double delta) {
        Assert.assertEquals(null, expected, actual, delta);
    }

    /**
     * Verifies that expected and actual are within delta, or are both NaN or
     * infinities of the same sign.
     */
    public static void assertEquals(String msg, double expected, double actual, double delta) {
        // check for NaN
        if(Double.isNaN(expected)){
            Assert.assertTrue("" + actual + " is not NaN.",
                Double.isNaN(actual));
        } else {
            Assert.assertEquals(msg, expected, actual, delta);
        }
    }

    /**
     * Verifies that the two arguments are exactly the same, either
     * both NaN or infinities of same sign, or identical floating point values.
     */
    public static void assertSame(double expected, double actual) {
     Assert.assertEquals(expected, actual, 0);
    }

    /**
     * Verifies that real and imaginary parts of the two complex arguments
     * are exactly the same.  Also ensures that NaN / infinite components match.
     */
    public static void assertSame(Complex expected, Complex actual) {
        assertSame(expected.getReal(), actual.getReal());
        assertSame(expected.getImaginary(), actual.getImaginary());
    }

    /**
     * Verifies that real and imaginary parts of the two complex arguments
     * differ by at most delta.  Also ensures that NaN / infinite components match.
     */
    public static void assertEquals(Complex expected, Complex actual, double delta) {
        Assert.assertEquals(expected.getReal(), actual.getReal(), delta);
        Assert.assertEquals(expected.getImaginary(), actual.getImaginary(), delta);
    }

    /**
     * Verifies that two double arrays have equal entries, up to tolerance
     */
    public static void assertEquals(double expected[], double observed[], double tolerance) {
        assertEquals("Array comparison failure", expected, observed, tolerance);
    }

    /**
     * Serializes an object to a bytes array and then recovers the object from the bytes array.
     * Returns the deserialized object.
     *
     * @param o  object to serialize and recover
     * @return  the recovered, deserialized object
     */
    public static Object serializeAndRecover(Object o) {
        try {
            // serialize the Object
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            ObjectOutputStream so = new ObjectOutputStream(bos);
            so.writeObject(o);

            // deserialize the Object
            ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
            ObjectInputStream si = new ObjectInputStream(bis);
            return si.readObject();
        } catch (IOException ioe) {
            return null;
        } catch (ClassNotFoundException cnfe) {
            return null;
        }
    }

    /**
     * Verifies that serialization preserves equals and hashCode.
     * Serializes the object, then recovers it and checks equals and hash code.
     *
     * @param object  the object to serialize and recover
     */
    public static void checkSerializedEquality(Object object) {
        Object object2 = serializeAndRecover(object);
        Assert.assertEquals("Equals check", object, object2);
        Assert.assertEquals("HashCode check", object.hashCode(), object2.hashCode());
    }

    /**
     * Verifies that the relative error in actual vs. expected is less than or
     * equal to relativeError.  If expected is infinite or NaN, actual must be
     * the same (NaN or infinity of the same sign).
     *
     * @param expected expected value
     * @param actual  observed value
     * @param relativeError  maximum allowable relative error
     */
    public static void assertRelativelyEquals(double expected, double actual,
            double relativeError) {
        assertRelativelyEquals(null, expected, actual, relativeError);
    }

    /**
     * Verifies that the relative error in actual vs. expected is less than or
     * equal to relativeError.  If expected is infinite or NaN, actual must be
     * the same (NaN or infinity of the same sign).
     *
     * @param msg  message to return with failure
     * @param expected expected value
     * @param actual  observed value
     * @param relativeError  maximum allowable relative error
     */
    public static void assertRelativelyEquals(String msg, double expected,
            double actual, double relativeError) {
        if (Double.isNaN(expected)) {
            Assert.assertTrue(msg, Double.isNaN(actual));
        } else if (Double.isNaN(actual)) {
            Assert.assertTrue(msg, Double.isNaN(expected));
        } else if (Double.isInfinite(actual) || Double.isInfinite(expected)) {
            Assert.assertEquals(expected, actual, relativeError);
        } else if (expected == 0.0) {
            Assert.assertEquals(msg, actual, expected, relativeError);
        } else {
            double absError = FastMath.abs(expected) * relativeError;
            Assert.assertEquals(msg, expected, actual, absError);
        }
    }

    /**
     * Fails iff values does not contain a number within epsilon of z.
     *
     * @param msg  message to return with failure
     * @param values complex array to search
     * @param z  value sought
     * @param epsilon  tolerance
     */
    public static void assertContains(String msg, Complex[] values,
                                      Complex z, double epsilon) {
        for (Complex value : values) {
            if (Precision.equals(value.getReal(), z.getReal(), epsilon) &&
                Precision.equals(value.getImaginary(), z.getImaginary(), epsilon)) {
                return;
            }
        }
        Assert.fail(msg + " Unable to find " + (new ComplexFormat()).format(z));
    }

    /**
     * Fails iff values does not contain a number within epsilon of z.
     *
     * @param values complex array to search
     * @param z  value sought
     * @param epsilon  tolerance
     */
    public static void assertContains(Complex[] values,
            Complex z, double epsilon) {
        assertContains(null, values, z, epsilon);
    }

    /**
     * Fails iff values does not contain a number within epsilon of x.
     *
     * @param msg  message to return with failure
     * @param values double array to search
     * @param x value sought
     * @param epsilon  tolerance
     */
    public static void assertContains(String msg, double[] values,
            double x, double epsilon) {
        for (double value : values) {
            if (Precision.equals(value, x, epsilon)) {
                return;
            }
        }
        Assert.fail(msg + " Unable to find " + x);
    }

    /**
     * Fails iff values does not contain a number within epsilon of x.
     *
     * @param values double array to search
     * @param x value sought
     * @param epsilon  tolerance
     */
    public static void assertContains(double[] values, double x,
            double epsilon) {
       assertContains(null, values, x, epsilon);
    }

    /**
     * Asserts that all entries of the specified vectors are equal to within a
     * positive {@code delta}.
     *
     * @param message the identifying message for the assertion error (can be
     * {@code null})
     * @param expected expected value
     * @param actual actual value
     * @param delta the maximum difference between the entries of the expected
     * and actual vectors for which both entries are still considered equal
     */
    public static void assertEquals(final String message,
        final double[] expected, final RealVector actual, final double delta) {
        final String msgAndSep = message.equals("") ? "" : message + ", ";
        Assert.assertEquals(msgAndSep + "dimension", expected.length,
            actual.getDimension());
        for (int i = 0; i < expected.length; i++) {
            Assert.assertEquals(msgAndSep + "entry #" + i, expected[i],
                actual.getEntry(i), delta);
        }
    }

    /**
     * Asserts that all entries of the specified vectors are equal to within a
     * positive {@code delta}.
     *
     * @param message the identifying message for the assertion error (can be
     * {@code null})
     * @param expected expected value
     * @param actual actual value
     * @param delta the maximum difference between the entries of the expected
     * and actual vectors for which both entries are still considered equal
     */
    public static void assertEquals(final String message,
        final RealVector expected, final RealVector actual, final double delta) {
        final String msgAndSep = message.equals("") ? "" : message + ", ";
        Assert.assertEquals(msgAndSep + "dimension", expected.getDimension(),
            actual.getDimension());
        final int dim = expected.getDimension();
        for (int i = 0; i < dim; i++) {
            Assert.assertEquals(msgAndSep + "entry #" + i,
                expected.getEntry(i), actual.getEntry(i), delta);
        }
    }

    /** verifies that two matrices are close (1-norm) */
    public static void assertEquals(String msg, RealMatrix expected, RealMatrix observed, double tolerance) {

        Assert.assertNotNull(msg + "\nObserved should not be null",observed);

        if (expected.getColumnDimension() != observed.getColumnDimension() ||
                expected.getRowDimension() != observed.getRowDimension()) {
            StringBuilder messageBuffer = new StringBuilder(msg);
            messageBuffer.append("\nObserved has incorrect dimensions.");
            messageBuffer.append("\nobserved is " + observed.getRowDimension() +
                    " x " + observed.getColumnDimension());
            messageBuffer.append("\nexpected " + expected.getRowDimension() +
                    " x " + expected.getColumnDimension());
            Assert.fail(messageBuffer.toString());
        }

        RealMatrix delta = expected.subtract(observed);
        if (delta.getNorm() >= tolerance) {
            StringBuilder messageBuffer = new StringBuilder(msg);
            messageBuffer.append("\nExpected: " + expected);
            messageBuffer.append("\nObserved: " + observed);
            messageBuffer.append("\nexpected - observed: " + delta);
            Assert.fail(messageBuffer.toString());
        }
    }

    /** verifies that two matrices are equal */
    public static void assertEquals(FieldMatrix<? extends FieldElement<?>> expected,
                                    FieldMatrix<? extends FieldElement<?>> observed) {

        Assert.assertNotNull("Observed should not be null",observed);

        if (expected.getColumnDimension() != observed.getColumnDimension() ||
                expected.getRowDimension() != observed.getRowDimension()) {
            StringBuilder messageBuffer = new StringBuilder();
            messageBuffer.append("Observed has incorrect dimensions.");
            messageBuffer.append("\nobserved is " + observed.getRowDimension() +
                    " x " + observed.getColumnDimension());
            messageBuffer.append("\nexpected " + expected.getRowDimension() +
                    " x " + expected.getColumnDimension());
            Assert.fail(messageBuffer.toString());
        }

        for (int i = 0; i < expected.getRowDimension(); ++i) {
            for (int j = 0; j < expected.getColumnDimension(); ++j) {
                FieldElement<?> eij = expected.getEntry(i, j);
                FieldElement<?> oij = observed.getEntry(i, j);
                Assert.assertEquals(eij, oij);
            }
        }
    }

    /** verifies that two arrays are close (sup norm) */
    public static void assertEquals(String msg, double[] expected, double[] observed, double tolerance) {
        StringBuilder out = new StringBuilder(msg);
        if (expected.length != observed.length) {
            out.append("\n Arrays not same length. \n");
            out.append("expected has length ");
            out.append(expected.length);
            out.append(" observed length = ");
            out.append(observed.length);
            Assert.fail(out.toString());
        }
        boolean failure = false;
        for (int i=0; i < expected.length; i++) {
            if (!Precision.equalsIncludingNaN(expected[i], observed[i], tolerance)) {
                failure = true;
                out.append("\n Elements at index ");
                out.append(i);
                out.append(" differ. ");
                out.append(" expected = ");
                out.append(expected[i]);
                out.append(" observed = ");
                out.append(observed[i]);
            }
        }
        if (failure) {
            Assert.fail(out.toString());
        }
    }

    /** verifies that two int arrays are equal */
    public static void assertEquals(int[] expected, int[] observed) {
        StringBuilder out = new StringBuilder();
        if (expected.length != observed.length) {
            out.append("\n Arrays not same length. \n");
            out.append("expected has length ");
            out.append(expected.length);
            out.append(" observed length = ");
            out.append(observed.length);
            Assert.fail(out.toString());
        }
        boolean failure = false;
        for (int i=0; i < expected.length; i++) {
            if (expected[i] != observed[i]) {
                failure = true;
                out.append("\n Elements at index ");
                out.append(i);
                out.append(" differ. ");
                out.append(" expected = ");
                out.append(expected[i]);
                out.append(" observed = ");
                out.append(observed[i]);
            }
        }
        if (failure) {
            Assert.fail(out.toString());
        }
    }

    /**
     * verifies that for i = 0,..., observed.length, observed[i] is within epsilon of one of the values in expected[i]
     * or observed[i] is NaN and expected[i] contains a NaN.
     */
    public static void assertContains(double[][] expected, double[] observed, double epsilon) {
        StringBuilder out = new StringBuilder();
        if (expected.length != observed.length) {
            out.append("\n Arrays not same length. \n");
            out.append("expected has length ");
            out.append(expected.length);
            out.append(" observed length = ");
            out.append(observed.length);
            Assert.fail(out.toString());
        }
        boolean failure = false;
        for (int i = 0; i < expected.length; i++) {
            boolean found = false;
            for (int j = 0; j < expected[i].length; j++) {
                if (Precision.equalsIncludingNaN(expected[i][j], observed[i], epsilon)) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                out.append("\n Observed element at index ");
                out.append(i);
                out.append(" is not among the expected values. ");
                out.append(" expected = " + Arrays.toString(expected[i]));
                out.append(" observed = ");
                out.append(observed[i]);
                failure = true;
            }
        }
        if (failure) {
            Assert.fail(out.toString());
        }
    }



    /** verifies that two int arrays are equal */
    public static void assertEquals(long[] expected, long[] observed) {
        StringBuilder out = new StringBuilder();
        if (expected.length != observed.length) {
            out.append("\n Arrays not same length. \n");
            out.append("expected has length ");
            out.append(expected.length);
            out.append(" observed length = ");
            out.append(observed.length);
            Assert.fail(out.toString());
        }
        boolean failure = false;
        for (int i=0; i < expected.length; i++) {
            if (expected[i] != observed[i]) {
                failure = true;
                out.append("\n Elements at index ");
                out.append(i);
                out.append(" differ. ");
                out.append(" expected = ");
                out.append(expected[i]);
                out.append(" observed = ");
                out.append(observed[i]);
            }
        }
        if (failure) {
            Assert.fail(out.toString());
        }
    }

    /** verifies that two arrays are equal */
    public static <T extends FieldElement<T>> void assertEquals(T[] m, T[] n) {
        if (m.length != n.length) {
            Assert.fail("vectors not same length");
        }
        for (int i = 0; i < m.length; i++) {
            Assert.assertEquals(m[i],n[i]);
        }
    }

    /**
     * Computes the sum of squared deviations of <values> from <target>
     * @param values array of deviates
     * @param target value to compute deviations from
     *
     * @return sum of squared deviations
     */
    public static double sumSquareDev(double[] values, double target) {
        double sumsq = 0d;
        for (int i = 0; i < values.length; i++) {
            final double dev = values[i] - target;
            sumsq += (dev * dev);
        }
        return sumsq;
    }

    /**
     * Asserts the null hypothesis for a ChiSquare test.  Fails and dumps arguments and test
     * statistics if the null hypothesis can be rejected with confidence 100 * (1 - alpha)%
     *
     * @param valueLabels labels for the values of the discrete distribution under test
     * @param expected expected counts
     * @param observed observed counts
     * @param alpha significance level of the test
     */
    public static void assertChiSquareAccept(String[] valueLabels, double[] expected, long[] observed, double alpha) {

        // Fail if we can reject null hypothesis that distributions are the same
        if (chiSquareTest(expected, observed) <= alpha) {
            StringBuilder msgBuffer = new StringBuilder();
            DecimalFormat df = new DecimalFormat("#.##");
            msgBuffer.append("Chisquare test failed");
            msgBuffer.append(" p-value = ");
            msgBuffer.append(chiSquareTest(expected, observed));
            msgBuffer.append(" chisquare statistic = ");
            msgBuffer.append(chiSquare(expected, observed));
            msgBuffer.append(". \n");
            msgBuffer.append("value\texpected\tobserved\n");
            for (int i = 0; i < expected.length; i++) {
                msgBuffer.append(valueLabels[i]);
                msgBuffer.append("\t");
                msgBuffer.append(df.format(expected[i]));
                msgBuffer.append("\t\t");
                msgBuffer.append(observed[i]);
                msgBuffer.append("\n");
            }
            msgBuffer.append("This test can fail randomly due to sampling error with probability ");
            msgBuffer.append(alpha);
            msgBuffer.append(".");
            Assert.fail(msgBuffer.toString());
        }
    }

    /**
     * Asserts the null hypothesis for a ChiSquare test.  Fails and dumps arguments and test
     * statistics if the null hypothesis can be rejected with confidence 100 * (1 - alpha)%
     *
     * @param values integer values whose observed and expected counts are being compared
     * @param expected expected counts
     * @param observed observed counts
     * @param alpha significance level of the test
     */
    public static void assertChiSquareAccept(int[] values, double[] expected, long[] observed, double alpha) {
        String[] labels = new String[values.length];
        for (int i = 0; i < values.length; i++) {
            labels[i] = Integer.toString(values[i]);
        }
        assertChiSquareAccept(labels, expected, observed, alpha);
    }

    /**
     * Asserts the null hypothesis for a ChiSquare test.  Fails and dumps arguments and test
     * statistics if the null hypothesis can be rejected with confidence 100 * (1 - alpha)%
     *
     * @param expected expected counts
     * @param observed observed counts
     * @param alpha significance level of the test
     */
    public static void assertChiSquareAccept(double[] expected, long[] observed, double alpha) {
        String[] labels = new String[expected.length];
        for (int i = 0; i < labels.length; i++) {
            labels[i] = Integer.toString(i + 1);
        }
        assertChiSquareAccept(labels, expected, observed, alpha);
    }

    /**
     * Asserts the null hypothesis that the sample follows the given distribution, using a G-test
     *
     * @param expectedDistribution distribution values are supposed to follow
     * @param values sample data
     * @param alpha significance level of the test
     */
    public static void assertGTest(final RealDistribution expectedDistribution, final double[] values, double alpha) {
        final int numBins = values.length / 30;
        final double[] breaks = new double[numBins];
        for (int b = 0; b < breaks.length; b++) {
            breaks[b] = expectedDistribution.inverseCumulativeProbability((double) b / numBins);
        }

        final long[] observed = new long[numBins];
        for (final double value : values) {
            int b = 0;
            do {
                b++;
            } while (b < numBins && value >= breaks[b]);

            observed[b - 1]++;
        }

        final double[] expected = new double[numBins];
        Arrays.fill(expected, (double) values.length / numBins);

        assertGTest(expected, observed, alpha);
    }

    /**
     * Asserts the null hypothesis that the observed counts follow the given distribution implied by expected,
     * using a G-test
     *
     * @param expected expected counts
     * @param observed observed counts
     * @param alpha significance level of the test
     */
    public static void assertGTest(final double[] expected, long[] observed, double alpha) {
        if (gTest(expected, observed) <  alpha) {
            StringBuilder msgBuffer = new StringBuilder();
            DecimalFormat df = new DecimalFormat("#.##");
            msgBuffer.append("G test failed");
            msgBuffer.append(" p-value = ");
            msgBuffer.append(gTest(expected, observed));
            msgBuffer.append(". \n");
            msgBuffer.append("value\texpected\tobserved\n");
            for (int i = 0; i < expected.length; i++) {
                msgBuffer.append(df.format(expected[i]));
                msgBuffer.append("\t\t");
                msgBuffer.append(observed[i]);
                msgBuffer.append("\n");
            }
            msgBuffer.append("This test can fail randomly due to sampling error with probability ");
            msgBuffer.append(alpha);
            msgBuffer.append(".");
            Assert.fail(msgBuffer.toString());
        }
    }

    /**
     * Computes the 25th, 50th and 75th percentiles of the given distribution and returns
     * these values in an array.
     */
    public static double[] getDistributionQuartiles(RealDistribution distribution) {
        double[] quantiles = new double[3];
        quantiles[0] = distribution.inverseCumulativeProbability(0.25d);
        quantiles[1] = distribution.inverseCumulativeProbability(0.5d);
        quantiles[2] = distribution.inverseCumulativeProbability(0.75d);
        return quantiles;
    }

    /**
     * Updates observed counts of values in quartiles.
     * counts[0] ↔ 1st quartile ... counts[3] ↔ top quartile
     */
    public static void updateCounts(double value, long[] counts, double[] quartiles) {
        if (value < quartiles[0]) {
            counts[0]++;
        } else if (value > quartiles[2]) {
            counts[3]++;
        } else if (value > quartiles[1]) {
            counts[2]++;
        } else {
            counts[1]++;
        }
    }

    /**
     * Eliminates points with zero mass from densityPoints and densityValues parallel
     * arrays.  Returns the number of positive mass points and collapses the arrays so
     * that the first <returned value> elements of the input arrays represent the positive
     * mass points.
     */
    public static int eliminateZeroMassPoints(int[] densityPoints, double[] densityValues) {
        int positiveMassCount = 0;
        for (int i = 0; i < densityValues.length; i++) {
            if (densityValues[i] > 0) {
                positiveMassCount++;
            }
        }
        if (positiveMassCount < densityValues.length) {
            int[] newPoints = new int[positiveMassCount];
            double[] newValues = new double[positiveMassCount];
            int j = 0;
            for (int i = 0; i < densityValues.length; i++) {
                if (densityValues[i] > 0) {
                    newPoints[j] = densityPoints[i];
                    newValues[j] = densityValues[i];
                    j++;
                }
            }
            System.arraycopy(newPoints,0,densityPoints,0,positiveMassCount);
            System.arraycopy(newValues,0,densityValues,0,positiveMassCount);
        }
        return positiveMassCount;
    }

    /*************************************************************************************
     * Stripped-down implementations of some basic statistics borrowed from hipparchus-stat.
     * NOTE: These implementations are NOT intended for reuse.  They are neither robust,
     * nor efficient; nor do they handle NaN, infinity or other corner cases in
     * a predictable way. They DO NOT CHECK PARAMETERS - the assumption is that incorrect
     * or meaningless results from bad parameters will trigger test failures in unit
     * tests using these methods.
     ************************************************************************************/

    /**
     * Returns p-value associated with null hypothesis that observed counts follow
     * expected distribution.  Will normalize inputs if necessary.
     *
     * @param expected expected counts
     * @param observed observed counts
     * @return p-value of Chi-square test
     */
    public static double chiSquareTest(final double[] expected, final long[] observed) {
            final org.hipparchus.distribution.continuous.ChiSquaredDistribution distribution =
                new ChiSquaredDistribution(expected.length - 1.0);
            return 1.0 - distribution.cumulativeProbability(chiSquare(expected, observed));
    }

    /**
     * Returns chi-square test statistic for expected and observed arrays. Rescales arrays
     * if necessary.
     *
     * @param expected expected counts
     * @param observed observed counts
     * @return chi-square statistic
     */
    public static double chiSquare(final double[] expected, final long[] observed) {
            double sumExpected = 0d;
            double sumObserved = 0d;
            for (int i = 0; i < observed.length; i++) {
                sumExpected += expected[i];
                sumObserved += observed[i];
            }
            double ratio = 1.0d;
            boolean rescale = false;
            if (FastMath.abs(sumExpected - sumObserved) > 10E-6) {
                ratio = sumObserved / sumExpected;
                rescale = true;
            }
            double sumSq = 0.0d;
            for (int i = 0; i < observed.length; i++) {
                if (rescale) {
                    final double dev = observed[i] - ratio * expected[i];
                    sumSq += dev * dev / (ratio * expected[i]);
                } else {
                    final double dev = observed[i] - expected[i];
                    sumSq += dev * dev / expected[i];
                }
            }
            return sumSq;

        }

    /**
     * Computes G-test statistic for expected, observed counts.
     * @param expected expected counts
     * @param observed observed counts
     * @return G statistic
     */
    private static double g(final double[] expected, final long[] observed) {
        double sumExpected = 0d;
        double sumObserved = 0d;
        for (int i = 0; i < observed.length; i++) {
            sumExpected += expected[i];
            sumObserved += observed[i];
        }
        double ratio = 1d;
        boolean rescale = false;
        if (FastMath.abs(sumExpected - sumObserved) > 10E-6) {
            ratio = sumObserved / sumExpected;
            rescale = true;
        }
        double sum = 0d;
        for (int i = 0; i < observed.length; i++) {
            final double dev = rescale ?
                    FastMath.log(observed[i] / (ratio * expected[i])) :
                        FastMath.log(observed[i] / expected[i]);
            sum += (observed[i]) * dev;
        }
        return 2d * sum;
    }

    /**
     * Computes p-value for G-test.
     *
     * @param expected expected counts
     * @param observed observed counts
     * @return p-value
     */
    private static double gTest(final double[] expected, final long[] observed) {
        final ChiSquaredDistribution distribution =
                new ChiSquaredDistribution(expected.length - 1.0);
        return 1.0 - distribution.cumulativeProbability(g(expected, observed));
    }

    /**
     * Computes the mean of the values in the array.
     *
     * @param values input values
     * @return arithmetic mean
     */
    public static double mean(final double[] values) {
        double sum = 0;
        for (double val : values) {
            sum += val;
        }
        return sum / values.length;
    }

    /**
     * Computes the (bias-adjusted) variance of the values in the input array.
     *
     * @param values input values
     * @return bias-adjusted variance
     */
    public static double variance(final double[] values) {
        final int length = values.length;
        final double mean = mean(values);
        double var = Double.NaN;
        if (length == 1) {
            var = 0.0;
        } else if (length > 1) {
            double accum = 0.0;
            double dev = 0.0;
            double accum2 = 0.0;
            for (int i = 0; i < length; i++) {
                dev = values[i] - mean;
                accum += dev * dev;
                accum2 += dev;
            }
            final double len = length;
            var = (accum - (accum2 * accum2 / len)) / (len - 1.0);
        }
        return var;
    }

    /**
     * Computes the standard deviation of the values in the input array.
     *
     * @param values input values
     * @return standard deviation
     */
    public static double standardDeviation(final double[] values) {
        return FastMath.sqrt(variance(values));
    }

    /**
     * Computes the median of the values in the input array.
     *
     * @param values input values
     * @return estimated median
     */
    public static double median(final double[] values) {
        final int len = values.length;
        final double[] sortedValues = Arrays.copyOf(values, len);
        Arrays.sort(sortedValues);
        if (len % 2 == 0) {
            return ((double)sortedValues[len/2] + (double)sortedValues[len/2 - 1])/2;
        } else {
            return (double) sortedValues[len/2];
        }
    }

    /**
     * Computes the covariance of the two input arrays.
     *
     * @param xArray first covariate
     * @param yArray second covariate
     * @return covariance
     */
    public static double covariance(final double[] xArray, final double[] yArray) {
        double result = 0d;
        final int length = xArray.length;
        final double xMean = mean(xArray);
        final double yMean = mean(yArray);
        for (int i = 0; i < length; i++) {
            final double xDev = xArray[i] - xMean;
            final double yDev = yArray[i] - yMean;
            result += (xDev * yDev - result) / (i + 1);
        }
        return result * ((double) length / (double)(length - 1));
    }

    /**
     * Computes a covariance matrix from a matrix whose columns represent covariates.
     *
     * @param matrix input matrix
     * @return covariance matrix
     */
    public static RealMatrix covarianceMatrix(RealMatrix matrix) {
        int dimension = matrix.getColumnDimension();
        final RealMatrix outMatrix = new BlockRealMatrix(dimension, dimension);
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < i; j++) {
                final double cov = covariance(matrix.getColumn(i), matrix.getColumn(j));
                outMatrix.setEntry(i, j, cov);
                outMatrix.setEntry(j, i, cov);
            }
            outMatrix.setEntry(i, i, variance(matrix.getColumn(i)));
        }
        return outMatrix;
    }

    public static double min(final double[] values) {
        double min = values[0];
        for (int i = 1; i < values.length; i++) {
            if (values[i] < min) {
                min = values[i];
            }
        }
        return min;
    }

    /**
     * Computes the maximum of the values in the input array.
     *
     * @param values input array
     * @return the maximum value
     */
    public static double max(final double[] values) {
        double max = values[0];
        for (int i = 1; i < values.length; i++) {
            if (values[i] > max) {
                max = values[i];
            }
        }
        return max;
    }

    /**
     * Unpacks a list of Doubles into a double[].
     *
     * @param values list of Double
     * @return double array
     */
    private static double[] unpack(List<Double> values) {
        int n = values.size();
        if (values == null || n == 0) {
            return new double[] {};
        }
        double[] out = new double[n];
        for (int i = 0; i < n; i++) {
            out[i] = values.get(i);
        }
        return out;
    }

    /**
     * Keeps track of the number of occurrences of distinct T instances
     * added via {@link #addValue(Object)}.
     *
     * @param <T> type of objects being tracked
     */
    public static class Frequency<T> {
        private Map<T, Integer> counts = new HashMap<>();
        public void addValue(T value) {
           Integer old = counts.put(value, 0);
           if (old != null) {
               counts.put(value, old++);
           }
        }
        public int getCount(T value) {
           Integer ret = counts.get(value);
           return ret == null ? 0 : ret;
        }
    }

    /**
     * Stripped down implementation of StreamingStatistics from o.h.stat.descriptive.
     * Actually holds all values in memory, so not suitable for very large streams of data.
     */
    public static class SimpleStatistics {
        private final List<Double> values = new ArrayList<>();
        public void addValue(double value) {
            values.add(value);
        }
        public double getMean() {
            return mean(unpack(values));
        }
        public double getStandardDeviation() {
            return standardDeviation(unpack(values));
        }
        public double getMin() {
            return min(unpack(values));
        }
        public double getMax() {
            return max(unpack(values));
        }
        public double getMedian() {
            return median(unpack(values));
        }
        public double getVariance() {
            return variance(unpack(values));
        }
        public long getN() {
            return values.size();
        }
    }

    /**
     * Stripped-down version of the bivariate regression class with the same name
     * in o.h.stat.regression.
     * Always estimates the model with an intercept term.
     */
    public static class SimpleRegression {
        private double sumX = 0d;
        private double sumXX = 0d;
        private double sumY = 0d;
        private double sumXY = 0d;
        private long n = 0;
        private double xbar = 0;
        private double ybar = 0;

        public void addData(double x, double y) {
            if (n == 0) {
                xbar = x;
                ybar = y;
            } else {
                final double fact1 = 1.0 + n;
                final double fact2 = n / (1.0 + n);
                final double dx = x - xbar;
                final double dy = y - ybar;
                sumXX += dx * dx * fact2;
                sumXY += dx * dy * fact2;
                xbar += dx / fact1;
                ybar += dy / fact1;
            }
            sumX += x;
            sumY += y;
            n++;
        }

        public double getSlope() {
            if (n < 2) {
                return Double.NaN; //not enough data
            }
            if (FastMath.abs(sumXX) < 10 * Double.MIN_VALUE) {
                return Double.NaN; //not enough variation in x
            }
            return sumXY / sumXX;
        }

        public double getIntercept() {
            return (sumY - getSlope() * sumX) / n;
        }
    }
}
