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
package org.hipparchus.analysis.interpolation;

import org.hipparchus.analysis.BivariateFunction;
import org.hipparchus.analysis.CalculusFieldBivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test case for the piecewise bicubic function.
 */
final class PiecewiseBicubicSplineInterpolatingFunctionTest {
    /**
     * Test preconditions.
     */
    @Test
    void testPreconditions() {
        double[] xval = new double[] { 3, 4, 5, 6.5, 7.5 };
        double[] yval = new double[] { -4, -3, -1, 2.5, 3.5 };
        double[][] zval = new double[xval.length][yval.length];

        @SuppressWarnings("unused")
        PiecewiseBicubicSplineInterpolatingFunction bcf = new PiecewiseBicubicSplineInterpolatingFunction(xval, yval, zval);

        try {
            new PiecewiseBicubicSplineInterpolatingFunction(null, yval, zval);
            fail("Failed to detect x null pointer");
        } catch (NullArgumentException iae) {
            // Expected.
        }

        try {
            new PiecewiseBicubicSplineInterpolatingFunction(xval, null, zval);
            fail("Failed to detect y null pointer");
        } catch (NullArgumentException iae) {
            // Expected.
        }

        try {
            new PiecewiseBicubicSplineInterpolatingFunction(xval, yval, null);
            fail("Failed to detect z null pointer");
        } catch (NullArgumentException iae) {
            // Expected.
        }

        try {
            final double[][] fnull = new double[1][];
            new PiecewiseBicubicSplineInterpolatingFunction(xval, yval, fnull);
            fail("Failed to detect z[0] null pointer");
        } catch (NullArgumentException iae) {
            // Expected.
        }

        try {
            final double[][] f = new double[1][1];
            new PiecewiseBicubicSplineInterpolatingFunction(new double[0], yval, f);
            fail("Failed to detect empty x pointer");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }

        try {
            final double[][] f = new double[1][1];
            new PiecewiseBicubicSplineInterpolatingFunction(xval, new double[0], f);
            fail("Failed to detect empty y pointer");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }


        try {
            final double[][] f = new double[1][0];
            new PiecewiseBicubicSplineInterpolatingFunction(xval, yval, f);
            fail("Failed to detect empty z[0] pointer");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }

        try {
            final double[] xval1 = { 0.0, 1.0, 2.0, 3.0 };
            new PiecewiseBicubicSplineInterpolatingFunction(xval1, yval, zval);
            fail("Failed to detect insufficient x data");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }

        try {
            final double[] yval1 = { 0.0, 1.0, 2.0, 3.0 };
            new PiecewiseBicubicSplineInterpolatingFunction(xval, yval1, zval);
            fail("Failed to detect insufficient y data");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }

        try {
            final double[][] zval1 = new double[4][4];
            new PiecewiseBicubicSplineInterpolatingFunction(xval, yval, zval1);
            fail("Failed to detect insufficient z data");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }

        try {
            final double[][] zval1 = new double[5][4];
            new PiecewiseBicubicSplineInterpolatingFunction(xval, yval, zval1);
            fail("Failed to detect insufficient z data");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }

        try {
            final double[] xval1 = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
            new PiecewiseBicubicSplineInterpolatingFunction(xval1, yval, zval);
            fail("Failed to detect data set array with different sizes.");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }

        try {
            final double[] yval1 = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
            new PiecewiseBicubicSplineInterpolatingFunction(xval, yval1, zval);
            fail("Failed to detect data set array with different sizes.");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }

        // X values not sorted.
        try {
            final double[] xval1 = { 0.0, 1.0, 0.5, 7.0, 3.5 };
            new PiecewiseBicubicSplineInterpolatingFunction(xval1, yval, zval);
            fail("Failed to detect unsorted x arguments.");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }

        // Y values not sorted.
        try {
            final double[] yval1 = { 0.0, 1.0, 1.5, 0.0, 3.0 };
            new PiecewiseBicubicSplineInterpolatingFunction(xval, yval1, zval);
            fail("Failed to detect unsorted y arguments.");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }
    }

    /**
     * Interpolating a plane.
     * <p>
     * z = 2 x - 3 y + 5
     */
    @Test
    void testPlane() {
        final int numberOfElements = 10;
        final double minimumX = -10;
        final double maximumX = 10;
        final double minimumY = -10;
        final double maximumY = 10;
        final int numberOfSamples = 100;

        final double interpolationTolerance = 7e-15;
        final double maxTolerance = 6e-14;

        // Function values
        BivariateFunction f = (x, y) -> 2 * x - 3 * y + 5;
        CalculusFieldBivariateFunction<Binary64> fT = (x, y) -> x.twice().subtract(y.multiply(3)).add(5);

        testInterpolation(minimumX,
                          maximumX,
                          minimumY,
                          maximumY,
                          numberOfElements,
                          numberOfSamples,
                          f, fT,
                          interpolationTolerance,
                          maxTolerance);
    }

    /**
     * Interpolating a paraboloid.
     * <p>
     * z = 2 x<sup>2</sup> - 3 y<sup>2</sup> + 4 x y - 5
     */
    @Test
    void testParabaloid() {
        final int numberOfElements = 10;
        final double minimumX = -10;
        final double maximumX = 10;
        final double minimumY = -10;
        final double maximumY = 10;
        final int numberOfSamples = 100;

        final double interpolationTolerance = 2e-14;
        final double maxTolerance = 6e-14;

        // Function values
        BivariateFunction f = (x, y) -> 2 * x * x - 3 * y * y + 4 * x * y - 5;
        CalculusFieldBivariateFunction<Binary64> fT = (x, y) -> x.square().twice().
                                                             subtract(y.square().multiply(3)).
                                                             add(x.multiply(y).multiply(4)).
                                                             subtract(5);

        testInterpolation(minimumX,
                          maximumX,
                          minimumY,
                          maximumY,
                          numberOfElements,
                          numberOfSamples,
                          f, fT,
                          interpolationTolerance,
                          maxTolerance);
    }

    /**
     * @param minimumX Lower bound of interpolation range along the x-coordinate.
     * @param maximumX Higher bound of interpolation range along the x-coordinate.
     * @param minimumY Lower bound of interpolation range along the y-coordinate.
     * @param maximumY Higher bound of interpolation range along the y-coordinate.
     * @param numberOfElements Number of data points (along each dimension).
     * @param numberOfSamples Number of test points.
     * @param f Function to test.
     * @param fT Binary64 version of the function to test
     * @param meanTolerance Allowed average error (mean error on all interpolated values).
     * @param maxTolerance Allowed error on each interpolated value.
     */
    private void testInterpolation(double minimumX,
                                   double maximumX,
                                   double minimumY,
                                   double maximumY,
                                   int numberOfElements,
                                   int numberOfSamples,
                                   BivariateFunction f,
                                   final CalculusFieldBivariateFunction<Binary64> fT,
                                   double meanTolerance,
                                   double maxTolerance) {
        double expected;
        double actual;
        Binary64 expected64;
        Binary64 actual64;
        double currentX;
        double currentY;
        Binary64 currentX64;
        Binary64 currentY64;
        final double deltaX = (maximumX - minimumX) / ((double) numberOfElements);
        final double deltaY = (maximumY - minimumY) / ((double) numberOfElements);
        final double[] xValues = new double[numberOfElements];
        final double[] yValues = new double[numberOfElements];
        final double[][] zValues = new double[numberOfElements][numberOfElements];

        for (int i = 0; i < numberOfElements; i++) {
            xValues[i] = minimumX + deltaX * (double) i;
            for (int j = 0; j < numberOfElements; j++) {
                yValues[j] = minimumY + deltaY * (double) j;
                zValues[i][j] = f.value(xValues[i], yValues[j]);
            }
        }

        final PiecewiseBicubicSplineInterpolatingFunction interpolation
            = new PiecewiseBicubicSplineInterpolatingFunction(xValues,
                                                              yValues,
                                                              zValues);

        for (int i = 0; i < numberOfElements; i++) {
            currentX = xValues[i];
            currentX64 = new Binary64(currentX);
            for (int j = 0; j < numberOfElements; j++) {
                currentY = yValues[j];
                currentY64 = new Binary64(currentY);
                assertTrue(Precision.equals(f.value(currentX, currentY),
                                                   interpolation.value(currentX, currentY)));
                assertTrue(Precision.equals(fT.value(currentX64, currentY64).getReal(),
                                                   interpolation.value(currentX64, currentY64).getReal()));
            }
        }

        final RandomDataGenerator gen = new RandomDataGenerator(1234567L);

        double sumError   = 0;
        double sumError64 = 0;
        for (int i = 0; i < numberOfSamples; i++) {
            currentX = gen.nextUniform(xValues[0], xValues[xValues.length - 1]);
            currentY = gen.nextUniform(yValues[0], yValues[yValues.length - 1]);
            currentX64 = new Binary64(currentX);
            currentY64 = new Binary64(currentY);
            expected = f.value(currentX, currentY);
            actual = interpolation.value(currentX, currentY);
            expected64 = fT.value(currentX64, currentY64);
            actual64 = interpolation.value(currentX64, currentY64);
            sumError   += FastMath.abs(actual - expected);
            sumError64 += FastMath.abs(actual64.subtract(expected64)).getReal();
            assertEquals(expected, actual, maxTolerance);
            assertEquals(expected64.getReal(), actual64.getReal(), maxTolerance);
        }

        final double meanError = sumError / numberOfSamples;
        assertEquals(0, meanError, meanTolerance);
        final double meanError64 = sumError64 / numberOfSamples;
        assertEquals(0, meanError64, meanTolerance);

    }

    @Test
    void testIsValidPoint() {
        // GIVEN
        final double[] x = new double[] { 0., 1., 2., 3., 4. };
        final double[] y = x.clone();
        final PiecewiseBicubicSplineInterpolatingFunction interpolatingFunction =
                new PiecewiseBicubicSplineInterpolatingFunction(x, y, new double[x.length][y.length]);
        // WHEN & THEN
        assertFalse(interpolatingFunction.isValidPoint(x[0] - 1, y[0]));
        assertFalse(interpolatingFunction.isValidPoint(x[x.length - 1] + 1, y[0]));
        assertFalse(interpolatingFunction.isValidPoint(x[0], y[0] - 1));
        assertFalse(interpolatingFunction.isValidPoint(x[0], y[y.length - 1] + 1));
    }

}
