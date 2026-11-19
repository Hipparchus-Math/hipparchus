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

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.analysis.BivariateFunction;
import org.hipparchus.analysis.CalculusFieldBivariateFunction;
import org.hipparchus.analysis.FieldBivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test case for the bicubic interpolator.
 */
final class BicubicInterpolatorTest {
    /**
     * Test preconditions.
     */
    @Test
    void testPreconditions() {
        double[] xval = new double[] {3, 4, 5, 6.5};
        double[] yval = new double[] {-4, -3, -1, 2.5};
        double[][] zval = new double[xval.length][yval.length];

        BivariateGridInterpolator interpolator = new BicubicInterpolator();

        @SuppressWarnings("unused")
        BivariateFunction p = interpolator.interpolate(xval, yval, zval);

        double[] wxval = new double[] {3, 2, 5, 6.5};
        try {
            p = interpolator.interpolate(wxval, yval, zval);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }

        double[] wyval = new double[] {-4, -3, -1, -1};
        try {
            p = interpolator.interpolate(xval, wyval, zval);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }

        double[][] wzval = new double[xval.length][yval.length + 1];
        try {
            p = interpolator.interpolate(xval, yval, wzval);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        wzval = new double[xval.length - 1][yval.length];
        try {
            p = interpolator.interpolate(xval, yval, wzval);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
    }

    /**
     * Interpolating a plane.
     * <p>
     * z = 2 x - 3 y + 5
     */
    @Test
    void testPlane() {
        BivariateFunction f = (x, y) -> 2 * x - 3 * y + 5;
        CalculusFieldBivariateFunction<Binary64> fT = new FieldBivariateFunction() {
            @Override
            public <T extends CalculusFieldElement<T>> T value(T x, T y) {
                return x.twice().add(y.multiply(-3)).add(5);
            }
        }.toCalculusFieldBivariateFunction(Binary64Field.getInstance());

        testInterpolation(3000,
                          1e-13,
                          f,
                          fT,
                          false);
    }

    /**
     * Interpolating a paraboloid.
     * <p>
     * z = 2 x<sup>2</sup> - 3 y<sup>2</sup> + 4 x y - 5
     */
    @Test
    void testParaboloid() {
        BivariateFunction f = (x, y) -> 2 * x * x - 3 * y * y + 4 * x * y - 5;
        CalculusFieldBivariateFunction<Binary64> fT = (x, y) -> x.multiply(x).twice()
                .add(y.multiply(y).multiply(-3))
                .add(x.multiply(y).multiply(4))
                .add(-5);

        testInterpolation(3000,
                          1e-12,
                          f,
                          fT,
                          false);
    }

    /**
     * @param numSamples Number of test samples.
     * @param tolerance Allowed tolerance on the interpolated value.
     * @param f Test function.
     * @param print Whether to print debugging output to the console.
     */
    private void testInterpolation(int numSamples,
                                   double tolerance,
                                   BivariateFunction f,
                                   CalculusFieldBivariateFunction<Binary64> fT,
                                   boolean print) {
        final int sz = 21;
        final double[] xval = new double[sz];
        final double[] yval = new double[sz];
        // Coordinate values
        final double delta = 1d / (sz - 1);
        for (int i = 0; i < sz; i++) {
            xval[i] = -1 + 15 * i * delta;
            yval[i] = -20 + 30 * i * delta;
        }

        final double[][] zval = new double[xval.length][yval.length];
        for (int i = 0; i < xval.length; i++) {
            for (int j = 0; j < yval.length; j++) {
                zval[i][j] = f.value(xval[i], yval[j]);
            }
        }

        final BicubicInterpolator interpolator = new BicubicInterpolator();
        final BicubicInterpolatingFunction p = interpolator.interpolate(xval, yval, zval);
        double x, y;

        final RandomDataGenerator gen = new RandomDataGenerator(1234567L);
        int count = 0;
        while (true) {
            x = gen.nextUniform(xval[0], xval[xval.length - 1]);
            y = gen.nextUniform(yval[0], yval[yval.length - 1]);
            if (!p.isValidPoint(x, y)) {
                if (print) {
                    System.out.println("# " + x + " " + y);
                }
                continue;
            }

            if (count++ > numSamples) {
                break;
            }
            final double expected = f.value(x, y);
            final double actual = p.value(x, y);

            if (print) {
                System.out.println(x + " " + y + " " + expected + " " + actual);
            }

            assertEquals(expected, actual, tolerance);

            final Binary64 x64 = new Binary64(x);
            final Binary64 y64 = new Binary64(y);
            final Binary64 expectedT = fT.value(x64, y64);
            final Binary64 actualT = p.value(x64, y64);

            assertEquals(expectedT.getReal(), actualT.getReal(), tolerance);
        }
    }
}
