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

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import org.hipparchus.analysis.MultivariateFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.UnitSphereRandomVectorGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.random.Well19937a;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test case for the {@link MicrosphereProjectionInterpolator
 * "microsphere projection"} interpolator.
 */
public final class MicrosphereProjectionInterpolatorTest {
    /**
     * Test of interpolator for a plane.
     * <p>
     * y = 2 x<sub>1</sub> - 3 x<sub>2</sub> + 5
     */
    @Test
    public void testLinearFunction2D() {
        MultivariateFunction f = new MultivariateFunction() {
            @Override
            public double value(double[] x) {
                if (x.length != 2) {
                    throw new IllegalArgumentException();
                }
                return 2 * x[0] - 3 * x[1] + 5;
            }
        };

        final double darkFraction = 0.5;
        final double darkThreshold = 1e-2;
        final double background = Double.NaN;
        final double exponent = 1.1;
        final boolean shareSphere = true;
        final double noInterpolationTolerance = Math.ulp(1d);
        final RandomGenerator random = new Well1024a(0x1c7a150c83a6d9dal);

        // N-dimensional interpolator.
        final MultivariateInterpolator interpolator
            = new MicrosphereProjectionInterpolator(new InterpolatingMicrosphere(2, 500,
                                                    darkFraction,
                                                    darkThreshold,
                                                    background,
                                                    new UnitSphereRandomVectorGenerator(2, random)),
                                                    exponent,
                                                    shareSphere,
                                                    noInterpolationTolerance);

        // 2D interpolator.
        final MultivariateInterpolator interpolator2D
            = new MicrosphereProjectionInterpolator(new InterpolatingMicrosphere2D(16,
                                                                                   darkFraction,
                                                                                   darkThreshold,
                                                                                   background),
                                                    exponent,
                                                    shareSphere,
                                                    noInterpolationTolerance);

        final double min = -1;
        final double max = 1;
        final double range = max - min;
        final int res = 5;
        final int n = res * res; // Number of sample points.
        final int dim = 2;
        double[][] x = new double[n][dim];
        double[] y = new double[n];
        int index = 0;
        for (int i = 0; i < res; i++) {
            final double x1Val = toCoordinate(min, range, res, i);
            for (int j = 0; j < res; j++) {
                final double x2Val = toCoordinate(min, range, res, j);
                x[index][0] = x1Val;
                x[index][1] = x2Val;
                y[index] = f.value(x[index]);
                ++index;
            }
        }

        final MultivariateFunction p = interpolator.interpolate(x, y);
        final MultivariateFunction p2D = interpolator2D.interpolate(x, y);

        double[] c = new double[dim];
        double expected, result, result2D;

        final int sampleIndex = 2;
        c[0] = x[sampleIndex][0];
        c[1] = x[sampleIndex][1];
        expected = f.value(c);
        result = p.value(c);
        result2D = p2D.value(c);
        Assert.assertEquals("on sample point (exact)", expected, result2D, FastMath.ulp(1d));
        Assert.assertEquals("on sample point (ND vs 2D)", result2D, result, FastMath.ulp(1d));

        // Interpolation.
        c[0] = 0.654321;
        c[1] = -0.345678;
        expected = f.value(c);
        result = p.value(c);
        result2D = p2D.value(c);
        Assert.assertEquals("interpolation (exact)", expected, result2D, 1e-1);
        Assert.assertEquals("interpolation (ND vs 2D)", result2D, result, 1e-1);

        // Extrapolation.
        c[0] = 0 - 1e-2;
        c[1] = 1 + 1e-2;
        expected = f.value(c);
        result = p.value(c);
        result2D = p2D.value(c);
        Assert.assertFalse(Double.isNaN(result));
        Assert.assertFalse(Double.isNaN(result2D));
        Assert.assertEquals("extrapolation (exact)", expected, result2D, 1e-1);
        Assert.assertEquals("extrapolation (ND vs 2D)", result2D, result, 1e-2);

        // Far away.
        c[0] = 20;
        c[1] = -30;
        result = p.value(c);
        Assert.assertTrue(result + " should be NaN", Double.isNaN(result));
        result2D = p2D.value(c);
        Assert.assertTrue(result2D + " should be NaN", Double.isNaN(result2D));
    }
    
    @Test
    public void testWrongDimensions() {
        checkWrongArguments(0, 1, 0.5, 0.0, 0.0,
                            LocalizedCoreFormats.NUMBER_TOO_SMALL_BOUND_EXCLUDED);
        checkWrongArguments(1, 0, 0.5, 0.0, 0.0,
                            LocalizedCoreFormats.NUMBER_TOO_SMALL_BOUND_EXCLUDED);
        checkWrongArguments(1, 1, 0.5, -1.0, 0.0,
                            LocalizedCoreFormats.NUMBER_TOO_SMALL);
    }

    private void checkWrongArguments(int dimension,
                                     int size,
                                     double maxDarkFraction,
                                     double darkThreshold,
                                     double background,
                                     LocalizedCoreFormats expected) {
        try {
            new InterpolatingMicrosphere(dimension, size, maxDarkFraction, darkThreshold, background, null);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(expected, miae.getSpecifier());
        }
    }

    @Test
    public void testCopy() {
        UnitSphereRandomVectorGenerator random =
                        new UnitSphereRandomVectorGenerator(3, new Well19937a(0x265318ael));
        InterpolatingMicrosphere original = new InterpolatingMicrosphere(3, 30, 0.5, 0.2, 0.0, random);
        InterpolatingMicrosphere copy     = original.copy();
        Assert.assertFalse(original == copy);
        Assert.assertEquals(original.getDimension(), copy.getDimension());
        Assert.assertEquals(original.getSize(),      copy.getSize());
    }

    @Test
    public void testSizeLimit() {
        UnitSphereRandomVectorGenerator random =
                        new UnitSphereRandomVectorGenerator(3, new Well19937a(0x453l));
        InterpolatingMicrosphere ims = new InterpolatingMicrosphere(3, 30, 0.5, 0.2, 0.0, random);
        try {
            Method add = InterpolatingMicrosphere.class.getDeclaredMethod("add", double[].class, Boolean.TYPE);
            add.setAccessible(true);
            try {
                add.invoke(ims, random.nextVector(), true);
                Assert.fail("an exception should have been thrown");
            } catch (InvocationTargetException ite) {
                MathIllegalStateException miae = (MathIllegalStateException) ite.getCause();
                Assert.assertEquals(LocalizedCoreFormats.MAX_COUNT_EXCEEDED, miae.getSpecifier());
            }
        } catch (NoSuchMethodException | IllegalAccessException | IllegalArgumentException e) {
            Assert.fail(e.getLocalizedMessage());
        }
    }

    @Test
    public void testInconsistentDimensions() {
        final int d1 = 5;
        final int d2 = 3;
        try {
            UnitSphereRandomVectorGenerator random =
                            new UnitSphereRandomVectorGenerator(d1, new Well19937a(0x1l));
            new InterpolatingMicrosphere(d2, 30, 0.5, 0.2, 0.0, random);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            Assert.assertEquals(d1, ((Integer) miae.getParts()[0]).intValue());
            Assert.assertEquals(d2, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    /**
     * @param min Minimum of the coordinate range.
     * @param range Extent of the coordinate interval.
     * @param res Number of pixels.
     * @param pixel Pixel index.
     */
    private static double toCoordinate(double min,
                                       double range,
                                       int res,
                                       int pixel) {
        return pixel * range / (res - 1) + min;
    }
}
