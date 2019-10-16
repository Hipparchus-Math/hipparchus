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
package org.hipparchus.analysis.interpolation;

import org.hipparchus.RealFieldElement;
import org.hipparchus.analysis.BivariateFunction;
import org.hipparchus.analysis.FieldBivariateFunction;
import org.hipparchus.analysis.RealFieldBivariateFunction;
import org.hipparchus.random.RandomVectorGenerator;
import org.hipparchus.random.SobolSequenceGenerator;
import org.hipparchus.util.Decimal64;
import org.hipparchus.util.Decimal64Field;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class BilinearInterpolatorTest {

    @Test
    public void testConstant() {

        double xMin = 0.0;
        double xMax = 7.0;
        int    nx   = 15;
        double[] xVal = createLinearGrid(xMin, xMax, nx);

        double yMin = -5.0;
        double yMax = +5.0;
        int    ny   = 11;
        double[] yVal = createLinearGrid(yMin, yMax, ny);

        BivariateFunction f = (x, y) -> 3.5;
        RealFieldBivariateFunction<Decimal64> fT = (x, y) -> new Decimal64(3.5);
        BilinearInterpolatingFunction bif = createInterpolatingFunction(xVal, yVal, f);

        Assert.assertEquals(xMin, bif.getXInf(), 1.0e-15);
        Assert.assertEquals(xMax, bif.getXSup(), 1.0e-15);
        Assert.assertEquals(yMin, bif.getYInf(), 1.0e-15);
        Assert.assertEquals(yMax, bif.getYSup(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, bif, f, fT, 1.0e-15);
        checkInterpolationRandom(new SobolSequenceGenerator(2), xMin, xMax, yMin, yMax, bif, f, fT, 1.0e-15);

    }

    @Test
    public void testLinear() {

        double xMin = -5.0;
        double xMax = +5.0;
        int    nx   = 11;
        double[] xVal = createLinearGrid(xMin, xMax, nx);

        double yMin = 0.0;
        double yMax = 7.0;
        int    ny   = 15;
        double[] yVal = createLinearGrid(yMin, yMax, ny);

        BivariateFunction f = (x, y) -> 2 * x - y;
        RealFieldBivariateFunction<Decimal64> fT = new FieldBivariateFunction() {
            @Override
            public <T extends RealFieldElement<T>> T value(T x, T y) {
                return x.multiply(2).subtract(y);
            }
        }.toRealFieldBivariateFunction(Decimal64Field.getInstance());
        BilinearInterpolatingFunction bif = createInterpolatingFunction(xVal, yVal, f);

        Assert.assertEquals(xMin, bif.getXInf(), 1.0e-15);
        Assert.assertEquals(xMax, bif.getXSup(), 1.0e-15);
        Assert.assertEquals(yMin, bif.getYInf(), 1.0e-15);
        Assert.assertEquals(yMax, bif.getYSup(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, bif, f, fT, 1.0e-15);
        checkInterpolationRandom(new SobolSequenceGenerator(2), xMin, xMax, yMin, yMax, bif, f, fT, 1.0e-15);

    }

    @Test
    public void testQuadratic() {

        double xMin = -5.0;
        double xMax = +5.0;
        int    nx   = 11;
        double[] xVal = createLinearGrid(xMin, xMax, nx);

        double yMin = 0.0;
        double yMax = 7.0;
        int    ny   = 15;
        double[] yVal = createLinearGrid(yMin, yMax, ny);

        BivariateFunction f = (x, y) -> (3 * x - 2) * (6 - 0.5 * y);
        RealFieldBivariateFunction<Decimal64> fT = (x, y) -> x.multiply(3).subtract(2).multiply(y.multiply(-0.5).add(6));
        BilinearInterpolatingFunction bif = createInterpolatingFunction(xVal, yVal, f);

        Assert.assertEquals(xMin, bif.getXInf(), 1.0e-15);
        Assert.assertEquals(xMax, bif.getXSup(), 1.0e-15);
        Assert.assertEquals(yMin, bif.getYInf(), 1.0e-15);
        Assert.assertEquals(yMax, bif.getYSup(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, bif, f, fT, 1.0e-15);
        checkInterpolationRandom(new SobolSequenceGenerator(2), xMin, xMax, yMin, yMax, bif, f, fT, 1.0e-15);

    }

    @Test
    public void testSinCos() {
        doTestSinCos(  10,   10, 1.8e-2);
        doTestSinCos( 100,  100, 1.5e-4);
        doTestSinCos(1000, 1000, 1.4e-6);
    }

    private void doTestSinCos(final int nx, final int ny, final double tol) {
        double xMin = -1.0;
        double xMax = +2.0;
        double[] xVal = createLinearGrid(xMin, xMax, nx);

        double yMin = 0.0;
        double yMax = 1.5;
        double[] yVal = createLinearGrid(yMin, yMax, ny);

        BivariateFunction f = (x, y) -> FastMath.sin(x) * FastMath.cos(y);
        RealFieldBivariateFunction<Decimal64> fT = (x, y) -> FastMath.sin(x).multiply(FastMath.cos(y));
        BilinearInterpolatingFunction bif = createInterpolatingFunction(xVal, yVal, f);

        Assert.assertEquals(xMin, bif.getXInf(), 1.0e-15);
        Assert.assertEquals(xMax, bif.getXSup(), 1.0e-15);
        Assert.assertEquals(yMin, bif.getYInf(), 1.0e-15);
        Assert.assertEquals(yMax, bif.getYSup(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, bif, f, fT, 1.0e-15);
        checkInterpolationRandom(new SobolSequenceGenerator(2), xMin, xMax, yMin, yMax, bif, f, fT, tol);

    }

    private double[] createLinearGrid(final double min, final double max, final int n) {
        final double[] grid = new double[n];
        for (int i = 0; i < n; ++i) {
            grid[i] = ((n - 1 - i) * min + i * max) / (n - 1);
        }
        return grid;
    }

    private BilinearInterpolatingFunction createInterpolatingFunction(double[] xVal, double[] yVal,
                                                                      BivariateFunction f) {
        final double[][] fVal = new double[xVal.length][yVal.length];
        for (int i = 0; i < xVal.length; ++i) {
            for (int j = 0; j < yVal.length; ++j) {
                fVal[i][j] = f.value(xVal[i], yVal[j]);
            }
        }
        return new BilinearInterpolator().interpolate(xVal, yVal, fVal);
    }

    private void checkInterpolationAtNodes(final double[] xVal,
                                           final double[] yVal,
                                           final BilinearInterpolatingFunction bif,
                                           final BivariateFunction f,
                                           final RealFieldBivariateFunction<Decimal64> fT,
                                           final double tol) {

        for (int i = 0; i < xVal.length; ++i) {
            for (int j = 0; j < yVal.length; ++j) {

                final double x = xVal[i];
                final double y = yVal[j];
                Assert.assertEquals(f.value(x, y), bif.value(x, y), tol);

                final Decimal64 x64 = new Decimal64(x);
                final Decimal64 y64 = new Decimal64(y);
                Assert.assertEquals(fT.value(x64, y64).getReal(), bif.value(x64, y64).getReal(), tol);

            }
        }
    }

    private void checkInterpolationRandom(final RandomVectorGenerator random,
                                          final double xMin, final double xMax,
                                          final double yMin, final double yMax,
                                          final BilinearInterpolatingFunction bif,
                                          final BivariateFunction f,
                                          final RealFieldBivariateFunction<Decimal64> fT,
                                          final double tol) {
        double maxError = 0.0;
        for (int i = 0; i < 10000; ++i) {

            final double[] v = random.nextVector();

            final double x = xMin + v[0] * (xMax - xMin);
            final double y = yMin + v[1] * (yMax - yMin);
            maxError = FastMath.max(maxError, FastMath.abs(f.value(x, y) - bif.value(x, y)));

            final Decimal64 x64 = new Decimal64(x);
            final Decimal64 y64 = new Decimal64(y);
            maxError = FastMath.max(maxError, FastMath.abs(fT.value(x64, y64).getReal()- bif.value(x64, y64).getReal()));
        }

        Assert.assertEquals(0.0, maxError, tol);

    }

}
