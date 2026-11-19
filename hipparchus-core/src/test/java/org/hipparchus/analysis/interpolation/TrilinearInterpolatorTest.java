/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.analysis.interpolation;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.analysis.*;
import org.hipparchus.random.RandomVectorGenerator;
import org.hipparchus.random.SobolSequenceGenerator;
import org.hipparchus.util.*;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class TrilinearInterpolatorTest {

    @Test
    void testConstant() {

        final double xMin = 0.0;
        final double xMax = 7.0;
        final int    nx   = 15;
        final double[] xVal = createLinearGrid(xMin, xMax, nx);

        final double yMin = -5.0;
        final double yMax = +5.0;
        final int    ny   = 11;
        final double[] yVal = createLinearGrid(yMin, yMax, ny);

        final double zMin = -2.0;
        final double zMax = 6.0;
        final int    nz   = 19;
        final double[] zVal = createLinearGrid(zMin, zMax, nz);

        final TrivariateFunction f = (x, y, z) -> 3.5;
        final CalculusFieldTrivariateFunction<Binary64> fT = (x, y, z) -> new Binary64(3.5);
        final TrilinearInterpolatingFunction bif = createInterpolatingFunction(xVal, yVal, zVal, f);

        assertEquals(xMin, bif.getXInf(), 1.0e-15);
        assertEquals(xMax, bif.getXSup(), 1.0e-15);
        assertEquals(yMin, bif.getYInf(), 1.0e-15);
        assertEquals(yMax, bif.getYSup(), 1.0e-15);
        assertEquals(zMin, bif.getZInf(), 1.0e-15);
        assertEquals(zMax, bif.getZSup(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, zVal, bif, f, fT, 1.0e-15);
        checkInterpolationRandom(new SobolSequenceGenerator(3), xMin, xMax, yMin, yMax, zMin, zMax,
                bif, f, fT, 1.0e-15);

    }

    @Test
    void testLinear() {

        final double xMin = -5.0;
        final double xMax = +5.0;
        final int    nx   = 11;
        final double[] xVal = createLinearGrid(xMin, xMax, nx);

        final double yMin = 0.0;
        final double yMax = 7.0;
        final int    ny   = 15;
        final double[] yVal = createLinearGrid(yMin, yMax, ny);

        final double zMin = -2.0;
        final double zMax = 6.0;
        final int    nz   = 19;
        final double[] zVal = createLinearGrid(zMin, zMax, nz);


        final TrivariateFunction f = (x, y, z) -> x + 2 * y - 3 * z;
        final CalculusFieldTrivariateFunction<Binary64> fT = new FieldTrivariateFunction() {
            @Override
            public <T extends CalculusFieldElement<T>> T value(T x, T y, T z) {
                return x.add(y.twice()).add(z.multiply(-3));
            }
        }.toCalculusFieldTrivariateFunction(Binary64Field.getInstance());
        TrilinearInterpolatingFunction bif = createInterpolatingFunction(xVal, yVal, zVal, f);

        assertEquals(xMin, bif.getXInf(), 1.0e-15);
        assertEquals(xMax, bif.getXSup(), 1.0e-15);
        assertEquals(yMin, bif.getYInf(), 1.0e-15);
        assertEquals(yMax, bif.getYSup(), 1.0e-15);
        assertEquals(zMin, bif.getZInf(), 1.0e-15);
        assertEquals(zMax, bif.getZSup(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, zVal, bif, f, fT, 5.0e-15);
        checkInterpolationRandom(new SobolSequenceGenerator(3), xMin, xMax, yMin, yMax, zMin, zMax,
                bif, f, fT, 5.0e-15);

    }

    @Test
    void testCubic() {

        final double xMin = -5.0;
        final double xMax = +5.0;
        final int    nx   = 11;
        final double[] xVal = createLinearGrid(xMin, xMax, nx);

        final double yMin = 0.0;
        final double yMax = 7.0;
        final int    ny   = 15;
        final double[] yVal = createLinearGrid(yMin, yMax, ny);

        final double zMin = -2.0;
        final double zMax = 6.0;
        final int    nz   = 19;
        final double[] zVal = createLinearGrid(zMin, zMax, nz);

        final TrivariateFunction f = (x, y, z) -> (3 * x - 2) * (6 - 0.5 * y) * (1 - 2 * z);
        final CalculusFieldTrivariateFunction<Binary64> fT = (x, y, z) -> x.multiply(3).subtract(2).multiply(y.multiply(-0.5).add(6))
                .multiply(z.multiply(-2).add(1));
        final TrilinearInterpolatingFunction bif = createInterpolatingFunction(xVal, yVal, zVal, f);

        assertEquals(xMin, bif.getXInf(), 1.0e-15);
        assertEquals(xMax, bif.getXSup(), 1.0e-15);
        assertEquals(yMin, bif.getYInf(), 1.0e-15);
        assertEquals(yMax, bif.getYSup(), 1.0e-15);
        assertEquals(zMin, bif.getZInf(), 1.0e-15);
        assertEquals(zMax, bif.getZSup(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, zVal, bif, f, fT, 3.0e-13);
        checkInterpolationRandom(new SobolSequenceGenerator(3), xMin, xMax, yMin, yMax, zMin, zMax,
                bif, f, fT, 3.0e-13);

    }

    @Test
    void testSinCos() {
        doTestSinCos(  10,   10, 10,1.2e-2);
        doTestSinCos( 100,  100, 100, 1.2e-4);
    }

    private void doTestSinCos(final int nx, final int ny, final int nz, final double tol) {
        final double xMin = -1.0;
        final double xMax = +2.0;
        final double[] xVal = createLinearGrid(xMin, xMax, nx);

        final double yMin = 0.0;
        final double yMax = 1.5;
        final double[] yVal = createLinearGrid(yMin, yMax, ny);

        final double zMin = -1.5;
        final double zMax = 1.0;
        final double[] zVal = createLinearGrid(zMin, zMax, nz);

        TrivariateFunction f = (x, y, z) -> FastMath.sin(x) * FastMath.cos(y) * FastMath.sin(0.5*z);
        CalculusFieldTrivariateFunction<Binary64> fT = (x, y, z) -> x.sin().multiply(y.cos()).multiply(z.half().sin());
        TrilinearInterpolatingFunction bif = createInterpolatingFunction(xVal, yVal, zVal, f);

        assertEquals(xMin, bif.getXInf(), 1.0e-15);
        assertEquals(xMax, bif.getXSup(), 1.0e-15);
        assertEquals(yMin, bif.getYInf(), 1.0e-15);
        assertEquals(yMax, bif.getYSup(), 1.0e-15);
        assertEquals(zMin, bif.getZInf(), 1.0e-15);
        assertEquals(zMax, bif.getZSup(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, zVal, bif, f, fT, 1.0e-15);
        checkInterpolationRandom(new SobolSequenceGenerator(3), xMin, xMax, yMin, yMax, zMin, zMax,
                bif, f, fT, tol);
    }

    private double[] createLinearGrid(final double min, final double max, final int n) {
        final double[] grid = new double[n];
        for (int i = 0; i < n; ++i) {
            grid[i] = ((n - 1 - i) * min + i * max) / (n - 1);
        }
        return grid;
    }

    private TrilinearInterpolatingFunction createInterpolatingFunction(double[] xVal, double[] yVal, double[] zVal,
                                                                      TrivariateFunction f) {
        final double[][][] fVal = new double[xVal.length][yVal.length][zVal.length];
        for (int i = 0; i < xVal.length; ++i) {
            for (int j = 0; j < yVal.length; ++j) {
                for (int k = 0; k < zVal.length; ++k) {
                    fVal[i][j][k] = f.value(xVal[i], yVal[j], zVal[k]);
                }
            }
        }
        return new TrilinearInterpolator().interpolate(xVal, yVal, zVal, fVal);
    }

    private void checkInterpolationAtNodes(final double[] xVal,
                                           final double[] yVal,
                                           final double[] zVal,
                                           final TrilinearInterpolatingFunction bif,
                                           final TrivariateFunction f,
                                           final CalculusFieldTrivariateFunction<Binary64> fT,
                                           final double tol) {

        for (final double x : xVal) {
            for (final double y : yVal) {
                for (final double z : zVal) {
                    assertEquals(f.value(x, y, z), bif.value(x, y, z), tol);

                    final Binary64 x64 = new Binary64(x);
                    final Binary64 y64 = new Binary64(y);
                    final Binary64 z64 = new Binary64(z);
                    assertEquals(fT.value(x64, y64, z64).getReal(), bif.value(x64, y64, z64).getReal(), tol);

                }
            }
        }
    }

    private void checkInterpolationRandom(final RandomVectorGenerator random,
                                          final double xMin, final double xMax,
                                          final double yMin, final double yMax,
                                          final double zMin, final double zMax,
                                          final TrilinearInterpolatingFunction bif,
                                          final TrivariateFunction f,
                                          final CalculusFieldTrivariateFunction<Binary64> fT,
                                          final double tol) {
        double maxError = 0.0;
        for (int i = 0; i < 10000; ++i) {

            final double[] v = random.nextVector();

            final double x = xMin + v[0] * (xMax - xMin);
            final double y = yMin + v[1] * (yMax - yMin);
            final double z = zMin + v[2] * (zMax - zMin);
            maxError = FastMath.max(maxError, FastMath.abs(f.value(x, y, z) - bif.value(x, y, z)));

            final Binary64 x64 = new Binary64(x);
            final Binary64 y64 = new Binary64(y);
            final Binary64 z64 = new Binary64(z);
            maxError = FastMath.max(maxError, FastMath.abs(fT.value(x64, y64, z64).getReal()- bif.value(x64, y64, z64).getReal()));
        }

        assertEquals(0.0, maxError, tol);

    }

}
