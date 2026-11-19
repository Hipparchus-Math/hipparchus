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

import org.hipparchus.analysis.CalculusFieldTrivariateFunction;
import org.hipparchus.analysis.differentiation.UnivariateDerivative1;
import org.hipparchus.random.RandomVectorGenerator;
import org.hipparchus.random.SobolSequenceGenerator;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test case for the field trilinear function.
 */
final class FieldTrilinearInterpolatingFunctionTest {

    @Test
    void testConstant() {

        final double xMin = 0.0;
        final double xMax = 7.0;
        final int    nx   = 15;
        final UnivariateDerivative1[] xVal = createLinearGrid(xMin, xMax, nx, 1.0);

        final double yMin = -5.0;
        final double yMax =  5.0;
        final int    ny   = 11;
        final UnivariateDerivative1[] yVal = createLinearGrid(yMin, yMax, ny, 2.0);

        final double zMin = -2.0;
        final double zMax = 6.0;
        final int    nz   = 19;
        final UnivariateDerivative1[] zVal = createLinearGrid(zMin, zMax, nz, 3.0);

        final CalculusFieldTrivariateFunction<UnivariateDerivative1> f  = (x, y, z) -> new UnivariateDerivative1(3.5, 1.0);
        final FieldTrilinearInterpolatingFunction<UnivariateDerivative1> bif = createInterpolatingFunction(xVal, yVal, zVal, f);

        assertEquals(xMin, bif.getXInf().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXInf().getFirstDerivative(), 1.0e-15);
        assertEquals(xMax, bif.getXSup().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXSup().getFirstDerivative(), 1.0e-15);
        assertEquals(yMin, bif.getYInf().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYInf().getFirstDerivative(), 1.0e-15);
        assertEquals(yMax, bif.getYSup().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYSup().getFirstDerivative(), 1.0e-15);
        assertEquals(zMin, bif.getZInf().getValue(),           1.0e-15);
        assertEquals(3.0,  bif.getZInf().getFirstDerivative(), 1.0e-15);
        assertEquals(zMax, bif.getZSup().getValue(),           1.0e-15);
        assertEquals(3.0,  bif.getZSup().getFirstDerivative(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, zVal, bif, f, 1.0e-15, 1.0e-15);
        checkInterpolationRandom(new SobolSequenceGenerator(3), xMin, xMax, yMin, yMax, zMin, zMax,
                bif, f, 1.0e-15, 2.0e-15);

    }

    @Test
    void testLinear() {

        final double xMin = -5.0;
        final double xMax =  5.0;
        final int    nx   = 11;
        final UnivariateDerivative1[] xVal = createLinearGrid(xMin, xMax, nx, 1.0);

        final double yMin = 0.0;
        final double yMax = 7.0;
        final int    ny   = 15;
        final UnivariateDerivative1[] yVal = createLinearGrid(yMin, yMax, ny, 2.0);

        final double zMin = -2.0;
        final double zMax = 6.0;
        final int    nz   = 19;
        final UnivariateDerivative1[] zVal = createLinearGrid(zMin, zMax, nz, 3.0);

        final CalculusFieldTrivariateFunction<UnivariateDerivative1> f = (x, y, z) -> x.add(y.twice()).add(z.multiply(-3));
        final FieldTrilinearInterpolatingFunction<UnivariateDerivative1> bif = createInterpolatingFunction(xVal, yVal, zVal, f);

        assertEquals(xMin, bif.getXInf().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXInf().getFirstDerivative(), 1.0e-15);
        assertEquals(xMax, bif.getXSup().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXSup().getFirstDerivative(), 1.0e-15);
        assertEquals(yMin, bif.getYInf().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYInf().getFirstDerivative(), 1.0e-15);
        assertEquals(yMax, bif.getYSup().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYSup().getFirstDerivative(), 1.0e-15);
        assertEquals(zMin, bif.getZInf().getValue(),           1.0e-15);
        assertEquals(3.0,  bif.getZInf().getFirstDerivative(), 1.0e-15);
        assertEquals(zMax, bif.getZSup().getValue(),           1.0e-15);
        assertEquals(3.0,  bif.getZSup().getFirstDerivative(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, zVal, bif, f, 1.0e-15, 1.0e-15);
        checkInterpolationRandom(new SobolSequenceGenerator(3), xMin, xMax, yMin, yMax, zMin, zMax,
                bif, f, 5.0e-15, 5.0e-14);
    }

    @Test
    void testCubic() {

        final double xMin = -5.0;
        final double xMax =  5.0;
        final int    nx   = 11;
        final UnivariateDerivative1[] xVal = createLinearGrid(xMin, xMax, nx, 1.0);

        final double yMin = 0.0;
        final double yMax = 7.0;
        final int    ny   = 15;
        final UnivariateDerivative1[] yVal = createLinearGrid(yMin, yMax, ny, 2.0);

        final double zMin = -2.0;
        final double zMax = 6.0;
        final int    nz   = 19;
        final UnivariateDerivative1[] zVal = createLinearGrid(zMin, zMax, nz, 3.0);

        CalculusFieldTrivariateFunction<UnivariateDerivative1> f = (x, y, z) ->
                x.multiply(3).subtract(2)
                        .multiply(y.multiply(-0.5).add(6))
                        .multiply(z.multiply(-2).add(1));
        FieldTrilinearInterpolatingFunction<UnivariateDerivative1> bif = createInterpolatingFunction(xVal, yVal, zVal, f);

        assertEquals(xMin, bif.getXInf().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXInf().getFirstDerivative(), 1.0e-15);
        assertEquals(xMax, bif.getXSup().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXSup().getFirstDerivative(), 1.0e-15);
        assertEquals(yMin, bif.getYInf().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYInf().getFirstDerivative(), 1.0e-15);
        assertEquals(yMax, bif.getYSup().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYSup().getFirstDerivative(), 1.0e-15);
        assertEquals(zMin, bif.getZInf().getValue(),           1.0e-15);
        assertEquals(3.0,  bif.getZInf().getFirstDerivative(), 1.0e-15);
        assertEquals(zMax, bif.getZSup().getValue(),           1.0e-15);
        assertEquals(3.0,  bif.getZSup().getFirstDerivative(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, zVal, bif, f, 1.0e-15, 1.0e-14);
        checkInterpolationRandom(new SobolSequenceGenerator(3), xMin, xMax, yMin, yMax, zMin, zMax,
                bif, f, 2.0e-13, 2.0e-12);

    }

    @Test
    void testSinCos() {
        doTestSinCos(  10,   10, 10, 1.2e-2, 0.1);
        doTestSinCos( 100,  100, 100, 1.2e-4, 8e-3);
    }

    private void doTestSinCos(final int nx, final int ny, final int nz, final double tol0, final double tol1) {
        final double xMin = -1.0;
        final double xMax =  2.0;
        final UnivariateDerivative1[] xVal = createLinearGrid(xMin, xMax, nx, 1.0);

        final double yMin = 0.0;
        final double yMax = 1.5;
        final UnivariateDerivative1[] yVal = createLinearGrid(yMin, yMax, ny, 2.0);

        final double zMin = -1.5;
        final double zMax = 1.0;
        final UnivariateDerivative1[] zVal = createLinearGrid(zMin, zMax, nz, 3.0);

        CalculusFieldTrivariateFunction<UnivariateDerivative1> f = (x, y, z) ->
                x.sin().multiply(y.cos()).multiply(z.half().sin());
        FieldTrilinearInterpolatingFunction<UnivariateDerivative1> bif = createInterpolatingFunction(xVal, yVal, zVal, f);

        assertEquals(xMin, bif.getXInf().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXInf().getFirstDerivative(), 1.0e-15);
        assertEquals(xMax, bif.getXSup().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXSup().getFirstDerivative(), 1.0e-15);
        assertEquals(yMin, bif.getYInf().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYInf().getFirstDerivative(), 1.0e-15);
        assertEquals(yMax, bif.getYSup().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYSup().getFirstDerivative(), 1.0e-15);
        assertEquals(zMin, bif.getZInf().getValue(),           1.0e-15);
        assertEquals(3.0,  bif.getZInf().getFirstDerivative(), 1.0e-15);
        assertEquals(zMax, bif.getZSup().getValue(),           1.0e-15);
        assertEquals(3.0,  bif.getZSup().getFirstDerivative(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, zVal, bif, f, 1.0e-15, 1.0e-15);
        checkInterpolationRandom(new SobolSequenceGenerator(3), xMin, xMax, yMin, yMax, zMin, zMax,
                bif, f, tol0, tol1);

    }

    private UnivariateDerivative1[] createLinearGrid(final double min, final double max, final int n,
                                                     final double derivative) {
        final UnivariateDerivative1[] grid = new UnivariateDerivative1[n];
        for (int i = 0; i < n; ++i) {
            grid[i] = new UnivariateDerivative1(((n - 1 - i) * min + i * max) / (n - 1), derivative);
        }
        return grid;
    }

    private FieldTrilinearInterpolatingFunction<UnivariateDerivative1>
    createInterpolatingFunction(UnivariateDerivative1[] xVal,
                                UnivariateDerivative1[] yVal,
                                UnivariateDerivative1[] zVal,
                                CalculusFieldTrivariateFunction<UnivariateDerivative1> f) {
        final UnivariateDerivative1[][][] fVal = new UnivariateDerivative1[xVal.length][yVal.length][zVal.length];
        for (int i = 0; i < xVal.length; ++i) {
            for (int j = 0; j < yVal.length; ++j) {
                for (int k = 0; k < zVal.length; ++k) {
                    fVal[i][j][k] = f.value(xVal[i], yVal[j], zVal[k]);
                }
            }
        }
        return new FieldTrilinearInterpolator<UnivariateDerivative1>().interpolate(xVal, yVal, zVal, fVal);
    }

    private void checkInterpolationAtNodes(final UnivariateDerivative1[] xVal,
                                           final UnivariateDerivative1[] yVal,
                                           final UnivariateDerivative1[] zVal,
                                           final FieldTrilinearInterpolatingFunction<UnivariateDerivative1> bif,
                                           final CalculusFieldTrivariateFunction<UnivariateDerivative1> f,
                                           final double tol0, final double tol1) {
        for (final UnivariateDerivative1 x : xVal) {
            for (final UnivariateDerivative1 y : yVal) {
                for (final UnivariateDerivative1 z : zVal) {
                    assertEquals(f.value(x, y, z).getValue(), bif.value(x, y, z).getValue(), tol0);
                    assertEquals(f.value(x, y, z).getFirstDerivative(), bif.value(x, y, z).getFirstDerivative(), tol1);
                }
            }
        }
    }

    private void checkInterpolationRandom(final RandomVectorGenerator random,
                                          final double xMin, final double xMax,
                                          final double yMin, final double yMax,
                                          final double zMin, final double zMax,
                                          final FieldTrilinearInterpolatingFunction<UnivariateDerivative1> bif,
                                          final CalculusFieldTrivariateFunction<UnivariateDerivative1> f,
                                          final double tol0, final double tol1) {
        double maxError0 = 0.0;
        double maxError1 = 0.0;
        for (int i = 0; i < 10000; ++i) {

            final double[] v = random.nextVector();

            final UnivariateDerivative1 x = new UnivariateDerivative1(xMin + v[0] * (xMax - xMin), 1.0);
            final UnivariateDerivative1 y = new UnivariateDerivative1(yMin + v[1] * (yMax - yMin), 1.0);
            final UnivariateDerivative1 z = new UnivariateDerivative1(zMin + v[2] * (zMax - zMin), 1.0);
            final UnivariateDerivative1 delta = f.value(x, y, z).subtract(bif.value(x, y, z));
            maxError0 = FastMath.max(maxError0, FastMath.abs(delta.getValue()));
            maxError1 = FastMath.max(maxError1, FastMath.abs(delta.getFirstDerivative()));

        }

        assertEquals(0.0, maxError0, tol0);
        assertEquals(0.0, maxError1, tol1);

    }

}
