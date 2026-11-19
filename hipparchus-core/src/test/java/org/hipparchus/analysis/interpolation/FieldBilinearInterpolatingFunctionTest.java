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

import org.hipparchus.analysis.CalculusFieldBivariateFunction;
import org.hipparchus.analysis.differentiation.UnivariateDerivative1;
import org.hipparchus.random.RandomVectorGenerator;
import org.hipparchus.random.SobolSequenceGenerator;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test case for the field bilinear function.
 */
final class FieldBilinearInterpolatingFunctionTest {

    @Test
    void testConstant() {

        double xMin = 0.0;
        double xMax = 7.0;
        int    nx   = 15;
        UnivariateDerivative1[] xVal = createLinearGrid(xMin, xMax, nx, 1.0);

        double yMin = -5.0;
        double yMax =  5.0;
        int    ny   = 11;
        UnivariateDerivative1[] yVal = createLinearGrid(yMin, yMax, ny, 2.0);

        CalculusFieldBivariateFunction<UnivariateDerivative1>     f  = (x, y) -> new UnivariateDerivative1(3.5, 1.0);
        FieldBilinearInterpolatingFunction<UnivariateDerivative1> bif = createInterpolatingFunction(xVal, yVal, f);

        assertEquals(xMin, bif.getXInf().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXInf().getFirstDerivative(), 1.0e-15);
        assertEquals(xMax, bif.getXSup().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXSup().getFirstDerivative(), 1.0e-15);
        assertEquals(yMin, bif.getYInf().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYInf().getFirstDerivative(), 1.0e-15);
        assertEquals(yMax, bif.getYSup().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYSup().getFirstDerivative(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, bif, f, 1.0e-15, 1.0e-15);
        checkInterpolationRandom(new SobolSequenceGenerator(2), xMin, xMax, yMin, yMax, bif, f, 1.0e-15, 1.0e-15);

    }

    @Test
    void testLinear() {

        double xMin = -5.0;
        double xMax =  5.0;
        int    nx   = 11;
        UnivariateDerivative1[] xVal = createLinearGrid(xMin, xMax, nx, 1.0);

        double yMin = 0.0;
        double yMax = 7.0;
        int    ny   = 15;
        UnivariateDerivative1[] yVal = createLinearGrid(yMin, yMax, ny, 2.0);

        CalculusFieldBivariateFunction<UnivariateDerivative1> f = (x, y) -> x.twice().subtract(y);
        FieldBilinearInterpolatingFunction<UnivariateDerivative1> bif = createInterpolatingFunction(xVal, yVal, f);

        assertEquals(xMin, bif.getXInf().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXInf().getFirstDerivative(), 1.0e-15);
        assertEquals(xMax, bif.getXSup().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXSup().getFirstDerivative(), 1.0e-15);
        assertEquals(yMin, bif.getYInf().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYInf().getFirstDerivative(), 1.0e-15);
        assertEquals(yMax, bif.getYSup().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYSup().getFirstDerivative(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, bif, f, 1.0e-15, 1.0e-15);
        checkInterpolationRandom(new SobolSequenceGenerator(2), xMin, xMax, yMin, yMax, bif, f, 1.0e-15, 1.0e-15);

    }

    @Test
    void testQuadratic() {

        double xMin = -5.0;
        double xMax =  5.0;
        int    nx   = 11;
        UnivariateDerivative1[] xVal = createLinearGrid(xMin, xMax, nx, 1.0);

        double yMin = 0.0;
        double yMax = 7.0;
        int    ny   = 15;
        UnivariateDerivative1[] yVal = createLinearGrid(yMin, yMax, ny, 2.0);

        CalculusFieldBivariateFunction<UnivariateDerivative1> f = (x, y) -> x.multiply(3).subtract(2).multiply(y.multiply(-0.5).add(6));
        FieldBilinearInterpolatingFunction<UnivariateDerivative1> bif = createInterpolatingFunction(xVal, yVal, f);

        assertEquals(xMin, bif.getXInf().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXInf().getFirstDerivative(), 1.0e-15);
        assertEquals(xMax, bif.getXSup().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXSup().getFirstDerivative(), 1.0e-15);
        assertEquals(yMin, bif.getYInf().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYInf().getFirstDerivative(), 1.0e-15);
        assertEquals(yMax, bif.getYSup().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYSup().getFirstDerivative(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, bif, f, 1.0e-15, 1.0e-15);
        checkInterpolationRandom(new SobolSequenceGenerator(2), xMin, xMax, yMin, yMax, bif, f, 1.0e-15, 1.0e-15);

    }

    @Test
    void testSinCos() {
        doTestSinCos(  10,   10, 1.8e-2, 8.3e-2);
        doTestSinCos( 100,  100, 1.5e-4, 7.6e-3);
        doTestSinCos(1000, 1000, 1.4e-6, 7.6e-4);
    }

    private void doTestSinCos(final int nx, final int ny, final double tol0, final double tol1) {
        double xMin = -1.0;
        double xMax =  2.0;
        UnivariateDerivative1[] xVal = createLinearGrid(xMin, xMax, nx, 1.0);

        double yMin = 0.0;
        double yMax = 1.5;
        UnivariateDerivative1[] yVal = createLinearGrid(yMin, yMax, ny, 2.0);

        CalculusFieldBivariateFunction<UnivariateDerivative1> f = (x, y) -> FastMath.sin(x).multiply(FastMath.cos(y));
        FieldBilinearInterpolatingFunction<UnivariateDerivative1> bif = createInterpolatingFunction(xVal, yVal, f);

        assertEquals(xMin, bif.getXInf().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXInf().getFirstDerivative(), 1.0e-15);
        assertEquals(xMax, bif.getXSup().getValue(),           1.0e-15);
        assertEquals(1.0,  bif.getXSup().getFirstDerivative(), 1.0e-15);
        assertEquals(yMin, bif.getYInf().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYInf().getFirstDerivative(), 1.0e-15);
        assertEquals(yMax, bif.getYSup().getValue(),           1.0e-15);
        assertEquals(2.0,  bif.getYSup().getFirstDerivative(), 1.0e-15);

        checkInterpolationAtNodes(xVal, yVal, bif, f, 1.0e-15, 1.0e-15);
        checkInterpolationRandom(new SobolSequenceGenerator(2), xMin, xMax, yMin, yMax, bif, f, tol0, tol1);

    }

    private UnivariateDerivative1[] createLinearGrid(final double min, final double max, final int n,
                                                     final double derivative) {
        final UnivariateDerivative1[] grid = new UnivariateDerivative1[n];
        for (int i = 0; i < n; ++i) {
            grid[i] = new UnivariateDerivative1(((n - 1 - i) * min + i * max) / (n - 1), derivative);
        }
        return grid;
    }

    private FieldBilinearInterpolatingFunction<UnivariateDerivative1>
    createInterpolatingFunction(UnivariateDerivative1[] xVal, UnivariateDerivative1[] yVal,
                                CalculusFieldBivariateFunction<UnivariateDerivative1> f) {
        final UnivariateDerivative1[][] fVal = new UnivariateDerivative1[xVal.length][yVal.length];
        for (int i = 0; i < xVal.length; ++i) {
            for (int j = 0; j < yVal.length; ++j) {
                fVal[i][j] = f.value(xVal[i], yVal[j]);
            }
        }
        return new FieldBilinearInterpolator<UnivariateDerivative1>().interpolate(xVal, yVal, fVal);
    }

    private void checkInterpolationAtNodes(final UnivariateDerivative1[] xVal,
                                           final UnivariateDerivative1[] yVal,
                                           final FieldBilinearInterpolatingFunction<UnivariateDerivative1> bif,
                                           final CalculusFieldBivariateFunction<UnivariateDerivative1> f,
                                           final double tol0, final double tol1) {
        for (final UnivariateDerivative1 x : xVal) {
            for (final UnivariateDerivative1 y : yVal) {
                assertEquals(f.value(x, y).getValue(), bif.value(x, y).getValue(), tol0);
                assertEquals(f.value(x, y).getFirstDerivative(), bif.value(x, y).getFirstDerivative(), tol1);
            }
        }
    }

    private void checkInterpolationRandom(final RandomVectorGenerator random,
                                          final double xMin, final double xMax,
                                          final double yMin, final double yMax,
                                          final FieldBilinearInterpolatingFunction<UnivariateDerivative1> bif,
                                          final CalculusFieldBivariateFunction<UnivariateDerivative1> f,
                                          final double tol0, final double tol1) {
        double maxError0 = 0.0;
        double maxError1 = 0.0;
        for (int i = 0; i < 10000; ++i) {

            final double[] v = random.nextVector();

            final UnivariateDerivative1 x = new UnivariateDerivative1(xMin + v[0] * (xMax - xMin), 1.0);
            final UnivariateDerivative1 y = new UnivariateDerivative1(yMin + v[1] * (yMax - yMin), 1.0);
            final UnivariateDerivative1 delta = f.value(x, y).subtract(bif.value(x, y));
            maxError0 = FastMath.max(maxError0, FastMath.abs(delta.getValue()));
            maxError1 = FastMath.max(maxError1, FastMath.abs(delta.getFirstDerivative()));

        }

        assertEquals(0.0, maxError0, tol0);
        assertEquals(0.0, maxError1, tol1);

    }

}
