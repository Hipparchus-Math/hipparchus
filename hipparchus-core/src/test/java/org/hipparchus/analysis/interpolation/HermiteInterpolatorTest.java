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

import java.util.Random;

import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.analysis.differentiation.Gradient;
import org.hipparchus.analysis.polynomials.PolynomialFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class HermiteInterpolatorTest {

    @Test
    public void testZero() {
        HermiteInterpolator interpolator = new HermiteInterpolator();
        interpolator.addSamplePoint(0.0, new double[] { 0.0 });
        DSFactory factory = new DSFactory(1, 1);
        for (double x = -10; x < 10; x += 1.0) {
            DerivativeStructure y = interpolator.value(factory.variable(0, x))[0];
            Assert.assertEquals(0.0, y.getValue(), 1.0e-15);
            Assert.assertEquals(0.0, y.getPartialDerivative(1), 1.0e-15);
        }
        checkPolynomial(new PolynomialFunction(new double[] { 0.0 }),
                        interpolator.getPolynomials()[0]);
    }

    @Test
    public void testQuadratic() {
        HermiteInterpolator interpolator = new HermiteInterpolator();
        interpolator.addSamplePoint(0.0, new double[] { 2.0 });
        interpolator.addSamplePoint(1.0, new double[] { 0.0 });
        interpolator.addSamplePoint(2.0, new double[] { 0.0 });
        DSFactory factory = new DSFactory(1, 1);
        for (double x = -10; x < 10; x += 1.0) {
            DerivativeStructure y = interpolator.value(factory.variable(0, x))[0];
            Assert.assertEquals((x - 1.0) * (x - 2.0), y.getValue(), 1.0e-15);
            Assert.assertEquals(2 * x - 3.0, y.getPartialDerivative(1), 1.0e-15);
        }
        checkPolynomial(new PolynomialFunction(new double[] { 2.0, -3.0, 1.0 }),
                        interpolator.getPolynomials()[0]);
    }

    @Test
    public void testMixedDerivatives() {
        HermiteInterpolator interpolator = new HermiteInterpolator();
        interpolator.addSamplePoint(0.0, new double[] { 1.0 }, new double[] { 2.0 });
        interpolator.addSamplePoint(1.0, new double[] { 4.0 });
        interpolator.addSamplePoint(2.0, new double[] { 5.0 }, new double[] { 2.0 });
        DSFactory factory = new DSFactory(1, 1);
        Assert.assertEquals(4, interpolator.getPolynomials()[0].degree());
        DerivativeStructure y0 = interpolator.value(factory.variable(0, 0.0))[0];
        Assert.assertEquals(1.0, y0.getValue(), 1.0e-15);
        Assert.assertEquals(2.0, y0.getPartialDerivative(1), 1.0e-15);
        double[][] d0 = interpolator.derivatives(0.0, 1);
        Assert.assertEquals(1.0, d0[0][0], 1.0e-15);
        Assert.assertEquals(2.0, d0[1][0], 1.0e-15);
        Assert.assertEquals(4.0, interpolator.value(1.0)[0], 1.0e-15);
        DerivativeStructure y2 = interpolator.value(factory.variable(0, 2.0))[0];
        Assert.assertEquals(5.0, y2.getValue(), 1.0e-15);
        Assert.assertEquals(2.0, y2.getPartialDerivative(1), 1.0e-15);
        double[][] d2 = interpolator.derivatives(2.0, 1);
        Assert.assertEquals(5.0, d2[0][0], 1.0e-15);
        Assert.assertEquals(2.0, d2[1][0], 1.0e-15);
        checkPolynomial(new PolynomialFunction(new double[] { 1.0, 2.0, 4.0, -4.0, 1.0 }),
                        interpolator.getPolynomials()[0]);
    }

    @Test
    public void testRandomPolynomialsValuesOnly() {

        Random random = new Random(0x42b1e7dbd361a932l);

        for (int i = 0; i < 100; ++i) {

            int maxDegree = 0;
            PolynomialFunction[] p = new PolynomialFunction[5];
            for (int k = 0; k < p.length; ++k) {
                int degree = random.nextInt(7);
                p[k] = randomPolynomial(degree, random);
                maxDegree = FastMath.max(maxDegree, degree);
            }

            HermiteInterpolator interpolator = new HermiteInterpolator();
            for (int j = 0; j < 1 + maxDegree; ++j) {
                double x = 0.1 * j;
                double[] values = new double[p.length];
                for (int k = 0; k < p.length; ++k) {
                    values[k] = p[k].value(x);
                }
                interpolator.addSamplePoint(x, values);
            }

            for (double x = 0; x < 2; x += 0.1) {
                double[] values = interpolator.value(x);
                Assert.assertEquals(p.length, values.length);
                for (int k = 0; k < p.length; ++k) {
                    Assert.assertEquals(p[k].value(x), values[k], 1.0e-8 * FastMath.abs(p[k].value(x)));
                }
            }

            PolynomialFunction[] result = interpolator.getPolynomials();
            for (int k = 0; k < p.length; ++k) {
                checkPolynomial(p[k], result[k]);
            }

        }
    }

    @Test
    public void testRandomPolynomialsFirstDerivative() {

        Random random = new Random(0x570803c982ca5d3bl);

        for (int i = 0; i < 100; ++i) {

            int maxDegree = 0;
            PolynomialFunction[] p      = new PolynomialFunction[5];
            PolynomialFunction[] pPrime = new PolynomialFunction[5];
            for (int k = 0; k < p.length; ++k) {
                int degree = random.nextInt(7);
                p[k]      = randomPolynomial(degree, random);
                pPrime[k] = p[k].polynomialDerivative();
                maxDegree = FastMath.max(maxDegree, degree);
            }

            HermiteInterpolator interpolator = new HermiteInterpolator();
            for (int j = 0; j < 1 + maxDegree / 2; ++j) {
                double x = 0.1 * j;
                double[] values      = new double[p.length];
                double[] derivatives = new double[p.length];
                for (int k = 0; k < p.length; ++k) {
                    values[k]      = p[k].value(x);
                    derivatives[k] = pPrime[k].value(x);
                }
                interpolator.addSamplePoint(x, values, derivatives);
            }

            DSFactory factory = new DSFactory(1, 1);
            for (double x = 0; x < 2; x += 0.1) {
                DerivativeStructure[] y = interpolator.value(factory.variable(0, x));
                Assert.assertEquals(p.length, y.length);
                for (int k = 0; k < p.length; ++k) {
                    Assert.assertEquals(p[k].value(x), y[k].getValue(), 1.0e-8 * FastMath.abs(p[k].value(x)));
                    Assert.assertEquals(pPrime[k].value(x), y[k].getPartialDerivative(1), 4.0e-8 * FastMath.abs(p[k].value(x)));
                }
            }

            PolynomialFunction[] result = interpolator.getPolynomials();
            for (int k = 0; k < p.length; ++k) {
                checkPolynomial(p[k], result[k]);
            }

        }
    }

    @Test
    public void testSine() {
        HermiteInterpolator interpolator = new HermiteInterpolator();
        for (double x = 0; x < FastMath.PI; x += 0.5) {
            interpolator.addSamplePoint(x, new double[] { FastMath.sin(x) });
        }
        DSFactory factory = new DSFactory(1, 2);
        for (double x = 0.1; x <= 2.9; x += 0.01) {
            DerivativeStructure y = interpolator.value(factory.variable(0, x))[0];
            Assert.assertEquals( FastMath.sin(x), y.getValue(), 3.5e-5);
            Assert.assertEquals( FastMath.cos(x), y.getPartialDerivative(1), 1.3e-4);
            Assert.assertEquals(-FastMath.sin(x), y.getPartialDerivative(2), 2.9e-3);
        }
    }

    @Test
    public void testSquareRoot() {
        HermiteInterpolator interpolator = new HermiteInterpolator();
        for (double x = 1.0; x < 3.6; x += 0.5) {
            interpolator.addSamplePoint(x, new double[] { FastMath.sqrt(x) });
        }
        DSFactory factory = new DSFactory(1, 1);
        for (double x = 1.1; x < 3.5; x += 0.01) {
            DerivativeStructure y = interpolator.value(factory.variable(0, x))[0];
            Assert.assertEquals(FastMath.sqrt(x), y.getValue(), 1.5e-4);
            Assert.assertEquals(0.5 / FastMath.sqrt(x), y.getPartialDerivative(1), 8.5e-4);
        }
    }

    @Test
    public void testWikipedia() {
        // this test corresponds to the example from Wikipedia page:
        // http://en.wikipedia.org/wiki/Hermite_interpolation
        HermiteInterpolator interpolator = new HermiteInterpolator();
        interpolator.addSamplePoint(-1, new double[] { 2 }, new double[] { -8 }, new double[] { 56 });
        interpolator.addSamplePoint( 0, new double[] { 1 }, new double[] {  0 }, new double[] {  0 });
        interpolator.addSamplePoint( 1, new double[] { 2 }, new double[] {  8 }, new double[] { 56 });
        DSFactory factory = new DSFactory(1, 1);
        for (double x = -1.0; x <= 1.0; x += 0.125) {
            DerivativeStructure y = interpolator.value(factory.variable(0, x))[0];
            double x2 = x * x;
            double x4 = x2 * x2;
            double x8 = x4 * x4;
            Assert.assertEquals(x8 + 1, y.getValue(), 1.0e-15);
            Assert.assertEquals(8 * x4 * x2 * x, y.getPartialDerivative(1), 1.0e-15);
        }
        checkPolynomial(new PolynomialFunction(new double[] { 1, 0, 0, 0, 0, 0, 0, 0, 1 }),
                        interpolator.getPolynomials()[0]);
    }

    @Test
    public void testWikipediaGradient() {
        // this test corresponds to the example from Wikipedia page:
        // http://en.wikipedia.org/wiki/Hermite_interpolation
        HermiteInterpolator interpolator = new HermiteInterpolator();
        interpolator.addSamplePoint(-1, new double[] { 2 }, new double[] { -8 }, new double[] { 56 });
        interpolator.addSamplePoint( 0, new double[] { 1 }, new double[] {  0 }, new double[] {  0 });
        interpolator.addSamplePoint( 1, new double[] { 2 }, new double[] {  8 }, new double[] { 56 });
        for (double x = -1.0; x <= 1.0; x += 0.125) {
            Gradient y = interpolator.value(Gradient.variable(1, 0, x))[0];
            double x2 = x * x;
            double x4 = x2 * x2;
            double x8 = x4 * x4;
            Assert.assertEquals(x8 + 1, y.getValue(), 1.0e-15);
            Assert.assertEquals(8 * x4 * x2 * x, y.getPartialDerivative(0), 1.0e-15);
        }
        checkPolynomial(new PolynomialFunction(new double[] { 1, 0, 0, 0, 0, 0, 0, 0, 1 }),
                        interpolator.getPolynomials()[0]);
    }

    @Test
    public void testOnePointParabola() {
        HermiteInterpolator interpolator = new HermiteInterpolator();
        interpolator.addSamplePoint(0, new double[] { 1 }, new double[] { 1 }, new double[] { 2 });
        DSFactory factory = new DSFactory(1, 1);
        for (double x = -1.0; x <= 1.0; x += 0.125) {
            DerivativeStructure y = interpolator.value(factory.variable(0, x))[0];
            Assert.assertEquals(1 + x * (1 + x), y.getValue(), 1.0e-15);
            Assert.assertEquals(1 + 2 * x, y.getPartialDerivative(1), 1.0e-15);
        }
        checkPolynomial(new PolynomialFunction(new double[] { 1, 1, 1 }),
                        interpolator.getPolynomials()[0]);
    }

    private PolynomialFunction randomPolynomial(int degree, Random random) {
        double[] coeff = new double[ 1 + degree];
        for (int j = 0; j < degree; ++j) {
            coeff[j] = random.nextDouble();
        }
        return new PolynomialFunction(coeff);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testEmptySample() {
        new HermiteInterpolator().value(0.0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testDuplicatedAbscissa() {
        HermiteInterpolator interpolator = new HermiteInterpolator();
        interpolator.addSamplePoint(1.0, new double[] { 0.0 });
        interpolator.addSamplePoint(1.0, new double[] { 1.0 });
    }

    private void checkPolynomial(PolynomialFunction expected, PolynomialFunction result) {
        Assert.assertTrue(result.degree() >= expected.degree());
        double[] cE = expected.getCoefficients();
        double[] cR = result.getCoefficients();
        for (int i = 0; i < cE.length; ++i) {
            Assert.assertEquals(cE[i], cR[i], 1.0e-8 * FastMath.abs(cE[i]));
        }
        for (int i = cE.length; i < cR.length; ++i) {
            Assert.assertEquals(0.0, cR[i], 1.0e-9);
        }
    }

}

