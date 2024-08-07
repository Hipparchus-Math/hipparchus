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

package org.hipparchus.analysis.differentiation;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.analysis.QuinticFunction;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.UnivariateMatrixFunction;
import org.hipparchus.analysis.UnivariateVectorFunction;
import org.hipparchus.analysis.function.Gaussian;
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test for class {@link FiniteDifferencesDifferentiator}.
 */
class FiniteDifferencesDifferentiatorTest {

    @Test
    void testWrongNumberOfPoints() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new FiniteDifferencesDifferentiator(1, 1.0);
        });
    }

    @Test
    void testWrongStepSize() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new FiniteDifferencesDifferentiator(3, 0.0);
        });
    }

    @Test
    void testSerialization() {
        FiniteDifferencesDifferentiator differentiator =
                new FiniteDifferencesDifferentiator(3, 1.0e-3);
        FiniteDifferencesDifferentiator recovered =
                (FiniteDifferencesDifferentiator) UnitTestUtils.serializeAndRecover(differentiator);
        assertEquals(differentiator.getNbPoints(), recovered.getNbPoints());
        assertEquals(differentiator.getStepSize(), recovered.getStepSize(), 1.0e-15);
    }

    @Test
    void testConstant() {
        FiniteDifferencesDifferentiator differentiator =
                new FiniteDifferencesDifferentiator(5, 0.01);
        UnivariateDifferentiableFunction f =
                differentiator.differentiate(new UnivariateFunction() {
                    @Override
                    public double value(double x) {
                        return 42.0;
                    }
                });
        DSFactory factory = new DSFactory(1, 2);
        for (double x = -10; x < 10; x += 0.1) {
            DerivativeStructure y = f.value(factory.variable(0, x));
            assertEquals(42.0, y.getValue(), 1.0e-15);
            assertEquals( 0.0, y.getPartialDerivative(1), 1.0e-15);
            assertEquals( 0.0, y.getPartialDerivative(2), 1.0e-15);
        }
    }

    @Test
    void testLinear() {
        FiniteDifferencesDifferentiator differentiator =
                new FiniteDifferencesDifferentiator(5, 0.01);
        UnivariateDifferentiableFunction f =
                differentiator.differentiate(new UnivariateFunction() {
                    @Override
                    public double value(double x) {
                        return 2 - 3 * x;
                    }
                });
        DSFactory factory = new DSFactory(1, 2);
        for (double x = -10; x < 10; x += 0.1) {
            DerivativeStructure y = f.value(factory.variable(0, x));
            assertEquals(2 - 3 * x, y.getValue(), 2.0e-15, "" + (2 - 3 * x - y.getValue()));
            assertEquals(-3.0, y.getPartialDerivative(1), 4.0e-13);
            assertEquals( 0.0, y.getPartialDerivative(2), 9.0e-11);
        }
    }

    @Test
    void testGaussian() {
        FiniteDifferencesDifferentiator differentiator =
                new FiniteDifferencesDifferentiator(9, 0.02);
        UnivariateDifferentiableFunction gaussian = new Gaussian(1.0, 2.0);
        UnivariateDifferentiableFunction f =
                differentiator.differentiate(gaussian);
        double[] expectedError = new double[] {
            6.939e-18, 1.284e-15, 2.477e-13, 1.168e-11, 2.840e-9, 7.971e-8
        };
        double[] maxError = new double[expectedError.length];
        DSFactory factory = new DSFactory(1, maxError.length - 1);
        for (double x = -10; x < 10; x += 0.1) {
            DerivativeStructure dsX  = factory.variable(0, x);
            DerivativeStructure yRef = gaussian.value(dsX);
            DerivativeStructure y    = f.value(dsX);
            assertEquals(f.value(dsX.getValue()), f.value(dsX).getValue(), 1.0e-15);
            for (int order = 0; order <= yRef.getOrder(); ++order) {
                maxError[order] = FastMath.max(maxError[order],
                                        FastMath.abs(yRef.getPartialDerivative(order) -
                                                     y.getPartialDerivative(order)));
            }
        }
        for (int i = 0; i < maxError.length; ++i) {
            assertEquals(expectedError[i], maxError[i], 0.01 * expectedError[i]);
        }
    }

    @Test
    void testStepSizeUnstability() {
        UnivariateDifferentiableFunction quintic = new QuinticFunction();
        UnivariateDifferentiableFunction goodStep =
                new FiniteDifferencesDifferentiator(7, 0.25).differentiate(quintic);
        UnivariateDifferentiableFunction badStep =
                new FiniteDifferencesDifferentiator(7, 1.0e-6).differentiate(quintic);
        double[] maxErrorGood = new double[7];
        double[] maxErrorBad  = new double[7];
        DSFactory factory = new DSFactory(1, maxErrorGood.length - 1);
        for (double x = -10; x < 10; x += 0.1) {
            DerivativeStructure dsX  = factory.variable(0, x);
            DerivativeStructure yRef  = quintic.value(dsX);
            DerivativeStructure yGood = goodStep.value(dsX);
            DerivativeStructure yBad  = badStep.value(dsX);
            for (int order = 0; order <= 6; ++order) {
                maxErrorGood[order] = FastMath.max(maxErrorGood[order],
                                                   FastMath.abs(yRef.getPartialDerivative(order) -
                                                                yGood.getPartialDerivative(order)));
                maxErrorBad[order]  = FastMath.max(maxErrorBad[order],
                                                   FastMath.abs(yRef.getPartialDerivative(order) -
                                                                yBad.getPartialDerivative(order)));
            }
        }

        // the 0.25 step size is good for finite differences in the quintic on this abscissa range for 7 points
        // the errors are fair
        final double[] expectedGood = new double[] {
            7.276e-12, 7.276e-11, 9.968e-10, 3.092e-9, 5.432e-8, 8.196e-8, 1.818e-6
        };

        // the 1.0e-6 step size is far too small for finite differences in the quintic on this abscissa range for 7 points
        // the errors are huge!
        final double[] expectedBad = new double[] {
            2.910e-11, 2.087e-5, 147.7, 3.820e7, 6.354e14, 6.548e19, 1.543e27
        };

        for (int i = 0; i < maxErrorGood.length; ++i) {
            assertEquals(expectedGood[i], maxErrorGood[i], 0.01 * expectedGood[i]);
            assertEquals(expectedBad[i],  maxErrorBad[i],  0.01 * expectedBad[i]);
        }

    }

    @Test
    void testWrongOrder() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            UnivariateDifferentiableFunction f =
                new FiniteDifferencesDifferentiator(3, 0.01).differentiate(new UnivariateFunction() {
                    @Override
                    public double value(double x) {
                        // this exception should not be thrown because wrong order
                        // should be detected before function call
                        throw MathRuntimeException.createInternalError();
                    }
                });
            f.value(new DSFactory(1, 3).variable(0, 1.0));
        });
    }

    @Test
    void testWrongOrderVector() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            UnivariateDifferentiableVectorFunction f =
                new FiniteDifferencesDifferentiator(3, 0.01).differentiate(new UnivariateVectorFunction() {
                    @Override
                    public double[] value(double x) {
                        // this exception should not be thrown because wrong order
                        // should be detected before function call
                        throw MathRuntimeException.createInternalError();
                    }
                });
            f.value(new DSFactory(1, 3).variable(0, 1.0));
        });
    }

    @Test
    void testWrongOrderMatrix() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            UnivariateDifferentiableMatrixFunction f =
                new FiniteDifferencesDifferentiator(3, 0.01).differentiate(new UnivariateMatrixFunction() {
                    @Override
                    public double[][] value(double x) {
                        // this exception should not be thrown because wrong order
                        // should be detected before function call
                        throw MathRuntimeException.createInternalError();
                    }
                });
            f.value(new DSFactory(1, 3).variable(0, 1.0));
        });
    }

    @Test
    void testTooLargeStep() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new FiniteDifferencesDifferentiator(3, 2.5, 0.0, 1.0);
        });
    }

    @Test
    void testBounds() {

        final double slope = 2.5;
        UnivariateFunction f = new UnivariateFunction() {
            @Override
            public double value(double x) {
                if (x < 0) {
                    throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL,
                                                           x, 0);
                } else if (x > 1) {
                    throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_LARGE,
                                                           x, 1);
                } else {
                    return slope * x;
                }
            }
        };

        UnivariateDifferentiableFunction missingBounds =
                new FiniteDifferencesDifferentiator(3, 0.1).differentiate(f);
        UnivariateDifferentiableFunction properlyBounded =
                new FiniteDifferencesDifferentiator(3, 0.1, 0.0, 1.0).differentiate(f);
        DSFactory factory = new DSFactory(1, 1);
        DerivativeStructure tLow  = factory.variable(0, 0.05);
        DerivativeStructure tHigh = factory.variable(0, 0.95);

        try {
            // here, we did not set the bounds, so the differences are evaluated out of domain
            // using f(-0.05), f(0.05), f(0.15)
            missingBounds.value(tLow);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException nse) {
            assertEquals(LocalizedCoreFormats.NUMBER_TOO_SMALL, nse.getSpecifier());
            assertEquals(-0.05, ((Double) nse.getParts()[0]).doubleValue(), 1.0e-10);
        } catch (Exception e) {
            fail("wrong exception caught: " + e.getClass().getName());
        }

        try {
            // here, we did not set the bounds, so the differences are evaluated out of domain
            // using f(0.85), f(0.95), f(1.05)
            missingBounds.value(tHigh);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException nle) {
            assertEquals(LocalizedCoreFormats.NUMBER_TOO_LARGE, nle.getSpecifier());
            assertEquals(1.05, ((Double) nle.getParts()[0]).doubleValue(), 1.0e-10);
        } catch (Exception e) {
            fail("wrong exception caught: " + e.getClass().getName());
        }

        // here, we did set the bounds, so evaluations are done within domain
        // using f(0.0), f(0.1), f(0.2)
        assertEquals(slope, properlyBounded.value(tLow).getPartialDerivative(1), 1.0e-10);

        // here, we did set the bounds, so evaluations are done within domain
        // using f(0.8), f(0.9), f(1.0)
        assertEquals(slope, properlyBounded.value(tHigh).getPartialDerivative(1), 1.0e-10);

    }

    @Test
    void testBoundedSqrt() {

        UnivariateFunctionDifferentiator differentiator =
                new FiniteDifferencesDifferentiator(9, 1.0 / 32, 0.0, Double.POSITIVE_INFINITY);
        UnivariateDifferentiableFunction sqrt = differentiator.differentiate(new UnivariateFunction() {
            @Override
            public double value(double x) {
                return FastMath.sqrt(x);
            }
        });

        // we are able to compute derivative near 0, but the accuracy is much poorer there
        DSFactory factory = new DSFactory(1, 1);
        DerivativeStructure t001 = factory.variable(0, 0.01);
        assertEquals(0.5 / FastMath.sqrt(t001.getValue()), sqrt.value(t001).getPartialDerivative(1), 1.6);
        DerivativeStructure t01 = factory.variable(0, 0.1);
        assertEquals(0.5 / FastMath.sqrt(t01.getValue()), sqrt.value(t01).getPartialDerivative(1), 7.0e-3);
        DerivativeStructure t03 = factory.variable(0, 0.3);
        assertEquals(0.5 / FastMath.sqrt(t03.getValue()), sqrt.value(t03).getPartialDerivative(1), 2.1e-7);

    }

    @Test
    void testVectorFunction() {

        FiniteDifferencesDifferentiator differentiator =
                new FiniteDifferencesDifferentiator(7, 0.01);
        UnivariateDifferentiableVectorFunction f =
                differentiator.differentiate(new UnivariateVectorFunction() {

            @Override
            public double[] value(double x) {
                return new double[] { FastMath.cos(x), FastMath.sin(x) };
            }

        });

        DSFactory factory = new DSFactory(1, 2);
        for (double x = -10; x < 10; x += 0.1) {
            DerivativeStructure dsX = factory.variable(0, x);
            DerivativeStructure[] y = f.value(dsX);
            double cos = FastMath.cos(x);
            double sin = FastMath.sin(x);
            double[] f1 = f.value(dsX.getValue());
            DerivativeStructure[] f2 = f.value(dsX);
            assertEquals(f1.length, f2.length);
            for (int i = 0; i < f1.length; ++i) {
                assertEquals(f1[i], f2[i].getValue(), 1.0e-15);
            }
            assertEquals( cos, y[0].getValue(), 7.0e-16);
            assertEquals( sin, y[1].getValue(), 7.0e-16);
            assertEquals(-sin, y[0].getPartialDerivative(1), 6.0e-14);
            assertEquals( cos, y[1].getPartialDerivative(1), 6.0e-14);
            assertEquals(-cos, y[0].getPartialDerivative(2), 2.0e-11);
            assertEquals(-sin, y[1].getPartialDerivative(2), 2.0e-11);
        }

    }

    @Test
    void testMatrixFunction() {

        FiniteDifferencesDifferentiator differentiator =
                new FiniteDifferencesDifferentiator(7, 0.01);
        UnivariateDifferentiableMatrixFunction f =
                differentiator.differentiate(new UnivariateMatrixFunction() {

            @Override
            public double[][] value(double x) {
                return new double[][] {
                    { FastMath.cos(x),  FastMath.sin(x)  },
                    { FastMath.cosh(x), FastMath.sinh(x) }
                };
            }

        });

        DSFactory factory = new DSFactory(1, 2);
        for (double x = -1; x < 1; x += 0.02) {
            DerivativeStructure dsX = factory.variable(0, x);
            DerivativeStructure[][] y = f.value(dsX);
            double cos = FastMath.cos(x);
            double sin = FastMath.sin(x);
            double cosh = FastMath.cosh(x);
            double sinh = FastMath.sinh(x);
            double[][] f1 = f.value(dsX.getValue());
            DerivativeStructure[][] f2 = f.value(dsX);
            assertEquals(f1.length, f2.length);
            for (int i = 0; i < f1.length; ++i) {
                assertEquals(f1[i].length, f2[i].length);
                for (int j = 0; j < f1[i].length; ++j) {
                    assertEquals(f1[i][j], f2[i][j].getValue(), 1.0e-15);
                }
            }
            assertEquals(cos,   y[0][0].getValue(), 7.0e-18);
            assertEquals(sin,   y[0][1].getValue(), 6.0e-17);
            assertEquals(cosh,  y[1][0].getValue(), 3.0e-16);
            assertEquals(sinh,  y[1][1].getValue(), 3.0e-16);
            assertEquals(-sin,  y[0][0].getPartialDerivative(1), 2.0e-14);
            assertEquals( cos,  y[0][1].getPartialDerivative(1), 2.0e-14);
            assertEquals( sinh, y[1][0].getPartialDerivative(1), 3.0e-14);
            assertEquals( cosh, y[1][1].getPartialDerivative(1), 3.0e-14);
            assertEquals(-cos,  y[0][0].getPartialDerivative(2), 3.0e-12);
            assertEquals(-sin,  y[0][1].getPartialDerivative(2), 3.0e-12);
            assertEquals( cosh, y[1][0].getPartialDerivative(2), 6.0e-12);
            assertEquals( sinh, y[1][1].getPartialDerivative(2), 6.0e-12);
        }

    }

    @Test
    void testSeveralFreeParameters() {
        FiniteDifferencesDifferentiator differentiator =
                new FiniteDifferencesDifferentiator(5, 0.001);
        UnivariateDifferentiableFunction sine = new Sin();
        UnivariateDifferentiableFunction f =
                differentiator.differentiate(sine);
        double[] expectedError = new double[] {
            6.696e-16, 1.371e-12, 2.007e-8, 1.754e-5
        };
        double[] maxError = new double[expectedError.length];
        DSFactory factory = new DSFactory(2, maxError.length - 1);
        for (double x = -2; x < 2; x += 0.1) {
           for (double y = -2; y < 2; y += 0.1) {
               DerivativeStructure dsX  = factory.variable(0, x);
               DerivativeStructure dsY  = factory.variable(1, y);
               DerivativeStructure dsT  = dsX.multiply(3).subtract(dsY.multiply(2));
               DerivativeStructure sRef = sine.value(dsT);
               DerivativeStructure s    = f.value(dsT);
               for (int xOrder = 0; xOrder <= sRef.getOrder(); ++xOrder) {
                   for (int yOrder = 0; yOrder <= sRef.getOrder(); ++yOrder) {
                       if (xOrder + yOrder <= sRef.getOrder()) {
                           maxError[xOrder +yOrder] = FastMath.max(maxError[xOrder + yOrder],
                                                                    FastMath.abs(sRef.getPartialDerivative(xOrder, yOrder) -
                                                                                 s.getPartialDerivative(xOrder, yOrder)));
                       }
                   }
               }
           }
       }
       for (int i = 0; i < maxError.length; ++i) {
           assertEquals(expectedError[i], maxError[i], 0.01 * expectedError[i]);
       }
    }

}
