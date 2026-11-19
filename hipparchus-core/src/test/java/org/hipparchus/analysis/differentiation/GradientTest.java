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
package org.hipparchus.analysis.differentiation;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.CalculusFieldElementAbstractTest;
import org.hipparchus.Field;
import org.hipparchus.UnitTestUtils;
import org.hipparchus.analysis.FieldUnivariateFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.MathArrays;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test for class {@link UnivariateDerivative}.
 */
class GradientTest extends CalculusFieldElementAbstractTest<Gradient> {

    @Override
    protected Gradient build(final double x) {
        // the function is really a two variables function : f(x) = g(x, 0) with g(x, y) = x + y / 1024
        return new Gradient(x, 1.0, FastMath.scalb(1.0, -10));
    }

    @Test
    void testNorm() {
        final Gradient x = new Gradient(-1, 2, 3);
        assertEquals(6., x.norm());
    }

    @Test
    void testGetGradient() {
        Gradient g = new Gradient(-0.5, 2.5, 10.0, -1.0);
        assertEquals(-0.5, g.getReal(), 1.0e-15);
        assertEquals(-0.5, g.getValue(), 1.0e-15);
        assertEquals(+2.5, g.getGradient()[0], 1.0e-15);
        assertEquals(10.0, g.getGradient()[1], 1.0e-15);
        assertEquals(-1.0, g.getGradient()[2], 1.0e-15);
        assertEquals(+2.5, g.getPartialDerivative(0), 1.0e-15);
        assertEquals(10.0, g.getPartialDerivative(1), 1.0e-15);
        assertEquals(-1.0, g.getPartialDerivative(2), 1.0e-15);
        assertEquals(3, g.getFreeParameters());
        try {
            g.getPartialDerivative(-1);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, miae.getSpecifier());
        }
        try {
            g.getPartialDerivative(+3);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, miae.getSpecifier());
        }
    }

    @Test
    void testConstant() {
        Gradient g = Gradient.constant(5, -4.5);
        assertEquals(5, g.getFreeParameters());
        assertEquals(-4.5, g.getValue(), 1.0e-15);
        for (int i = 0 ; i < g.getFreeParameters(); ++i) {
            assertEquals(0.0, g.getPartialDerivative(i), 1.0e-15);
        }
    }

    @Test
    void testVariable() {
        Gradient g = Gradient.variable(5, 1, -4.5);
        assertEquals(5, g.getFreeParameters());
        assertEquals(-4.5, g.getValue(), 1.0e-15);
        for (int i = 0 ; i < g.getFreeParameters(); ++i) {
            assertEquals(i == 1 ? 1.0 : 0.0, g.getPartialDerivative(i), 1.0e-15);
        }
    }

    @Test
    void testStackVariable() {
        // GIVEN
        final Gradient gradient = new Gradient(1, 2, 3);
        // WHEN
        final Gradient gradientWithMoreVariable = gradient.stackVariable();
        // THEN
        Assertions.assertEquals(gradient.getValue(), gradientWithMoreVariable.getValue());
        Assertions.assertEquals(gradient.getFreeParameters() + 1, gradientWithMoreVariable.getFreeParameters());
        Assertions.assertEquals(0., gradientWithMoreVariable.getGradient()[gradient.getFreeParameters()]);
        Assertions.assertArrayEquals(gradient.getGradient(), Arrays.copyOfRange(gradientWithMoreVariable.getGradient(),
                0, gradient.getFreeParameters()));
    }

    @Test
    void testDoublePow() {
        assertSame(build(3).getField().getZero(), Gradient.pow(0.0, build(1.5)));
        Gradient g = Gradient.pow(2.0, build(1.5));
        DSFactory factory = new DSFactory(2, 1);
        DerivativeStructure ds = factory.constant(2.0).pow(factory.build(1.5, 1.0, FastMath.scalb(1.0, -10)));
        assertEquals(ds.getValue(), g.getValue(), 1.0e-15);
        final int[] indices = new int[ds.getFreeParameters()];
        for (int i = 0; i < g.getFreeParameters(); ++i) {
            indices[i] = 1;
            assertEquals(ds.getPartialDerivative(indices), g.getPartialDerivative(i), 1.0e-15);
            indices[i] = 0;
        }
    }

    @Test
    void testTaylor() {
        assertEquals(2.75, new Gradient(2, 1, 0.125).taylor(0.5, 2.0), 1.0e-15);
    }

    @Test
    void testOrder() {
        assertEquals(1, new Gradient(2,  1, 0.125).getOrder());
    }

    @Test
    void testGetPartialDerivative() {
        final Gradient g = new Gradient(2,  1, 0.125);
        assertEquals(2.0,   g.getPartialDerivative(0, 0), 1.0e-15); // f(x,y)
        assertEquals(1.0,   g.getPartialDerivative(1, 0), 1.0e-15); // ∂f/∂x
        assertEquals(0.125, g.getPartialDerivative(0, 1), 1.0e-15); // ∂f/∂y
    }

    @Test
    void testGetPartialDerivativeErrors() {
        final Gradient g = new Gradient(2,  1, 0.125);
        try {
            g.getPartialDerivative(0, 0, 0);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(2, ((Integer) miae.getParts()[1]).intValue());
        }
        try {
            g.getPartialDerivative(0, 5);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DERIVATION_ORDER_NOT_ALLOWED, miae.getSpecifier());
            assertEquals(5, ((Integer) miae.getParts()[0]).intValue());
        }
        try {
            g.getPartialDerivative(1, 1);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DERIVATION_ORDER_NOT_ALLOWED, miae.getSpecifier());
            assertEquals(1, ((Integer) miae.getParts()[0]).intValue());
        }
    }

    @Test
    void testHashcode() {
        assertEquals(1608501298, new Gradient(2, 1, -0.25).hashCode());
    }

    @Test
    void testEquals() {
        Gradient g = new Gradient(12, -34, 56);
        assertEquals(g, g);
        assertNotEquals("", g);
        assertEquals(new Gradient(12, -34, 56), g);
        assertNotEquals(new Gradient(21, -34, 56), g);
        assertNotEquals(new Gradient(12, -43, 56), g);
        assertNotEquals(new Gradient(12, -34, 65), g);
        assertNotEquals(new Gradient(21, -43, 65), g);
    }

    @Test
    void testRunTimeClass() {
        Field<Gradient> field = build(0.0).getField();
        assertEquals(Gradient.class, field.getRuntimeClass());
    }

    @Test
    void testConversion() {
        Gradient gA = new Gradient(-0.5, 2.5, 4.5);
        DerivativeStructure ds = gA.toDerivativeStructure();
        assertEquals(2, ds.getFreeParameters());
        assertEquals(1, ds.getOrder());
        assertEquals(-0.5, ds.getValue(), 1.0e-15);
        assertEquals(-0.5, ds.getPartialDerivative(0, 0), 1.0e-15);
        assertEquals( 2.5, ds.getPartialDerivative(1, 0), 1.0e-15);
        assertEquals( 4.5, ds.getPartialDerivative(0, 1), 1.0e-15);
        Gradient gB = new Gradient(ds);
        assertNotSame(gA, gB);
        assertEquals(gA, gB);
        try {
            new Gradient(new DSFactory(1, 2).variable(0, 1.0));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
        }
    }

    @Test
    public void testNewInstance() {
        Gradient g = build(5.25);
        assertEquals(5.25, g.getValue(), 1.0e-15);
        assertEquals(1.0,  g.getPartialDerivative(0), 1.0e-15);
        assertEquals(0.0009765625,  g.getPartialDerivative(1), 1.0e-15);
        Gradient newInstance = g.newInstance(7.5);
        assertEquals(7.5, newInstance.getValue(), 1.0e-15);
        assertEquals(0.0, newInstance.getPartialDerivative(0), 1.0e-15);
        assertEquals(0.0, newInstance.getPartialDerivative(1), 1.0e-15);
    }

    protected void checkAgainstDS(final double x,
                                  final FieldUnivariateFunction f) {
        final Gradient xG = build(x);
        final Gradient yG = f.value(xG);
        final DerivativeStructure yDS = f.value(xG.toDerivativeStructure());
        assertEquals(yDS.getFreeParameters(),
                                yG.getFreeParameters());

        if (Double.isNaN(yDS.getValue())) {
            assertEquals(yDS.getValue(), yG.getValue());
        } else {
            assertEquals(yDS.getValue(), yG.getValue(),
                                    1.0e-15 * FastMath.abs(yDS.getValue()));
        }
        final int[] indices = new int[yDS.getFreeParameters()];
        for (int i = 0; i < yG.getFreeParameters(); ++i) {
            indices[i] = 1;
            if (Double.isNaN(yDS.getPartialDerivative(indices))) {
                assertEquals(yDS.getPartialDerivative(indices),
                                        yG.getPartialDerivative(i));
            } else {
                assertEquals(yDS.getPartialDerivative(indices),
                                        yG.getPartialDerivative(i),
                                        4.0e-14 * FastMath.abs(
                                                        yDS.getPartialDerivative(
                                                                        indices)));

            }
            indices[i] = 0;
        }
    }

    @Test
    void testArithmeticVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   final S y = x.add(3).multiply(x).subtract(5).multiply(0.5);
                                   return y.negate().divide(4).divide(x).add(y).subtract(x).twice().reciprocal();
                               }
                           });
        }
    }

    @Test
    void testRemainderDoubleVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.remainder(0.5);
                               }
                           });
        }
    }

    @Test
    void testRemainderGVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                              public <S extends CalculusFieldElement<S>> S value(S x) {
                                  return x.remainder(x.divide(0.7));
                              }
                           });
        }
    }

    @Test
    void testAbsVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.abs();
                               }
                           });
        }
    }

    @Test
    void testHypotVsDS() {
        for (double x = -3.25; x < 3.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.cos().multiply(5).hypot(x.sin().twice());
                               }
                           });
        }
    }

    @Test
    void testAtan2VsDS() {
        for (double x = -3.25; x < 3.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.cos().multiply(5).atan2(x.sin().twice());
                               }
                           });
        }
    }

    @Test
    void testPowersVsDS() {
        for (double x = -3.25; x < 3.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   final FieldSinCos<S> sc = x.sinCos();
                                   return x.pow(3.2).add(x.pow(2)).subtract(sc.cos().abs().pow(sc.sin()));
                               }
                           });
        }
    }

    @Test
    void testRootsVsDS() {
        for (double x = 0.001; x < 3.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.rootN(5);//x.sqrt().add(x.cbrt()).subtract(x.rootN(5));
                               }
                           });
        }
    }

    @Test
    void testExpsLogsVsDS() {
        for (double x = 2.5; x < 3.25; x+= 0.125) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.exp().add(x.multiply(0.5).expm1()).log().log10().log1p();
                               }
                           });
        }
    }

    @Test
    void testTrigonometryVsDS() {
        for (double x = -3.25; x < 3.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.cos().multiply(x.sin()).atan().divide(12).asin().multiply(0.1).acos().tan();
                               }
                           });
        }
    }

    @Test
    void testHyperbolicVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.cosh().multiply(x.sinh()).multiply(12).abs().acosh().asinh().divide(7).tanh().multiply(0.1).atanh();
                               }
                           });
        }
    }

    @Test
    void testConvertersVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.multiply(5).toDegrees().subtract(x).toRadians();
                               }
                           });
        }
    }

    @Test
    void testLinearCombination2D2FVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.linearCombination(1.0, x.multiply(0.9),
                                                              2.0, x.multiply(0.8));
                               }
                           });
        }
    }

    @Test
    void testLinearCombination2F2FVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.linearCombination(x.add(1), x.multiply(0.9),
                                                              x.add(2), x.multiply(0.8));
                               }
                           });
        }
    }

    @Test
    void testLinearCombination3D3FVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.linearCombination(1.0, x.multiply(0.9),
                                                              2.0, x.multiply(0.8),
                                                              3.0, x.multiply(0.7));
                               }
                           });
        }
    }

    @Test
    void testLinearCombination3F3FVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.linearCombination(x.add(1), x.multiply(0.9),
                                                              x.add(2), x.multiply(0.8),
                                                              x.add(3), x.multiply(0.7));
                               }
                           });
        }
    }

    @Test
    void testLinearCombination4D4FVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.linearCombination(1.0, x.multiply(0.9),
                                                              2.0, x.multiply(0.8),
                                                              3.0, x.multiply(0.7),
                                                              4.0, x.multiply(0.6));
                               }
                           });
        }
    }

    @Test
    void testLinearCombination4F4FVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.linearCombination(x.add(1), x.multiply(0.9),
                                                              x.add(2), x.multiply(0.8),
                                                              x.add(3), x.multiply(0.7),
                                                              x.add(4), x.multiply(0.6));
                               }
                           });
        }
    }

    @Test
    void testLinearCombinationnDnFVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   final S[] b = MathArrays.buildArray(x.getField(), 4);
                                   b[0] = x.add(0.9);
                                   b[1] = x.add(0.8);
                                   b[2] = x.add(0.7);
                                   b[3] = x.add(0.6);
                                   return x.linearCombination(new double[] { 1, 2, 3, 4 }, b);
                               }
                           });
        }
    }

    @Test
    void testLinearCombinationnFnFVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   final S[] a = MathArrays.buildArray(x.getField(), 4);
                                   a[0] = x.add(1);
                                   a[1] = x.add(2);
                                   a[2] = x.add(3);
                                   a[3] = x.add(4);
                                   final S[] b = MathArrays.buildArray(x.getField(), 4);
                                   b[0] = x.add(0.9);
                                   b[1] = x.add(0.8);
                                   b[2] = x.add(0.7);
                                   b[3] = x.add(0.6);
                                   return x.linearCombination(a, b);
                               }
                           });
        }
    }

    @Test
    void testSerialization() {
        Gradient a = build(1.3);
        Gradient b = (Gradient) UnitTestUtils.serializeAndRecover(a);
        assertEquals(a, b);
        assertNotSame(a, b);
    }

    @Test
    void testZero() {
        Gradient zero = build(17.0).getField().getZero();
        assertEquals(0.0, zero.getValue(), 1.0e-15);
        for (int i = 0; i < zero.getFreeParameters(); ++i) {
            assertEquals(0.0, zero.getPartialDerivative(i), 1.0e-15);
        }
    }

    @Test
    void testOne() {
        Gradient one = build(17.0).getField().getOne();
        assertEquals(1.0, one.getValue(), 1.0e-15);
        for (int i = 0; i < one.getFreeParameters(); ++i) {
            assertEquals(0.0, one.getPartialDerivative(i), 1.0e-15);
        }
    }


    @Test
    public void testIsSmallDimensionMismatch() {
        try {
            new Gradient(0.0, 0.0, 0.0).isSmall(new Gradient(0.0, 0.0), 1.0e-15);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            Assertions.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, e.getSpecifier());
        }
    }

    @Test
    public void testIsSmall() {
        Assertions.assertTrue(new Gradient(1.0e-13, 2.5e-13).isSmall(new Gradient(1.0, 1.0), 1.0e-12));
        Assertions.assertFalse(new Gradient(1.0e-11, 2.5e-13).isSmall(new Gradient(1.0, 1.0), 1.0e-12));
        Assertions.assertFalse(new Gradient(1.0e-13, 2.5e-11).isSmall(new Gradient(1.0, 1.0), 1.0e-12));
    }

}
