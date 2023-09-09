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

import org.hipparchus.Field;
import org.hipparchus.CalculusFieldElement;
import org.hipparchus.CalculusFieldElementAbstractTest;
import org.hipparchus.UnitTestUtils;
import org.hipparchus.analysis.FieldUnivariateFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.MathArrays;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test for class {@link UnivariateDerivative}.
 */
public class GradientTest extends CalculusFieldElementAbstractTest<Gradient> {

    @Override
    protected Gradient build(final double x) {
        // the function is really a two variables function : f(x) = g(x, 0) with g(x, y) = x + y / 1024
        return new Gradient(x, 1.0, FastMath.scalb(1.0, -10));
    }

    @Test
    public void testGetGradient() {
        Gradient g = new Gradient(-0.5, 2.5, 10.0, -1.0);
        Assert.assertEquals(-0.5, g.getReal(), 1.0e-15);
        Assert.assertEquals(-0.5, g.getValue(), 1.0e-15);
        Assert.assertEquals(+2.5, g.getGradient()[0], 1.0e-15);
        Assert.assertEquals(10.0, g.getGradient()[1], 1.0e-15);
        Assert.assertEquals(-1.0, g.getGradient()[2], 1.0e-15);
        Assert.assertEquals(+2.5, g.getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(10.0, g.getPartialDerivative(1), 1.0e-15);
        Assert.assertEquals(-1.0, g.getPartialDerivative(2), 1.0e-15);
        Assert.assertEquals(3, g.getFreeParameters());
        try {
            g.getPartialDerivative(-1);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, miae.getSpecifier());
        }
        try {
            g.getPartialDerivative(+3);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, miae.getSpecifier());
        }
    }

    @Test
    public void testConstant() {
        Gradient g = Gradient.constant(5, -4.5);
        Assert.assertEquals(5, g.getFreeParameters());
        Assert.assertEquals(-4.5, g.getValue(), 1.0e-15);
        for (int i = 0 ; i < g.getFreeParameters(); ++i) {
            Assert.assertEquals(0.0, g.getPartialDerivative(i), 1.0e-15);
        }
    }

    @Test
    public void testVariable() {
        Gradient g = Gradient.variable(5, 1, -4.5);
        Assert.assertEquals(5, g.getFreeParameters());
        Assert.assertEquals(-4.5, g.getValue(), 1.0e-15);
        for (int i = 0 ; i < g.getFreeParameters(); ++i) {
            Assert.assertEquals(i == 1 ? 1.0 : 0.0, g.getPartialDerivative(i), 1.0e-15);
        }
    }

    @Test
    public void testDoublePow() {
        Assert.assertSame(build(3).getField().getZero(), Gradient.pow(0.0, build(1.5)));
        Gradient g = Gradient.pow(2.0, build(1.5));
        DSFactory factory = new DSFactory(2, 1);
        DerivativeStructure ds = factory.constant(2.0).pow(factory.build(1.5, 1.0, FastMath.scalb(1.0, -10)));
        Assert.assertEquals(ds.getValue(), g.getValue(), 1.0e-15);
        final int[] indices = new int[ds.getFreeParameters()];
        for (int i = 0; i < g.getFreeParameters(); ++i) {
            indices[i] = 1;
            Assert.assertEquals(ds.getPartialDerivative(indices), g.getPartialDerivative(i), 1.0e-15);
            indices[i] = 0;
        }
    }

    @Test
    public void testTaylor() {
        Assert.assertEquals(2.75, new Gradient(2, 1, 0.125).taylor(0.5, 2.0), 1.0e-15);
    }

    @Test
    public void testOrder() {
        Assert.assertEquals(1, new Gradient(2,  1, 0.125).getOrder());
    }

    @Test
    public void testGetPartialDerivative() {
        final Gradient g = new Gradient(2,  1, 0.125);
        Assert.assertEquals(2.0,   g.getPartialDerivative(0, 0), 1.0e-15); // f(x,y)
        Assert.assertEquals(1.0,   g.getPartialDerivative(1, 0), 1.0e-15); // ∂f/∂x
        Assert.assertEquals(0.125, g.getPartialDerivative(0, 1), 1.0e-15); // ∂f/∂y
    }

    @Test
    public void testGetPartialDerivativeErrors() {
        final Gradient g = new Gradient(2,  1, 0.125);
        try {
            g.getPartialDerivative(0, 0, 0);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            Assert.assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            Assert.assertEquals(2, ((Integer) miae.getParts()[1]).intValue());
        }
        try {
            g.getPartialDerivative(0, 5);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DERIVATION_ORDER_NOT_ALLOWED, miae.getSpecifier());
            Assert.assertEquals(5, ((Integer) miae.getParts()[0]).intValue());
        }
        try {
            g.getPartialDerivative(1, 1);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DERIVATION_ORDER_NOT_ALLOWED, miae.getSpecifier());
            Assert.assertEquals(1, ((Integer) miae.getParts()[0]).intValue());
        }
    }

    @Test
    public void testHashcode() {
        Assert.assertEquals(1608501298, new Gradient(2, 1, -0.25).hashCode());
    }

    @Test
    public void testEquals() {
        Gradient g = new Gradient(12, -34, 56);
        Assert.assertEquals(g, g);
        Assert.assertNotEquals(g, "");
        Assert.assertEquals(g, new Gradient(12, -34, 56));
        Assert.assertNotEquals(g, new Gradient(21, -34, 56));
        Assert.assertNotEquals(g, new Gradient(12, -43, 56));
        Assert.assertNotEquals(g, new Gradient(12, -34, 65));
        Assert.assertNotEquals(g, new Gradient(21, -43, 65));
    }

    @Test
    public void testRunTimeClass() {
        Field<Gradient> field = build(0.0).getField();
        Assert.assertEquals(Gradient.class, field.getRuntimeClass());
    }

    @Test
    public void testConversion() {
        Gradient gA = new Gradient(-0.5, 2.5, 4.5);
        DerivativeStructure ds = gA.toDerivativeStructure();
        Assert.assertEquals(2, ds.getFreeParameters());
        Assert.assertEquals(1, ds.getOrder());
        Assert.assertEquals(-0.5, ds.getValue(), 1.0e-15);
        Assert.assertEquals(-0.5, ds.getPartialDerivative(0, 0), 1.0e-15);
        Assert.assertEquals( 2.5, ds.getPartialDerivative(1, 0), 1.0e-15);
        Assert.assertEquals( 4.5, ds.getPartialDerivative(0, 1), 1.0e-15);
        Gradient gB = new Gradient(ds);
        Assert.assertNotSame(gA, gB);
        Assert.assertEquals(gA, gB);
        try {
            new Gradient(new DSFactory(1, 2).variable(0, 1.0));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
        }
    }

    @Test
    public void testNewInstance() {
        Gradient g = build(5.25);
        Assert.assertEquals(5.25, g.getValue(), 1.0e-15);
        Assert.assertEquals(1.0,  g.getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(0.0009765625,  g.getPartialDerivative(1), 1.0e-15);
        Gradient newInstance = g.newInstance(7.5);
        Assert.assertEquals(7.5, newInstance.getValue(), 1.0e-15);
        Assert.assertEquals(0.0, newInstance.getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(0.0, newInstance.getPartialDerivative(1), 1.0e-15);
    }

    protected void checkAgainstDS(final double x, final FieldUnivariateFunction f) {
        final Gradient xG = build(x);
        final Gradient yG = f.value(xG);
        final DerivativeStructure yDS = f.value(xG.toDerivativeStructure());
        Assert.assertEquals(yDS.getFreeParameters(), yG.getFreeParameters());
        Assert.assertEquals(yDS.getValue(), yG.getValue(), 1.0e-15 * FastMath.abs(yDS.getValue()));
        final int[] indices = new int[yDS.getFreeParameters()];
        for (int i = 0; i < yG.getFreeParameters(); ++i) {
            indices[i] = 1;
            Assert.assertEquals(yDS.getPartialDerivative(indices),
                                yG.getPartialDerivative(i),
                                4.0e-14* FastMath.abs(yDS.getPartialDerivative(indices)));
            indices[i] = 0;
        }
    }

    @Test
    public void testArithmeticVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   final S y = x.add(3).multiply(x).subtract(5).multiply(0.5);
                                   return y.negate().divide(4).divide(x).add(y).subtract(x).multiply(2).reciprocal();
                               }
                           });
        }
    }

    @Test
    public void testRemainderDoubleVsDS() {
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
    public void testRemainderGVsDS() {
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
    public void testAbsVsDS() {
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
    public void testHypotVsDS() {
        for (double x = -3.25; x < 3.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.cos().multiply(5).hypot(x.sin().multiply(2));
                               }
                           });
        }
    }

    @Test
    public void testAtan2VsDS() {
        for (double x = -3.25; x < 3.25; x+= 0.5) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.cos().multiply(5).atan2(x.sin().multiply(2));
                               }
                           });
        }
    }

    @Test
    public void testPowersVsDS() {
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
    public void testRootsVsDS() {
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
    public void testExpsLogsVsDS() {
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
    public void testTrigonometryVsDS() {
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
    public void testHyperbolicVsDS() {
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
    public void testConvertersVsDS() {
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
    public void testLinearCombination2D2FVsDS() {
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
    public void testLinearCombination2F2FVsDS() {
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
    public void testLinearCombination3D3FVsDS() {
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
    public void testLinearCombination3F3FVsDS() {
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
    public void testLinearCombination4D4FVsDS() {
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
    public void testLinearCombination4F4FVsDS() {
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
    public void testLinearCombinationnDnFVsDS() {
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
    public void testLinearCombinationnFnFVsDS() {
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
    public void testSerialization() {
        Gradient a = build(1.3);
        Gradient b = (Gradient) UnitTestUtils.serializeAndRecover(a);
        Assert.assertEquals(a, b);
        Assert.assertNotSame(a, b);
    }

    @Test
    public void testZero() {
        Gradient zero = build(17.0).getField().getZero();
        Assert.assertEquals(0.0, zero.getValue(), 1.0e-15);
        for (int i = 0; i < zero.getFreeParameters(); ++i) {
            Assert.assertEquals(0.0, zero.getPartialDerivative(i), 1.0e-15);
        }
    }

    @Test
    public void testOne() {
        Gradient one = build(17.0).getField().getOne();
        Assert.assertEquals(1.0, one.getValue(), 1.0e-15);
        for (int i = 0; i < one.getFreeParameters(); ++i) {
            Assert.assertEquals(0.0, one.getPartialDerivative(i), 1.0e-15);
        }
    }

}
