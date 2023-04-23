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
import org.hipparchus.analysis.FieldUnivariateFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.MathArrays;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test for class {@link FieldGradiant}.
 */
public abstract class FieldGradientAbstractTest<T extends CalculusFieldElement<T>>
    extends CalculusFieldElementAbstractTest<FieldGradient<T>> {

    protected abstract Field<T> getValueField();

    protected FieldGradient<T> build(final double x) {
        // the function is really a two variables function : f(x) = g(x, 0) with g(x, y) = x + y / 1024
        return build(x, 1.0, FastMath.scalb(1.0, -10));
    }

    protected FieldGradient<T> build(final double x, final double... derivatives) {
        final Field<T> valueField = getValueField();
        final T[] gradient = MathArrays.buildArray(valueField, derivatives.length);
        for (int i = 0; i < gradient.length;++i) {
            gradient[i] = valueField.getZero().newInstance(derivatives[i]);
        }
        return new FieldGradient<>(valueField.getZero().newInstance(x), gradient);
    }

    protected T buildScalar(double value) {
        return getValueField().getZero().newInstance(value);
    }

    @Test
    public void testFieldAdd() {
        check(build(1.0).add(buildScalar(5.0)), 6.0, 1.0);
    }

    @Test
    public void testFieldSubtract() {
        check(build(1.0).subtract(buildScalar(5.0)), -4.0, 1.0);
    }

    @Test
    public void testFieldMultiply() {
        check(build(1.0).multiply(buildScalar(5.0)), 5.0, 5.0);
    }

    @Test
    public void testFieldDivide() {
        check(build(1.0).divide(buildScalar(5.0)), 0.2, 0.2);
    }

    @Test
    public void testgetGradient() {
        FieldGradient<T> g = build(-0.5, 2.5, 10.0, -1.0);
        Assert.assertEquals(-0.5, g.getReal(), 1.0e-15);
        Assert.assertEquals(-0.5, g.getValue().getReal(), 1.0e-15);
        Assert.assertEquals(+2.5, g.getGradient()[0].getReal(), 1.0e-15);
        Assert.assertEquals(10.0, g.getGradient()[1].getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, g.getGradient()[2].getReal(), 1.0e-15);
        Assert.assertEquals(+2.5, g.getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(10.0, g.getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, g.getPartialDerivative(2).getReal(), 1.0e-15);
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
        FieldGradient<T> g = FieldGradient.constant(5, getValueField().getZero().newInstance(-4.5));
        Assert.assertEquals(5, g.getFreeParameters());
        Assert.assertEquals(getValueField(), g.getValue().getField());
        Assert.assertEquals(-4.5, g.getValue().getReal(), 1.0e-15);
        for (int i = 0 ; i < g.getFreeParameters(); ++i) {
            Assert.assertEquals(0.0, g.getPartialDerivative(i).getReal(), 1.0e-15);
        }
    }

    @Test
    public void testVariable() {
        FieldGradient<T> g = FieldGradient.variable(5, 1, getValueField().getZero().newInstance(-4.5));
        Assert.assertEquals(5, g.getFreeParameters());
        Assert.assertEquals(getValueField(), g.getValue().getField());
        Assert.assertEquals(-4.5, g.getValue().getReal(), 1.0e-15);
        for (int i = 0 ; i < g.getFreeParameters(); ++i) {
            Assert.assertEquals(i == 1 ? 1.0 : 0.0, g.getPartialDerivative(i).getReal(), 1.0e-15);
        }
    }

    @Test
    public void testDoublePow() {
        Assert.assertSame(build(3).getField().getZero(), FieldGradient.pow(0.0, build(1.5)));
        FieldGradient<T> g = FieldGradient.pow(2.0, build(1.5));
        FDSFactory<T> factory = new FDSFactory<>(getValueField(), 2, 1);
        FieldDerivativeStructure<T> ds = factory.constant(2.0).pow(factory.build(1.5, 1.0, FastMath.scalb(1.0, -10)));
        Assert.assertEquals(ds.getValue().getReal(), g.getValue().getReal(), 1.0e-15);
        final int[] indices = new int[ds.getFreeParameters()];
        for (int i = 0; i < g.getFreeParameters(); ++i) {
            indices[i] = 1;
            Assert.assertEquals(ds.getPartialDerivative(indices).getReal(), g.getPartialDerivative(i).getReal(), 1.0e-15);
            indices[i] = 0;
        }
    }

    @Test
    public void testTaylor() {
        Assert.assertEquals(2.75, build(2, 1, 0.125).taylor(0.5, 2.0).getReal(), 1.0e-15);
    }

    @Test
    public void testOrder() {
        Assert.assertEquals(1, build(2).getOrder());
    }

    @Test
    public void testGetPartialDerivative() {
        final FieldGradient<T> g = build(2);
        Assert.assertEquals(2.0,        g.getPartialDerivative(0, 0).getReal(), 1.0e-15); // f(x,y)
        Assert.assertEquals(1.0,        g.getPartialDerivative(1, 0).getReal(), 1.0e-15); // ∂f/∂x
        Assert.assertEquals(1.0 / 1024, g.getPartialDerivative(0, 1).getReal(), 1.0e-15); // ∂f/∂y
    }

    @Test
    public void testGetPartialDerivativeErrors() {
        final FieldGradient<T> g = build(2);
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
        Assert.assertEquals(1608501298, build(2, 1, -0.25).hashCode());
    }

    @Test
    public void testEquals() {
        FieldGradient<T> g = build(12, -34, 56);
        Assert.assertEquals(g, g);
        Assert.assertNotEquals(g, "");
        Assert.assertEquals(g, build(12, -34, 56));
        Assert.assertNotEquals(g, build(21, -34, 56));
        Assert.assertNotEquals(g, build(12, -43, 56));
        Assert.assertNotEquals(g, build(12, -34, 65));
        Assert.assertNotEquals(g, build(21, -43, 65));
    }

    @Test
    public void testRunTimeClass() {
        Field<FieldGradient<T>> field = build(0.0).getField();
        Assert.assertEquals(FieldGradient.class, field.getRuntimeClass());
    }

    @Test
    public void testConversion() {
        FieldGradient<T> gA = build(-0.5, 2.5, 4.5);
        FieldDerivativeStructure<T> ds = gA.toDerivativeStructure();
        Assert.assertEquals(2, ds.getFreeParameters());
        Assert.assertEquals(1, ds.getOrder());
        Assert.assertEquals(-0.5, ds.getValue().getReal(), 1.0e-15);
        Assert.assertEquals(-0.5, ds.getPartialDerivative(0, 0).getReal(), 1.0e-15);
        Assert.assertEquals( 2.5, ds.getPartialDerivative(1, 0).getReal(), 1.0e-15);
        Assert.assertEquals( 4.5, ds.getPartialDerivative(0, 1).getReal(), 1.0e-15);
        FieldGradient<T> gB = new FieldGradient<>(ds);
        Assert.assertNotSame(gA, gB);
        Assert.assertEquals(gA, gB);
        try {
            new FieldGradient<>(new FDSFactory<>(getValueField(), 1, 2).variable(0, 1.0));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
        }
    }

    @Test
    public void testNewInstance() {
        FieldGradient<T> g = build(5.25);
        Assert.assertEquals(5.25, g.getValue().getReal(), 1.0e-15);
        Assert.assertEquals(1.0,  g.getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(0.0009765625,  g.getPartialDerivative(1).getReal(), 1.0e-15);
        FieldGradient<T> newInstance = g.newInstance(7.5);
        Assert.assertEquals(7.5, newInstance.getValue().getReal(), 1.0e-15);
        Assert.assertEquals(0.0, newInstance.getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(0.0, newInstance.getPartialDerivative(1).getReal(), 1.0e-15);
    }

    protected void checkAgainstDS(final double x, final FieldUnivariateFunction f) {
        final FieldGradient<T> xG = build(x);
        final FieldGradient<T> yG = f.value(xG);
        final FieldDerivativeStructure<T> yDS = f.value(xG.toDerivativeStructure());
        Assert.assertEquals(yDS.getFreeParameters(), yG.getFreeParameters());
        Assert.assertEquals(yDS.getValue().getReal(), yG.getValue().getReal(), 1.0e-15 * FastMath.abs(yDS.getValue().getReal()));
        final int[] indices = new int[yDS.getFreeParameters()];
        for (int i = 0; i < yG.getFreeParameters(); ++i) {
            indices[i] = 1;
            Assert.assertEquals(yDS.getPartialDerivative(indices).getReal(),
                                yG.getPartialDerivative(i).getReal(),
                                4.0e-14* FastMath.abs(yDS.getPartialDerivative(indices).getReal()));
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
    public void testCopySignField() {

        FieldGradient<T> minusOne = build(-1.0);
        Assert.assertEquals(+1.0, minusOne.copySign(buildScalar(+1.0)).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(buildScalar(-1.0)).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(buildScalar(+0.0)).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(buildScalar(-0.0)).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(buildScalar(Double.NaN)).getReal(), 1.0e-15);

        FieldGradient<T> plusOne = build(1.0);
        Assert.assertEquals(+1.0, plusOne.copySign(buildScalar(+1.0)).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, plusOne.copySign(buildScalar(-1.0)).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(buildScalar(+0.0)).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, plusOne.copySign(buildScalar(-0.0)).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(buildScalar(Double.NaN)).getReal(), 1.0e-15);

    }

    @Test
    public void testRemainderField() {
        double epsilon = 2.0e-15;
        for (double x = -1.7; x < 2; x += 0.2) {
            FieldGradient<T> dsX = build(x);
            for (double y = -1.7; y < 2; y += 0.2) {
                FieldGradient<T> remainder = dsX.remainder(buildScalar(y));
                FieldGradient<T> ref = dsX.subtract(x - FastMath.IEEEremainder(x, y));
                FieldGradient<T> zero = remainder.subtract(ref);
                Assert.assertEquals(0, zero.getPartialDerivative(0).getReal(), epsilon);
            }
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
    public void testScalbVsDS() {
        for (int n = -4; n < 4; ++n) {
            final int theN = n;
            for (double x = -1.25; x < 1.25; x+= 0.5) {
                checkAgainstDS(x,
                               new FieldUnivariateFunction() {
                                   public <S extends CalculusFieldElement<S>> S value(S x) {
                                       return x.scalb(theN);
                                   }
                               });
            }
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
    public void testLinearCombinationField() {
        final T[] a = MathArrays.buildArray(getValueField(), 3);
        a[0] = buildScalar(-1321008684645961.0 / 268435456.0);
        a[1] = buildScalar(-5774608829631843.0 / 268435456.0);
        a[2] = buildScalar(-7645843051051357.0 / 8589934592.0);
        final FieldGradient<T>[] b = MathArrays.buildArray(FieldGradientField.getField(getValueField(), 1), 3);
        b[0] = build(-5712344449280879.0 / 2097152.0);
        b[1] = build(-4550117129121957.0 / 2097152.0);
        b[2] = build(8846951984510141.0 / 131072.0);

        final FieldGradient<T> abSumInline = b[0].linearCombination(a[0], b[0],
                                                                    a[1], b[1],
                                                                    a[2], b[2]);
        final FieldGradient<T> abSumArray = b[0].linearCombination(a, b);
        Assert.assertEquals(abSumInline.getReal(), abSumArray.getReal(), 3.0e-8);
        Assert.assertEquals(-1.8551294182586248737720779899, abSumInline.getReal(), 5.0e-8);
        Assert.assertEquals(abSumInline.getPartialDerivative(0).getReal(), abSumArray.getPartialDerivative(0).getReal(), 3.0e-8);
    }

    @Test
    public void testZero() {
        FieldGradient<T> zero = build(17.0).getField().getZero();
        Assert.assertEquals(0.0, zero.getValue().getReal(), 1.0e-15);
        for (int i = 0; i < zero.getFreeParameters(); ++i) {
            Assert.assertEquals(0.0, zero.getPartialDerivative(i).getReal(), 1.0e-15);
        }
    }

    @Test
    public void testOne() {
        FieldGradient<T> one = build(17.0).getField().getOne();
        Assert.assertEquals(1.0, one.getValue().getReal(), 1.0e-15);
        for (int i = 0; i < one.getFreeParameters(); ++i) {
            Assert.assertEquals(0.0, one.getPartialDerivative(i).getReal(), 1.0e-15);
        }
    }

    private void check(FieldGradient<T> g, double value, double... derivatives) {

        // check value
        Assert.assertEquals(value, g.getReal(), 1.0e-15);

        // check derivatives
        for (int i = 0; i < derivatives.length; ++i) {
            Assert.assertEquals(derivatives[i], g.getPartialDerivative(i).getReal(), 1.0e-15);
        }

    }
}
