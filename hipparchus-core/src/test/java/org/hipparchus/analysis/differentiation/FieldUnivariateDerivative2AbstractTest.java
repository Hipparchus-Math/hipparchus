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
 * Test for class {@link FieldUnivariateDerivative2}.
 */
public abstract class FieldUnivariateDerivative2AbstractTest<T extends CalculusFieldElement<T>>
    extends CalculusFieldElementAbstractTest<FieldUnivariateDerivative2<T>> {

    protected abstract Field<T> getValueField();

    protected FieldUnivariateDerivative2<T> build(double value) {
        final Field<T> valueField = getValueField();
        return new FieldUnivariateDerivative2<>(valueField.getZero().newInstance(value),
                                                valueField.getOne(),
                                                valueField.getZero());
    }

    protected T buildScalar(double value) {
        return getValueField().getZero().newInstance(value);
    }

    private int getMaxOrder() {
        return 2;
    }

    protected FieldUnivariateDerivative2<T> build(final double f0, final double f1, final double f2) {
        T prototype = build(0.0).getValue();
        return new FieldUnivariateDerivative2<>(prototype.newInstance(f0),
                                                prototype.newInstance(f1),
                                                prototype.newInstance(f2));
    }

    @Test
    public void testFieldAdd() {
        check(build(1.0, 1.0, 1.0).add(buildScalar(5.0)), 6.0, 1.0, 1.0);
    }

    @Test
    public void testFieldSubtract() {
        check(build(1.0, 1.0, 1.0).subtract(buildScalar(5.0)), -4.0, 1.0, 1.0);
    }

    @Test
    public void testFieldMultiply() {
        check(build(1.0, 1.0, 1.0).multiply(buildScalar(5.0)), 5.0, 5.0, 5.0);
    }

    @Test
    public void testFieldDivide() {
        check(build(1.0, 1.0, 1.0).divide(buildScalar(5.0)), 0.2, 0.2, 0.2);
    }

    @Test
    public void testOrder() {
        Assert.assertEquals(getMaxOrder(), build(0).getOrder());
    }

    @Test
    public void testNewInstance() {
        FieldUnivariateDerivative2<T> ud = build(5.25);
        Assert.assertEquals(5.25, ud.getValue().getReal(), 1.0e-15);
        Assert.assertEquals(1.0,  ud.getDerivative(1).getReal(), 1.0e-15);
        FieldUnivariateDerivative2<T> newInstance = ud.newInstance(7.5);
        Assert.assertEquals(7.5, newInstance.getValue().getReal(), 1.0e-15);
        Assert.assertEquals(0.0, newInstance.getDerivative(1).getReal(), 1.0e-15);
    }

    @Test
    public void testGetPartialDerivative() {
        try {
            build(3.0).getPartialDerivative(0, 1);
            Assert.fail("an exception should have been thrown");
        } catch( MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            Assert.assertEquals(2, ((Integer) miae.getParts()[0]).intValue());
            Assert.assertEquals(1, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    public void testGetDerivative() {
        FieldUnivariateDerivative2<T> x  = build(3.0);
        FieldUnivariateDerivative2<T> ud = x.multiply(x);
        try {
            ud.getDerivative(-1);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DERIVATION_ORDER_NOT_ALLOWED, miae.getSpecifier());
        }
        Assert.assertEquals(9.0, ud.getValue().getReal(), 1.0e-15);
        Assert.assertEquals(9.0, ud.getDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(6.0, ud.getDerivative(1).getReal(), 1.0e-15);
        for (int n = 2; n <= getMaxOrder(); ++n) {
            Assert.assertEquals(n == 2 ? 2.0 : 0.0, ud.getDerivative(n).getReal(), 1.0e-15);
        }
        try {
            ud.getDerivative(getMaxOrder() + 1);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DERIVATION_ORDER_NOT_ALLOWED, miae.getSpecifier());
        }
    }

    @Test
    public void testGetFreeParameters() {
        Assert.assertEquals(1, build(3.0).getFreeParameters());
    }

    protected void checkAgainstDS(final double x, final FieldUnivariateFunction f) {
        final FieldUnivariateDerivative2<T> xUD = build(x);
        final FieldUnivariateDerivative2<T> yUD = f.value(xUD);
        final FieldDerivativeStructure<T> yDS = f.value(xUD.toDerivativeStructure());
        for (int i = 0; i <= yUD.getOrder(); ++i) {
            Assert.assertEquals(yDS.getPartialDerivative(i).getReal(),
                                yUD.getDerivative(i).getReal(),
                                4.0e-14* FastMath.abs(yDS.getPartialDerivative(i).getReal()));
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

        FieldUnivariateDerivative2<T> minusOne = build(-1.0);
        Assert.assertEquals(+1.0, minusOne.copySign(buildScalar(+1.0)).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(buildScalar(-1.0)).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(buildScalar(+0.0)).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(buildScalar(-0.0)).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(buildScalar(Double.NaN)).getReal(), 1.0e-15);

        FieldUnivariateDerivative2<T> plusOne = build(1.0);
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
            FieldUnivariateDerivative2<T> dsX = build(x);
            for (double y = -1.7; y < 2; y += 0.2) {
                FieldUnivariateDerivative2<T> remainder = dsX.remainder(buildScalar(y));
                FieldUnivariateDerivative2<T> ref = dsX.subtract(x - FastMath.IEEEremainder(x, y));
                FieldUnivariateDerivative2<T> zero = remainder.subtract(ref);
                Assert.assertEquals(0, zero.getFirstDerivative().getReal(), epsilon);
                Assert.assertEquals(0, zero.getSecondDerivative().getReal(), epsilon);
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
    public void testRemainderUdVsDS() {
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
    public void testUlpVsDS() {
        for (double x = -1.25; x < 1.25; x+= 0.0001) {
            checkAgainstDS(x,
                           new FieldUnivariateFunction() {
                               public <S extends CalculusFieldElement<S>> S value(S x) {
                                   return x.ulp();
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
    public void testSqrtVsDS() {
        for (double x = 0.001; x < 3.25; x+= 0.5) {
            checkAgainstDS(x,
                    new FieldUnivariateFunction() {
                        public <S extends CalculusFieldElement<S>> S value(S x) {
                            return x.sqrt();
                        }
                    });
        }
    }

    @Test
    public void testCbrtVsDS() {
        for (double x = 0.001; x < 3.25; x+= 0.5) {
            checkAgainstDS(x,
                    new FieldUnivariateFunction() {
                        public <S extends CalculusFieldElement<S>> S value(S x) {
                            return x.cbrt();
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
        final FieldUnivariateDerivative2<T>[] b = MathArrays.buildArray(FieldUnivariateDerivative2Field.getUnivariateDerivative2Field(getValueField()), 3);
        b[0] = build(-5712344449280879.0 / 2097152.0);
        b[1] = build(-4550117129121957.0 / 2097152.0);
        b[2] = build(8846951984510141.0 / 131072.0);

        final FieldUnivariateDerivative2<T> abSumInline = b[0].linearCombination(a[0], b[0],
                                                                                 a[1], b[1],
                                                                                 a[2], b[2]);
        final FieldUnivariateDerivative2<T> abSumArray = b[0].linearCombination(a, b);
        Assert.assertEquals(abSumInline.getReal(), abSumArray.getReal(), 3.0e-8);
        Assert.assertEquals(-1.8551294182586248737720779899, abSumInline.getReal(), 5.0e-8);
        Assert.assertEquals(abSumInline.getFirstDerivative().getReal(), abSumArray.getFirstDerivative().getReal(), 3.0e-8);
        Assert.assertEquals(abSumInline.getSecondDerivative().getReal(), abSumArray.getSecondDerivative().getReal(), 3.0e-8);
    }

    @Test
    public void testGetFirstAndSecondDerivative() {
        FieldUnivariateDerivative2<T> ud2 = build(-0.5, 2.5, 4.5);
        Assert.assertEquals(-0.5, ud2.getReal(), 1.0e-15);
        Assert.assertEquals(-0.5, ud2.getValue().getReal(), 1.0e-15);
        Assert.assertEquals(+2.5, ud2.getFirstDerivative().getReal(), 1.0e-15);
        Assert.assertEquals(+4.5, ud2.getSecondDerivative().getReal(), 1.0e-15);
    }

    @Test
    public void testConversion() {
        FieldUnivariateDerivative2<T> udA = build(-0.5, 2.5, 4.5);
        FieldDerivativeStructure<T> ds = udA.toDerivativeStructure();
        Assert.assertEquals(1, ds.getFreeParameters());
        Assert.assertEquals(2, ds.getOrder());
        Assert.assertEquals(-0.5, ds.getValue().getReal(), 1.0e-15);
        Assert.assertEquals(-0.5, ds.getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals( 2.5, ds.getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals( 4.5, ds.getPartialDerivative(2).getReal(), 1.0e-15);
        FieldUnivariateDerivative2<T> udB = new FieldUnivariateDerivative2<>(ds);
        Assert.assertNotSame(udA, udB);
        Assert.assertEquals(udA, udB);
        try {
            new FieldUnivariateDerivative2<>(new FDSFactory<>(getValueField(), 2, 2).variable(0, 1.0));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
        }
        try {
            new FieldUnivariateDerivative2<>(new FDSFactory<>(getValueField(), 1, 1).variable(0, 1.0));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
        }
    }

    @Test
    public void testDoublePow() {
        Assert.assertSame(build(3).getField().getZero(), FieldUnivariateDerivative2.pow(0.0, build(1.5)));
        FieldUnivariateDerivative2<T> ud = FieldUnivariateDerivative2.pow(2.0, build(1.5));
        DSFactory factory = new DSFactory(1, 2);
        DerivativeStructure ds = factory.constant(2.0).pow(factory.variable(0, 1.5));
        Assert.assertEquals(ds.getValue(), ud.getValue().getReal(), 1.0e-15);
        Assert.assertEquals(ds.getPartialDerivative(1), ud.getFirstDerivative().getReal(), 1.0e-15);
        Assert.assertEquals(ds.getPartialDerivative(2), ud.getSecondDerivative().getReal(), 1.0e-15);
    }

    @Test
    public void testTaylor() {
        Assert.assertEquals(-0.125, build(1, -3, 4).taylor(0.75).getReal(), 1.0e-15);
        Assert.assertEquals(-0.125, build(1, -3, 4).taylor(getValueField().getZero().newInstance(0.75)).getReal(), 1.0e-15);
    }

    @Test
    public void testHashcode() {
        Assert.assertEquals(-1025507011, build(2, 1, -1).hashCode());
    }

    @Test
    public void testEquals() {
        FieldUnivariateDerivative2<T> ud2 = build(12, -34, 56);
        Assert.assertEquals(ud2, ud2);
        Assert.assertNotEquals(ud2, "");
        Assert.assertEquals(ud2, build(12, -34, 56));
        Assert.assertNotEquals(ud2, build(21, -34, 56));
        Assert.assertNotEquals(ud2, build(12, -43, 56));
        Assert.assertNotEquals(ud2, build(12, -34, 65));
        Assert.assertNotEquals(ud2, build(21, -43, 65));
    }

    @Test
    public void testRunTimeClass() {
        Field<FieldUnivariateDerivative2<T>> field = build(0.0).getField();
        Assert.assertEquals(FieldUnivariateDerivative2.class, field.getRuntimeClass());
    }

    private void check(FieldUnivariateDerivative2<T> ud2, double value,
                       double derivative1, double derivative2) {

        // check value
        Assert.assertEquals(value, ud2.getReal(), 1.0e-15);

        // check derivatives
        Assert.assertEquals(derivative1, ud2.getFirstDerivative().getReal(), 1.0e-15);
        Assert.assertEquals(derivative2, ud2.getSecondDerivative().getReal(), 1.0e-15);

    }

}
