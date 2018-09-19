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

package org.hipparchus.analysis.differentiation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.hipparchus.ExtendedFieldElementAbstractTest;
import org.hipparchus.Field;
import org.hipparchus.analysis.polynomials.PolynomialFunction;
import org.hipparchus.dfp.Dfp;
import org.hipparchus.dfp.DfpField;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.ArithmeticUtils;
import org.hipparchus.util.CombinatoricsUtils;
import org.hipparchus.util.Decimal64;
import org.hipparchus.util.Decimal64Field;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.Precision;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test for class {@link FieldDerivativeStructure}.
 */
public class FieldDerivativeStructureTest extends ExtendedFieldElementAbstractTest<FieldDerivativeStructure<Decimal64>> {

    @Override
    protected FieldDerivativeStructure<Decimal64> build(final double x) {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 2, 1);
        return factory.variable(0, x);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testWrongFieldVariableIndex() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, 1);
        factory.variable(3, new Decimal64(1.0));
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testWrongPrimitiveVariableIndex() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, 1);
        factory.variable(3, 1.0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testMissingOrders() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, 1);
        factory.variable(0, 1.0).getPartialDerivative(0, 1);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testWrongDimensionField() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, 1);
        factory.build(Decimal64.ONE, Decimal64.ONE, Decimal64.ONE, Decimal64.ONE, Decimal64.ONE);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testWrongDimensionPrimitive() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, 1);
        factory.build(1.0, 1.0, 1.0, 1.0, 1.0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testTooLargeOrder() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, 1);
        factory.variable(0, 1.0).getPartialDerivative(1, 1, 2);
    }

    @Test
    public void testVariableWithoutDerivativeField() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, 0);
        FieldDerivativeStructure<Decimal64> v = factory.variable(0, new Decimal64(1.0));
        Assert.assertEquals(1.0, v.getReal(), 1.0e-15);
    }

    @Test
    public void testVariableWithoutDerivativePrimitive() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, 0);
        FieldDerivativeStructure<Decimal64> v = factory.variable(0, 1.0);
        Assert.assertEquals(1.0, v.getReal(), 1.0e-15);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testVariableWithoutDerivative1() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, 0);
        FieldDerivativeStructure<Decimal64> v = factory.variable(0, 1.0);
        Assert.assertEquals(1.0, v.getPartialDerivative(1).getReal(), 1.0e-15);
    }

    @Test
    public void testVariable() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            checkF0F1(factory.variable(0, 1.0),
                      1.0, 1.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0),
                      2.0, 0.0, 1.0, 0.0);
            checkF0F1(factory.variable(2, 3.0),
                      3.0, 0.0, 0.0, 1.0);
        }
    }

    @Test
    public void testConstant() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            checkF0F1(factory.constant(FastMath.PI),
                      FastMath.PI, 0.0, 0.0, 0.0);
        }
    }

    @Test
    public void testFieldAdd() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).add(new Decimal64(5)), 6.0, 1.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).add(new Decimal64(5)), 7.0, 0.0, 1.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).add(new Decimal64(5)), 8.0, 0.0, 0.0, 1.0);
        }
    }

    @Test
    public void testPrimitiveAdd() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).add(5), 6.0, 1.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).add(5), 7.0, 0.0, 1.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).add(5), 8.0, 0.0, 0.0, 1.0);
        }
    }

    @Test
    public void testAdd() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            FieldDerivativeStructure<Decimal64> x = factory.variable(0, 1.0);
            FieldDerivativeStructure<Decimal64> y = factory.variable(1, 2.0);
            FieldDerivativeStructure<Decimal64> z = factory.variable(2, 3.0);
            FieldDerivativeStructure<Decimal64> xyz = x.add(y.add(z));
            checkF0F1(xyz, x.getValue().getReal() + y.getValue().getReal() + z.getValue().getReal(), 1.0, 1.0, 1.0);
        }
    }

    @Test
    public void testFieldSubtract() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).subtract(new Decimal64(5)), -4.0, 1.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).subtract(new Decimal64(5)), -3.0, 0.0, 1.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).subtract(new Decimal64(5)), -2.0, 0.0, 0.0, 1.0);
        }
    }

    @Test
    public void testPrimitiveSubtract() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).subtract(5), -4.0, 1.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).subtract(5), -3.0, 0.0, 1.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).subtract(5), -2.0, 0.0, 0.0, 1.0);
        }
    }

    @Test
    public void testSubtract() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            FieldDerivativeStructure<Decimal64> x = factory.variable(0, 1.0);
            FieldDerivativeStructure<Decimal64> y = factory.variable(1, 2.0);
            FieldDerivativeStructure<Decimal64> z = factory.variable(2, 3.0);
            FieldDerivativeStructure<Decimal64> xyz = x.subtract(y.subtract(z));
            checkF0F1(xyz, x.getReal() - (y.getReal() - z.getReal()), 1.0, -1.0, 1.0);
        }
    }

    @Test
    public void testFieldMultiply() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).multiply(new Decimal64(5)),  5.0, 5.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).multiply(new Decimal64(5)), 10.0, 0.0, 5.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).multiply(new Decimal64(5)), 15.0, 0.0, 0.0, 5.0);
        }
    }

    @Test
    public void testPrimitiveMultiply() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).multiply(5),  5.0, 5.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).multiply(5), 10.0, 0.0, 5.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).multiply(5), 15.0, 0.0, 0.0, 5.0);
        }
    }

    @Test
    public void testMultiply() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            FieldDerivativeStructure<Decimal64> x = factory.variable(0, 1.0);
            FieldDerivativeStructure<Decimal64> y = factory.variable(1, 2.0);
            FieldDerivativeStructure<Decimal64> z = factory.variable(2, 3.0);
            FieldDerivativeStructure<Decimal64> xyz = x.multiply(y.multiply(z));
            for (int i = 0; i <= maxOrder; ++i) {
                for (int j = 0; j <= maxOrder; ++j) {
                    for (int k = 0; k <= maxOrder; ++k) {
                        if (i + j + k <= maxOrder) {
                            Assert.assertEquals((i == 0 ? x.getReal() : (i == 1 ? 1.0 : 0.0)) *
                                                (j == 0 ? y.getReal() : (j == 1 ? 1.0 : 0.0)) *
                                                (k == 0 ? z.getReal() : (k == 1 ? 1.0 : 0.0)),
                                                xyz.getPartialDerivative(i, j, k).getReal(),
                                                1.0e-15);
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testFieldDivide() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).divide(new Decimal64(2)),  0.5, 0.5, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).divide(new Decimal64(2)),  1.0, 0.0, 0.5, 0.0);
            checkF0F1(factory.variable(2, 3.0).divide(new Decimal64(2)),  1.5, 0.0, 0.0, 0.5);
        }
    }

    @Test
    public void testPrimitiveDivide() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).divide(2),  0.5, 0.5, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).divide(2),  1.0, 0.0, 0.5, 0.0);
            checkF0F1(factory.variable(2, 3.0).divide(2),  1.5, 0.0, 0.0, 0.5);
        }
    }

    @Test
    public void testNegate() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).negate(), -1.0, -1.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).negate(), -2.0, 0.0, -1.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).negate(), -3.0, 0.0, 0.0, -1.0);
        }
    }

    @Test
    public void testReciprocal() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, 6);
        for (double x = 0.1; x < 1.2; x += 0.1) {
            FieldDerivativeStructure<Decimal64> r = factory.variable(0, x).reciprocal();
            Assert.assertEquals(1 / x, r.getReal(), 1.0e-15);
            for (int i = 1; i < r.getOrder(); ++i) {
                double expected = ArithmeticUtils.pow(-1, i) * CombinatoricsUtils.factorial(i) /
                                  FastMath.pow(x, i + 1);
                Assert.assertEquals(expected, r.getPartialDerivative(i).getReal(), 1.0e-15 * FastMath.abs(expected));
            }
        }
    }

    @Test
    public void testPow() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            for (int n = 0; n < 10; ++n) {

                FieldDerivativeStructure<Decimal64> x = factory.variable(0, 1.0);
                FieldDerivativeStructure<Decimal64> y = factory.variable(1, 2.0);
                FieldDerivativeStructure<Decimal64> z = factory.variable(2, 3.0);
                List<FieldDerivativeStructure<Decimal64>> list = Arrays.asList(x, y, z,
                                                                               x.add(y).add(z),
                                                                               x.multiply(y).multiply(z));

                if (n == 0) {
                    for (FieldDerivativeStructure<Decimal64> ds : list) {
                        checkEquals(ds.getField().getOne(), ds.pow(n), 1.0e-15);
                    }
                } else if (n == 1) {
                    for (FieldDerivativeStructure<Decimal64> ds : list) {
                        checkEquals(ds, ds.pow(n), 1.0e-15);
                    }
                } else {
                    for (FieldDerivativeStructure<Decimal64> ds : list) {
                        FieldDerivativeStructure<Decimal64> p = ds.getField().getOne();
                        for (int i = 0; i < n; ++i) {
                            p = p.multiply(ds);
                        }
                        checkEquals(p, ds.pow(n), 1.0e-15);
                    }
                }
            }
        }
    }

    @Test
    public void testPowDoubleDS() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {

            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            FieldDerivativeStructure<Decimal64> x = factory.variable(0, 0.1);
            FieldDerivativeStructure<Decimal64> y = factory.variable(1, 0.2);
            FieldDerivativeStructure<Decimal64> z = factory.variable(2, 0.3);
            List<FieldDerivativeStructure<Decimal64>> list = Arrays.asList(x, y, z,
                                                                           x.add(y).add(z),
                                                                           x.multiply(y).multiply(z));

            for (FieldDerivativeStructure<Decimal64> ds : list) {
                // the special case a = 0 is included here
                for (double a : new double[] { 0.0, 0.1, 1.0, 2.0, 5.0 }) {
                    FieldDerivativeStructure<Decimal64> reference = (a == 0) ?
                                                    x.getField().getZero() :
                                                    factory.constant(a).pow(ds);
                    FieldDerivativeStructure<Decimal64> result = FieldDerivativeStructure.pow(a, ds);
                    checkEquals(reference, result, 1.0e-15);
                }

            }

            // negative base: -1^x can be evaluated for integers only, so value is sometimes OK, derivatives are always NaN
            FieldDerivativeStructure<Decimal64> negEvenInteger = FieldDerivativeStructure.pow(-2.0, factory.variable(0, 2.0));
            Assert.assertEquals(4.0, negEvenInteger.getReal(), 1.0e-15);
            Assert.assertTrue(Double.isNaN(negEvenInteger.getPartialDerivative(1, 0, 0).getReal()));
            FieldDerivativeStructure<Decimal64> negOddInteger = FieldDerivativeStructure.pow(-2.0, factory.variable(0, 3.0));
            Assert.assertEquals(-8.0, negOddInteger.getReal(), 1.0e-15);
            Assert.assertTrue(Double.isNaN(negOddInteger.getPartialDerivative(1, 0, 0).getReal()));
            FieldDerivativeStructure<Decimal64> negNonInteger = FieldDerivativeStructure.pow(-2.0, factory.variable(0, 2.001));
            Assert.assertTrue(Double.isNaN(negNonInteger.getReal()));
            Assert.assertTrue(Double.isNaN(negNonInteger.getPartialDerivative(1, 0, 0).getReal()));

            FieldDerivativeStructure<Decimal64> zeroNeg = FieldDerivativeStructure.pow(0.0, factory.variable(0, -1.0));
            Assert.assertTrue(Double.isNaN(zeroNeg.getReal()));
            Assert.assertTrue(Double.isNaN(zeroNeg.getPartialDerivative(1, 0, 0).getReal()));
            FieldDerivativeStructure<Decimal64> posNeg = FieldDerivativeStructure.pow(2.0, factory.variable(0, -2.0));
            Assert.assertEquals(1.0 / 4.0, posNeg.getReal(), 1.0e-15);
            Assert.assertEquals(FastMath.log(2.0) / 4.0, posNeg.getPartialDerivative(1, 0, 0).getReal(), 1.0e-15);

            // very special case: a = 0 and power = 0
            FieldDerivativeStructure<Decimal64> zeroZero = FieldDerivativeStructure.pow(0.0, factory.variable(0, 0.0));

            // this should be OK for simple first derivative with one variable only ...
            Assert.assertEquals(1.0, zeroZero.getReal(), 1.0e-15);
            Assert.assertEquals(Double.NEGATIVE_INFINITY, zeroZero.getPartialDerivative(1, 0, 0).getReal(), 1.0e-15);

            // the following checks show a LIMITATION of the current implementation
            // we have no way to tell x is a pure linear variable x = 0
            // we only say: "x is a structure with value = 0.0,
            // first derivative with respect to x = 1.0, and all other derivatives
            // (first order with respect to y and z and higher derivatives) all 0.0.
            // We have function f(x) = a^x and x = 0 so we compute:
            // f(0) = 1, f'(0) = ln(a), f''(0) = ln(a)^2. The limit of these values
            // when a converges to 0 implies all derivatives keep switching between
            // +infinity and -infinity.
            //
            // Function composition rule for first derivatives is:
            // d[f(g(x,y,z))]/dy = f'(g(x,y,z)) * dg(x,y,z)/dy
            // so given that in our case x represents g and does not depend
            // on y or z, we have dg(x,y,z)/dy = 0
            // applying the composition rules gives:
            // d[f(g(x,y,z))]/dy = f'(g(x,y,z)) * dg(x,y,z)/dy
            //                 = -infinity * 0
            //                 = NaN
            // if we knew x is really the x variable and not the identity
            // function applied to x, we would not have computed f'(g(x,y,z)) * dg(x,y,z)/dy
            // and we would have found that the result was 0 and not NaN
            Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(0, 1, 0).getReal()));
            Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(0, 0, 1).getReal()));

            // Function composition rule for second derivatives is:
            // d2[f(g(x))]/dx2 = f''(g(x)) * [g'(x)]^2 + f'(g(x)) * g''(x)
            // when function f is the a^x root and x = 0 we have:
            // f(0) = 1, f'(0) = ln(a), f''(0) = ln(a)^2 which for a = 0 implies
            // all derivatives keep switching between +infinity and -infinity
            // so given that in our case x represents g, we have g(x) = 0,
            // g'(x) = 1 and g''(x) = 0
            // applying the composition rules gives:
            // d2[f(g(x))]/dx2 = f''(g(x)) * [g'(x)]^2 + f'(g(x)) * g''(x)
            //                 = +infinity * 1^2 + -infinity * 0
            //                 = +infinity + NaN
            //                 = NaN
            // if we knew x is really the x variable and not the identity
            // function applied to x, we would not have computed f'(g(x)) * g''(x)
            // and we would have found that the result was +infinity and not NaN
            if (maxOrder > 1) {
                Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(2, 0, 0).getReal()));
                Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(0, 2, 0).getReal()));
                Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(0, 0, 2).getReal()));
                Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(1, 1, 0).getReal()));
                Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(0, 1, 1).getReal()));
                Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(1, 1, 0).getReal()));
            }

            // very special case: 0^0 where the power is a primitive
            FieldDerivativeStructure<Decimal64> zeroDsZeroDouble = factory.variable(0, 0.0).pow(0.0);
            boolean first = true;
            for (final Decimal64 d : zeroDsZeroDouble.getAllDerivatives()) {
                if (first) {
                    Assert.assertEquals(1.0, d.getReal(), Precision.EPSILON);
                    first = false;
                } else {
                    Assert.assertEquals(0.0, d.getReal(), Precision.SAFE_MIN);
                }
            }
            FieldDerivativeStructure<Decimal64> zeroDsZeroInt = factory.variable(0, 0.0).pow(0);
            first = true;
            for (final Decimal64 d : zeroDsZeroInt.getAllDerivatives()) {
                if (first) {
                    Assert.assertEquals(1.0, d.getReal(), Precision.EPSILON);
                    first = false;
                } else {
                    Assert.assertEquals(0.0, d.getReal(), Precision.SAFE_MIN);
                }
            }

            // 0^p with p smaller than 1.0
            FieldDerivativeStructure<Decimal64> u = factory.variable(1, -0.0).pow(0.25);
            for (int i0 = 0; i0 <= maxOrder; ++i0) {
                for (int i1 = 0; i1 <= maxOrder; ++i1) {
                    for (int i2 = 0; i2 <= maxOrder; ++i2) {
                        if (i0 + i1 + i2 <= maxOrder) {
                            Assert.assertEquals(0.0, u.getPartialDerivative(i0, i1, i2).getReal(), 1.0e-10);
                        }
                    }
                }
            }
        }

    }

    @Test
    public void testExpression() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, 5);
        double epsilon = 2.5e-13;
        for (double x = 0; x < 2; x += 0.2) {
            FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
            for (double y = 0; y < 2; y += 0.2) {
                FieldDerivativeStructure<Decimal64> dsY = factory.variable(1, y);
                for (double z = 0; z >- 2; z -= 0.2) {
                    FieldDerivativeStructure<Decimal64> dsZ = factory.variable(2, z);

                    // f(x, y, z) = x + 5 x y - 2 z + (8 z x - y)^3
                    FieldDerivativeStructure<Decimal64> ds =
                            dsX.linearCombination(1, dsX,
                                                    5, dsX.multiply(dsY),
                                                    -2, dsZ,
                                                    1, dsX.linearCombination(8, dsZ.multiply(dsX),
                                                                               -1, dsY).pow(3));
                    FieldDerivativeStructure<Decimal64> dsOther =
                                    dsX.linearCombination(1, dsX,
                                                    5, dsX.multiply(dsY),
                                                    -2, dsZ).add(dsX.linearCombination   (8, dsZ.multiply(dsX),
                                                                                         -1, dsY).pow(3));
                    double f = x + 5 * x * y - 2 * z + FastMath.pow(8 * z * x - y, 3);
                    Assert.assertEquals(f, ds.getReal(),
                                        FastMath.abs(epsilon * f));
                    Assert.assertEquals(f, dsOther.getReal(),
                                        FastMath.abs(epsilon * f));

                    // df/dx = 1 + 5 y + 24 (8 z x - y)^2 z
                    double dfdx = 1 + 5 * y + 24 * z * FastMath.pow(8 * z * x - y, 2);
                    Assert.assertEquals(dfdx, ds.getPartialDerivative(1, 0, 0).getReal(),
                                        FastMath.abs(epsilon * dfdx));
                    Assert.assertEquals(dfdx, dsOther.getPartialDerivative(1, 0, 0).getReal(),
                                        FastMath.abs(epsilon * dfdx));

                    // df/dxdy = 5 + 48 z*(y - 8 z x)
                    double dfdxdy = 5 + 48 * z * (y - 8 * z * x);
                    Assert.assertEquals(dfdxdy, ds.getPartialDerivative(1, 1, 0).getReal(),
                                        FastMath.abs(epsilon * dfdxdy));
                    Assert.assertEquals(dfdxdy, dsOther.getPartialDerivative(1, 1, 0).getReal(),
                                        FastMath.abs(epsilon * dfdxdy));

                    // df/dxdydz = 48 (y - 16 z x)
                    double dfdxdydz = 48 * (y - 16 * z * x);
                    Assert.assertEquals(dfdxdydz, ds.getPartialDerivative(1, 1, 1).getReal(),
                                        FastMath.abs(epsilon * dfdxdydz));
                    Assert.assertEquals(dfdxdydz, dsOther.getPartialDerivative(1, 1, 1).getReal(),
                                        FastMath.abs(epsilon * dfdxdydz));

                }

            }
        }
    }

    @Test
    public void testCompositionOneVariableX() {
        double epsilon = 1.0e-13;
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.1) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                for (double y = 0.1; y < 1.2; y += 0.1) {
                    FieldDerivativeStructure<Decimal64> dsY = factory.constant(y);
                    FieldDerivativeStructure<Decimal64> f = dsX.divide(dsY).sqrt();
                    double f0 = FastMath.sqrt(x / y);
                    Assert.assertEquals(f0, f.getReal(), FastMath.abs(epsilon * f0));
                    if (f.getOrder() > 0) {
                        double f1 = 1 / (2 * FastMath.sqrt(x * y));
                        Assert.assertEquals(f1, f.getPartialDerivative(1).getReal(), FastMath.abs(epsilon * f1));
                        if (f.getOrder() > 1) {
                            double f2 = -f1 / (2 * x);
                            Assert.assertEquals(f2, f.getPartialDerivative(2).getReal(), FastMath.abs(epsilon * f2));
                            if (f.getOrder() > 2) {
                                double f3 = (f0 + x / (2 * y * f0)) / (4 * x * x * x);
                                Assert.assertEquals(f3, f.getPartialDerivative(3).getReal(), FastMath.abs(epsilon * f3));
                            }
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testTrigo() {
        double epsilon = 2.0e-12;
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.1) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                for (double y = 0.1; y < 1.2; y += 0.1) {
                    FieldDerivativeStructure<Decimal64> dsY = factory.variable(1, y);
                    for (double z = 0.1; z < 1.2; z += 0.1) {
                        FieldDerivativeStructure<Decimal64> dsZ = factory.variable(2, z);
                        FieldDerivativeStructure<Decimal64> f = dsX.divide(dsY.cos().add(dsZ.tan())).sin();
                        double a = FastMath.cos(y) + FastMath.tan(z);
                        double f0 = FastMath.sin(x / a);
                        Assert.assertEquals(f0, f.getReal(), FastMath.abs(epsilon * f0));
                        if (f.getOrder() > 0) {
                            double dfdx = FastMath.cos(x / a) / a;
                            Assert.assertEquals(dfdx, f.getPartialDerivative(1, 0, 0).getReal(), FastMath.abs(epsilon * dfdx));
                            double dfdy =  x * FastMath.sin(y) * dfdx / a;
                            Assert.assertEquals(dfdy, f.getPartialDerivative(0, 1, 0).getReal(), FastMath.abs(epsilon * dfdy));
                            double cz = FastMath.cos(z);
                            double cz2 = cz * cz;
                            double dfdz = -x * dfdx / (a * cz2);
                            Assert.assertEquals(dfdz, f.getPartialDerivative(0, 0, 1).getReal(), FastMath.abs(epsilon * dfdz));
                            if (f.getOrder() > 1) {
                                double df2dx2 = -(f0 / (a * a));
                                Assert.assertEquals(df2dx2, f.getPartialDerivative(2, 0, 0).getReal(), FastMath.abs(epsilon * df2dx2));
                                double df2dy2 = x * FastMath.cos(y) * dfdx / a -
                                                x * x * FastMath.sin(y) * FastMath.sin(y) * f0 / (a * a * a * a) +
                                                2 * FastMath.sin(y) * dfdy / a;
                                Assert.assertEquals(df2dy2, f.getPartialDerivative(0, 2, 0).getReal(), FastMath.abs(epsilon * df2dy2));
                                double c4 = cz2 * cz2;
                                double df2dz2 = x * (2 * a * (1 - a * cz * FastMath.sin(z)) * dfdx - x * f0 / a ) / (a * a * a * c4);
                                Assert.assertEquals(df2dz2, f.getPartialDerivative(0, 0, 2).getReal(), FastMath.abs(epsilon * df2dz2));
                                double df2dxdy = dfdy / x  - x * FastMath.sin(y) * f0 / (a * a * a);
                                Assert.assertEquals(df2dxdy, f.getPartialDerivative(1, 1, 0).getReal(), FastMath.abs(epsilon * df2dxdy));
                            }
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testSqrtDefinition() {
        double[] epsilon = new double[] { 5.0e-16, 5.0e-16, 2.0e-15, 5.0e-14, 2.0e-12 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> sqrt1 = dsX.pow(0.5);
                FieldDerivativeStructure<Decimal64> sqrt2 = dsX.sqrt();
                FieldDerivativeStructure<Decimal64> zero = sqrt1.subtract(sqrt2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testRootNSingularity() {
        for (int n = 2; n < 10; ++n) {
            for (int maxOrder = 0; maxOrder < 12; ++maxOrder) {
                final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
                FieldDerivativeStructure<Decimal64> dsZero = factory.variable(0, 0.0);
                FieldDerivativeStructure<Decimal64> rootN  = dsZero.rootN(n);
                Assert.assertEquals(0.0, rootN.getReal(), 1.0e-20);
                if (maxOrder > 0) {
                    Assert.assertTrue(Double.isInfinite(rootN.getPartialDerivative(1).getReal()));
                    Assert.assertTrue(rootN.getPartialDerivative(1).getReal() > 0);
                    for (int order = 2; order <= maxOrder; ++order) {
                        // the following checks shows a LIMITATION of the current implementation
                        // we have no way to tell dsZero is a pure linear variable x = 0
                        // we only say: "dsZero is a structure with value = 0.0,
                        // first derivative = 1.0, second and higher derivatives = 0.0".
                        // Function composition rule for second derivatives is:
                        // d2[f(g(x))]/dx2 = f''(g(x)) * [g'(x)]^2 + f'(g(x)) * g''(x)
                        // when function f is the nth root and x = 0 we have:
                        // f(0) = 0, f'(0) = +infinity, f''(0) = -infinity (and higher
                        // derivatives keep switching between +infinity and -infinity)
                        // so given that in our case dsZero represents g, we have g(x) = 0,
                        // g'(x) = 1 and g''(x) = 0
                        // applying the composition rules gives:
                        // d2[f(g(x))]/dx2 = f''(g(x)) * [g'(x)]^2 + f'(g(x)) * g''(x)
                        //                 = -infinity * 1^2 + +infinity * 0
                        //                 = -infinity + NaN
                        //                 = NaN
                        // if we knew dsZero is really the x variable and not the identity
                        // function applied to x, we would not have computed f'(g(x)) * g''(x)
                        // and we would have found that the result was -infinity and not NaN
                        Assert.assertTrue(Double.isNaN(rootN.getPartialDerivative(order).getReal()));
                    }
                }

                // the following shows that the limitation explained above is NOT a bug...
                // if we set up the higher order derivatives for g appropriately, we do
                // compute the higher order derivatives of the composition correctly
                double[] gDerivatives = new double[ 1 + maxOrder];
                gDerivatives[0] = 0.0;
                for (int k = 1; k <= maxOrder; ++k) {
                    gDerivatives[k] = FastMath.pow(-1.0, k + 1);
                }
                FieldDerivativeStructure<Decimal64> correctRoot = factory.build(gDerivatives).rootN(n);
                Assert.assertEquals(0.0, correctRoot.getReal(), 1.0e-20);
                if (maxOrder > 0) {
                    Assert.assertTrue(Double.isInfinite(correctRoot.getPartialDerivative(1).getReal()));
                    Assert.assertTrue(correctRoot.getPartialDerivative(1).getReal() > 0);
                    for (int order = 2; order <= maxOrder; ++order) {
                        Assert.assertTrue(Double.isInfinite(correctRoot.getPartialDerivative(order).getReal()));
                        if ((order % 2) == 0) {
                            Assert.assertTrue(correctRoot.getPartialDerivative(order).getReal() < 0);
                        } else {
                            Assert.assertTrue(correctRoot.getPartialDerivative(order).getReal() > 0);
                        }
                    }
                }

            }

        }

    }

    @Test
    public void testSqrtPow2() {
        double[] epsilon = new double[] { 1.0e-16, 3.0e-16, 2.0e-15, 6.0e-14, 6.0e-12 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> rebuiltX = dsX.multiply(dsX).sqrt();
                FieldDerivativeStructure<Decimal64> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testCbrtDefinition() {
        double[] epsilon = new double[] { 4.0e-16, 9.0e-16, 6.0e-15, 2.0e-13, 4.0e-12 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> cbrt1 = dsX.pow(1.0 / 3.0);
                FieldDerivativeStructure<Decimal64> cbrt2 = dsX.cbrt();
                FieldDerivativeStructure<Decimal64> zero = cbrt1.subtract(cbrt2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testCbrtPow3() {
        double[] epsilon = new double[] { 1.0e-16, 5.0e-16, 8.0e-15, 4.0e-13, 3.0e-11 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> rebuiltX = dsX.multiply(dsX.multiply(dsX)).cbrt();
                FieldDerivativeStructure<Decimal64> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testPowReciprocalPow() {
        double[] epsilon = new double[] { 2.0e-15, 2.0e-14, 3.0e-13, 8.0e-12, 3.0e-10 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 2, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.01) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                for (double y = 0.1; y < 1.2; y += 0.01) {
                    FieldDerivativeStructure<Decimal64> dsY = factory.variable(1, y);
                    FieldDerivativeStructure<Decimal64> rebuiltX = dsX.pow(dsY).pow(dsY.reciprocal());
                    FieldDerivativeStructure<Decimal64> zero = rebuiltX.subtract(dsX);
                    for (int n = 0; n <= maxOrder; ++n) {
                        for (int m = 0; m <= maxOrder; ++m) {
                            if (n + m <= maxOrder) {
                                Assert.assertEquals(0.0, zero.getPartialDerivative(n, m).getReal(), epsilon[n + m]);
                            }
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testHypotDefinition() {
        double epsilon = 1.0e-20;
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 2, maxOrder);
            for (double x = -1.7; x < 2; x += 0.2) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                for (double y = -1.7; y < 2; y += 0.2) {
                    FieldDerivativeStructure<Decimal64> dsY = factory.variable(1, y);
                    FieldDerivativeStructure<Decimal64> hypot = FieldDerivativeStructure.hypot(dsY, dsX);
                    FieldDerivativeStructure<Decimal64> ref = dsX.multiply(dsX).add(dsY.multiply(dsY)).sqrt();
                    FieldDerivativeStructure<Decimal64> zero = hypot.subtract(ref);
                    for (int n = 0; n <= maxOrder; ++n) {
                        for (int m = 0; m <= maxOrder; ++m) {
                            if (n + m <= maxOrder) {
                                Assert.assertEquals(0, zero.getPartialDerivative(n, m).getReal(), epsilon);
                            }
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testHypotNoOverflow() {

        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 2, 5);
        FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, +3.0e250);
        FieldDerivativeStructure<Decimal64> dsY = factory.variable(1, -4.0e250);
        FieldDerivativeStructure<Decimal64> hypot = FieldDerivativeStructure.hypot(dsX, dsY);
        Assert.assertEquals(5.0e250, hypot.getReal(), 1.0e235);
        Assert.assertEquals(dsX.getReal() / hypot.getReal(), hypot.getPartialDerivative(1, 0).getReal(), 1.0e-10);
        Assert.assertEquals(dsY.getReal() / hypot.getReal(), hypot.getPartialDerivative(0, 1).getReal(), 1.0e-10);

        FieldDerivativeStructure<Decimal64> sqrt  = dsX.multiply(dsX).add(dsY.multiply(dsY)).sqrt();
        Assert.assertTrue(Double.isInfinite(sqrt.getReal()));

    }

    @Test
    public void testHypotNeglectible() {

        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 2, 5);
        FieldDerivativeStructure<Decimal64> dsSmall = factory.variable(0, +3.0e-10);
        FieldDerivativeStructure<Decimal64> dsLarge = factory.variable(1, -4.0e25);

        Assert.assertEquals(dsLarge.abs().getReal(),
                            FieldDerivativeStructure.hypot(dsSmall, dsLarge).getReal(),
                            1.0e-10);
        Assert.assertEquals(0,
                            FieldDerivativeStructure.hypot(dsSmall, dsLarge).getPartialDerivative(1, 0).getReal(),
                            1.0e-10);
        Assert.assertEquals(-1,
                            FieldDerivativeStructure.hypot(dsSmall, dsLarge).getPartialDerivative(0, 1).getReal(),
                            1.0e-10);

        Assert.assertEquals(dsLarge.abs().getReal(),
                            FieldDerivativeStructure.hypot(dsLarge, dsSmall).getReal(),
                            1.0e-10);
        Assert.assertEquals(0,
                            FieldDerivativeStructure.hypot(dsLarge, dsSmall).getPartialDerivative(1, 0).getReal(),
                            1.0e-10);
        Assert.assertEquals(-1,
                            FieldDerivativeStructure.hypot(dsLarge, dsSmall).getPartialDerivative(0, 1).getReal(),
                            1.0e-10);

    }

    @Test
    public void testHypotSpecial() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 2, 5);
        Assert.assertTrue(Double.isNaN(FieldDerivativeStructure.hypot(factory.variable(0, Double.NaN),
                                                                 factory.variable(0, +3.0e250)).getReal()));
        Assert.assertTrue(Double.isNaN(FieldDerivativeStructure.hypot(factory.variable(0, +3.0e250),
                                                                 factory.variable(0, Double.NaN)).getReal()));
        Assert.assertTrue(Double.isInfinite(FieldDerivativeStructure.hypot(factory.variable(0, Double.POSITIVE_INFINITY),
                                                                      factory.variable(0, +3.0e250)).getReal()));
        Assert.assertTrue(Double.isInfinite(FieldDerivativeStructure.hypot(factory.variable(0, +3.0e250),
                                                                      factory.variable(0, Double.POSITIVE_INFINITY)).getReal()));
    }

    @Test
    public void testFieldRemainder() {
        double epsilon = 1.0e-15;
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 2, maxOrder);
            for (double x = -1.7; x < 2; x += 0.2) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                for (double y = -1.7; y < 2; y += 0.2) {
                    FieldDerivativeStructure<Decimal64> remainder = dsX.remainder(new Decimal64(y));
                    FieldDerivativeStructure<Decimal64> ref = dsX.subtract(x - FastMath.IEEEremainder(x, y));
                    FieldDerivativeStructure<Decimal64> zero = remainder.subtract(ref);
                    for (int n = 0; n <= maxOrder; ++n) {
                        for (int m = 0; m <= maxOrder; ++m) {
                            if (n + m <= maxOrder) {
                                Assert.assertEquals(0, zero.getPartialDerivative(n, m).getReal(), epsilon);
                            }
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testPrimitiveRemainder() {
        double epsilon = 1.0e-15;
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 2, maxOrder);
            for (double x = -1.7; x < 2; x += 0.2) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                for (double y = -1.7; y < 2; y += 0.2) {
                    FieldDerivativeStructure<Decimal64> remainder = dsX.remainder(y);
                    FieldDerivativeStructure<Decimal64> ref = dsX.subtract(x - FastMath.IEEEremainder(x, y));
                    FieldDerivativeStructure<Decimal64> zero = remainder.subtract(ref);
                    for (int n = 0; n <= maxOrder; ++n) {
                        for (int m = 0; m <= maxOrder; ++m) {
                            if (n + m <= maxOrder) {
                                Assert.assertEquals(0, zero.getPartialDerivative(n, m).getReal(), epsilon);
                            }
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testRemainder() {
        double epsilon = 2.0e-15;
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 2, maxOrder);
            for (double x = -1.7; x < 2; x += 0.2) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                for (double y = -1.7; y < 2; y += 0.2) {
                    FieldDerivativeStructure<Decimal64> dsY = factory.variable(1, y);
                    FieldDerivativeStructure<Decimal64> remainder = dsX.remainder(dsY);
                    FieldDerivativeStructure<Decimal64> ref = dsX.subtract(dsY.multiply((x - FastMath.IEEEremainder(x, y)) / y));
                    FieldDerivativeStructure<Decimal64> zero = remainder.subtract(ref);
                    for (int n = 0; n <= maxOrder; ++n) {
                        for (int m = 0; m <= maxOrder; ++m) {
                            if (n + m <= maxOrder) {
                                Assert.assertEquals(0, zero.getPartialDerivative(n, m).getReal(), epsilon);
                            }
                        }
                    }
                }
            }
        }
    }

    @Override
    @Test
    public void testExp() {
        double[] epsilon = new double[] { 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                double refExp = FastMath.exp(x);
                FieldDerivativeStructure<Decimal64> exp = factory.variable(0, x).exp();
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(refExp, exp.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testExpm1Definition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> expm11 = dsX.expm1();
                FieldDerivativeStructure<Decimal64> expm12 = dsX.exp().subtract(dsX.getField().getOne());
                FieldDerivativeStructure<Decimal64> zero = expm11.subtract(expm12);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n).getReal(), epsilon);
                }
            }
        }
    }

    @Override
    @Test
    public void testLog() {
        double[] epsilon = new double[] { 1.0e-16, 1.0e-16, 3.0e-14, 7.0e-13, 3.0e-11 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> log = factory.variable(0, x).log();
                Assert.assertEquals(FastMath.log(x), log.getReal(), epsilon[0]);
                for (int n = 1; n <= maxOrder; ++n) {
                    double refDer = -CombinatoricsUtils.factorial(n - 1) / FastMath.pow(-x, n);
                    Assert.assertEquals(refDer, log.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testLog1pDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> log1p1 = dsX.log1p();
                FieldDerivativeStructure<Decimal64> log1p2 = dsX.add(dsX.getField().getOne()).log();
                FieldDerivativeStructure<Decimal64> zero = log1p1.subtract(log1p2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n).getReal(), epsilon);
                }
            }
        }
    }

    @Test
    public void testLog10Definition() {
        double[] epsilon = new double[] { 3.0e-16, 9.0e-16, 8.0e-15, 3.0e-13, 8.0e-12 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> log101 = dsX.log10();
                FieldDerivativeStructure<Decimal64> log102 = dsX.log().divide(FastMath.log(10.0));
                FieldDerivativeStructure<Decimal64> zero = log101.subtract(log102);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testLogExp() {
        double[] epsilon = new double[] { 2.0e-16, 2.0e-16, 3.0e-16, 2.0e-15, 6.0e-15 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> rebuiltX = dsX.exp().log();
                FieldDerivativeStructure<Decimal64> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testLog1pExpm1() {
        double[] epsilon = new double[] { 6.0e-17, 3.0e-16, 5.0e-16, 9.0e-16, 6.0e-15 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> rebuiltX = dsX.expm1().log1p();
                FieldDerivativeStructure<Decimal64> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testLog10Power() {
        double[] epsilon = new double[] { 3.0e-16, 3.0e-16, 9.0e-16, 6.0e-15, 7.0e-14 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> rebuiltX = factory.constant(10.0).pow(dsX).log10();
                FieldDerivativeStructure<Decimal64> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testSinCosSeparated() {
        double epsilon = 5.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> sin = dsX.sin();
                FieldDerivativeStructure<Decimal64> cos = dsX.cos();
                double s = FastMath.sin(x);
                double c = FastMath.cos(x);
                for (int n = 0; n <= maxOrder; ++n) {
                    switch (n % 4) {
                    case 0 :
                        Assert.assertEquals( s, sin.getPartialDerivative(n).getReal(), epsilon);
                        Assert.assertEquals( c, cos.getPartialDerivative(n).getReal(), epsilon);
                        break;
                    case 1 :
                        Assert.assertEquals( c, sin.getPartialDerivative(n).getReal(), epsilon);
                        Assert.assertEquals(-s, cos.getPartialDerivative(n).getReal(), epsilon);
                        break;
                    case 2 :
                        Assert.assertEquals(-s, sin.getPartialDerivative(n).getReal(), epsilon);
                        Assert.assertEquals(-c, cos.getPartialDerivative(n).getReal(), epsilon);
                        break;
                    default :
                        Assert.assertEquals(-c, sin.getPartialDerivative(n).getReal(), epsilon);
                        Assert.assertEquals( s, cos.getPartialDerivative(n).getReal(), epsilon);
                        break;
                    }
                }
            }
        }
    }

    @Test
    public void testSinCosCombined() {
        double epsilon = 5.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldSinCos<FieldDerivativeStructure<Decimal64>> sinCos = dsX.sinCos();
                double s = FastMath.sin(x);
                double c = FastMath.cos(x);
                for (int n = 0; n <= maxOrder; ++n) {
                    switch (n % 4) {
                    case 0 :
                        Assert.assertEquals( s, sinCos.sin().getPartialDerivative(n).getReal(), epsilon);
                        Assert.assertEquals( c, sinCos.cos().getPartialDerivative(n).getReal(), epsilon);
                        break;
                    case 1 :
                        Assert.assertEquals( c, sinCos.sin().getPartialDerivative(n).getReal(), epsilon);
                        Assert.assertEquals(-s, sinCos.cos().getPartialDerivative(n).getReal(), epsilon);
                        break;
                    case 2 :
                        Assert.assertEquals(-s, sinCos.sin().getPartialDerivative(n).getReal(), epsilon);
                        Assert.assertEquals(-c, sinCos.cos().getPartialDerivative(n).getReal(), epsilon);
                        break;
                    default :
                        Assert.assertEquals(-c, sinCos.sin().getPartialDerivative(n).getReal(), epsilon);
                        Assert.assertEquals( s, sinCos.cos().getPartialDerivative(n).getReal(), epsilon);
                        break;
                    }
                }
            }
        }
    }

    @Test
    public void testSinAsin() {
        double[] epsilon = new double[] { 3.0e-16, 5.0e-16, 3.0e-15, 2.0e-14, 4.0e-13 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> rebuiltX = dsX.sin().asin();
                FieldDerivativeStructure<Decimal64> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testCosAcos() {
        double[] epsilon = new double[] { 6.0e-16, 6.0e-15, 2.0e-13, 4.0e-12, 2.0e-10 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> rebuiltX = dsX.cos().acos();
                FieldDerivativeStructure<Decimal64> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testTanAtan() {
        double[] epsilon = new double[] { 6.0e-17, 2.0e-16, 2.0e-15, 4.0e-14, 2.0e-12 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> rebuiltX = dsX.tan().atan();
                FieldDerivativeStructure<Decimal64> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testTangentDefinition() {
        double[] epsilon = new double[] { 5.0e-16, 2.0e-15, 3.0e-14, 5.0e-13, 2.0e-11 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> tan1 = dsX.sin().divide(dsX.cos());
                FieldDerivativeStructure<Decimal64> tan2 = dsX.tan();
                FieldDerivativeStructure<Decimal64> zero = tan1.subtract(tan2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Override
    @Test
    public void testAtan2() {
        double[] epsilon = new double[] { 5.0e-16, 3.0e-15, 2.2e-14, 1.0e-12, 8.0e-11 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 2, maxOrder);
            for (double x = -1.7; x < 2; x += 0.2) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                for (double y = -1.7; y < 2; y += 0.2) {
                    FieldDerivativeStructure<Decimal64> dsY = factory.variable(1, y);
                    FieldDerivativeStructure<Decimal64> atan2 = FieldDerivativeStructure.atan2(dsY, dsX);
                    FieldDerivativeStructure<Decimal64> ref = dsY.divide(dsX).atan();
                    if (x < 0) {
                        ref = (y < 0) ? ref.subtract(FastMath.PI) : ref.add(FastMath.PI);
                    }
                    FieldDerivativeStructure<Decimal64> zero = atan2.subtract(ref);
                    for (int n = 0; n <= maxOrder; ++n) {
                        for (int m = 0; m <= maxOrder; ++m) {
                            if (n + m <= maxOrder) {
                                Assert.assertEquals(0, zero.getPartialDerivative(n, m).getReal(), epsilon[n + m]);
                            }
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testAtan2SpecialCases() {

        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 2, 2);
        FieldDerivativeStructure<Decimal64> pp =
                FieldDerivativeStructure.atan2(factory.variable(1, new Decimal64(+0.0)), factory.variable(1, new Decimal64(+0.0)));
        Assert.assertEquals(0, pp.getReal(), 1.0e-15);
        Assert.assertEquals(+1, FastMath.copySign(1, pp.getReal()), 1.0e-15);

        FieldDerivativeStructure<Decimal64> pn =
                FieldDerivativeStructure.atan2(factory.variable(1, new Decimal64(+0.0)), factory.variable(1, new Decimal64(-0.0)));
        Assert.assertEquals(FastMath.PI, pn.getReal(), 1.0e-15);

        FieldDerivativeStructure<Decimal64> np =
                FieldDerivativeStructure.atan2(factory.variable(1, new Decimal64(-0.0)), factory.variable(1, new Decimal64(+0.0)));
        Assert.assertEquals(0, np.getReal(), 1.0e-15);
        Assert.assertEquals(-1, FastMath.copySign(1, np.getReal()), 1.0e-15);

        FieldDerivativeStructure<Decimal64> nn =
                FieldDerivativeStructure.atan2(factory.variable(1, new Decimal64(-0.0)), factory.variable(1, new Decimal64(-0.0)));
        Assert.assertEquals(-FastMath.PI, nn.getReal(), 1.0e-15);

    }

    @Test
    public void testSinhDefinition() {
        double[] epsilon = new double[] { 3.0e-16, 3.0e-16, 5.0e-16, 2.0e-15, 6.0e-15 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> sinh1 = dsX.exp().subtract(dsX.exp().reciprocal()).multiply(0.5);
                FieldDerivativeStructure<Decimal64> sinh2 = dsX.sinh();
                FieldDerivativeStructure<Decimal64> zero = sinh1.subtract(sinh2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testCoshDefinition() {
        double[] epsilon = new double[] { 3.0e-16, 3.0e-16, 5.0e-16, 2.0e-15, 6.0e-15 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> cosh1 = dsX.exp().add(dsX.exp().reciprocal()).multiply(0.5);
                FieldDerivativeStructure<Decimal64> cosh2 = dsX.cosh();
                FieldDerivativeStructure<Decimal64> zero = cosh1.subtract(cosh2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testTanhDefinition() {
        double[] epsilon = new double[] { 3.0e-16, 5.0e-16, 7.0e-16, 3.0e-15, 2.0e-14 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> tanh1 = dsX.exp().subtract(dsX.exp().reciprocal()).divide(dsX.exp().add(dsX.exp().reciprocal()));
                FieldDerivativeStructure<Decimal64> tanh2 = dsX.tanh();
                FieldDerivativeStructure<Decimal64> zero = tanh1.subtract(tanh2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testSinhAsinh() {
        double[] epsilon = new double[] { 3.0e-16, 3.0e-16, 4.0e-16, 7.0e-16, 3.0e-15, 8.0e-15 };
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> rebuiltX = dsX.sinh().asinh();
                FieldDerivativeStructure<Decimal64> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testCoshAcosh() {
        double[] epsilon = new double[] { 2.0e-15, 1.0e-14, 2.0e-13, 6.0e-12, 3.0e-10, 2.0e-8 };
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> rebuiltX = dsX.cosh().acosh();
                FieldDerivativeStructure<Decimal64> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testTanhAtanh() {
        double[] epsilon = new double[] { 3.0e-16, 2.0e-16, 7.0e-16, 4.0e-15, 3.0e-14, 4.0e-13 };
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> rebuiltX = dsX.tanh().atanh();
                FieldDerivativeStructure<Decimal64> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testCompositionOneVariableY() {
        double epsilon = 1.0e-13;
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.1) {
                FieldDerivativeStructure<Decimal64> dsX = factory.constant(x);
                for (double y = 0.1; y < 1.2; y += 0.1) {
                    FieldDerivativeStructure<Decimal64> dsY = factory.variable(0, y);
                    FieldDerivativeStructure<Decimal64> f = dsX.divide(dsY).sqrt();
                    double f0 = FastMath.sqrt(x / y);
                    Assert.assertEquals(f0, f.getReal(), FastMath.abs(epsilon * f0));
                    if (f.getOrder() > 0) {
                        double f1 = -x / (2 * y * y * f0);
                        Assert.assertEquals(f1, f.getPartialDerivative(1).getReal(), FastMath.abs(epsilon * f1));
                        if (f.getOrder() > 1) {
                            double f2 = (f0 - x / (4 * y * f0)) / (y * y);
                            Assert.assertEquals(f2, f.getPartialDerivative(2).getReal(), FastMath.abs(epsilon * f2));
                            if (f.getOrder() > 2) {
                                double f3 = (x / (8 * y * f0) - 2 * f0) / (y * y * y);
                                Assert.assertEquals(f3, f.getPartialDerivative(3).getReal(), FastMath.abs(epsilon * f3));
                            }
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testTaylorPrimitivePolynomial() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, 4);
        for (double x = 0; x < 1.2; x += 0.1) {
            FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
            for (double y = 0; y < 1.2; y += 0.2) {
                FieldDerivativeStructure<Decimal64> dsY = factory.variable(1, y);
                for (double z = 0; z < 1.2; z += 0.2) {
                    FieldDerivativeStructure<Decimal64> dsZ = factory.variable(2, z);
                    FieldDerivativeStructure<Decimal64> f = dsX.multiply(dsY).add(dsZ).multiply(dsX).multiply(dsY);
                    for (double dx = -0.2; dx < 0.2; dx += 0.2) {
                        for (double dy = -0.2; dy < 0.2; dy += 0.1) {
                            for (double dz = -0.2; dz < 0.2; dz += 0.1) {
                                double ref = (x + dx) * (y + dy) * ((x + dx) * (y + dy) + (z + dz));
                                Assert.assertEquals(ref, f.taylor(dx, dy, dz).getReal(), 2.0e-15);
                            }
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testTaylorFieldPolynomial() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, 4);
        for (double x = 0; x < 1.2; x += 0.1) {
            FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
            for (double y = 0; y < 1.2; y += 0.2) {
                FieldDerivativeStructure<Decimal64> dsY = factory.variable(1, y);
                for (double z = 0; z < 1.2; z += 0.2) {
                    FieldDerivativeStructure<Decimal64> dsZ = factory.variable(2, z);
                    FieldDerivativeStructure<Decimal64> f = dsX.multiply(dsY).add(dsZ).multiply(dsX).multiply(dsY);
                    for (double dx = -0.2; dx < 0.2; dx += 0.2) {
                        Decimal64 dx64 = new Decimal64(dx);
                        for (double dy = -0.2; dy < 0.2; dy += 0.1) {
                            Decimal64 dy64 = new Decimal64(dy);
                            for (double dz = -0.2; dz < 0.2; dz += 0.1) {
                                Decimal64 dz64 = new Decimal64(dz);
                                double ref = (x + dx) * (y + dy) * ((x + dx) * (y + dy) + (z + dz));
                                Assert.assertEquals(ref, f.taylor(dx64, dy64, dz64).getReal(), 2.0e-15);
                            }
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testTaylorAtan2() {
        double[] expected = new double[] { 0.214, 0.0241, 0.00422, 6.48e-4, 8.04e-5 };
        double x0 =  0.1;
        double y0 = -0.3;
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 2, maxOrder);
            FieldDerivativeStructure<Decimal64> dsX   = factory.variable(0, x0);
            FieldDerivativeStructure<Decimal64> dsY   = factory.variable(1, y0);
            FieldDerivativeStructure<Decimal64> atan2 = FieldDerivativeStructure.atan2(dsY, dsX);
            double maxError = 0;
            for (double dx = -0.05; dx < 0.05; dx += 0.001) {
                for (double dy = -0.05; dy < 0.05; dy += 0.001) {
                    double ref = FastMath.atan2(y0 + dy, x0 + dx);
                    maxError = FastMath.max(maxError, FastMath.abs(ref - atan2.taylor(dx, dy).getReal()));
                }
            }
            Assert.assertEquals(0.0, expected[maxOrder] - maxError, 0.01 * expected[maxOrder]);
        }
    }

    @Override
    @Test
    public void testAbs() {

        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, 1);
        FieldDerivativeStructure<Decimal64> minusOne = factory.variable(0, -1.0);
        Assert.assertEquals(+1.0, minusOne.abs().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.abs().getPartialDerivative(1).getReal(), 1.0e-15);

        FieldDerivativeStructure<Decimal64> plusOne = factory.variable(0, +1.0);
        Assert.assertEquals(+1.0, plusOne.abs().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.abs().getPartialDerivative(1).getReal(), 1.0e-15);

        FieldDerivativeStructure<Decimal64> minusZero = factory.variable(0, new Decimal64(-0.0));
        Assert.assertEquals(+0.0, minusZero.abs().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusZero.abs().getPartialDerivative(1).getReal(), 1.0e-15);

        FieldDerivativeStructure<Decimal64> plusZero = factory.variable(0, new Decimal64(+0.0));
        Assert.assertEquals(+0.0, plusZero.abs().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusZero.abs().getPartialDerivative(1).getReal(), 1.0e-15);

    }

    @Override
    @Test
    public void testSignum() {

        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, 1);
        FieldDerivativeStructure<Decimal64> minusOne = factory.variable(0, -1.0);
        Assert.assertEquals(-1.0, minusOne.signum().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals( 0.0, minusOne.signum().getPartialDerivative(1).getReal(), 1.0e-15);

        FieldDerivativeStructure<Decimal64> plusOne = factory.variable(0, +1.0);
        Assert.assertEquals(+1.0, plusOne.signum().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals( 0.0, plusOne.signum().getPartialDerivative(1).getReal(), 1.0e-15);

        FieldDerivativeStructure<Decimal64> minusZero = factory.variable(0, new Decimal64(-0.0));
        Assert.assertEquals(-0.0, minusZero.signum().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertTrue(Double.doubleToLongBits(minusZero.signum().getReal()) < 0);
        Assert.assertEquals( 0.0, minusZero.signum().getPartialDerivative(1).getReal(), 1.0e-15);

        FieldDerivativeStructure<Decimal64> plusZero = factory.variable(0, new Decimal64(+0.0));
        Assert.assertEquals(+0.0, plusZero.signum().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertTrue(Double.doubleToLongBits(plusZero.signum().getReal()) == 0);
        Assert.assertEquals( 0.0, plusZero.signum().getPartialDerivative(1).getReal(), 1.0e-15);

    }

    @Test
    public void testCeilFloorRintLong() {

        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, 1);
        FieldDerivativeStructure<Decimal64> x = factory.variable(0, -1.5);
        Assert.assertEquals(-1.5, x.getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, x.getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, x.ceil().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+0.0, x.ceil().getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-2.0, x.floor().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+0.0, x.floor().getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-2.0, x.rint().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+0.0, x.rint().getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-2.0, x.subtract(x.getField().getOne()).rint().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1l, x.round());

    }

    @Test
    public void testCopySign() {

        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, 1);
        FieldDerivativeStructure<Decimal64> minusOne = factory.variable(0, -1.0);
        Assert.assertEquals(+1.0, minusOne.copySign(+1.0).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(+1.0).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(-1.0).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(-1.0).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(+0.0).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(+0.0).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(-0.0).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(-0.0).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(Double.NaN).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(Double.NaN).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(new Decimal64(+1.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(new Decimal64(+1.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(new Decimal64(-1.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(new Decimal64(-1.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(new Decimal64(+0.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(new Decimal64(+0.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(new Decimal64(-0.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(new Decimal64(-0.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(new Decimal64(Double.NaN)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(new Decimal64(Double.NaN)).getPartialDerivative(1).getReal(), 1.0e-15);

        FieldDerivativeStructure<Decimal64> plusOne = factory.variable(0, +1.0);
        Assert.assertEquals(+1.0, plusOne.copySign(+1.0).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(+1.0).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, plusOne.copySign(-1.0).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, plusOne.copySign(-1.0).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(+0.0).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(+0.0).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, plusOne.copySign(-0.0).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, plusOne.copySign(-0.0).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(Double.NaN).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(Double.NaN).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(new Decimal64(+1.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(new Decimal64(+1.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, plusOne.copySign(new Decimal64(-1.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, plusOne.copySign(new Decimal64(-1.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(new Decimal64(+0.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(new Decimal64(+0.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, plusOne.copySign(new Decimal64(-0.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, plusOne.copySign(new Decimal64(-0.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(new Decimal64(Double.NaN)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(new Decimal64(Double.NaN)).getPartialDerivative(1).getReal(), 1.0e-15);

    }

    @Test
    public void testToDegreesDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                Assert.assertEquals(FastMath.toDegrees(x), dsX.toDegrees().getReal(), epsilon * FastMath.toDegrees(x));
                for (int n = 1; n <= maxOrder; ++n) {
                    if (n == 1) {
                        Assert.assertEquals(180 / FastMath.PI, dsX.toDegrees().getPartialDerivative(1).getReal(), epsilon);
                    } else {
                        Assert.assertEquals(0.0, dsX.toDegrees().getPartialDerivative(n).getReal(), epsilon);
                    }
                }
            }
        }
    }

    @Test
    public void testToRadiansDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                Assert.assertEquals(FastMath.toRadians(x), dsX.toRadians().getReal(), epsilon);
                for (int n = 1; n <= maxOrder; ++n) {
                    if (n == 1) {
                        Assert.assertEquals(FastMath.PI / 180, dsX.toRadians().getPartialDerivative(1).getReal(), epsilon);
                    } else {
                        Assert.assertEquals(0.0, dsX.toRadians().getPartialDerivative(n).getReal(), epsilon);
                    }
                }
            }
        }
    }

    @Test
    public void testDegRad() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> rebuiltX = dsX.toDegrees().toRadians();
                FieldDerivativeStructure<Decimal64> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon);
                }
            }
        }
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testComposeMismatchedDimensions() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, 3);
        factory.variable(0, 1.2).compose(new double[3]);
    }

    @Test
    public void testComposeField() {
        double[] epsilon = new double[] { 1.0e-20, 5.0e-14, 2.0e-13, 3.0e-13, 2.0e-13, 1.0e-20 };
        PolynomialFunction poly =
                new PolynomialFunction(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 });
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            PolynomialFunction[] p = new PolynomialFunction[maxOrder + 1];
            p[0] = poly;
            for (int i = 1; i <= maxOrder; ++i) {
                p[i] = p[i - 1].polynomialDerivative();
            }
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> dsY1 = dsX.getField().getZero();
                for (int i = poly.degree(); i >= 0; --i) {
                    dsY1 = dsY1.multiply(dsX).add(poly.getCoefficients()[i]);
                }
                Decimal64[] f = new Decimal64[maxOrder + 1];
                for (int i = 0; i < f.length; ++i) {
                    f[i] = new Decimal64(p[i].value(x));
                }
                FieldDerivativeStructure<Decimal64> dsY2 = dsX.compose(f);
                FieldDerivativeStructure<Decimal64> zero = dsY1.subtract(dsY2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testComposePrimitive() {
        double[] epsilon = new double[] { 1.0e-20, 5.0e-14, 2.0e-13, 3.0e-13, 2.0e-13, 1.0e-20 };
        PolynomialFunction poly =
                new PolynomialFunction(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 });
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, maxOrder);
            PolynomialFunction[] p = new PolynomialFunction[maxOrder + 1];
            p[0] = poly;
            for (int i = 1; i <= maxOrder; ++i) {
                p[i] = p[i - 1].polynomialDerivative();
            }
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<Decimal64> dsX = factory.variable(0, x);
                FieldDerivativeStructure<Decimal64> dsY1 = dsX.getField().getZero();
                for (int i = poly.degree(); i >= 0; --i) {
                    dsY1 = dsY1.multiply(dsX).add(poly.getCoefficients()[i]);
                }
                double[] f = new double[maxOrder + 1];
                for (int i = 0; i < f.length; ++i) {
                    f[i] = p[i].value(x);
                }
                FieldDerivativeStructure<Decimal64> dsY2 = dsX.compose(f);
                FieldDerivativeStructure<Decimal64> zero = dsY1.subtract(dsY2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testField() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, maxOrder);
            FieldDerivativeStructure<Decimal64> x = factory.variable(0, 1.0);
            checkF0F1(x.getField().getZero(), 0.0, 0.0, 0.0, 0.0);
            checkF0F1(x.getField().getOne(), 1.0, 0.0, 0.0, 0.0);
            Assert.assertEquals(maxOrder, x.getField().getZero().getOrder());
            Assert.assertEquals(3, x.getField().getZero().getFreeParameters());
            Assert.assertEquals(FieldDerivativeStructure.class, x.getField().getRuntimeClass());
        }
    }

    @Test
    public void testOneParameterConstructor() {
        double x = 1.2;
        double cos = FastMath.cos(x);
        double sin = FastMath.sin(x);
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 1, 4);
        FieldDerivativeStructure<Decimal64> yRef = factory.variable(0, x).cos();
        try {
            factory.build(0.0, 0.0);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            // expected
        } catch (Exception e) {
            Assert.fail("wrong exceptionc caught " + e.getClass().getName());
        }
        double[] derivatives = new double[] { cos, -sin, -cos, sin, cos };
        FieldDerivativeStructure<Decimal64> y = factory.build(derivatives);
        checkEquals(yRef, y, 1.0e-15);
        Decimal64[] all = y.getAllDerivatives();
        Assert.assertEquals(derivatives.length, all.length);
        for (int i = 0; i < all.length; ++i) {
            Assert.assertEquals(derivatives[i], all[i].doubleValue(), 1.0e-15);
        }
    }

    @Test
    public void testOneOrderConstructor() {
        double x =  1.2;
        double y =  2.4;
        double z = 12.5;
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, 1);
        FieldDerivativeStructure<Decimal64> xRef = factory.variable(0, x);
        FieldDerivativeStructure<Decimal64> yRef = factory.variable(1, y);
        FieldDerivativeStructure<Decimal64> zRef = factory.variable(2, z);
        try {
            factory.build(x + y - z, 1.0, 1.0);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            // expected
        } catch (Exception e) {
            Assert.fail("wrong exceptionc caught " + e.getClass().getName());
        }
        double[] derivatives = new double[] { x + y - z, 1.0, 1.0, -1.0 };
        FieldDerivativeStructure<Decimal64> t = factory.build(derivatives);
        checkEquals(xRef.add(yRef.subtract(zRef)), t, 1.0e-15);
        Decimal64[] all = xRef.add(yRef.subtract(zRef)).getAllDerivatives();
        Assert.assertEquals(derivatives.length, all.length);
        for (int i = 0; i < all.length; ++i) {
            Assert.assertEquals(derivatives[i], all[i].doubleValue(), 1.0e-15);
        }
    }

    @Test
    public void testLinearCombination1DSDS() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 6, 1);
        final FieldDerivativeStructure<Decimal64>[] a = MathArrays.buildArray(factory.getDerivativeField(), 3);
        a[0] = factory.variable(0, -1321008684645961.0 / 268435456.0);
        a[1] = factory.variable(1, -5774608829631843.0 / 268435456.0);
        a[2] = factory.variable(2, -7645843051051357.0 / 8589934592.0);
        final FieldDerivativeStructure<Decimal64>[] b = MathArrays.buildArray(factory.getDerivativeField(), 3);
        b[0] = factory.variable(3, -5712344449280879.0 / 2097152.0);
        b[1] = factory.variable(4, -4550117129121957.0 / 2097152.0);
        b[2] = factory.variable(5, 8846951984510141.0 / 131072.0);

        final FieldDerivativeStructure<Decimal64> abSumInline = a[0].linearCombination(a[0], b[0], a[1], b[1], a[2], b[2]);
        final FieldDerivativeStructure<Decimal64> abSumArray = a[0].linearCombination(a, b);

        Assert.assertEquals(abSumInline.getReal(), abSumArray.getReal(), 0);
        Assert.assertEquals(-1.8551294182586248737720779899, abSumInline.getReal(), 1.0e-15);
        Assert.assertEquals(b[0].getReal(), abSumInline.getPartialDerivative(1, 0, 0, 0, 0, 0).getReal(), 1.0e-15);
        Assert.assertEquals(b[1].getReal(), abSumInline.getPartialDerivative(0, 1, 0, 0, 0, 0).getReal(), 1.0e-15);
        Assert.assertEquals(b[2].getReal(), abSumInline.getPartialDerivative(0, 0, 1, 0, 0, 0).getReal(), 1.0e-15);
        Assert.assertEquals(a[0].getReal(), abSumInline.getPartialDerivative(0, 0, 0, 1, 0, 0).getReal(), 1.0e-15);
        Assert.assertEquals(a[1].getReal(), abSumInline.getPartialDerivative(0, 0, 0, 0, 1, 0).getReal(), 1.0e-15);
        Assert.assertEquals(a[2].getReal(), abSumInline.getPartialDerivative(0, 0, 0, 0, 0, 1).getReal(), 1.0e-15);

    }

    @Test
    public void testLinearCombination1FieldDS() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, 1);
        final Decimal64[] a = new Decimal64[] {
            new Decimal64(-1321008684645961.0 / 268435456.0),
            new Decimal64(-5774608829631843.0 / 268435456.0),
            new Decimal64(-7645843051051357.0 / 8589934592.0)
        };
        final FieldDerivativeStructure<Decimal64>[] b = MathArrays.buildArray(factory.getDerivativeField(), 3);
        b[0] = factory.variable(0, -5712344449280879.0 / 2097152.0);
        b[1] = factory.variable(1, -4550117129121957.0 / 2097152.0);
        b[2] = factory.variable(2, 8846951984510141.0 / 131072.0);

        final FieldDerivativeStructure<Decimal64> abSumInline = b[0].linearCombination(a[0], b[0],
                                                                                       a[1], b[1],
                                                                                       a[2], b[2]);
        final FieldDerivativeStructure<Decimal64> abSumArray = b[0].linearCombination(a, b);

        Assert.assertEquals(abSumInline.getReal(), abSumArray.getReal(), 0);
        Assert.assertEquals(-1.8551294182586248737720779899, abSumInline.getReal(), 1.0e-15);
        Assert.assertEquals(a[0].getReal(), abSumInline.getPartialDerivative(1, 0, 0).getReal(), 1.0e-15);
        Assert.assertEquals(a[1].getReal(), abSumInline.getPartialDerivative(0, 1, 0).getReal(), 1.0e-15);
        Assert.assertEquals(a[2].getReal(), abSumInline.getPartialDerivative(0, 0, 1).getReal(), 1.0e-15);

    }

    @Test
    public void testLinearCombination1DoubleDS() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, 1);
        final double[] a = new double[] {
            -1321008684645961.0 / 268435456.0,
            -5774608829631843.0 / 268435456.0,
            -7645843051051357.0 / 8589934592.0
        };
        final FieldDerivativeStructure<Decimal64>[] b = MathArrays.buildArray(factory.getDerivativeField(), 3);
        b[0] = factory.variable(0, -5712344449280879.0 / 2097152.0);
        b[1] = factory.variable(1, -4550117129121957.0 / 2097152.0);
        b[2] = factory.variable(2, 8846951984510141.0 / 131072.0);

        final FieldDerivativeStructure<Decimal64> abSumInline = b[0].linearCombination(a[0], b[0],
                                                                       a[1], b[1],
                                                                       a[2], b[2]);
        final FieldDerivativeStructure<Decimal64> abSumArray = b[0].linearCombination(a, b);

        Assert.assertEquals(abSumInline.getReal(), abSumArray.getReal(), 0);
        Assert.assertEquals(-1.8551294182586248737720779899, abSumInline.getReal(), 1.0e-15);
        Assert.assertEquals(a[0], abSumInline.getPartialDerivative(1, 0, 0).getReal(), 1.0e-15);
        Assert.assertEquals(a[1], abSumInline.getPartialDerivative(0, 1, 0).getReal(), 1.0e-15);
        Assert.assertEquals(a[2], abSumInline.getPartialDerivative(0, 0, 1).getReal(), 1.0e-15);

    }

    @Test
    public void testLinearCombination2DSDS() {

        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 4, 1);
        final FieldDerivativeStructure<Decimal64>[] u = MathArrays.buildArray(factory.getDerivativeField(), 4);
        final FieldDerivativeStructure<Decimal64>[] v = MathArrays.buildArray(factory.getDerivativeField(), 4);

        // we compare accurate versus naive dot product implementations
        // on regular vectors (i.e. not extreme cases like in the previous test)
        Well1024a random = new Well1024a(0xc6af886975069f11l);

        for (int i = 0; i < 10000; ++i) {
            for (int j = 0; j < u.length; ++j) {
                u[j] = factory.variable(j, 1e17 * random.nextDouble());
                v[j] = factory.constant(1e17 * random.nextDouble());
            }

            FieldDerivativeStructure<Decimal64> lin = u[0].linearCombination(u[0], v[0], u[1], v[1]);
            double ref = u[0].getReal() * v[0].getReal() +
                         u[1].getReal() * v[1].getReal();
            Assert.assertEquals(ref, lin.getReal(), 1.0e-15 * FastMath.abs(ref));
            Assert.assertEquals(v[0].getReal(), lin.getPartialDerivative(1, 0, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[0].getReal()));
            Assert.assertEquals(v[1].getReal(), lin.getPartialDerivative(0, 1, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[1].getReal()));

            lin = u[0].linearCombination(u[0], v[0], u[1], v[1], u[2], v[2]);
            ref = u[0].getReal() * v[0].getReal() +
                  u[1].getReal() * v[1].getReal() +
                  u[2].getReal() * v[2].getReal();
            Assert.assertEquals(ref, lin.getReal(), 1.0e-15 * FastMath.abs(ref));
            Assert.assertEquals(v[0].getReal(), lin.getPartialDerivative(1, 0, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[0].getReal()));
            Assert.assertEquals(v[1].getReal(), lin.getPartialDerivative(0, 1, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[1].getReal()));
            Assert.assertEquals(v[2].getReal(), lin.getPartialDerivative(0, 0, 1, 0).getReal(), 1.0e-15 * FastMath.abs(v[2].getReal()));

            lin = u[0].linearCombination(u[0], v[0], u[1], v[1], u[2], v[2], u[3], v[3]);
            ref = u[0].getReal() * v[0].getReal() +
                  u[1].getReal() * v[1].getReal() +
                  u[2].getReal() * v[2].getReal() +
                  u[3].getReal() * v[3].getReal();
            Assert.assertEquals(ref, lin.getReal(), 1.0e-15 * FastMath.abs(ref));
            Assert.assertEquals(v[0].getReal(), lin.getPartialDerivative(1, 0, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[0].getReal()));
            Assert.assertEquals(v[1].getReal(), lin.getPartialDerivative(0, 1, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[1].getReal()));
            Assert.assertEquals(v[2].getReal(), lin.getPartialDerivative(0, 0, 1, 0).getReal(), 1.0e-15 * FastMath.abs(v[2].getReal()));
            Assert.assertEquals(v[3].getReal(), lin.getPartialDerivative(0, 0, 0, 1).getReal(), 1.0e-15 * FastMath.abs(v[3].getReal()));

        }
    }

    @Test
    public void testLinearCombination2DoubleDS() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 4, 1);
        final double[] u = new double[4];
        final FieldDerivativeStructure<Decimal64>[] v = MathArrays.buildArray(factory.getDerivativeField(), 4);
        // we compare accurate versus naive dot product implementations
        // on regular vectors (i.e. not extreme cases like in the previous test)
        Well1024a random = new Well1024a(0xc6af886975069f11l);

        for (int i = 0; i < 10000; ++i) {
            for (int j = 0; j < u.length; ++j) {
                u[j] = 1e17 * random.nextDouble();
                v[j] = factory.variable(j, 1e17 * random.nextDouble());
            }

            FieldDerivativeStructure<Decimal64> lin = v[0].linearCombination(u[0], v[0], u[1], v[1]);
            double ref = u[0] * v[0].getReal() +
                         u[1] * v[1].getReal();
            Assert.assertEquals(ref, lin.getReal(), 1.0e-15 * FastMath.abs(ref));
            Assert.assertEquals(u[0], lin.getPartialDerivative(1, 0, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[0].getReal()));
            Assert.assertEquals(u[1], lin.getPartialDerivative(0, 1, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[1].getReal()));

            lin = v[0].linearCombination(u[0], v[0], u[1], v[1], u[2], v[2]);
            ref = u[0] * v[0].getReal() +
                  u[1] * v[1].getReal() +
                  u[2] * v[2].getReal();
            Assert.assertEquals(ref, lin.getReal(), 1.0e-15 * FastMath.abs(ref));
            Assert.assertEquals(u[0], lin.getPartialDerivative(1, 0, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[0].getReal()));
            Assert.assertEquals(u[1], lin.getPartialDerivative(0, 1, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[1].getReal()));
            Assert.assertEquals(u[2], lin.getPartialDerivative(0, 0, 1, 0).getReal(), 1.0e-15 * FastMath.abs(v[2].getReal()));

            lin = v[0].linearCombination(u[0], v[0], u[1], v[1], u[2], v[2], u[3], v[3]);
            ref = u[0] * v[0].getReal() +
                  u[1] * v[1].getReal() +
                  u[2] * v[2].getReal() +
                  u[3] * v[3].getReal();
            Assert.assertEquals(ref, lin.getReal(), 1.0e-15 * FastMath.abs(ref));
            Assert.assertEquals(u[0], lin.getPartialDerivative(1, 0, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[0].getReal()));
            Assert.assertEquals(u[1], lin.getPartialDerivative(0, 1, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[1].getReal()));
            Assert.assertEquals(u[2], lin.getPartialDerivative(0, 0, 1, 0).getReal(), 1.0e-15 * FastMath.abs(v[2].getReal()));
            Assert.assertEquals(u[3], lin.getPartialDerivative(0, 0, 0, 1).getReal(), 1.0e-15 * FastMath.abs(v[3].getReal()));

        }
    }

    @Test
    public void testLinearCombination2FieldDS() {
        final FDSFactory<Decimal64> factory = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 4, 1);
        final Decimal64[] u = new Decimal64[4];
        final FieldDerivativeStructure<Decimal64>[] v = MathArrays.buildArray(factory.getDerivativeField(), 4);
        // we compare accurate versus naive dot product implementations
        // on regular vectors (i.e. not extreme cases like in the previous test)
        Well1024a random = new Well1024a(0xc6af886975069f11l);

        for (int i = 0; i < 10000; ++i) {
            for (int j = 0; j < u.length; ++j) {
                u[j] = new Decimal64(1e17 * random.nextDouble());
                v[j] = factory.variable(j, 1e17 * random.nextDouble());
            }

            FieldDerivativeStructure<Decimal64> lin = v[0].linearCombination(u[0], v[0], u[1], v[1]);
            double ref = u[0].getReal() * v[0].getReal() +
                         u[1].getReal() * v[1].getReal();
            Assert.assertEquals(ref, lin.getReal(), 1.0e-15 * FastMath.abs(ref));
            Assert.assertEquals(u[0].getReal(), lin.getPartialDerivative(1, 0, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[0].getReal()));
            Assert.assertEquals(u[1].getReal(), lin.getPartialDerivative(0, 1, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[1].getReal()));

            lin = v[0].linearCombination(u[0], v[0], u[1], v[1], u[2], v[2]);
            ref = u[0].getReal() * v[0].getReal() +
                  u[1].getReal() * v[1].getReal() +
                  u[2].getReal() * v[2].getReal();
            Assert.assertEquals(ref, lin.getReal(), 1.0e-15 * FastMath.abs(ref));
            Assert.assertEquals(u[0].getReal(), lin.getPartialDerivative(1, 0, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[0].getReal()));
            Assert.assertEquals(u[1].getReal(), lin.getPartialDerivative(0, 1, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[1].getReal()));
            Assert.assertEquals(u[2].getReal(), lin.getPartialDerivative(0, 0, 1, 0).getReal(), 1.0e-15 * FastMath.abs(v[2].getReal()));

            lin = v[0].linearCombination(u[0], v[0], u[1], v[1], u[2], v[2], u[3], v[3]);
            ref = u[0].getReal() * v[0].getReal() +
                  u[1].getReal() * v[1].getReal() +
                  u[2].getReal() * v[2].getReal() +
                  u[3].getReal() * v[3].getReal();
            Assert.assertEquals(ref, lin.getReal(), 1.0e-15 * FastMath.abs(ref));
            Assert.assertEquals(u[0].getReal(), lin.getPartialDerivative(1, 0, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[0].getReal()));
            Assert.assertEquals(u[1].getReal(), lin.getPartialDerivative(0, 1, 0, 0).getReal(), 1.0e-15 * FastMath.abs(v[1].getReal()));
            Assert.assertEquals(u[2].getReal(), lin.getPartialDerivative(0, 0, 1, 0).getReal(), 1.0e-15 * FastMath.abs(v[2].getReal()));
            Assert.assertEquals(u[3].getReal(), lin.getPartialDerivative(0, 0, 0, 1).getReal(), 1.0e-15 * FastMath.abs(v[3].getReal()));

        }
    }

    @Test
    public void testZero() {
        FDSFactory<Decimal64> factory = new FDSFactory<>(Decimal64Field.getInstance(), 3, 2);
        FieldDerivativeStructure<Decimal64> zero = factory.constant(17).getField().getZero();
        Decimal64[] a = zero.getAllDerivatives();
        Assert.assertEquals(10, a.length);
        for (int i = 0; i < a.length; ++i) {
            Assert.assertEquals(Decimal64.ZERO, a[i]);
        }
    }

    @Test
    public void testOne() {
        FDSFactory<Decimal64> factory = new FDSFactory<>(Decimal64Field.getInstance(), 3, 2);
        FieldDerivativeStructure<Decimal64> one = factory.constant(17).getField().getOne();
        Decimal64[] a = one.getAllDerivatives();
        Assert.assertEquals(10, a.length);
        for (int i = 0; i < a.length; ++i) {
            Assert.assertEquals(i == 0 ? Decimal64.ONE : Decimal64.ZERO, a[i]);
        }
    }

    @Test
    public void testMap() {
        List<int[]> pairs = new ArrayList<>();
        for (int parameters = 1; parameters < 5; ++parameters) {
            for (int order = 0; order < 3; ++order) {
                pairs.add(new int[] { parameters, order });
            }
        }
        Map<Field<?>, Integer> map = new HashMap<>();
        for (int i = 0; i < 1000; ++i) {
            // create a brand new factory for each derivative
            int parameters = pairs.get(i % pairs.size())[0];
            int order      = pairs.get(i % pairs.size())[1];
            FDSFactory<Decimal64> factory = new FDSFactory<>(Decimal64Field.getInstance(), parameters, order);
            map.put(factory.constant(new Decimal64(17)).getField(), 0);
        }

        // despite we have created numerous factories,
        // there should be only one field for each pair parameters/order
        Assert.assertEquals(pairs.size(), map.size());
        @SuppressWarnings("unchecked")
        Field<FieldDerivativeStructure<Decimal64>> first = (Field<FieldDerivativeStructure<Decimal64>>) map.entrySet().iterator().next().getKey();
        Assert.assertTrue(first.equals(first));
        Assert.assertFalse(first.equals(Decimal64Field.getInstance()));

        // even at same parameters and differentiation orders, different values generate different fields
        FieldDerivativeStructure<Decimal64> zero64 = new FDSFactory<Decimal64>(Decimal64Field.getInstance(), 3, 2).build();
        FieldDerivativeStructure<Dfp> zeroDFP = new FDSFactory<Dfp>(new DfpField(15), 3, 2).build();
        Assert.assertEquals(zero64.getFreeParameters(), zeroDFP.getFreeParameters());
        Assert.assertEquals(zero64.getOrder(), zeroDFP.getOrder());
        Assert.assertFalse(zero64.getField().equals(zeroDFP.getField()));
    }

    @Test
    public void testRunTimeClass() {
        FDSFactory<Decimal64> factory = new FDSFactory<>(Decimal64Field.getInstance(), 3, 2);
        Field<FieldDerivativeStructure<Decimal64>> field = factory.getDerivativeField();
        Assert.assertEquals(FieldDerivativeStructure.class, field.getRuntimeClass());
        Assert.assertEquals(Decimal64Field.getInstance(), factory.getValueField());
        Assert.assertEquals("org.hipparchus.analysis.differentiation.FDSFactory$DerivativeField",
                            factory.getDerivativeField().getClass().getName());
    }

    private void checkF0F1(FieldDerivativeStructure<Decimal64> ds, double value, double...derivatives) {

        // check dimension
        Assert.assertEquals(derivatives.length, ds.getFreeParameters());

        // check value, directly and also as 0th order derivative
        Assert.assertEquals(value, ds.getReal(), 1.0e-15);
        Assert.assertEquals(value, ds.getPartialDerivative(new int[ds.getFreeParameters()]).getReal(), 1.0e-15);

        // check first order derivatives
        for (int i = 0; i < derivatives.length; ++i) {
            int[] orders = new int[derivatives.length];
            orders[i] = 1;
            Assert.assertEquals(derivatives[i], ds.getPartialDerivative(orders).getReal(), 1.0e-15);
        }

    }

    private void checkEquals(FieldDerivativeStructure<Decimal64> ds1, FieldDerivativeStructure<Decimal64> ds2, double epsilon) {

        // check dimension
        Assert.assertEquals(ds1.getFreeParameters(), ds2.getFreeParameters());
        Assert.assertEquals(ds1.getOrder(), ds2.getOrder());

        int[] derivatives = new int[ds1.getFreeParameters()];
        int sum = 0;
        while (true) {

            if (sum <= ds1.getOrder()) {
                Assert.assertEquals(ds1.getPartialDerivative(derivatives).getReal(),
                                    ds2.getPartialDerivative(derivatives).getReal(),
                                    epsilon);
            }

            boolean increment = true;
            sum = 0;
            for (int i = derivatives.length - 1; i >= 0; --i) {
                if (increment) {
                    if (derivatives[i] == ds1.getOrder()) {
                        derivatives[i] = 0;
                    } else {
                        derivatives[i]++;
                        increment = false;
                    }
                }
                sum += derivatives[i];
            }
            if (increment) {
                return;
            }

        }

    }

}
