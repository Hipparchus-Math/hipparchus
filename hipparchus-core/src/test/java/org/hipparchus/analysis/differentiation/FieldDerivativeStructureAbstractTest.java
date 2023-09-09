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

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.CalculusFieldElementAbstractTest;
import org.hipparchus.Field;
import org.hipparchus.analysis.CalculusFieldMultivariateFunction;
import org.hipparchus.analysis.CalculusFieldMultivariateVectorFunction;
import org.hipparchus.analysis.polynomials.FieldPolynomialFunction;
import org.hipparchus.analysis.polynomials.PolynomialFunction;
import org.hipparchus.dfp.Dfp;
import org.hipparchus.dfp.DfpField;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.random.Well19937a;
import org.hipparchus.util.ArithmeticUtils;
import org.hipparchus.util.CombinatoricsUtils;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.FieldSinhCosh;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.Precision;
import org.junit.Assert;
import org.junit.Test;

/**
 * Abstract test for class {@link FieldDerivativeStructure}.
 */
public abstract class FieldDerivativeStructureAbstractTest<T extends CalculusFieldElement<T>>
    extends CalculusFieldElementAbstractTest<FieldDerivativeStructure<T>> {

    protected abstract Field<T> getField();

    protected T buildScalar(double value) {
        return getField().getZero().newInstance(value);
    }

    protected FDSFactory<T> buildFactory(int parameters, int order) {
        return new FDSFactory<>(getField(), parameters, order);
    }

    @Override
    protected FieldDerivativeStructure<T> build(final double x) {
        return buildFactory(2, 1).variable(0, x);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testWrongFieldVariableIndex() {
        buildFactory(3, 1).variable(3, buildScalar(1.0));
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testWrongPrimitiveVariableIndex() {
        final FDSFactory<T> factory = buildFactory(3, 1);
        factory.variable(3, 1.0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testMissingOrders() {
        final FDSFactory<T> factory = buildFactory(3, 1);
        factory.variable(0, 1.0).getPartialDerivative(0, 1);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testWrongDimensionField() {
        final FDSFactory<T> factory = buildFactory(3, 1);
        factory.build(buildScalar(1.0), buildScalar(1.0), buildScalar(1.0), buildScalar(1.0), buildScalar(1.0));
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testWrongDimensionPrimitive() {
        final FDSFactory<T> factory = buildFactory(3, 1);
        factory.build(1.0, 1.0, 1.0, 1.0, 1.0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testTooLargeOrder() {
        final FDSFactory<T> factory = buildFactory(3, 1);
        factory.variable(0, 1.0).getPartialDerivative(1, 1, 2);
    }

    @Test
    public void testVariableWithoutDerivativeField() {
        final FDSFactory<T> factory = buildFactory(1, 0);
        FieldDerivativeStructure<T> v = factory.variable(0, buildScalar(1.0));
        Assert.assertEquals(1.0, v.getReal(), 1.0e-15);
    }

    @Test
    public void testVariableWithoutDerivativePrimitive() {
        final FDSFactory<T> factory = buildFactory(1, 0);
        FieldDerivativeStructure<T> v = factory.variable(0, 1.0);
        Assert.assertEquals(1.0, v.getReal(), 1.0e-15);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testVariableWithoutDerivative1() {
        final FDSFactory<T> factory = buildFactory(1, 0);
        FieldDerivativeStructure<T> v = factory.variable(0, 1.0);
        Assert.assertEquals(1.0, v.getPartialDerivative(1).getReal(), 1.0e-15);
    }

    @Test
    public void testVariable() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
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
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            checkF0F1(factory.constant(FastMath.PI),
                      FastMath.PI, 0.0, 0.0, 0.0);
        }
    }

    @Test
    public void testFieldAdd() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).add(buildScalar(5)), 6.0, 1.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).add(buildScalar(5)), 7.0, 0.0, 1.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).add(buildScalar(5)), 8.0, 0.0, 0.0, 1.0);
        }
    }

    @Test
    public void testPrimitiveAdd() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).add(5), 6.0, 1.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).add(5), 7.0, 0.0, 1.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).add(5), 8.0, 0.0, 0.0, 1.0);
        }
    }

    @Test
    public void testAdd() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            FieldDerivativeStructure<T> x = factory.variable(0, 1.0);
            FieldDerivativeStructure<T> y = factory.variable(1, 2.0);
            FieldDerivativeStructure<T> z = factory.variable(2, 3.0);
            FieldDerivativeStructure<T> xyz = x.add(y.add(z));
            checkF0F1(xyz, x.getValue().getReal() + y.getValue().getReal() + z.getValue().getReal(), 1.0, 1.0, 1.0);
        }
    }

    @Test
    public void testFieldSubtract() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).subtract(buildScalar(5)), -4.0, 1.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).subtract(buildScalar(5)), -3.0, 0.0, 1.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).subtract(buildScalar(5)), -2.0, 0.0, 0.0, 1.0);
        }
    }

    @Test
    public void testPrimitiveSubtract() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).subtract(5), -4.0, 1.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).subtract(5), -3.0, 0.0, 1.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).subtract(5), -2.0, 0.0, 0.0, 1.0);
        }
    }

    @Test
    public void testSubtract() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            FieldDerivativeStructure<T> x = factory.variable(0, 1.0);
            FieldDerivativeStructure<T> y = factory.variable(1, 2.0);
            FieldDerivativeStructure<T> z = factory.variable(2, 3.0);
            FieldDerivativeStructure<T> xyz = x.subtract(y.subtract(z));
            checkF0F1(xyz, x.getReal() - (y.getReal() - z.getReal()), 1.0, -1.0, 1.0);
        }
    }

    @Test
    public void testFieldMultiply() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).multiply(buildScalar(5)),  5.0, 5.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).multiply(buildScalar(5)), 10.0, 0.0, 5.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).multiply(buildScalar(5)), 15.0, 0.0, 0.0, 5.0);
        }
    }

    @Test
    public void testPrimitiveMultiply() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).multiply(5),  5.0, 5.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).multiply(5), 10.0, 0.0, 5.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).multiply(5), 15.0, 0.0, 0.0, 5.0);
        }
    }

    @Test
    public void testMultiply() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            FieldDerivativeStructure<T> x = factory.variable(0, 1.0);
            FieldDerivativeStructure<T> y = factory.variable(1, 2.0);
            FieldDerivativeStructure<T> z = factory.variable(2, 3.0);
            FieldDerivativeStructure<T> xyz = x.multiply(y.multiply(z));
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
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).divide(buildScalar(2)),  0.5, 0.5, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).divide(buildScalar(2)),  1.0, 0.0, 0.5, 0.0);
            checkF0F1(factory.variable(2, 3.0).divide(buildScalar(2)),  1.5, 0.0, 0.0, 0.5);
        }
    }

    @Test
    public void testPrimitiveDivide() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).divide(2),  0.5, 0.5, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).divide(2),  1.0, 0.0, 0.5, 0.0);
            checkF0F1(factory.variable(2, 3.0).divide(2),  1.5, 0.0, 0.0, 0.5);
        }
    }

    @Test
    public void testNegate() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).negate(), -1.0, -1.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).negate(), -2.0, 0.0, -1.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).negate(), -3.0, 0.0, 0.0, -1.0);
        }
    }

    @Test
    public void testReciprocal() {
        final FDSFactory<T> factory = buildFactory(1, 6);
        for (double x = 0.1; x < 1.2; x += 0.1) {
            FieldDerivativeStructure<T> r = factory.variable(0, x).reciprocal();
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
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            for (int n = 0; n < 10; ++n) {

                FieldDerivativeStructure<T> x = factory.variable(0, 1.0);
                FieldDerivativeStructure<T> y = factory.variable(1, 2.0);
                FieldDerivativeStructure<T> z = factory.variable(2, 3.0);
                List<FieldDerivativeStructure<T>> list = Arrays.asList(x, y, z,
                                                                               x.add(y).add(z),
                                                                               x.multiply(y).multiply(z));

                if (n == 0) {
                    for (FieldDerivativeStructure<T> ds : list) {
                        checkEquals(ds.getField().getOne(), ds.pow(n), 1.0e-15);
                    }
                } else if (n == 1) {
                    for (FieldDerivativeStructure<T> ds : list) {
                        checkEquals(ds, ds.pow(n), 1.0e-15);
                    }
                } else {
                    for (FieldDerivativeStructure<T> ds : list) {
                        FieldDerivativeStructure<T> p = ds.getField().getOne();
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

            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            FieldDerivativeStructure<T> x = factory.variable(0, 0.1);
            FieldDerivativeStructure<T> y = factory.variable(1, 0.2);
            FieldDerivativeStructure<T> z = factory.variable(2, 0.3);
            List<FieldDerivativeStructure<T>> list = Arrays.asList(x, y, z,
                                                                           x.add(y).add(z),
                                                                           x.multiply(y).multiply(z));

            for (FieldDerivativeStructure<T> ds : list) {
                // the special case a = 0 is included here
                for (double a : new double[] { 0.0, 0.1, 1.0, 2.0, 5.0 }) {
                    FieldDerivativeStructure<T> reference = (a == 0) ?
                                                    x.getField().getZero() :
                                                    factory.constant(a).pow(ds);
                    FieldDerivativeStructure<T> result = FieldDerivativeStructure.pow(a, ds);
                    checkEquals(reference, result, 2.0e-14 * FastMath.abs(reference.getReal()));
                }

            }

            // negative base: -1^x can be evaluated for integers only, so value is sometimes OK, derivatives are always NaN
            FieldDerivativeStructure<T> negEvenInteger = FieldDerivativeStructure.pow(-2.0, factory.variable(0, 2.0));
            Assert.assertEquals(4.0, negEvenInteger.getReal(), 1.0e-15);
            Assert.assertTrue(Double.isNaN(negEvenInteger.getPartialDerivative(1, 0, 0).getReal()));
            FieldDerivativeStructure<T> negOddInteger = FieldDerivativeStructure.pow(-2.0, factory.variable(0, 3.0));
            Assert.assertEquals(-8.0, negOddInteger.getReal(), 1.0e-15);
            Assert.assertTrue(Double.isNaN(negOddInteger.getPartialDerivative(1, 0, 0).getReal()));
            FieldDerivativeStructure<T> negNonInteger = FieldDerivativeStructure.pow(-2.0, factory.variable(0, 2.001));
            Assert.assertTrue(Double.isNaN(negNonInteger.getReal()));
            Assert.assertTrue(Double.isNaN(negNonInteger.getPartialDerivative(1, 0, 0).getReal()));

            FieldDerivativeStructure<T> zeroNeg = FieldDerivativeStructure.pow(0.0, factory.variable(0, -1.0));
            Assert.assertTrue(Double.isNaN(zeroNeg.getReal()));
            Assert.assertTrue(Double.isNaN(zeroNeg.getPartialDerivative(1, 0, 0).getReal()));
            FieldDerivativeStructure<T> posNeg = FieldDerivativeStructure.pow(2.0, factory.variable(0, -2.0));
            Assert.assertEquals(1.0 / 4.0, posNeg.getReal(), 1.0e-15);
            Assert.assertEquals(FastMath.log(2.0) / 4.0, posNeg.getPartialDerivative(1, 0, 0).getReal(), 1.0e-15);

            // very special case: a = 0 and power = 0
            FieldDerivativeStructure<T> zeroZero = FieldDerivativeStructure.pow(0.0, factory.variable(0, 0.0));

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
            FieldDerivativeStructure<T> zeroDsZeroDouble = factory.variable(0, 0.0).pow(0.0);
            boolean first = true;
            for (final T d : zeroDsZeroDouble.getAllDerivatives()) {
                if (first) {
                    Assert.assertEquals(1.0, d.getReal(), Precision.EPSILON);
                    first = false;
                } else {
                    Assert.assertEquals(0.0, d.getReal(), Precision.SAFE_MIN);
                }
            }
            FieldDerivativeStructure<T> zeroDsZeroInt = factory.variable(0, 0.0).pow(0);
            first = true;
            for (final T d : zeroDsZeroInt.getAllDerivatives()) {
                if (first) {
                    Assert.assertEquals(1.0, d.getReal(), Precision.EPSILON);
                    first = false;
                } else {
                    Assert.assertEquals(0.0, d.getReal(), Precision.SAFE_MIN);
                }
            }

            // 0^p with p smaller than 1.0
            FieldDerivativeStructure<T> u = factory.variable(1, -0.0).pow(0.25);
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
        final FDSFactory<T> factory = buildFactory(3, 5);
        double epsilon = 2.5e-13;
        for (double x = 0; x < 2; x += 0.2) {
            FieldDerivativeStructure<T> dsX = factory.variable(0, x);
            for (double y = 0; y < 2; y += 0.2) {
                FieldDerivativeStructure<T> dsY = factory.variable(1, y);
                for (double z = 0; z >- 2; z -= 0.2) {
                    FieldDerivativeStructure<T> dsZ = factory.variable(2, z);

                    // f(x, y, z) = x + 5 x y - 2 z + (8 z x - y)^3
                    FieldDerivativeStructure<T> ds =
                            dsX.linearCombination(1, dsX,
                                                    5, dsX.multiply(dsY),
                                                    -2, dsZ,
                                                    1, dsX.linearCombination(8, dsZ.multiply(dsX),
                                                                               -1, dsY).pow(3));
                    FieldDerivativeStructure<T> dsOther =
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.1) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                for (double y = 0.1; y < 1.2; y += 0.1) {
                    FieldDerivativeStructure<T> dsY = factory.constant(y);
                    FieldDerivativeStructure<T> f = dsX.divide(dsY).sqrt();
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
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.1) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                for (double y = 0.1; y < 1.2; y += 0.1) {
                    FieldDerivativeStructure<T> dsY = factory.variable(1, y);
                    for (double z = 0.1; z < 1.2; z += 0.1) {
                        FieldDerivativeStructure<T> dsZ = factory.variable(2, z);
                        FieldDerivativeStructure<T> f = dsX.divide(dsY.cos().add(dsZ.tan())).sin();
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
        double[] epsilon = new double[] { 5.0e-16, 5.0e-16, 2.7e-15, 5.7e-14, 2.0e-12 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> sqrt1 = dsX.pow(0.5);
                FieldDerivativeStructure<T> sqrt2 = dsX.sqrt();
                FieldDerivativeStructure<T> zero = sqrt1.subtract(sqrt2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testRootNSingularity() {
        doTestRootNSingularity(true);
    }

    protected void doTestRootNSingularity(final boolean signedInfinities) {
        for (int n = 2; n < 10; ++n) {
            for (int maxOrder = 0; maxOrder < 12; ++maxOrder) {
                final FDSFactory<T> factory = buildFactory(1, maxOrder);
                FieldDerivativeStructure<T> dsZero = factory.variable(0, 0.0);
                FieldDerivativeStructure<T> rootN  = dsZero.rootN(n);
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
                        final double d = rootN.getPartialDerivative(order).getReal();
                        Assert.assertTrue(Double.isNaN(d) || Double.isInfinite(d));
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
                FieldDerivativeStructure<T> correctRoot = factory.build(gDerivatives).rootN(n);
                Assert.assertEquals(0.0, correctRoot.getReal(), 1.0e-20);
                if (maxOrder > 0) {
                    Assert.assertTrue(Double.isInfinite(correctRoot.getPartialDerivative(1).getReal()));
                    Assert.assertTrue(correctRoot.getPartialDerivative(1).getReal() > 0);
                    for (int order = 2; order <= maxOrder; ++order) {
                        Assert.assertTrue(Double.isInfinite(correctRoot.getPartialDerivative(order).getReal()));
                        if (signedInfinities) {
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

    }

    @Test
    public void testSqrtPow2() {
        double[] epsilon = new double[] { 1.0e-16, 3.0e-16, 2.0e-15, 6.0e-14, 6.0e-12 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> rebuiltX = dsX.multiply(dsX).sqrt();
                FieldDerivativeStructure<T> zero = rebuiltX.subtract(dsX);
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> cbrt1 = dsX.pow(1.0 / 3.0);
                FieldDerivativeStructure<T> cbrt2 = dsX.cbrt();
                FieldDerivativeStructure<T> zero = cbrt1.subtract(cbrt2);
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> rebuiltX = dsX.multiply(dsX.multiply(dsX)).cbrt();
                FieldDerivativeStructure<T> zero = rebuiltX.subtract(dsX);
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
            final FDSFactory<T> factory = buildFactory(2, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.01) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                for (double y = 0.1; y < 1.2; y += 0.01) {
                    FieldDerivativeStructure<T> dsY = factory.variable(1, y);
                    FieldDerivativeStructure<T> rebuiltX = dsX.pow(dsY).pow(dsY.reciprocal());
                    FieldDerivativeStructure<T> zero = rebuiltX.subtract(dsX);
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
            final FDSFactory<T> factory = buildFactory(2, maxOrder);
            for (double x = -1.7; x < 2; x += 0.2) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                for (double y = -1.7; y < 2; y += 0.2) {
                    FieldDerivativeStructure<T> dsY = factory.variable(1, y);
                    FieldDerivativeStructure<T> hypot = FieldDerivativeStructure.hypot(dsY, dsX);
                    FieldDerivativeStructure<T> ref = dsX.multiply(dsX).add(dsY.multiply(dsY)).sqrt();
                    FieldDerivativeStructure<T> zero = hypot.subtract(ref);
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
    public abstract void testHypotNoOverflow();

    protected void doTestHypotNoOverflow(int tenPower) {

        final FDSFactory<T> factory = buildFactory(2, 5);
        FieldDerivativeStructure<T> dsX = factory.variable(0, +3.0);
        FieldDerivativeStructure<T> dsY = factory.variable(1, -4.0);
        T scaling = factory.getValueField().getOne();
        for (int i = 0; i < tenPower; ++i) {
            scaling = scaling.multiply(10);
        }
        dsX = dsX.multiply(scaling);
        dsY = dsY.multiply(scaling);
        FieldDerivativeStructure<T> hypot = FieldDerivativeStructure.hypot(dsX, dsY);
        FieldDerivativeStructure<T> scaledDownHypot = hypot;
        scaledDownHypot = scaledDownHypot.divide(scaling);
        Assert.assertEquals(5.0, scaledDownHypot.getReal(), 5.0e-15);
        Assert.assertEquals(dsX.divide(hypot).getReal(), scaledDownHypot.getPartialDerivative(1, 0).getReal(), 1.0e-10);
        Assert.assertEquals(dsY.divide(hypot).getReal(), scaledDownHypot.getPartialDerivative(0, 1).getReal(), 1.0e-10);

        FieldDerivativeStructure<T> sqrt  = dsX.multiply(dsX).add(dsY.multiply(dsY)).sqrt();
        Assert.assertTrue(sqrt.getValue().isInfinite() || sqrt.getValue().isNaN());

    }

    @Test
    public void testHypotNeglectible() {

        final FDSFactory<T> factory = buildFactory(2, 5);
        FieldDerivativeStructure<T> dsSmall = factory.variable(0, +3.0e-10);
        FieldDerivativeStructure<T> dsLarge = factory.variable(1, -4.0e25);

        Assert.assertEquals(dsLarge.norm(),
                            FieldDerivativeStructure.hypot(dsSmall, dsLarge).getReal(),
                            1.0e-10);
        Assert.assertEquals(0,
                            FieldDerivativeStructure.hypot(dsSmall, dsLarge).getPartialDerivative(1, 0).getReal(),
                            1.0e-10);
        Assert.assertEquals(-1,
                            FieldDerivativeStructure.hypot(dsSmall, dsLarge).getPartialDerivative(0, 1).getReal(),
                            1.0e-10);

        Assert.assertEquals(dsLarge.norm(),
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
        final FDSFactory<T> factory = buildFactory(2, 5);
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
            final FDSFactory<T> factory = buildFactory(2, maxOrder);
            for (double x = -1.7; x < 2; x += 0.2) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                for (double y = -1.7; y < 2; y += 0.2) {
                    FieldDerivativeStructure<T> remainder = dsX.remainder(buildScalar(y));
                    FieldDerivativeStructure<T> ref = dsX.subtract(x - FastMath.IEEEremainder(x, y));
                    FieldDerivativeStructure<T> zero = remainder.subtract(ref);
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
            final FDSFactory<T> factory = buildFactory(2, maxOrder);
            for (double x = -1.7; x < 2; x += 0.2) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                for (double y = -1.7; y < 2; y += 0.2) {
                    FieldDerivativeStructure<T> remainder = dsX.remainder(y);
                    FieldDerivativeStructure<T> ref = dsX.subtract(x - FastMath.IEEEremainder(x, y));
                    FieldDerivativeStructure<T> zero = remainder.subtract(ref);
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
            final FDSFactory<T> factory = buildFactory(2, maxOrder);
            for (double x = -1.7; x < 2; x += 0.2) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                for (double y = -1.7; y < 2; y += 0.2) {
                    FieldDerivativeStructure<T> dsY = factory.variable(1, y);
                    FieldDerivativeStructure<T> remainder = dsX.remainder(dsY);
                    FieldDerivativeStructure<T> ref = dsX.subtract(dsY.multiply((x - FastMath.IEEEremainder(x, y)) / y));
                    FieldDerivativeStructure<T> zero = remainder.subtract(ref);
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                double refExp = FastMath.exp(x);
                FieldDerivativeStructure<T> exp = factory.variable(0, x).exp();
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> expm11 = dsX.expm1();
                FieldDerivativeStructure<T> expm12 = dsX.exp().subtract(dsX.getField().getOne());
                FieldDerivativeStructure<T> zero = expm11.subtract(expm12);
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> log = factory.variable(0, x).log();
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> log1p1 = dsX.log1p();
                FieldDerivativeStructure<T> log1p2 = dsX.add(dsX.getField().getOne()).log();
                FieldDerivativeStructure<T> zero = log1p1.subtract(log1p2);
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> log101 = dsX.log10();
                FieldDerivativeStructure<T> log102 = dsX.log().divide(FastMath.log(10.0));
                FieldDerivativeStructure<T> zero = log101.subtract(log102);
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> rebuiltX = dsX.exp().log();
                FieldDerivativeStructure<T> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testLog1pExpm1() {
        double[] epsilon = new double[] { 3.0e-16, 3.0e-16, 5.0e-16, 9.0e-16, 6.0e-15 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> rebuiltX = dsX.expm1().log1p();
                FieldDerivativeStructure<T> zero = rebuiltX.subtract(dsX);
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> rebuiltX = factory.constant(10.0).pow(dsX).log10();
                FieldDerivativeStructure<T> zero = rebuiltX.subtract(dsX);
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> sin = dsX.sin();
                FieldDerivativeStructure<T> cos = dsX.cos();
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldSinCos<FieldDerivativeStructure<T>> sinCos = dsX.sinCos();
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> rebuiltX = dsX.sin().asin();
                FieldDerivativeStructure<T> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testCosAcos() {
        double[] epsilon = new double[] { 7.0e-16, 6.0e-15, 2.0e-13, 4.0e-12, 2.0e-10 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> rebuiltX = dsX.cos().acos();
                FieldDerivativeStructure<T> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testTanAtan() {
        double[] epsilon = new double[] { 3.0e-16, 2.0e-16, 2.0e-15, 4.0e-14, 2.0e-12 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> rebuiltX = dsX.tan().atan();
                FieldDerivativeStructure<T> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testTangentDefinition() {
        double[] epsilon = new double[] { 9.0e-16, 4.0e-15, 4.0e-14, 5.0e-13, 2.0e-11 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> tan1 = dsX.sin().divide(dsX.cos());
                FieldDerivativeStructure<T> tan2 = dsX.tan();
                FieldDerivativeStructure<T> zero = tan1.subtract(tan2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Override
    @Test
    public void testAtan2() {
        double[] epsilon = new double[] { 5.0e-16, 3.0e-15, 2.9e-14, 1.0e-12, 8.0e-11 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(2, maxOrder);
            for (double x = -1.7; x < 2; x += 0.2) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                for (double y = -1.7; y < 2; y += 0.2) {
                    FieldDerivativeStructure<T> dsY = factory.variable(1, y);
                    FieldDerivativeStructure<T> atan2 = FieldDerivativeStructure.atan2(dsY, dsX);
                    FieldDerivativeStructure<T> ref = dsY.divide(dsX).atan();
                    if (x < 0) {
                        ref = (y < 0) ? ref.subtract(FastMath.PI) : ref.add(FastMath.PI);
                    }
                    FieldDerivativeStructure<T> zero = atan2.subtract(ref);
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
    public void testAtan2SpecialCasesDerivatives() {

        final FDSFactory<T> factory = buildFactory(2, 2);
        FieldDerivativeStructure<T> pp =
                FieldDerivativeStructure.atan2(factory.variable(1, buildScalar(+0.0)), factory.variable(1, buildScalar(+0.0)));
        Assert.assertEquals(0, pp.getReal(), 1.0e-15);
        Assert.assertEquals(+1, FastMath.copySign(1, pp.getReal()), 1.0e-15);

        FieldDerivativeStructure<T> pn =
                FieldDerivativeStructure.atan2(factory.variable(1, buildScalar(+0.0)), factory.variable(1, buildScalar(-0.0)));
        Assert.assertEquals(FastMath.PI, pn.getReal(), 1.0e-15);

        FieldDerivativeStructure<T> np =
                FieldDerivativeStructure.atan2(factory.variable(1, buildScalar(-0.0)), factory.variable(1, buildScalar(+0.0)));
        Assert.assertEquals(0, np.getReal(), 1.0e-15);
        Assert.assertEquals(-1, FastMath.copySign(1, np.getReal()), 1.0e-15);

        FieldDerivativeStructure<T> nn =
                FieldDerivativeStructure.atan2(factory.variable(1, buildScalar(-0.0)), factory.variable(1, buildScalar(-0.0)));
        Assert.assertEquals(-FastMath.PI, nn.getReal(), 1.0e-15);

    }

    @Test
    public void testSinhCoshCombined() {
        double epsilon = 5.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldSinhCosh<FieldDerivativeStructure<T>> sinhCosh = dsX.sinhCosh();
                double sh = FastMath.sinh(x);
                double ch = FastMath.cosh(x);
                for (int n = 0; n <= maxOrder; ++n) {
                    if (n % 2 == 0) {
                        Assert.assertEquals(sh, sinhCosh.sinh().getPartialDerivative(n).getReal(), epsilon);
                        Assert.assertEquals(ch, sinhCosh.cosh().getPartialDerivative(n).getReal(), epsilon);
                    } else {
                        Assert.assertEquals(ch, sinhCosh.sinh().getPartialDerivative(n).getReal(), epsilon);
                        Assert.assertEquals(sh, sinhCosh.cosh().getPartialDerivative(n).getReal(), epsilon);
                    }
                }
            }
        }
    }

    @Test
    public void testSinhDefinition() {
        double[] epsilon = new double[] { 3.0e-16, 3.0e-16, 5.0e-16, 2.0e-15, 6.0e-15 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> sinh1 = dsX.exp().subtract(dsX.exp().reciprocal()).multiply(0.5);
                FieldDerivativeStructure<T> sinh2 = dsX.sinh();
                FieldDerivativeStructure<T> zero = sinh1.subtract(sinh2);
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> cosh1 = dsX.exp().add(dsX.exp().reciprocal()).multiply(0.5);
                FieldDerivativeStructure<T> cosh2 = dsX.cosh();
                FieldDerivativeStructure<T> zero = cosh1.subtract(cosh2);
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> tanh1 = dsX.exp().subtract(dsX.exp().reciprocal()).divide(dsX.exp().add(dsX.exp().reciprocal()));
                FieldDerivativeStructure<T> tanh2 = dsX.tanh();
                FieldDerivativeStructure<T> zero = tanh1.subtract(tanh2);
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> rebuiltX = dsX.sinh().asinh();
                FieldDerivativeStructure<T> zero = rebuiltX.subtract(dsX);
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> rebuiltX = dsX.cosh().acosh();
                FieldDerivativeStructure<T> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testTanhAtanh() {
        double[] epsilon = new double[] { 5.0e-16, 2.0e-16, 7.0e-16, 4.0e-15, 3.0e-14, 4.0e-13 };
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> rebuiltX = dsX.tanh().atanh();
                FieldDerivativeStructure<T> zero = rebuiltX.subtract(dsX);
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.1) {
                FieldDerivativeStructure<T> dsX = factory.constant(x);
                for (double y = 0.1; y < 1.2; y += 0.1) {
                    FieldDerivativeStructure<T> dsY = factory.variable(0, y);
                    FieldDerivativeStructure<T> f = dsX.divide(dsY).sqrt();
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
        final FDSFactory<T> factory = buildFactory(3, 4);
        for (double x = 0; x < 1.2; x += 0.1) {
            FieldDerivativeStructure<T> dsX = factory.variable(0, x);
            for (double y = 0; y < 1.2; y += 0.2) {
                FieldDerivativeStructure<T> dsY = factory.variable(1, y);
                for (double z = 0; z < 1.2; z += 0.2) {
                    FieldDerivativeStructure<T> dsZ = factory.variable(2, z);
                    FieldDerivativeStructure<T> f = dsX.multiply(dsY).add(dsZ).multiply(dsX).multiply(dsY);
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
        final FDSFactory<T> factory = buildFactory(3, 4);
        for (double x = 0; x < 1.2; x += 0.1) {
            FieldDerivativeStructure<T> dsX = factory.variable(0, x);
            for (double y = 0; y < 1.2; y += 0.2) {
                FieldDerivativeStructure<T> dsY = factory.variable(1, y);
                for (double z = 0; z < 1.2; z += 0.2) {
                    FieldDerivativeStructure<T> dsZ = factory.variable(2, z);
                    FieldDerivativeStructure<T> f = dsX.multiply(dsY).add(dsZ).multiply(dsX).multiply(dsY);
                    for (double dx = -0.2; dx < 0.2; dx += 0.2) {
                        T dxF = buildScalar(dx);
                        for (double dy = -0.2; dy < 0.2; dy += 0.1) {
                            T dyF = buildScalar(dy);
                            for (double dz = -0.2; dz < 0.2; dz += 0.1) {
                                T dzF = buildScalar(dz);
                                double ref = (x + dx) * (y + dy) * ((x + dx) * (y + dy) + (z + dz));
                                Assert.assertEquals(ref, f.taylor(dxF, dyF, dzF).getReal(), 2.0e-15);
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
            final FDSFactory<T> factory = buildFactory(2, maxOrder);
            FieldDerivativeStructure<T> dsX   = factory.variable(0, x0);
            FieldDerivativeStructure<T> dsY   = factory.variable(1, y0);
            FieldDerivativeStructure<T> atan2 = FieldDerivativeStructure.atan2(dsY, dsX);
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

    @Test
    public void testNorm() {

        final FDSFactory<T> factory = buildFactory(1, 1);
        FieldDerivativeStructure<T> minusOne = factory.variable(0, -1.0);
        Assert.assertEquals(+1.0, minusOne.abs().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.abs().getPartialDerivative(1).getReal(), 1.0e-15);

        FieldDerivativeStructure<T> plusOne = factory.variable(0, +1.0);
        Assert.assertEquals(+1.0, plusOne.abs().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.abs().getPartialDerivative(1).getReal(), 1.0e-15);

        FieldDerivativeStructure<T> minusZero = factory.variable(0, buildScalar(-0.0));
        Assert.assertEquals(+0.0, minusZero.abs().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusZero.abs().getPartialDerivative(1).getReal(), 1.0e-15);

        FieldDerivativeStructure<T> plusZero = factory.variable(0, buildScalar(+0.0));
        Assert.assertEquals(+0.0, plusZero.abs().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusZero.abs().getPartialDerivative(1).getReal(), 1.0e-15);

    }

    @Override
    @Test
    public void testSign() {

        final FDSFactory<T> factory = buildFactory(1, 1);
        FieldDerivativeStructure<T> minusOne = factory.variable(0, -1.0);
        Assert.assertEquals(-1.0, minusOne.sign().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals( 0.0, minusOne.sign().getPartialDerivative(1).getReal(), 1.0e-15);

        FieldDerivativeStructure<T> plusOne = factory.variable(0, +1.0);
        Assert.assertEquals(+1.0, plusOne.sign().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals( 0.0, plusOne.sign().getPartialDerivative(1).getReal(), 1.0e-15);

        FieldDerivativeStructure<T> minusZero = factory.variable(0, buildScalar(-0.0));
        Assert.assertEquals(-0.0, minusZero.sign().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertTrue(Double.doubleToLongBits(minusZero.sign().getReal()) < 0);
        Assert.assertEquals( 0.0, minusZero.sign().getPartialDerivative(1).getReal(), 1.0e-15);

        FieldDerivativeStructure<T> plusZero = factory.variable(0, buildScalar(+0.0));
        Assert.assertEquals(+0.0, plusZero.sign().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertTrue(Double.doubleToLongBits(plusZero.sign().getReal()) == 0);
        Assert.assertEquals( 0.0, plusZero.sign().getPartialDerivative(1).getReal(), 1.0e-15);

    }

    @Test
    public void testCeilFloorRintLong() {

        final FDSFactory<T> factory = buildFactory(1, 1);
        FieldDerivativeStructure<T> x = factory.variable(0, -1.5);
        Assert.assertEquals(-1.5, x.getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, x.getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, x.ceil().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+0.0, x.ceil().getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-2.0, x.floor().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+0.0, x.floor().getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-2.0, x.rint().getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+0.0, x.rint().getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-2.0, x.subtract(x.getField().getOne()).rint().getPartialDerivative(0).getReal(), 1.0e-15);

    }

    @Test
    public void testCopySign() {

        final FDSFactory<T> factory = buildFactory(1, 1);
        FieldDerivativeStructure<T> minusOne = factory.variable(0, -1.0);
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
        Assert.assertEquals(+1.0, minusOne.copySign(buildScalar(+1.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(buildScalar(+1.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(buildScalar(-1.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(buildScalar(-1.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(buildScalar(+0.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(buildScalar(+0.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(buildScalar(-0.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(buildScalar(-0.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, minusOne.copySign(buildScalar(Double.NaN)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, minusOne.copySign(buildScalar(Double.NaN)).getPartialDerivative(1).getReal(), 1.0e-15);

        FieldDerivativeStructure<T> plusOne = factory.variable(0, +1.0);
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
        Assert.assertEquals(+1.0, plusOne.copySign(buildScalar(+1.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(buildScalar(+1.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, plusOne.copySign(buildScalar(-1.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, plusOne.copySign(buildScalar(-1.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(buildScalar(+0.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(buildScalar(+0.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, plusOne.copySign(buildScalar(-0.0)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, plusOne.copySign(buildScalar(-0.0)).getPartialDerivative(1).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(buildScalar(Double.NaN)).getPartialDerivative(0).getReal(), 1.0e-15);
        Assert.assertEquals(+1.0, plusOne.copySign(buildScalar(Double.NaN)).getPartialDerivative(1).getReal(), 1.0e-15);

    }

    @Test
    public void testToDegreesDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
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
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> rebuiltX = dsX.toDegrees().toRadians();
                FieldDerivativeStructure<T> zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n).getReal(), epsilon);
                }
            }
        }
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testComposeMismatchedDimensions() {
        final FDSFactory<T> factory = buildFactory(1, 3);
        factory.variable(0, 1.2).compose(new double[3]);
    }

    @Test
    public abstract void testComposeField();

    protected void doTestComposeField(final double[] epsilon) {
        double[] maxError = new double[epsilon.length];
        for (int maxOrder = 0; maxOrder < epsilon.length; ++maxOrder) {
            @SuppressWarnings("unchecked")
            FieldPolynomialFunction<T>[] p = (FieldPolynomialFunction<T>[]) Array.newInstance(FieldPolynomialFunction.class,
                                                                                              maxOrder + 1);
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            T[] coefficients = MathArrays.buildArray(factory.getValueField(), epsilon.length);
            for (int i = 0; i < coefficients.length; ++i) {
                coefficients[i] = factory.getValueField().getZero().newInstance(i + 1);
            }
            p[0] = new FieldPolynomialFunction<>(coefficients);
            for (int i = 1; i <= maxOrder; ++i) {
                p[i] = p[i - 1].polynomialDerivative();
            }
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> dsY1 = dsX.getField().getZero();
                for (int i = p[0].degree(); i >= 0; --i) {
                    dsY1 = dsY1.multiply(dsX).add(p[0].getCoefficients()[i]);
                }
                T[] f = MathArrays.buildArray(getField(), maxOrder + 1);
                for (int i = 0; i < f.length; ++i) {
                    f[i] = p[i].value(x);
                }
                FieldDerivativeStructure<T> dsY2 = dsX.compose(f);
                FieldDerivativeStructure<T> zero = dsY1.subtract(dsY2);
                for (int n = 0; n <= maxOrder; ++n) {
                    maxError[n] = FastMath.max(maxError[n], FastMath.abs(zero.getPartialDerivative(n).getReal()));
                }
            }
        }
        for (int n = 0; n < maxError.length; ++n) {
            Assert.assertEquals(0.0, maxError[n], epsilon[n]);
        }
    }

    @Test
    public abstract void testComposePrimitive();

    protected void doTestComposePrimitive(final double[] epsilon) {
        PolynomialFunction poly =
                new PolynomialFunction(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 });
        double[] maxError = new double[epsilon.length];
        for (int maxOrder = 0; maxOrder < epsilon.length; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(1, maxOrder);
            PolynomialFunction[] p = new PolynomialFunction[maxOrder + 1];
            p[0] = poly;
            for (int i = 1; i <= maxOrder; ++i) {
                p[i] = p[i - 1].polynomialDerivative();
            }
            for (double x = 0.1; x < 1.2; x += 0.001) {
                FieldDerivativeStructure<T> dsX = factory.variable(0, x);
                FieldDerivativeStructure<T> dsY1 = dsX.getField().getZero();
                for (int i = poly.degree(); i >= 0; --i) {
                    dsY1 = dsY1.multiply(dsX).add(poly.getCoefficients()[i]);
                }
                double[] f = new double[maxOrder + 1];
                for (int i = 0; i < f.length; ++i) {
                    f[i] = p[i].value(x);
                }
                FieldDerivativeStructure<T> dsY2 = dsX.compose(f);
                FieldDerivativeStructure<T> zero = dsY1.subtract(dsY2);
                for (int n = 0; n <= maxOrder; ++n) {
                    maxError[n] = FastMath.max(maxError[n], FastMath.abs(zero.getPartialDerivative(n).getReal()));
                }
            }
        }
        for (int n = 0; n < maxError.length; ++n) {
            Assert.assertEquals(0.0, maxError[n], epsilon[n]);
        }
    }

    @Test
    public void testIntegration() {
        // check that first-order integration on two variables does not depend on sequence of operations
        final RandomGenerator random = new Well19937a(0x87bb96d6e11557bdl);
        final FDSFactory<T> factory = buildFactory(3, 7);
        final int size = factory.getCompiler().getSize();
        for (int count = 0; count < 100; ++count) {
            final double[] data = new double[size];
            for (int i = 0; i < size; i++) {
                data[i] = random.nextDouble();
            }
            final FieldDerivativeStructure<T> f       = factory.build(data);
            final FieldDerivativeStructure<T> i2fIxIy = f.integrate(0, 1).integrate(1, 1);
            final FieldDerivativeStructure<T> i2fIyIx = f.integrate(1, 1).integrate(0, 1);
            checkEquals(i2fIxIy, i2fIyIx, 0.);
        }
    }

    @Test
    public void testIntegrationGreaterThanOrder() {
        // check that integration to a too high order generates zero
        // as integration constants are set to zero
        final RandomGenerator random = new Well19937a(0x4744a847b11e4c6fl);
        final FDSFactory<T> factory = buildFactory(3, 7);
        final int size = factory.getCompiler().getSize();
        for (int count = 0; count < 100; ++count) {
            final double[] data = new double[size];
            for (int i = 0; i < size; i++) {
                data[i] = random.nextDouble();
            }
            final FieldDerivativeStructure<T> f = factory.build(data);
            for (int index = 0; index < factory.getCompiler().getFreeParameters(); ++index) {
                final FieldDerivativeStructure<T> integ = f.integrate(index, factory.getCompiler().getOrder() + 1);
                checkEquals(factory.constant(0), integ, 0.);
            }
        }
    }

    @Test
    public void testIntegrationNoOp() {
        // check that integration of order 0 is no-op
        final RandomGenerator random = new Well19937a(0x75a35152f30f644bl);
        final FDSFactory<T> factory = buildFactory(3, 7);
        final int size = factory.getCompiler().getSize();
        for (int count = 0; count < 100; ++count) {
            final double[] data = new double[size];
            for (int i = 0; i < size; i++) {
                data[i] = random.nextDouble();
            }
            final FieldDerivativeStructure<T> f = factory.build(data);
            for (int index = 0; index < factory.getCompiler().getFreeParameters(); ++index) {
                final FieldDerivativeStructure<T> integ = f.integrate(index, 0);
                checkEquals(f, integ, 0.);
            }
        }
    }

    @Test
    public void testDifferentiationNoOp() {
        // check that differentiation of order 0 is no-op
        final RandomGenerator random = new Well19937a(0x3b6ae4c2f1282949l);
        final FDSFactory<T> factory = buildFactory(3, 7);
        final int size = factory.getCompiler().getSize();
        for (int count = 0; count < 100; ++count) {
            final double[] data = new double[size];
            for (int i = 0; i < size; i++) {
                data[i] = random.nextDouble();
            }
            final FieldDerivativeStructure<T> f = factory.build(data);
            for (int index = 0; index < factory.getCompiler().getFreeParameters(); ++index) {
                final FieldDerivativeStructure<T> integ = f.differentiate(index, 0);
                checkEquals(f, integ, 0.);
            }
        }
    }

    @Test
    public void testIntegrationDifferentiation() {
        // check that integration and differentiation for univariate functions are each other inverse except for constant
        // term and highest order one
        final RandomGenerator random = new Well19937a(0x67fe66c05e5ee222l);
        final FDSFactory<T> factory = buildFactory(1, 25);
        final int size = factory.getCompiler().getSize();
        for (int count = 0; count < 100; ++count) {
            final double[] data = new double[size];
            for (int i = 1; i < size - 1; i++) {
                data[i] = random.nextDouble();
            }
            final int indexVar = 0;
            final FieldDerivativeStructure<T> f = factory.build(data);
            final FieldDerivativeStructure<T> f2 = f.integrate(indexVar, 1).differentiate(indexVar, 1);
            final FieldDerivativeStructure<T> f3 = f.differentiate(indexVar, 1).integrate(indexVar, 1);
            checkEquals(f2, f, 0.);
            checkEquals(f2, f3, 0.);
            // check special case when non-positive integration order actually returns differentiation
            final FieldDerivativeStructure<T> df = f.integrate(indexVar, -1);
            final FieldDerivativeStructure<T> df2 = f.differentiate(indexVar, 1);
            checkEquals(df, df2, 0.);
            // check special case when non-positive differentiation order actually returns integration
            final FieldDerivativeStructure<T> fi  = f.differentiate(indexVar, -1);
            final FieldDerivativeStructure<T> fi2 = f.integrate(indexVar, 1);
            checkEquals(fi, fi2, 0.);
        }
    }

    @Test
    public void testDifferentiation1() {
        // check differentiation operator with result obtained manually
        final int freeParam = 3;
        final int order = 5;
        final FDSFactory<T> factory = buildFactory(freeParam, order);
        final FieldDerivativeStructure<T> f = factory.variable(0, 1.0);
        final int[] orders = new int[freeParam];
        orders[0] = 2;
        orders[1] = 1;
        orders[2] = 1;
        final T value = factory.getValueField().getZero().newInstance(10.);
        f.setDerivativeComponent(factory.getCompiler().getPartialDerivativeIndex(orders), value);
        final FieldDerivativeStructure<T> dfDx = f.differentiate(0, 1);
        orders[0] -= 1;
        Assert.assertEquals(1., dfDx.getPartialDerivative(new int[freeParam]).getReal(), 0.);
        Assert.assertEquals(value.getReal(), dfDx.getPartialDerivative(orders).getReal(), 0.);
        checkEquals(factory.constant(0.0), f.differentiate(0, order + 1), 0.);
    }

    @Test
    public void testDifferentiation2() {
        // check that first-order differentiation twice is same as second-order differentiation
        final RandomGenerator random = new Well19937a(0xec293aaee352de94l);
        final FDSFactory<T> factory = buildFactory(5, 4);
        final int size = factory.getCompiler().getSize();
        for (int count = 0; count < 100; ++count) {
            final double[] data = new double[size];
            for (int i = 0; i < size; i++) {
                data[i] = random.nextDouble();
            }
            final FieldDerivativeStructure<T> f = factory.build(data);
            final FieldDerivativeStructure<T> d2fDx2 = f.differentiate(0, 1).differentiate(0, 1);
            final FieldDerivativeStructure<T> d2fDx2Bis = f.differentiate(0, 2);
            checkEquals(d2fDx2, d2fDx2Bis, 0.);
        }
    }

    @Test
    public void testDifferentiation3() {
        // check that first-order differentiation on two variables does not depend on sequence of operations
        final RandomGenerator random = new Well19937a(0x35409ecc1348e46cl);
        final FDSFactory<T> factory = buildFactory(3, 7);
        final int size = factory.getCompiler().getSize();
        for (int count = 0; count < 100; ++count) {
            final double[] data = new double[size];
            for (int i = 0; i < size; i++) {
                data[i] = random.nextDouble();
            }
            final FieldDerivativeStructure<T> f = factory.build(data);
            final FieldDerivativeStructure<T> d2fDxDy = f.differentiate(0, 1).differentiate(1, 1);
            final FieldDerivativeStructure<T> d2fDyDx = f.differentiate(1, 1).differentiate(0, 1);
            checkEquals(d2fDxDy, d2fDyDx, 0.);
        }
    }

    @Test
    public void testField() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            final FDSFactory<T> factory = buildFactory(3, maxOrder);
            FieldDerivativeStructure<T> x = factory.variable(0, 1.0);
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
        final FDSFactory<T> factory = buildFactory(1, 4);
        FieldDerivativeStructure<T> yRef = factory.variable(0, x).cos();
        try {
            factory.build(0.0, 0.0);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            // expected
        } catch (Exception e) {
            Assert.fail("wrong exceptionc caught " + e.getClass().getName());
        }
        double[] derivatives = new double[] { cos, -sin, -cos, sin, cos };
        FieldDerivativeStructure<T> y = factory.build(derivatives);
        checkEquals(yRef, y, 1.0e-15);
        T[] all = y.getAllDerivatives();
        Assert.assertEquals(derivatives.length, all.length);
        for (int i = 0; i < all.length; ++i) {
            Assert.assertEquals(derivatives[i], all[i].getReal(), 1.0e-15);
        }
    }

    @Test
    public void testOneOrderConstructor() {
        double x =  1.2;
        double y =  2.4;
        double z = 12.5;
        final FDSFactory<T> factory = buildFactory(3, 1);
        FieldDerivativeStructure<T> xRef = factory.variable(0, x);
        FieldDerivativeStructure<T> yRef = factory.variable(1, y);
        FieldDerivativeStructure<T> zRef = factory.variable(2, z);
        try {
            factory.build(x + y - z, 1.0, 1.0);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            // expected
        } catch (Exception e) {
            Assert.fail("wrong exceptionc caught " + e.getClass().getName());
        }
        double[] derivatives = new double[] { x + y - z, 1.0, 1.0, -1.0 };
        FieldDerivativeStructure<T> t = factory.build(derivatives);
        checkEquals(xRef.add(yRef.subtract(zRef)), t, 1.0e-15);
        T[] all = xRef.add(yRef.subtract(zRef)).getAllDerivatives();
        Assert.assertEquals(derivatives.length, all.length);
        for (int i = 0; i < all.length; ++i) {
            Assert.assertEquals(derivatives[i], all[i].getReal(), 1.0e-15);
        }
    }

    @Test
    public void testLinearCombination1DSDS() {
        doTestLinearCombination1DSDS(1.0e-15);
    }

    protected void doTestLinearCombination1DSDS(final double tol) {
        final FDSFactory<T> factory = buildFactory(6, 1);
        final FieldDerivativeStructure<T>[] a = MathArrays.buildArray(factory.getDerivativeField(), 3);
        a[0] = factory.variable(0, -1321008684645961.0 / 268435456.0);
        a[1] = factory.variable(1, -5774608829631843.0 / 268435456.0);
        a[2] = factory.variable(2, -7645843051051357.0 / 8589934592.0);
        final FieldDerivativeStructure<T>[] b = MathArrays.buildArray(factory.getDerivativeField(), 3);
        b[0] = factory.variable(3, -5712344449280879.0 / 2097152.0);
        b[1] = factory.variable(4, -4550117129121957.0 / 2097152.0);
        b[2] = factory.variable(5, 8846951984510141.0 / 131072.0);

        final FieldDerivativeStructure<T> abSumInline = a[0].linearCombination(a[0], b[0], a[1], b[1], a[2], b[2]);
        final FieldDerivativeStructure<T> abSumArray = a[0].linearCombination(a, b);

        Assert.assertEquals(abSumInline.getReal(), abSumArray.getReal(), 0);
        Assert.assertEquals(-1.8551294182586248737720779899, abSumInline.getReal(), tol);
        Assert.assertEquals(b[0].getReal(), abSumInline.getPartialDerivative(1, 0, 0, 0, 0, 0).getReal(), 1.0e-15);
        Assert.assertEquals(b[1].getReal(), abSumInline.getPartialDerivative(0, 1, 0, 0, 0, 0).getReal(), 1.0e-15);
        Assert.assertEquals(b[2].getReal(), abSumInline.getPartialDerivative(0, 0, 1, 0, 0, 0).getReal(), 1.0e-15);
        Assert.assertEquals(a[0].getReal(), abSumInline.getPartialDerivative(0, 0, 0, 1, 0, 0).getReal(), 1.0e-15);
        Assert.assertEquals(a[1].getReal(), abSumInline.getPartialDerivative(0, 0, 0, 0, 1, 0).getReal(), 1.0e-15);
        Assert.assertEquals(a[2].getReal(), abSumInline.getPartialDerivative(0, 0, 0, 0, 0, 1).getReal(), 1.0e-15);

    }

    @Test
    public void testLinearCombination1FieldDS() {
        doTestLinearCombination1FieldDS(1.0e-15);
    }

    protected void doTestLinearCombination1FieldDS(final double tol) {
        final FDSFactory<T> factory = buildFactory(3, 1);
        final T[] a = MathArrays.buildArray(getField(), 3);
        a[0] = buildScalar(-1321008684645961.0 / 268435456.0);
        a[1] = buildScalar(-5774608829631843.0 / 268435456.0);
        a[2] = buildScalar(-7645843051051357.0 / 8589934592.0);
        final FieldDerivativeStructure<T>[] b = MathArrays.buildArray(factory.getDerivativeField(), 3);
        b[0] = factory.variable(0, -5712344449280879.0 / 2097152.0);
        b[1] = factory.variable(1, -4550117129121957.0 / 2097152.0);
        b[2] = factory.variable(2, 8846951984510141.0 / 131072.0);

        final FieldDerivativeStructure<T> abSumInline = b[0].linearCombination(a[0], b[0],
                                                                               a[1], b[1],
                                                                               a[2], b[2]);
        final FieldDerivativeStructure<T> abSumArray = b[0].linearCombination(a, b);

        Assert.assertEquals(abSumInline.getReal(), abSumArray.getReal(), 0);
        Assert.assertEquals(-1.8551294182586248737720779899, abSumInline.getReal(), tol);
        Assert.assertEquals(a[0].getReal(), abSumInline.getPartialDerivative(1, 0, 0).getReal(), 1.0e-15);
        Assert.assertEquals(a[1].getReal(), abSumInline.getPartialDerivative(0, 1, 0).getReal(), 1.0e-15);
        Assert.assertEquals(a[2].getReal(), abSumInline.getPartialDerivative(0, 0, 1).getReal(), 1.0e-15);

    }

    @Test
    public void testLinearCombination1DoubleDS() {
        doTestLinearCombination1DoubleDS(1.0e-15);
    }

    protected void doTestLinearCombination1DoubleDS(final double tol) {
        final FDSFactory<T> factory = buildFactory(3, 1);
        final double[] a = new double[] {
            -1321008684645961.0 / 268435456.0,
            -5774608829631843.0 / 268435456.0,
            -7645843051051357.0 / 8589934592.0
        };
        final FieldDerivativeStructure<T>[] b = MathArrays.buildArray(factory.getDerivativeField(), 3);
        b[0] = factory.variable(0, -5712344449280879.0 / 2097152.0);
        b[1] = factory.variable(1, -4550117129121957.0 / 2097152.0);
        b[2] = factory.variable(2, 8846951984510141.0 / 131072.0);

        final FieldDerivativeStructure<T> abSumInline = b[0].linearCombination(a[0], b[0],
                                                                       a[1], b[1],
                                                                       a[2], b[2]);
        final FieldDerivativeStructure<T> abSumArray = b[0].linearCombination(a, b);

        Assert.assertEquals(abSumInline.getReal(), abSumArray.getReal(), 0);
        Assert.assertEquals(-1.8551294182586248737720779899, abSumInline.getReal(), tol);
        Assert.assertEquals(a[0], abSumInline.getPartialDerivative(1, 0, 0).getReal(), 1.0e-15);
        Assert.assertEquals(a[1], abSumInline.getPartialDerivative(0, 1, 0).getReal(), 1.0e-15);
        Assert.assertEquals(a[2], abSumInline.getPartialDerivative(0, 0, 1).getReal(), 1.0e-15);

    }

    @Test
    public void testLinearCombination2DSDS() {

        final FDSFactory<T> factory = buildFactory(4, 1);
        final FieldDerivativeStructure<T>[] u = MathArrays.buildArray(factory.getDerivativeField(), 4);
        final FieldDerivativeStructure<T>[] v = MathArrays.buildArray(factory.getDerivativeField(), 4);

        // we compare accurate versus naive dot product implementations
        // on regular vectors (i.e. not extreme cases like in the previous test)
        Well1024a random = new Well1024a(0xc6af886975069f11l);

        for (int i = 0; i < 10000; ++i) {
            for (int j = 0; j < u.length; ++j) {
                u[j] = factory.variable(j, 1e17 * random.nextDouble());
                v[j] = factory.constant(1e17 * random.nextDouble());
            }

            FieldDerivativeStructure<T> lin = u[0].linearCombination(u[0], v[0], u[1], v[1]);
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
        final FDSFactory<T> factory = buildFactory(4, 1);
        final double[] u = new double[4];
        final FieldDerivativeStructure<T>[] v = MathArrays.buildArray(factory.getDerivativeField(), 4);
        // we compare accurate versus naive dot product implementations
        // on regular vectors (i.e. not extreme cases like in the previous test)
        Well1024a random = new Well1024a(0xc6af886975069f11l);

        for (int i = 0; i < 10000; ++i) {
            for (int j = 0; j < u.length; ++j) {
                u[j] = 1e17 * random.nextDouble();
                v[j] = factory.variable(j, 1e17 * random.nextDouble());
            }

            FieldDerivativeStructure<T> lin = v[0].linearCombination(u[0], v[0], u[1], v[1]);
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
        final FDSFactory<T> factory = buildFactory(4, 1);
        final T[] u = MathArrays.buildArray(getField(), 4);
        final FieldDerivativeStructure<T>[] v = MathArrays.buildArray(factory.getDerivativeField(), 4);
        // we compare accurate versus naive dot product implementations
        // on regular vectors (i.e. not extreme cases like in the previous test)
        Well1024a random = new Well1024a(0xc6af886975069f11l);

        for (int i = 0; i < 10000; ++i) {
            for (int j = 0; j < u.length; ++j) {
                u[j] = buildScalar(1e17 * random.nextDouble());
                v[j] = factory.variable(j, 1e17 * random.nextDouble());
            }

            FieldDerivativeStructure<T> lin = v[0].linearCombination(u[0], v[0], u[1], v[1]);
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
        FDSFactory<T> factory = buildFactory(3, 2);
        FieldDerivativeStructure<T> zero = factory.constant(17).getField().getZero();
        T[] a = zero.getAllDerivatives();
        Assert.assertEquals(10, a.length);
        for (int i = 0; i < a.length; ++i) {
            Assert.assertEquals(buildScalar(0.0), a[i]);
        }
    }

    @Test
    public void testOne() {
        FDSFactory<T> factory = buildFactory(3, 2);
        FieldDerivativeStructure<T> one = factory.constant(17).getField().getOne();
        T[] a = one.getAllDerivatives();
        Assert.assertEquals(10, a.length);
        for (int i = 0; i < a.length; ++i) {
            Assert.assertEquals(i == 0 ? buildScalar(1.0) : buildScalar(0.0), a[i]);
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
            FDSFactory<T> factory = buildFactory(parameters, order);
            map.put(factory.constant(buildScalar(17)).getField(), 0);
        }

        // despite we have created numerous factories,
        // there should be only one field for each pair parameters/order
        Assert.assertEquals(pairs.size(), map.size());
        @SuppressWarnings("unchecked")
        Field<FieldDerivativeStructure<T>> first = (Field<FieldDerivativeStructure<T>>) map.entrySet().iterator().next().getKey();
        Assert.assertTrue(first.equals(first));
        Assert.assertFalse(first.equals(getField()));

        // even at same parameters and differentiation orders, different values generate different fields
        FieldDerivativeStructure<T> zero64 = buildFactory(3, 2).build();
        FieldDerivativeStructure<Dfp> zeroDFP = new FDSFactory<Dfp>(new DfpField(15), 3, 2).build();
        Assert.assertEquals(zero64.getFreeParameters(), zeroDFP.getFreeParameters());
        Assert.assertEquals(zero64.getOrder(), zeroDFP.getOrder());
        Assert.assertFalse(zero64.getField().equals(zeroDFP.getField()));
    }

    @SuppressWarnings("unchecked")
    @Test
    public void testRebaseConditions() {
        final FDSFactory<T> f32 = buildFactory(3, 2);
        final FDSFactory<T> f22 = buildFactory(2, 2);
        final FDSFactory<T> f31 = buildFactory(3, 1);
        try {
            f32.variable(0, 0).rebase(f22.variable(0, 0), f22.variable(1, 1.0));
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            Assert.assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            Assert.assertEquals(2, ((Integer) miae.getParts()[1]).intValue());
        }
        try {
            f32.variable(0, 0).rebase(f31.variable(0, 0), f31.variable(1, 1.0), f31.variable(2, 2.0));
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            Assert.assertEquals(2, ((Integer) miae.getParts()[0]).intValue());
            Assert.assertEquals(1, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @SuppressWarnings("unchecked")
    @Test
    public void testRebaseNoVariables() {
        final FieldDerivativeStructure<T> x = buildFactory(0, 2).constant(1.0);
        Assert.assertSame(x, x.rebase());
    }

    @Test
    public void testRebaseValueMoreIntermediateThanBase() {
        doTestRebaseValue(createBaseVariables(buildFactory(2, 4), 1.5, -2.0),
                          q -> {
                              final FieldDerivativeStructure<T>[] a = MathArrays.buildArray(q[0].getFactory().getDerivativeField(), 3);
                              a[0] = q[0].add(q[1].multiply(3));
                              a[1] = q[0].log();
                              a[2] = q[1].divide(q[0].sin());
                              return a;
                          },
                          buildFactory(3, 4),
                          p -> p[0].add(p[1].divide(p[2])),
                          1.0e-15);
    }

    @Test
    public void testRebaseValueLessIntermediateThanBase() {
        doTestRebaseValue(createBaseVariables(buildFactory(3, 4), 1.5, -2.0, 0.5),
                          q -> {
                              final FieldDerivativeStructure<T>[] a = MathArrays.buildArray(q[0].getFactory().getDerivativeField(), 2);
                              a[0] = q[0].add(q[1].multiply(3));
                              a[1] = q[0].add(q[1]).subtract(q[2]);
                              return a;
                          },
                          buildFactory(2, 4),
                          p -> p[0].multiply(p[1]),
                          1.0e-15);
    }

    @Test
    public void testRebaseValueEqualIntermediateAndBase() {
        doTestRebaseValue(createBaseVariables(buildFactory(2, 4), 1.5, -2.0),
                          q -> {
                              final FieldDerivativeStructure<T>[] a = MathArrays.buildArray(q[0].getFactory().getDerivativeField(), 2);
                              a[0] = q[0].add(q[1].multiply(3));
                              a[1] = q[0].add(q[1]);
                              return a;
                          },
                          buildFactory(2, 4),
                          p -> p[0].multiply(p[1]),
                          1.0e-15);
    }

    private void doTestRebaseValue(final FieldDerivativeStructure<T>[] q,
                                   final CalculusFieldMultivariateVectorFunction<FieldDerivativeStructure<T>> qToP,
                                   final FDSFactory<T> factoryP,
                                   final CalculusFieldMultivariateFunction<FieldDerivativeStructure<T>> f,
                                   final double tol) {

        // intermediate variables as functions of base variables
        final FieldDerivativeStructure<T>[] pBase = qToP.value(q);
        
        // reference function
        final FieldDerivativeStructure<T> ref = f.value(pBase);

        // intermediate variables as independent variables
        final FieldDerivativeStructure<T>[] pIntermediate = creatIntermediateVariables(factoryP, pBase);

        // function of the intermediate variables
        final FieldDerivativeStructure<T> fI = f.value(pIntermediate);

        // function rebased to base variables
        final FieldDerivativeStructure<T> rebased = fI.rebase(pBase);

        Assert.assertEquals(q[0].getFreeParameters(),                   ref.getFreeParameters());
        Assert.assertEquals(q[0].getOrder(),                            ref.getOrder());
        Assert.assertEquals(factoryP.getCompiler().getFreeParameters(), fI.getFreeParameters());
        Assert.assertEquals(factoryP.getCompiler().getOrder(),          fI.getOrder());
        Assert.assertEquals(ref.getFreeParameters(),                    rebased.getFreeParameters());
        Assert.assertEquals(ref.getOrder(),                             rebased.getOrder());

        checkEquals(ref, rebased, tol);

    }

    final FieldDerivativeStructure<T>[] createBaseVariables(final FDSFactory<T> factory, double... q) {
        final FieldDerivativeStructure<T>[] qDS = MathArrays.buildArray(factory.getDerivativeField(), q.length);
        for (int i = 0; i < q.length; ++i) {
            qDS[i] = factory.variable(i, q[i]);
        }
        return qDS;
    }

    final FieldDerivativeStructure<T>[] creatIntermediateVariables(final FDSFactory<T> factory,
                                                                   @SuppressWarnings("unchecked") FieldDerivativeStructure<T>... pBase) {
        final FieldDerivativeStructure<T>[] pIntermediate = MathArrays.buildArray(factory.getDerivativeField(), pBase.length);
        for (int i = 0; i < pBase.length; ++i) {
            pIntermediate[i] = factory.variable(i, pBase[i].getValue());
        }
        return pIntermediate;
    }

    @Test
    public void testRunTimeClass() {
        FDSFactory<T> factory = buildFactory(3, 2);
        Field<FieldDerivativeStructure<T>> field = factory.getDerivativeField();
        Assert.assertEquals(FieldDerivativeStructure.class, field.getRuntimeClass());
        Assert.assertEquals(getField(), factory.getValueField());
        Assert.assertEquals("org.hipparchus.analysis.differentiation.FDSFactory$DerivativeField",
                            factory.getDerivativeField().getClass().getName());
    }

    private void checkF0F1(FieldDerivativeStructure<T> ds, double value, double...derivatives) {

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

    public static <T extends CalculusFieldElement<T>> void checkEquals(FieldDerivativeStructure<T> ds1,
                                                                       FieldDerivativeStructure<T> ds2,
                                                                       double epsilon) {

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
