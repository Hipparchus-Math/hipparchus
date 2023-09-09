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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

import org.hipparchus.CalculusFieldElementAbstractTest;
import org.hipparchus.Field;
import org.hipparchus.UnitTestUtils;
import org.hipparchus.analysis.CalculusFieldMultivariateFunction;
import org.hipparchus.analysis.CalculusFieldMultivariateVectorFunction;
import org.hipparchus.analysis.polynomials.PolynomialFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.random.Well19937a;
import org.hipparchus.util.ArithmeticUtils;
import org.hipparchus.util.CombinatoricsUtils;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.FieldSinhCosh;
import org.hipparchus.util.Precision;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test for class {@link DerivativeStructure}.
 */
public class DerivativeStructureTest extends CalculusFieldElementAbstractTest<DerivativeStructure> {

    @Override
    protected DerivativeStructure build(final double x) {
        return new DSFactory(2, 1).variable(0, x);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testWrongVariableIndex() {
        new DSFactory(3, 1).variable(3, 1.0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testMissingOrders() {
        new DSFactory(3, 1).variable(0, 1.0).getPartialDerivative(0, 1);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testTooLargeOrder() {
        new DSFactory(3, 1).variable(0, 1.0).getPartialDerivative(1, 1, 2);
    }

    @Test
    public void testVariableWithoutDerivative0() {
        DerivativeStructure v = new DSFactory(1, 0).variable(0, 1.0);
        Assert.assertEquals(1.0, v.getValue(), 1.0e-15);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testVariableWithoutDerivative1() {
        DerivativeStructure v = new DSFactory(1, 0).variable(0, 1.0);
        Assert.assertEquals(1.0, v.getPartialDerivative(1), 1.0e-15);
    }

    @Test
    public void testVariable() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(3, maxOrder);
            checkF0F1(factory.variable(0, 1.0), 1.0, 1.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0), 2.0, 0.0, 1.0, 0.0);
            checkF0F1(factory.variable(2, 3.0), 3.0, 0.0, 0.0, 1.0);
        }
    }

    @Test
    public void testConstant() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(3, maxOrder);
            checkF0F1(factory.constant(FastMath.PI), FastMath.PI, 0.0, 0.0, 0.0);
        }
    }

    @Test
    public void testPrimitiveAdd() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).add(5), 6.0, 1.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).add(5), 7.0, 0.0, 1.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).add(5), 8.0, 0.0, 0.0, 1.0);
        }
    }

    @Test
    public void testAdd() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(3, maxOrder);
            DerivativeStructure x = factory.variable(0, 1.0);
            DerivativeStructure y = factory.variable(1, 2.0);
            DerivativeStructure z = factory.variable(2, 3.0);
            DerivativeStructure xyz = x.add(y.add(z));
            checkF0F1(xyz, x.getValue() + y.getValue() + z.getValue(), 1.0, 1.0, 1.0);
        }
    }

    @Test
    public void testPrimitiveSubtract() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).subtract(5), -4.0, 1.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).subtract(5), -3.0, 0.0, 1.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).subtract(5), -2.0, 0.0, 0.0, 1.0);
        }
    }

    @Test
    public void testSubtract() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(3, maxOrder);
            DerivativeStructure x = factory.variable(0, 1.0);
            DerivativeStructure y = factory.variable(1, 2.0);
            DerivativeStructure z = factory.variable(2, 3.0);
            DerivativeStructure xyz = x.subtract(y.subtract(z));
            checkF0F1(xyz, x.getValue() - (y.getValue() - z.getValue()), 1.0, -1.0, 1.0);
        }
    }

    @Test
    public void testPrimitiveMultiply() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).multiply(5),  5.0, 5.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).multiply(5), 10.0, 0.0, 5.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).multiply(5), 15.0, 0.0, 0.0, 5.0);
        }
    }

    @Test
    public void testMultiply() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(3, maxOrder);
            DerivativeStructure x = factory.variable(0, 1.0);
            DerivativeStructure y = factory.variable(1, 2.0);
            DerivativeStructure z = factory.variable(2, 3.0);
            DerivativeStructure xyz = x.multiply(y.multiply(z));
            for (int i = 0; i <= maxOrder; ++i) {
                for (int j = 0; j <= maxOrder; ++j) {
                    for (int k = 0; k <= maxOrder; ++k) {
                        if (i + j + k <= maxOrder) {
                            Assert.assertEquals((i == 0 ? x.getValue() : (i == 1 ? 1.0 : 0.0)) *
                                                (j == 0 ? y.getValue() : (j == 1 ? 1.0 : 0.0)) *
                                                (k == 0 ? z.getValue() : (k == 1 ? 1.0 : 0.0)),
                                                xyz.getPartialDerivative(i, j, k),
                                                1.0e-15);
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testNegate() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(3, maxOrder);
            checkF0F1(factory.variable(0, 1.0).negate(), -1.0, -1.0, 0.0, 0.0);
            checkF0F1(factory.variable(1, 2.0).negate(), -2.0, 0.0, -1.0, 0.0);
            checkF0F1(factory.variable(2, 3.0).negate(), -3.0, 0.0, 0.0, -1.0);
        }
    }

    @Test
    public void testReciprocal() {
        for (double x = 0.1; x < 1.2; x += 0.1) {
            DerivativeStructure r = new DSFactory(1, 6).variable(0, x).reciprocal();
            Assert.assertEquals(1 / x, r.getValue(), 1.0e-15);
            for (int i = 1; i < r.getOrder(); ++i) {
                double expected = ArithmeticUtils.pow(-1, i) * CombinatoricsUtils.factorial(i) /
                                  FastMath.pow(x, i + 1);
                Assert.assertEquals(expected, r.getPartialDerivative(i), 1.0e-15 * FastMath.abs(expected));
            }
        }
    }

    @Test
    public void testPow() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(3, maxOrder);
            for (int n = 0; n < 10; ++n) {

                DerivativeStructure x = factory.variable(0, 1.0);
                DerivativeStructure y = factory.variable(1, 2.0);
                DerivativeStructure z = factory.variable(2, 3.0);
                List<DerivativeStructure> list = Arrays.asList(x, y, z,
                                                               x.add(y).add(z),
                                                               x.multiply(y).multiply(z));

                if (n == 0) {
                    for (DerivativeStructure ds : list) {
                        checkEquals(ds.getField().getOne(), FastMath.pow(ds, n), 1.0e-15);
                    }
                } else if (n == 1) {
                    for (DerivativeStructure ds : list) {
                        checkEquals(ds, FastMath.pow(ds, n), 1.0e-15);
                    }
                } else {
                    for (DerivativeStructure ds : list) {
                        DerivativeStructure p = ds.getField().getOne();
                        for (int i = 0; i < n; ++i) {
                            p = p.multiply(ds);
                        }
                        checkEquals(p, FastMath.pow(ds, n), 1.0e-15);
                    }
                }
            }
        }
    }

    @Test
    public void testPowDoubleDS() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {

            DSFactory factory = new DSFactory(3, maxOrder);
            DerivativeStructure x = factory.variable(0, 0.1);
            DerivativeStructure y = factory.variable(1, 0.2);
            DerivativeStructure z = factory.variable(2, 0.3);
            List<DerivativeStructure> list = Arrays.asList(x, y, z,
                                                           x.add(y).add(z),
                                                           x.multiply(y).multiply(z));

            for (DerivativeStructure ds : list) {
                // the special case a = 0 is included here
                for (double a : new double[] { 0.0, 0.1, 1.0, 2.0, 5.0 }) {
                    DerivativeStructure reference = (a == 0) ?
                                                    x.getField().getZero() :
                                                    FastMath.pow(new DSFactory(3, maxOrder).constant(a), ds);
                    DerivativeStructure result = DerivativeStructure.pow(a, ds);
                    checkEquals(reference, result, 1.0e-15);
                }

            }

            // negative base: -1^x can be evaluated for integers only, so value is sometimes OK, derivatives are always NaN
            DerivativeStructure negEvenInteger = DerivativeStructure.pow(-2.0, factory.variable(0, 2.0));
            Assert.assertEquals(4.0, negEvenInteger.getValue(), 1.0e-15);
            Assert.assertTrue(Double.isNaN(negEvenInteger.getPartialDerivative(1, 0, 0)));
            DerivativeStructure negOddInteger = DerivativeStructure.pow(-2.0, factory.variable(0, 3.0));
            Assert.assertEquals(-8.0, negOddInteger.getValue(), 1.0e-15);
            Assert.assertTrue(Double.isNaN(negOddInteger.getPartialDerivative(1, 0, 0)));
            DerivativeStructure negNonInteger = DerivativeStructure.pow(-2.0, factory.variable(0, 2.001));
            Assert.assertTrue(Double.isNaN(negNonInteger.getValue()));
            Assert.assertTrue(Double.isNaN(negNonInteger.getPartialDerivative(1, 0, 0)));

            DerivativeStructure zeroNeg = DerivativeStructure.pow(0.0, factory.variable(0, -1.0));
            Assert.assertTrue(Double.isNaN(zeroNeg.getValue()));
            Assert.assertTrue(Double.isNaN(zeroNeg.getPartialDerivative(1, 0, 0)));
            DerivativeStructure posNeg = DerivativeStructure.pow(2.0, factory.variable(0, -2.0));
            Assert.assertEquals(1.0 / 4.0, posNeg.getValue(), 1.0e-15);
            Assert.assertEquals(FastMath.log(2.0) / 4.0, posNeg.getPartialDerivative(1, 0, 0), 1.0e-15);

            // very special case: a = 0 and power = 0
            DerivativeStructure zeroZero = DerivativeStructure.pow(0.0, factory.variable(0, 0.0));

            // this should be OK for simple first derivative with one variable only ...
            Assert.assertEquals(1.0, zeroZero.getValue(), 1.0e-15);
            Assert.assertEquals(Double.NEGATIVE_INFINITY, zeroZero.getPartialDerivative(1, 0, 0), 1.0e-15);

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
            Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(0, 1, 0)));
            Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(0, 0, 1)));

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
                Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(2, 0, 0)));
                Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(0, 2, 0)));
                Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(0, 0, 2)));
                Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(1, 1, 0)));
                Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(0, 1, 1)));
                Assert.assertTrue(Double.isNaN(zeroZero.getPartialDerivative(1, 1, 0)));
            }

            // very special case: 0^0 where the power is a primitive
            DerivativeStructure zeroDsZeroDouble = factory.variable(0, 0.0).pow(0.0);
            boolean first = true;
            for (final double d : zeroDsZeroDouble.getAllDerivatives()) {
                if (first) {
                    Assert.assertEquals(1.0, d, Precision.EPSILON);
                    first = false;
                } else {
                    Assert.assertEquals(0.0, d, Precision.SAFE_MIN);
                }
            }
            DerivativeStructure zeroDsZeroInt = factory.variable(0, 0.0).pow(0);
            first = true;
            for (final double d : zeroDsZeroInt.getAllDerivatives()) {
                if (first) {
                    Assert.assertEquals(1.0, d, Precision.EPSILON);
                    first = false;
                } else {
                    Assert.assertEquals(0.0, d, Precision.SAFE_MIN);
                }
            }

            // 0^p with p smaller than 1.0
            DerivativeStructure u = factory.variable(1, -0.0).pow(0.25);
            for (int i0 = 0; i0 <= maxOrder; ++i0) {
                for (int i1 = 0; i1 <= maxOrder; ++i1) {
                    for (int i2 = 0; i2 <= maxOrder; ++i2) {
                        if (i0 + i1 + i2 <= maxOrder) {
                            Assert.assertEquals(0.0, u.getPartialDerivative(i0, i1, i2), 1.0e-10);
                        }
                    }
                }
            }
        }

    }

    @Test
    public void testScalb() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(3, maxOrder);
            DerivativeStructure x = factory.variable(0, 1.0);
            DerivativeStructure y = factory.variable(1, 2.0);
            DerivativeStructure z = factory.variable(2, 3.0);
            DerivativeStructure xyz = x.multiply(y.multiply(z));
            double s = 0.125;
            for (int n = -3; n <= 3; ++n) {
                DerivativeStructure scaled = xyz.scalb(n);
                for (int i = 0; i <= maxOrder; ++i) {
                    for (int j = 0; j <= maxOrder; ++j) {
                        for (int k = 0; k <= maxOrder; ++k) {
                            if (i + j + k <= maxOrder) {
                                Assert.assertEquals((i == 0 ? x.getValue() : (i == 1 ? 1.0 : 0.0)) *
                                                    (j == 0 ? y.getValue() : (j == 1 ? 1.0 : 0.0)) *
                                                    (k == 0 ? z.getValue() : (k == 1 ? 1.0 : 0.0)) *
                                                    s,
                                                    scaled.getPartialDerivative(i, j, k),
                                                    1.0e-15);
                            }
                        }
                    }
                }
                s *= 2;
            }
        }
    }

    @Test
    public void testUlp() {
        final RandomGenerator random = new Well19937a(0x85d201920b5be954l);
        for (int k = 0; k < 10000; ++k) {
            int maxOrder = 1 + random.nextInt(5);
            DSFactory factory = new DSFactory(3, maxOrder);
            DerivativeStructure x = factory.variable(0, FastMath.scalb(2 * random.nextDouble() - 1, random.nextInt(600) - 300));
            DerivativeStructure y = factory.variable(1, FastMath.scalb(2 * random.nextDouble() - 1, random.nextInt(600) - 300));
            DerivativeStructure z = factory.variable(2, FastMath.scalb(2 * random.nextDouble() - 1, random.nextInt(600) - 300));
            DerivativeStructure xyz = x.multiply(y.multiply(z));
            DerivativeStructure ulp = xyz.ulp();
            boolean first = true;
            for (double d : ulp.getAllDerivatives()) {
                Assert.assertEquals(first ? FastMath.ulp(xyz.getValue()) : 0.0, d, 1.0e-15 * FastMath.ulp(xyz.getValue()));
                first = false;
            }
        }
    }

    @Test
    public void testExpression() {
        DSFactory factory = new DSFactory(3, 5);
        double epsilon = 2.5e-13;
        for (double x = 0; x < 2; x += 0.2) {
            DerivativeStructure dsX = factory.variable(0, x);
            for (double y = 0; y < 2; y += 0.2) {
                DerivativeStructure dsY = factory.variable(1, y);
                for (double z = 0; z >- 2; z -= 0.2) {
                    DerivativeStructure dsZ = factory.variable(2, z);

                    // f(x, y, z) = x + 5 x y - 2 z + (8 z x - y)^3
                    DerivativeStructure ds =
                                    dsX.linearCombination(1, dsX,
                                                          5, dsX.multiply(dsY),
                                                         -2, dsZ,
                                                          1, dsX.linearCombination(8, dsZ.multiply(dsX),
                                                                                  -1, dsY).pow(3));
                    DerivativeStructure dsOther =
                                    dsX.linearCombination(1, dsX,
                                                          5, dsX.multiply(dsY),
                                                         -2, dsZ).add(dsX.linearCombination(8, dsZ.multiply(dsX),
                                                                                           -1, dsY).pow(3));
                    double f = x + 5 * x * y - 2 * z + FastMath.pow(8 * z * x - y, 3);
                    Assert.assertEquals(f, ds.getValue(),
                                        FastMath.abs(epsilon * f));
                    Assert.assertEquals(f, dsOther.getValue(),
                                        FastMath.abs(epsilon * f));

                    // df/dx = 1 + 5 y + 24 (8 z x - y)^2 z
                    double dfdx = 1 + 5 * y + 24 * z * FastMath.pow(8 * z * x - y, 2);
                    Assert.assertEquals(dfdx, ds.getPartialDerivative(1, 0, 0),
                                        FastMath.abs(epsilon * dfdx));
                    Assert.assertEquals(dfdx, dsOther.getPartialDerivative(1, 0, 0),
                                        FastMath.abs(epsilon * dfdx));

                    // df/dxdy = 5 + 48 z*(y - 8 z x)
                    double dfdxdy = 5 + 48 * z * (y - 8 * z * x);
                    Assert.assertEquals(dfdxdy, ds.getPartialDerivative(1, 1, 0),
                                        FastMath.abs(epsilon * dfdxdy));
                    Assert.assertEquals(dfdxdy, dsOther.getPartialDerivative(1, 1, 0),
                                        FastMath.abs(epsilon * dfdxdy));

                    // df/dxdydz = 48 (y - 16 z x)
                    double dfdxdydz = 48 * (y - 16 * z * x);
                    Assert.assertEquals(dfdxdydz, ds.getPartialDerivative(1, 1, 1),
                                        FastMath.abs(epsilon * dfdxdydz));
                    Assert.assertEquals(dfdxdydz, dsOther.getPartialDerivative(1, 1, 1),
                                        FastMath.abs(epsilon * dfdxdydz));

                }

            }
        }
    }

    @Test
    public void testCompositionOneVariableX() {
        double epsilon = 1.0e-13;
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.1) {
                DerivativeStructure dsX = factory.variable(0, x);
                for (double y = 0.1; y < 1.2; y += 0.1) {
                    DerivativeStructure dsY = factory.constant(y);
                    DerivativeStructure f = dsX.divide(dsY).sqrt();
                    double f0 = FastMath.sqrt(x / y);
                    Assert.assertEquals(f0, f.getValue(), FastMath.abs(epsilon * f0));
                    if (f.getOrder() > 0) {
                        double f1 = 1 / (2 * FastMath.sqrt(x * y));
                        Assert.assertEquals(f1, f.getPartialDerivative(1), FastMath.abs(epsilon * f1));
                        if (f.getOrder() > 1) {
                            double f2 = -f1 / (2 * x);
                            Assert.assertEquals(f2, f.getPartialDerivative(2), FastMath.abs(epsilon * f2));
                            if (f.getOrder() > 2) {
                                double f3 = (f0 + x / (2 * y * f0)) / (4 * x * x * x);
                                Assert.assertEquals(f3, f.getPartialDerivative(3), FastMath.abs(epsilon * f3));
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
            DSFactory factory = new DSFactory(3, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.1) {
                DerivativeStructure dsX = factory.variable(0, x);
                for (double y = 0.1; y < 1.2; y += 0.1) {
                    DerivativeStructure dsY = factory.variable(1, y);
                    for (double z = 0.1; z < 1.2; z += 0.1) {
                        DerivativeStructure dsZ = factory.variable(2, z);
                        DerivativeStructure f = FastMath.sin(dsX.divide(FastMath.cos(dsY).add(FastMath.tan(dsZ))));
                        double a = FastMath.cos(y) + FastMath.tan(z);
                        double f0 = FastMath.sin(x / a);
                        Assert.assertEquals(f0, f.getValue(), FastMath.abs(epsilon * f0));
                        if (f.getOrder() > 0) {
                            double dfdx = FastMath.cos(x / a) / a;
                            Assert.assertEquals(dfdx, f.getPartialDerivative(1, 0, 0), FastMath.abs(epsilon * dfdx));
                            double dfdy =  x * FastMath.sin(y) * dfdx / a;
                            Assert.assertEquals(dfdy, f.getPartialDerivative(0, 1, 0), FastMath.abs(epsilon * dfdy));
                            double cz = FastMath.cos(z);
                            double cz2 = cz * cz;
                            double dfdz = -x * dfdx / (a * cz2);
                            Assert.assertEquals(dfdz, f.getPartialDerivative(0, 0, 1), FastMath.abs(epsilon * dfdz));
                            if (f.getOrder() > 1) {
                                double df2dx2 = -(f0 / (a * a));
                                Assert.assertEquals(df2dx2, f.getPartialDerivative(2, 0, 0), FastMath.abs(epsilon * df2dx2));
                                double df2dy2 = x * FastMath.cos(y) * dfdx / a -
                                                x * x * FastMath.sin(y) * FastMath.sin(y) * f0 / (a * a * a * a) +
                                                2 * FastMath.sin(y) * dfdy / a;
                                Assert.assertEquals(df2dy2, f.getPartialDerivative(0, 2, 0), FastMath.abs(epsilon * df2dy2));
                                double c4 = cz2 * cz2;
                                double df2dz2 = x * (2 * a * (1 - a * cz * FastMath.sin(z)) * dfdx - x * f0 / a ) / (a * a * a * c4);
                                Assert.assertEquals(df2dz2, f.getPartialDerivative(0, 0, 2), FastMath.abs(epsilon * df2dz2));
                                double df2dxdy = dfdy / x  - x * FastMath.sin(y) * f0 / (a * a * a);
                                Assert.assertEquals(df2dxdy, f.getPartialDerivative(1, 1, 0), FastMath.abs(epsilon * df2dxdy));
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
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure sqrt1 = dsX.pow(0.5);
                DerivativeStructure sqrt2 = FastMath.sqrt(dsX);
                DerivativeStructure zero = sqrt1.subtract(sqrt2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testRootNSingularity() {
        for (int n = 2; n < 10; ++n) {
            for (int maxOrder = 0; maxOrder < 12; ++maxOrder) {
                DSFactory factory = new DSFactory(1, maxOrder);
                DerivativeStructure dsZero = factory.variable(0, 0.0);
                DerivativeStructure rootN  = dsZero.rootN(n);
                Assert.assertEquals(0.0, rootN.getValue(), 1.0e-20);
                if (maxOrder > 0) {
                    Assert.assertTrue(Double.isInfinite(rootN.getPartialDerivative(1)));
                    Assert.assertTrue(rootN.getPartialDerivative(1) > 0);
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
                        Assert.assertTrue(Double.isNaN(rootN.getPartialDerivative(order)));
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
                DerivativeStructure correctRoot = factory.build(gDerivatives).rootN(n);
                Assert.assertEquals(0.0, correctRoot.getValue(), 1.0e-20);
                if (maxOrder > 0) {
                    Assert.assertTrue(Double.isInfinite(correctRoot.getPartialDerivative(1)));
                    Assert.assertTrue(correctRoot.getPartialDerivative(1) > 0);
                    for (int order = 2; order <= maxOrder; ++order) {
                        Assert.assertTrue(Double.isInfinite(correctRoot.getPartialDerivative(order)));
                        if ((order % 2) == 0) {
                            Assert.assertTrue(correctRoot.getPartialDerivative(order) < 0);
                        } else {
                            Assert.assertTrue(correctRoot.getPartialDerivative(order) > 0);
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
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure rebuiltX = dsX.multiply(dsX).sqrt();
                DerivativeStructure zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testCbrtDefinition() {
        double[] epsilon = new double[] { 4.0e-16, 9.0e-16, 6.0e-15, 2.0e-13, 4.0e-12 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure cbrt1 = dsX.pow(1.0 / 3.0);
                DerivativeStructure cbrt2 = FastMath.cbrt(dsX);
                DerivativeStructure zero = cbrt1.subtract(cbrt2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testCbrtPow3() {
        double[] epsilon = new double[] { 1.0e-16, 5.0e-16, 8.0e-15, 3.0e-13, 4.0e-11 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure rebuiltX = dsX.multiply(dsX.multiply(dsX)).cbrt();
                DerivativeStructure zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testPowReciprocalPow() {
        double[] epsilon = new double[] { 2.0e-15, 2.0e-14, 3.0e-13, 8.0e-12, 3.0e-10 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(2, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.01) {
                DerivativeStructure dsX = factory.variable(0, x);
                for (double y = 0.1; y < 1.2; y += 0.01) {
                    DerivativeStructure dsY = factory.variable(1, y);
                    DerivativeStructure rebuiltX = dsX.pow(dsY).pow(dsY.reciprocal());
                    DerivativeStructure zero = rebuiltX.subtract(dsX);
                    for (int n = 0; n <= maxOrder; ++n) {
                        for (int m = 0; m <= maxOrder; ++m) {
                            if (n + m <= maxOrder) {
                                Assert.assertEquals(0.0, zero.getPartialDerivative(n, m), epsilon[n + m]);
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
            DSFactory factory = new DSFactory(2, maxOrder);
            for (double x = -1.7; x < 2; x += 0.2) {
                DerivativeStructure dsX = factory.variable(0, x);
                for (double y = -1.7; y < 2; y += 0.2) {
                    DerivativeStructure dsY = factory.variable(1, y);
                    DerivativeStructure hypot = FastMath.hypot(dsY, dsX);
                    DerivativeStructure ref = dsX.multiply(dsX).add(dsY.multiply(dsY)).sqrt();
                    DerivativeStructure zero = hypot.subtract(ref);
                    for (int n = 0; n <= maxOrder; ++n) {
                        for (int m = 0; m <= maxOrder; ++m) {
                            if (n + m <= maxOrder) {
                                Assert.assertEquals(0, zero.getPartialDerivative(n, m), epsilon);
                            }
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testHypotNoOverflow() {

        DSFactory factory = new DSFactory(2, 5);
        DerivativeStructure dsX = factory.variable(0, +3.0e250);
        DerivativeStructure dsY = factory.variable(1, -4.0e250);
        DerivativeStructure hypot = FastMath.hypot(dsX, dsY);
        Assert.assertEquals(5.0e250, hypot.getValue(), 1.0e235);
        Assert.assertEquals(dsX.getValue() / hypot.getValue(), hypot.getPartialDerivative(1, 0), 1.0e-10);
        Assert.assertEquals(dsY.getValue() / hypot.getValue(), hypot.getPartialDerivative(0, 1), 1.0e-10);

        DerivativeStructure sqrt  = dsX.multiply(dsX).add(dsY.multiply(dsY)).sqrt();
        Assert.assertTrue(Double.isInfinite(sqrt.getValue()));

    }

    @Test
    public void testHypotNeglectible() {

        DSFactory factory = new DSFactory(2, 5);
        DerivativeStructure dsSmall = factory.variable(0, +3.0e-10);
        DerivativeStructure dsLarge = factory.variable(1, -4.0e25);

        Assert.assertEquals(dsLarge.abs().getValue(),
                            DerivativeStructure.hypot(dsSmall, dsLarge).getValue(),
                            1.0e-10);
        Assert.assertEquals(0,
                            DerivativeStructure.hypot(dsSmall, dsLarge).getPartialDerivative(1, 0),
                            1.0e-10);
        Assert.assertEquals(-1,
                            DerivativeStructure.hypot(dsSmall, dsLarge).getPartialDerivative(0, 1),
                            1.0e-10);

        Assert.assertEquals(dsLarge.abs().getValue(),
                            DerivativeStructure.hypot(dsLarge, dsSmall).getValue(),
                            1.0e-10);
        Assert.assertEquals(0,
                            DerivativeStructure.hypot(dsLarge, dsSmall).getPartialDerivative(1, 0),
                            1.0e-10);
        Assert.assertEquals(-1,
                            DerivativeStructure.hypot(dsLarge, dsSmall).getPartialDerivative(0, 1),
                            1.0e-10);

    }

    @Test
    public void testHypotSpecial() {
        DSFactory factory = new DSFactory(2, 5);
        Assert.assertTrue(Double.isNaN(DerivativeStructure.hypot(factory.variable(0, Double.NaN),
                                                                 factory.variable(0, +3.0e250)).getValue()));
        Assert.assertTrue(Double.isNaN(DerivativeStructure.hypot(factory.variable(0, +3.0e250),
                                                                 factory.variable(0, Double.NaN)).getValue()));
        Assert.assertTrue(Double.isInfinite(DerivativeStructure.hypot(factory.variable(0, Double.POSITIVE_INFINITY),
                                                                      factory.variable(0, +3.0e250)).getValue()));
        Assert.assertTrue(Double.isInfinite(DerivativeStructure.hypot(factory.variable(0, +3.0e250),
                                                                      factory.variable(0, Double.POSITIVE_INFINITY)).getValue()));
    }

    @Test
    public void testPrimitiveRemainder() {
        double epsilon = 1.0e-15;
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(2, maxOrder);
            for (double x = -1.7; x < 2; x += 0.2) {
                DerivativeStructure dsX = factory.variable(0, x);
                for (double y = -1.7; y < 2; y += 0.2) {
                    DerivativeStructure remainder = FastMath.IEEEremainder(dsX, y);
                    DerivativeStructure ref = dsX.subtract(x - FastMath.IEEEremainder(x, y));
                    DerivativeStructure zero = remainder.subtract(ref);
                    for (int n = 0; n <= maxOrder; ++n) {
                        for (int m = 0; m <= maxOrder; ++m) {
                            if (n + m <= maxOrder) {
                                Assert.assertEquals(0, zero.getPartialDerivative(n, m), epsilon);
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
            DSFactory factory = new DSFactory(2, maxOrder);
            for (double x = -1.7; x < 2; x += 0.2) {
                DerivativeStructure dsX = factory.variable(0, x);
                for (double y = -1.7; y < 2; y += 0.2) {
                    DerivativeStructure dsY = factory.variable(1, y);
                    DerivativeStructure remainder = FastMath.IEEEremainder(dsX, dsY);
                    DerivativeStructure ref = dsX.subtract(dsY.multiply((x - FastMath.IEEEremainder(x, y)) / y));
                    DerivativeStructure zero = remainder.subtract(ref);
                    for (int n = 0; n <= maxOrder; ++n) {
                        for (int m = 0; m <= maxOrder; ++m) {
                            if (n + m <= maxOrder) {
                                Assert.assertEquals(0, zero.getPartialDerivative(n, m), epsilon);
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
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                double refExp = FastMath.exp(x);
                DerivativeStructure exp = FastMath.exp(factory.variable(0, x));
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(refExp, exp.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testExpm1Definition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure expm11 = FastMath.expm1(dsX);
                DerivativeStructure expm12 = dsX.exp().subtract(dsX.getField().getOne());
                DerivativeStructure zero = expm11.subtract(expm12);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n), epsilon);
                }
            }
        }
    }

    @Override
    @Test
    public void testLog() {
        double[] epsilon = new double[] { 1.0e-16, 1.0e-16, 3.0e-14, 7.0e-13, 3.0e-11 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure log = FastMath.log(factory.variable(0, x));
                Assert.assertEquals(FastMath.log(x), log.getValue(), epsilon[0]);
                for (int n = 1; n <= maxOrder; ++n) {
                    double refDer = -CombinatoricsUtils.factorial(n - 1) / FastMath.pow(-x, n);
                    Assert.assertEquals(refDer, log.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testLog1pDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DSFactory factory = new DSFactory(1, maxOrder);
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure log1p1 = FastMath.log1p(dsX);
                DerivativeStructure log1p2 = FastMath.log(dsX.add(dsX.getField().getOne()));
                DerivativeStructure zero = log1p1.subtract(log1p2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n), epsilon);
                }
            }
        }
    }

    @Test
    public void testLog10Definition() {
        double[] epsilon = new double[] { 3.0e-16, 9.0e-16, 8.0e-15, 3.0e-13, 8.0e-12 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure log101 = FastMath.log10(dsX);
                DerivativeStructure log102 = dsX.log().divide(FastMath.log(10.0));
                DerivativeStructure zero = log101.subtract(log102);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testLogExp() {
        double[] epsilon = new double[] { 2.0e-16, 2.0e-16, 3.0e-16, 2.0e-15, 6.0e-15 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure rebuiltX = dsX.exp().log();
                DerivativeStructure zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testLog1pExpm1() {
        double[] epsilon = new double[] { 6.0e-17, 3.0e-16, 5.0e-16, 9.0e-16, 6.0e-15 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure rebuiltX = dsX.expm1().log1p();
                DerivativeStructure zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testLog10Power() {
        double[] epsilon = new double[] { 3.0e-16, 3.0e-16, 9.0e-16, 6.0e-15, 6.0e-14 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure rebuiltX = factory.constant(10.0).pow(dsX).log10();
                DerivativeStructure zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testSinCosSeparated() {
        double epsilon = 5.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure sin = FastMath.sin(dsX);
                DerivativeStructure cos = FastMath.cos(dsX);
                double s = FastMath.sin(x);
                double c = FastMath.cos(x);
                for (int n = 0; n <= maxOrder; ++n) {
                    switch (n % 4) {
                    case 0 :
                        Assert.assertEquals( s, sin.getPartialDerivative(n), epsilon);
                        Assert.assertEquals( c, cos.getPartialDerivative(n), epsilon);
                        break;
                    case 1 :
                        Assert.assertEquals( c, sin.getPartialDerivative(n), epsilon);
                        Assert.assertEquals(-s, cos.getPartialDerivative(n), epsilon);
                        break;
                    case 2 :
                        Assert.assertEquals(-s, sin.getPartialDerivative(n), epsilon);
                        Assert.assertEquals(-c, cos.getPartialDerivative(n), epsilon);
                        break;
                    default :
                        Assert.assertEquals(-c, sin.getPartialDerivative(n), epsilon);
                        Assert.assertEquals( s, cos.getPartialDerivative(n), epsilon);
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
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                FieldSinCos<DerivativeStructure> sinCos = FastMath.sinCos(dsX);
                double s = FastMath.sin(x);
                double c = FastMath.cos(x);
                for (int n = 0; n <= maxOrder; ++n) {
                    switch (n % 4) {
                    case 0 :
                        Assert.assertEquals( s, sinCos.sin().getPartialDerivative(n), epsilon);
                        Assert.assertEquals( c, sinCos.cos().getPartialDerivative(n), epsilon);
                        break;
                    case 1 :
                        Assert.assertEquals( c, sinCos.sin().getPartialDerivative(n), epsilon);
                        Assert.assertEquals(-s, sinCos.cos().getPartialDerivative(n), epsilon);
                        break;
                    case 2 :
                        Assert.assertEquals(-s, sinCos.sin().getPartialDerivative(n), epsilon);
                        Assert.assertEquals(-c, sinCos.cos().getPartialDerivative(n), epsilon);
                        break;
                    default :
                        Assert.assertEquals(-c, sinCos.sin().getPartialDerivative(n), epsilon);
                        Assert.assertEquals( s, sinCos.cos().getPartialDerivative(n), epsilon);
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
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure rebuiltX = FastMath.asin(FastMath.sin(dsX));
                DerivativeStructure zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testCosAcos() {
        double[] epsilon = new double[] { 6.0e-16, 6.0e-15, 2.0e-13, 4.0e-12, 2.0e-10 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure rebuiltX = FastMath.acos(FastMath.cos(dsX));
                DerivativeStructure zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testTanAtan() {
        double[] epsilon = new double[] { 6.0e-17, 2.0e-16, 2.0e-15, 4.0e-14, 2.0e-12 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure rebuiltX = FastMath.atan(FastMath.tan(dsX));
                DerivativeStructure zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testTangentDefinition() {
        double[] epsilon = new double[] { 5.0e-16, 2.7e-15, 3.0e-14, 5.0e-13, 2.0e-11 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure tan1 = dsX.sin().divide(dsX.cos());
                DerivativeStructure tan2 = dsX.tan();
                DerivativeStructure zero = tan1.subtract(tan2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Override
    @Test
    public void testAtan2() {
        double[] epsilon = new double[] { 5.0e-16, 3.0e-15, 2.9e-14, 1.0e-12, 8.0e-11 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(2, maxOrder);
            for (double x = -1.7; x < 2; x += 0.2) {
                DerivativeStructure dsX = factory.variable(0, x);
                for (double y = -1.7; y < 2; y += 0.2) {
                    DerivativeStructure dsY = factory.variable(1, y);
                    DerivativeStructure atan2 = FastMath.atan2(dsY, dsX);
                    DerivativeStructure ref = dsY.divide(dsX).atan();
                    if (x < 0) {
                        ref = (y < 0) ? ref.subtract(FastMath.PI) : ref.add(FastMath.PI);
                    }
                    DerivativeStructure zero = atan2.subtract(ref);
                    for (int n = 0; n <= maxOrder; ++n) {
                        for (int m = 0; m <= maxOrder; ++m) {
                            if (n + m <= maxOrder) {
                                Assert.assertEquals(0, zero.getPartialDerivative(n, m), epsilon[n + m]);
                            }
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testAtan2SpecialCasesDerivative() {

        DSFactory factory = new DSFactory(2, 2);
        DerivativeStructure pp =
                DerivativeStructure.atan2(factory.variable(1, +0.0),
                                          factory.variable(1, +0.0));
        Assert.assertEquals(0, pp.getValue(), 1.0e-15);
        Assert.assertEquals(+1, FastMath.copySign(1, pp.getValue()), 1.0e-15);

        DerivativeStructure pn =
                DerivativeStructure.atan2(factory.variable(1, +0.0),
                                          factory.variable(1, -0.0));
        Assert.assertEquals(FastMath.PI, pn.getValue(), 1.0e-15);

        DerivativeStructure np =
                DerivativeStructure.atan2(factory.variable(1, -0.0),
                                          factory.variable(1, +0.0));
        Assert.assertEquals(0, np.getValue(), 1.0e-15);
        Assert.assertEquals(-1, FastMath.copySign(1, np.getValue()), 1.0e-15);

        DerivativeStructure nn =
                DerivativeStructure.atan2(factory.variable(1, -0.0),
                                          factory.variable(1, -0.0));
        Assert.assertEquals(-FastMath.PI, nn.getValue(), 1.0e-15);

    }

    @Test
    public void testSinhCoshCombined() {
        double epsilon = 5.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                FieldSinhCosh<DerivativeStructure> sinhCosh = FastMath.sinhCosh(dsX);
                double sh = FastMath.sinh(x);
                double ch = FastMath.cosh(x);
                for (int n = 0; n <= maxOrder; ++n) {
                    if (n % 2 == 0) {
                        Assert.assertEquals(sh, sinhCosh.sinh().getPartialDerivative(n), epsilon);
                        Assert.assertEquals(ch, sinhCosh.cosh().getPartialDerivative(n), epsilon);
                    } else {
                        Assert.assertEquals(ch, sinhCosh.sinh().getPartialDerivative(n), epsilon);
                        Assert.assertEquals(sh, sinhCosh.cosh().getPartialDerivative(n), epsilon);
                    }
                }
            }
        }
    }

    @Test
    public void testSinhDefinition() {
        double[] epsilon = new double[] { 3.0e-16, 3.0e-16, 5.0e-16, 2.0e-15, 6.0e-15 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure sinh1 = dsX.exp().subtract(dsX.exp().reciprocal()).multiply(0.5);
                DerivativeStructure sinh2 = FastMath.sinh(dsX);
                DerivativeStructure zero = sinh1.subtract(sinh2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testCoshDefinition() {
        double[] epsilon = new double[] { 3.0e-16, 3.0e-16, 5.0e-16, 2.0e-15, 6.0e-15 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure cosh1 = dsX.exp().add(dsX.exp().reciprocal()).multiply(0.5);
                DerivativeStructure cosh2 = FastMath.cosh(dsX);
                DerivativeStructure zero = cosh1.subtract(cosh2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testTanhDefinition() {
        double[] epsilon = new double[] { 3.0e-16, 5.0e-16, 7.0e-16, 3.0e-15, 2.0e-14 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure tanh1 = dsX.exp().subtract(dsX.exp().reciprocal()).divide(dsX.exp().add(dsX.exp().reciprocal()));
                DerivativeStructure tanh2 = FastMath.tanh(dsX);
                DerivativeStructure zero = tanh1.subtract(tanh2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testSinhAsinh() {
        double[] epsilon = new double[] { 3.0e-16, 3.0e-16, 4.0e-16, 7.0e-16, 3.0e-15, 8.0e-15 };
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure rebuiltX = FastMath.asinh(dsX.sinh());
                DerivativeStructure zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testCoshAcosh() {
        double[] epsilon = new double[] { 2.0e-15, 1.0e-14, 2.0e-13, 6.0e-12, 3.0e-10, 2.0e-8 };
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure rebuiltX = FastMath.acosh(dsX.cosh());
                DerivativeStructure zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testTanhAtanh() {
        double[] epsilon = new double[] { 3.0e-16, 2.0e-16, 7.0e-16, 4.0e-15, 3.0e-14, 4.0e-13 };
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure rebuiltX = FastMath.atanh(dsX.tanh());
                DerivativeStructure zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testCompositionOneVariableY() {
        double epsilon = 1.0e-13;
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.1) {
                DerivativeStructure dsX = factory.constant(x);
                for (double y = 0.1; y < 1.2; y += 0.1) {
                    DerivativeStructure dsY = factory.variable(0, y);
                    DerivativeStructure f = dsX.divide(dsY).sqrt();
                    double f0 = FastMath.sqrt(x / y);
                    Assert.assertEquals(f0, f.getValue(), FastMath.abs(epsilon * f0));
                    if (f.getOrder() > 0) {
                        double f1 = -x / (2 * y * y * f0);
                        Assert.assertEquals(f1, f.getPartialDerivative(1), FastMath.abs(epsilon * f1));
                        if (f.getOrder() > 1) {
                            double f2 = (f0 - x / (4 * y * f0)) / (y * y);
                            Assert.assertEquals(f2, f.getPartialDerivative(2), FastMath.abs(epsilon * f2));
                            if (f.getOrder() > 2) {
                                double f3 = (x / (8 * y * f0) - 2 * f0) / (y * y * y);
                                Assert.assertEquals(f3, f.getPartialDerivative(3), FastMath.abs(epsilon * f3));
                            }
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testTaylorPolynomial() {
        DSFactory factory = new DSFactory(3, 4);
        for (double x = 0; x < 1.2; x += 0.1) {
            DerivativeStructure dsX = factory.variable(0, x);
            for (double y = 0; y < 1.2; y += 0.2) {
                DerivativeStructure dsY = factory.variable(1, y);
                for (double z = 0; z < 1.2; z += 0.2) {
                    DerivativeStructure dsZ = factory.variable(2, z);
                    DerivativeStructure f = dsX.multiply(dsY).add(dsZ).multiply(dsX).multiply(dsY);
                    for (double dx = -0.2; dx < 0.2; dx += 0.2) {
                        for (double dy = -0.2; dy < 0.2; dy += 0.1) {
                            for (double dz = -0.2; dz < 0.2; dz += 0.1) {
                                double ref = (x + dx) * (y + dy) * ((x + dx) * (y + dy) + (z + dz));
                                Assert.assertEquals(ref, f.taylor(dx, dy, dz), 2.0e-15);
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
            DSFactory factory = new DSFactory(2, maxOrder);
            DerivativeStructure dsX   = factory.variable(0, x0);
            DerivativeStructure dsY   = factory.variable(1, y0);
            DerivativeStructure atan2 = DerivativeStructure.atan2(dsY, dsX);
            double maxError = 0;
            for (double dx = -0.05; dx < 0.05; dx += 0.001) {
                for (double dy = -0.05; dy < 0.05; dy += 0.001) {
                    double ref = FastMath.atan2(y0 + dy, x0 + dx);
                    maxError = FastMath.max(maxError, FastMath.abs(ref - atan2.taylor(dx, dy)));
                }
            }
            Assert.assertEquals(0.0, expected[maxOrder] - maxError, 0.01 * expected[maxOrder]);
        }
    }

    @Test
    public void testAbs() {

        DSFactory factory = new DSFactory(1, 1);
        DerivativeStructure minusOne = factory.variable(0, -1.0);
        Assert.assertEquals(+1.0, FastMath.abs(minusOne).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(-1.0, FastMath.abs(minusOne).getPartialDerivative(1), 1.0e-15);

        DerivativeStructure plusOne = factory.variable(0, +1.0);
        Assert.assertEquals(+1.0, FastMath.abs(plusOne).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(+1.0, FastMath.abs(plusOne).getPartialDerivative(1), 1.0e-15);

        DerivativeStructure minusZero = factory.variable(0, -0.0);
        Assert.assertEquals(+0.0, FastMath.abs(minusZero).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(-1.0, FastMath.abs(minusZero).getPartialDerivative(1), 1.0e-15);

        DerivativeStructure plusZero = factory.variable(0, +0.0);
        Assert.assertEquals(+0.0, FastMath.abs(plusZero).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(+1.0, FastMath.abs(plusZero).getPartialDerivative(1), 1.0e-15);

    }

    @Override
    @Test
    public void testSign() {

        DSFactory factory = new DSFactory(1, 1);
        DerivativeStructure minusOne = factory.variable(0, -1.0);
        Assert.assertEquals(-1.0, FastMath.sign(minusOne).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals( 0.0, FastMath.sign(minusOne).getPartialDerivative(1), 1.0e-15);

        DerivativeStructure plusOne = factory.variable(0, +1.0);
        Assert.assertEquals(+1.0, FastMath.sign(plusOne).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals( 0.0, FastMath.sign(plusOne).getPartialDerivative(1), 1.0e-15);

        DerivativeStructure minusZero = factory.variable(0, -0.0);
        Assert.assertEquals(-0.0, FastMath.sign(minusZero).getPartialDerivative(0), 1.0e-15);
        Assert.assertTrue(Double.doubleToLongBits(FastMath.sign(minusZero).getValue()) < 0);
        Assert.assertEquals( 0.0, FastMath.sign(minusZero).getPartialDerivative(1), 1.0e-15);

        DerivativeStructure plusZero = factory.variable(0, +0.0);
        Assert.assertEquals(+0.0, FastMath.sign(plusZero).getPartialDerivative(0), 1.0e-15);
        Assert.assertTrue(Double.doubleToLongBits(FastMath.sign(plusZero).getValue()) == 0);
        Assert.assertEquals( 0.0, FastMath.sign(plusZero).getPartialDerivative(1), 1.0e-15);

    }

    @Test
    public void testCeilFloorRintLong() {

        DSFactory factory = new DSFactory(1, 1);
        DerivativeStructure x = factory.variable(0, -1.5);
        Assert.assertEquals(-1.5, x.getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(+1.0, x.getPartialDerivative(1), 1.0e-15);
        Assert.assertEquals(-1.0, FastMath.ceil(x).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(+0.0, FastMath.ceil(x).getPartialDerivative(1), 1.0e-15);
        Assert.assertEquals(-2.0, FastMath.floor(x).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(+0.0, FastMath.floor(x).getPartialDerivative(1), 1.0e-15);
        Assert.assertEquals(-2.0, FastMath.rint(x).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(+0.0, FastMath.rint(x).getPartialDerivative(1), 1.0e-15);
        Assert.assertEquals(-2.0, x.subtract(x.getField().getOne()).rint().getPartialDerivative(0), 1.0e-15);

    }

    @Test
    public void testCopySign() {

        DSFactory factory = new DSFactory(1, 1);
        DerivativeStructure minusOne = factory.variable(0, -1.0);
        Assert.assertEquals(+1.0, FastMath.copySign(minusOne, +1.0).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(-1.0, FastMath.copySign(minusOne, +1.0).getPartialDerivative(1), 1.0e-15);
        Assert.assertEquals(-1.0, FastMath.copySign(minusOne, -1.0).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(+1.0, FastMath.copySign(minusOne, -1.0).getPartialDerivative(1), 1.0e-15);
        Assert.assertEquals(+1.0, FastMath.copySign(minusOne, +0.0).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(-1.0, FastMath.copySign(minusOne, +0.0).getPartialDerivative(1), 1.0e-15);
        Assert.assertEquals(-1.0, FastMath.copySign(minusOne, -0.0).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(+1.0, FastMath.copySign(minusOne, -0.0).getPartialDerivative(1), 1.0e-15);
        Assert.assertEquals(+1.0, FastMath.copySign(minusOne, Double.NaN).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(-1.0, FastMath.copySign(minusOne, Double.NaN).getPartialDerivative(1), 1.0e-15);

        DerivativeStructure plusOne = factory.variable(0, +1.0);
        Assert.assertEquals(+1.0, FastMath.copySign(plusOne, factory.constant(+1.0)).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(+1.0, FastMath.copySign(plusOne, factory.constant(+1.0)).getPartialDerivative(1), 1.0e-15);
        Assert.assertEquals(-1.0, FastMath.copySign(plusOne, factory.constant(-1.0)).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(-1.0, FastMath.copySign(plusOne, factory.constant(-1.0)).getPartialDerivative(1), 1.0e-15);
        Assert.assertEquals(+1.0, FastMath.copySign(plusOne, factory.constant(+0.0)).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(+1.0, FastMath.copySign(plusOne, factory.constant(+0.0)).getPartialDerivative(1), 1.0e-15);
        Assert.assertEquals(-1.0, FastMath.copySign(plusOne, factory.constant(-0.0)).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(-1.0, FastMath.copySign(plusOne, factory.constant(-0.0)).getPartialDerivative(1), 1.0e-15);
        Assert.assertEquals(+1.0, FastMath.copySign(plusOne, factory.constant(Double.NaN)).getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals(+1.0, FastMath.copySign(plusOne, factory.constant(Double.NaN)).getPartialDerivative(1), 1.0e-15);

    }

    @Test
    public void testToDegreesDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                Assert.assertEquals(FastMath.toDegrees(x), dsX.toDegrees().getValue(), epsilon);
                for (int n = 1; n <= maxOrder; ++n) {
                    if (n == 1) {
                        Assert.assertEquals(180 / FastMath.PI, dsX.toDegrees().getPartialDerivative(1), epsilon);
                    } else {
                        Assert.assertEquals(0.0, dsX.toDegrees().getPartialDerivative(n), epsilon);
                    }
                }
            }
        }
    }

    @Test
    public void testToRadiansDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                Assert.assertEquals(FastMath.toRadians(x), dsX.toRadians().getValue(), epsilon);
                for (int n = 1; n <= maxOrder; ++n) {
                    if (n == 1) {
                        Assert.assertEquals(FastMath.PI / 180, dsX.toRadians().getPartialDerivative(1), epsilon);
                    } else {
                        Assert.assertEquals(0.0, dsX.toRadians().getPartialDerivative(n), epsilon);
                    }
                }
            }
        }
    }

    @Test
    public void testDegRad() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure rebuiltX = dsX.toDegrees().toRadians();
                DerivativeStructure zero = rebuiltX.subtract(dsX);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n), epsilon);
                }
            }
        }
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testComposeMismatchedDimensions() {
        new DSFactory(1, 3).variable(0, 1.2).compose(new double[3]);
    }

    @Test
    public void testCompose() {
        double[] epsilon = new double[] { 1.0e-20, 5.0e-14, 2.0e-13, 3.0e-13, 2.0e-13, 1.0e-20 };
        PolynomialFunction poly =
                new PolynomialFunction(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 });
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            PolynomialFunction[] p = new PolynomialFunction[maxOrder + 1];
            p[0] = poly;
            for (int i = 1; i <= maxOrder; ++i) {
                p[i] = p[i - 1].polynomialDerivative();
            }
            for (double x = 0.1; x < 1.2; x += 0.001) {
                DerivativeStructure dsX = factory.variable(0, x);
                DerivativeStructure dsY1 = dsX.getField().getZero();
                for (int i = poly.degree(); i >= 0; --i) {
                    dsY1 = dsY1.multiply(dsX).add(poly.getCoefficients()[i]);
                }
                double[] f = new double[maxOrder + 1];
                for (int i = 0; i < f.length; ++i) {
                    f[i] = p[i].value(x);
                }
                DerivativeStructure dsY2 = dsX.compose(f);
                DerivativeStructure zero = dsY1.subtract(dsY2);
                for (int n = 0; n <= maxOrder; ++n) {
                    Assert.assertEquals(0.0, zero.getPartialDerivative(n), epsilon[n]);
                }
            }
        }
    }

    @Test
    public void testIntegration() {
        // check that first-order integration on two variables does not depend on sequence of operations
        final RandomGenerator random = new Well19937a(0x87bb96d6e11557bdl);
        final DSFactory factory = new DSFactory(3, 7);
        final int size = factory.getCompiler().getSize();
        for (int count = 0; count < 100; ++count) {
            final double[] data = new double[size];
            for (int i = 0; i < size; i++) {
                data[i] = random.nextDouble();
            }
            final DerivativeStructure f       = factory.build(data);
            final DerivativeStructure i2fIxIy = f.integrate(0, 1).integrate(1, 1);
            final DerivativeStructure i2fIyIx = f.integrate(1, 1).integrate(0, 1);
            checkEquals(i2fIxIy, i2fIyIx, 0.);
        }
    }

    @Test
    public void testIntegrationGreaterThanOrder() {
        // check that integration to a too high order generates zero
        // as integration constants are set to zero
        final RandomGenerator random = new Well19937a(0x4744a847b11e4c6fl);
        final DSFactory factory = new DSFactory(3, 7);
        final int size = factory.getCompiler().getSize();
        for (int count = 0; count < 100; ++count) {
            final double[] data = new double[size];
            for (int i = 0; i < size; i++) {
                data[i] = random.nextDouble();
            }
            final DerivativeStructure f = factory.build(data);
            for (int index = 0; index < factory.getCompiler().getFreeParameters(); ++index) {
                final DerivativeStructure integ = f.integrate(index, factory.getCompiler().getOrder() + 1);
                checkEquals(factory.constant(0), integ, 0.);
            }
        }
    }

    @Test
    public void testIntegrationNoOp() {
        // check that integration of order 0 is no-op
        final RandomGenerator random = new Well19937a(0x75a35152f30f644bl);
        final DSFactory factory = new DSFactory(3, 7);
        final int size = factory.getCompiler().getSize();
        for (int count = 0; count < 100; ++count) {
            final double[] data = new double[size];
            for (int i = 0; i < size; i++) {
                data[i] = random.nextDouble();
            }
            final DerivativeStructure f = factory.build(data);
            for (int index = 0; index < factory.getCompiler().getFreeParameters(); ++index) {
                final DerivativeStructure integ = f.integrate(index, 0);
                checkEquals(f, integ, 0.);
            }
        }
    }

    @Test
    public void testDifferentiationNoOp() {
        // check that differentiation of order 0 is no-op
        final RandomGenerator random = new Well19937a(0x3b6ae4c2f1282949l);
        final DSFactory factory = new DSFactory(3, 7);
        final int size = factory.getCompiler().getSize();
        for (int count = 0; count < 100; ++count) {
            final double[] data = new double[size];
            for (int i = 0; i < size; i++) {
                data[i] = random.nextDouble();
            }
            final DerivativeStructure f = factory.build(data);
            for (int index = 0; index < factory.getCompiler().getFreeParameters(); ++index) {
                final DerivativeStructure integ = f.differentiate(index, 0);
                checkEquals(f, integ, 0.);
            }
        }
    }

    @Test
    public void testIntegrationDifferentiation() {
        // check that integration and differentiation for univariate functions are each other inverse except for constant
        // term and highest order one
        final RandomGenerator random = new Well19937a(0x67fe66c05e5ee222l);
        final DSFactory factory = new DSFactory(1, 25);
        final int size = factory.getCompiler().getSize();
        for (int count = 0; count < 100; ++count) {
            final double[] data = new double[size];
            for (int i = 1; i < size - 1; i++) {
                data[i] = random.nextDouble();
            }
            final int indexVar = 0;
            final DerivativeStructure f = factory.build(data);
            final DerivativeStructure f2 = f.integrate(indexVar, 1).differentiate(indexVar, 1);
            final DerivativeStructure f3 = f.differentiate(indexVar, 1).integrate(indexVar, 1);
            checkEquals(f2, f, 0.);
            checkEquals(f2, f3, 0.);
            // check special case when non-positive integration order actually returns differentiation
            final DerivativeStructure df = f.integrate(indexVar, -1);
            final DerivativeStructure df2 = f.differentiate(indexVar, 1);
            checkEquals(df, df2, 0.);
            // check special case when non-positive differentiation order actually returns integration
            final DerivativeStructure fi  = f.differentiate(indexVar, -1);
            final DerivativeStructure fi2 = f.integrate(indexVar, 1);
            checkEquals(fi, fi2, 0.);
        }
    }

    @Test
    public void testDifferentiation1() {
        // check differentiation operator with result obtained manually
        final int freeParam = 3;
        final int order = 5;
        final DSFactory factory = new DSFactory(freeParam, order);
        final DerivativeStructure f = factory.variable(0, 1.0);
        final int[] orders = new int[freeParam];
        orders[0] = 2;
        orders[1] = 1;
        orders[2] = 1;
        final double value = 10.;
        f.setDerivativeComponent(factory.getCompiler().getPartialDerivativeIndex(orders), value);
        final DerivativeStructure dfDx = f.differentiate(0, 1);
        orders[0] -= 1;
        Assert.assertEquals(1., dfDx.getPartialDerivative(new int[freeParam]), 0.);
        Assert.assertEquals(value, dfDx.getPartialDerivative(orders), 0.);
        checkEquals(factory.constant(0), f.differentiate(0, order + 1), 0.);
    }

    @Test
    public void testDifferentiation2() {
        // check that first-order differentiation twice is same as second-order differentiation
        final RandomGenerator random = new Well19937a(0xec293aaee352de94l);
        final DSFactory factory = new DSFactory(5, 4);
        final int size = factory.getCompiler().getSize();
        for (int count = 0; count < 100; ++count) {
            final double[] data = new double[size];
            for (int i = 0; i < size; i++) {
                data[i] = random.nextDouble();
            }
            final DerivativeStructure f = factory.build(data);
            final DerivativeStructure d2fDx2 = f.differentiate(0, 1).differentiate(0, 1);
            final DerivativeStructure d2fDx2Bis = f.differentiate(0, 2);
            checkEquals(d2fDx2, d2fDx2Bis, 0.);
        }
    }

    @Test
    public void testDifferentiation3() {
        // check that first-order differentiation on two variables does not depend on sequence of operations
        final RandomGenerator random = new Well19937a(0x35409ecc1348e46cl);
        final DSFactory factory = new DSFactory(3, 7);
        final int size = factory.getCompiler().getSize();
        for (int count = 0; count < 100; ++count) {
            final double[] data = new double[size];
            for (int i = 0; i < size; i++) {
                data[i] = random.nextDouble();
            }
            final DerivativeStructure f = factory.build(data);
            final DerivativeStructure d2fDxDy = f.differentiate(0, 1).differentiate(1, 1);
            final DerivativeStructure d2fDyDx = f.differentiate(1, 1).differentiate(0, 1);
            checkEquals(d2fDxDy, d2fDyDx, 0.);
        }
    }

    @Test
    public void testField() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            DSFactory factory = new DSFactory(3, maxOrder);
            DerivativeStructure x = factory.variable(0, 1.0);
            checkF0F1(x.getField().getZero(), 0.0, 0.0, 0.0, 0.0);
            checkF0F1(x.getField().getOne(), 1.0, 0.0, 0.0, 0.0);
            Assert.assertEquals(maxOrder, x.getField().getZero().getOrder());
            Assert.assertEquals(3, x.getField().getZero().getFreeParameters());
            Assert.assertEquals(DerivativeStructure.class, x.getField().getRuntimeClass());
        }
    }

    @Test
    public void testOneParameterConstructor() {
        double x = 1.2;
        double cos = FastMath.cos(x);
        double sin = FastMath.sin(x);
        DSFactory factory = new DSFactory(1, 4);
        DerivativeStructure yRef = factory.variable(0, x).cos();
        try {
            new DSFactory(1, 4).build(0.0, 0.0);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            // expected
        } catch (Exception e) {
            Assert.fail("wrong exceptionc caught " + e.getClass().getName());
        }
        double[] derivatives = new double[] { cos, -sin, -cos, sin, cos };
        DerivativeStructure y = factory.build(derivatives);
        checkEquals(yRef, y, 1.0e-15);
        UnitTestUtils.assertEquals(derivatives, y.getAllDerivatives(), 1.0e-15);
    }

    @Test
    public void testOneOrderConstructor() {
        DSFactory factory = new DSFactory(3, 1);
        double x =  1.2;
        double y =  2.4;
        double z = 12.5;
        DerivativeStructure xRef = factory.variable(0, x);
        DerivativeStructure yRef = factory.variable(1, y);
        DerivativeStructure zRef = factory.variable(2, z);
        try {
            new DSFactory(3, 1).build(x + y - z, 1.0, 1.0);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            // expected
        } catch (Exception e) {
            Assert.fail("wrong exceptionc caught " + e.getClass().getName());
        }
        double[] derivatives = new double[] { x + y - z, 1.0, 1.0, -1.0 };
        DerivativeStructure t = factory.build(derivatives);
        checkEquals(xRef.add(yRef.subtract(zRef)), t, 1.0e-15);
        UnitTestUtils.assertEquals(derivatives, xRef.add(yRef.subtract(zRef)).getAllDerivatives(), 1.0e-15);
    }

    @Test
    public void testLinearCombination1DSDS() {
        DSFactory factory = new DSFactory(6, 1);
        final DerivativeStructure[] a = new DerivativeStructure[] {
            factory.variable(0, -1321008684645961.0 / 268435456.0),
            factory.variable(1, -5774608829631843.0 / 268435456.0),
            factory.variable(2, -7645843051051357.0 / 8589934592.0)
        };
        final DerivativeStructure[] b = new DerivativeStructure[] {
            factory.variable(3, -5712344449280879.0 / 2097152.0),
            factory.variable(4, -4550117129121957.0 / 2097152.0),
            factory.variable(5, 8846951984510141.0 / 131072.0)
        };

        final DerivativeStructure abSumInline = a[0].linearCombination(a[0], b[0], a[1], b[1], a[2], b[2]);
        final DerivativeStructure abSumArray = a[0].linearCombination(a, b);

        Assert.assertEquals(abSumInline.getValue(), abSumArray.getValue(), 0);
        Assert.assertEquals(-1.8551294182586248737720779899, abSumInline.getValue(), 1.0e-15);
        Assert.assertEquals(b[0].getValue(), abSumInline.getPartialDerivative(1, 0, 0, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(b[1].getValue(), abSumInline.getPartialDerivative(0, 1, 0, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(b[2].getValue(), abSumInline.getPartialDerivative(0, 0, 1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(a[0].getValue(), abSumInline.getPartialDerivative(0, 0, 0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(a[1].getValue(), abSumInline.getPartialDerivative(0, 0, 0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(a[2].getValue(), abSumInline.getPartialDerivative(0, 0, 0, 0, 0, 1), 1.0e-15);

    }

    @Test
    public void testLinearCombination1DoubleDS() {
        DSFactory factory = new DSFactory(3, 1);
        final double[] a = new double[] {
            -1321008684645961.0 / 268435456.0,
            -5774608829631843.0 / 268435456.0,
            -7645843051051357.0 / 8589934592.0
        };
        final DerivativeStructure[] b = new DerivativeStructure[] {
            factory.variable(0, -5712344449280879.0 / 2097152.0),
            factory.variable(1, -4550117129121957.0 / 2097152.0),
            factory.variable(2, 8846951984510141.0 / 131072.0)
        };

        final DerivativeStructure abSumInline = b[0].linearCombination(a[0], b[0],
                                                                       a[1], b[1],
                                                                       a[2], b[2]);
        final DerivativeStructure abSumArray = b[0].linearCombination(a, b);

        Assert.assertEquals(abSumInline.getValue(), abSumArray.getValue(), 0);
        Assert.assertEquals(-1.8551294182586248737720779899, abSumInline.getValue(), 1.0e-15);
        Assert.assertEquals(a[0], abSumInline.getPartialDerivative(1, 0, 0), 1.0e-15);
        Assert.assertEquals(a[1], abSumInline.getPartialDerivative(0, 1, 0), 1.0e-15);
        Assert.assertEquals(a[2], abSumInline.getPartialDerivative(0, 0, 1), 1.0e-15);

    }

    @Test
    public void testLinearCombination2DSDS() {
        // we compare accurate versus naive dot product implementations
        // on regular vectors (i.e. not extreme cases like in the previous test)
        Well1024a random = new Well1024a(0xc6af886975069f11l);

        DSFactory factory = new DSFactory(4, 1);
        for (int i = 0; i < 10000; ++i) {
            final DerivativeStructure[] u = new DerivativeStructure[factory.getCompiler().getFreeParameters()];
            final DerivativeStructure[] v = new DerivativeStructure[factory.getCompiler().getFreeParameters()];
            for (int j = 0; j < u.length; ++j) {
                u[j] = factory.variable(j, 1e17 * random.nextDouble());
                v[j] = factory.constant(1e17 * random.nextDouble());
            }

            DerivativeStructure lin = u[0].linearCombination(u[0], v[0], u[1], v[1]);
            double ref = u[0].getValue() * v[0].getValue() +
                         u[1].getValue() * v[1].getValue();
            Assert.assertEquals(ref, lin.getValue(), 1.0e-15 * FastMath.abs(ref));
            Assert.assertEquals(v[0].getValue(), lin.getPartialDerivative(1, 0, 0, 0), 1.0e-15 * FastMath.abs(v[0].getValue()));
            Assert.assertEquals(v[1].getValue(), lin.getPartialDerivative(0, 1, 0, 0), 1.0e-15 * FastMath.abs(v[1].getValue()));

            lin = u[0].linearCombination(u[0], v[0], u[1], v[1], u[2], v[2]);
            ref = u[0].getValue() * v[0].getValue() +
                  u[1].getValue() * v[1].getValue() +
                  u[2].getValue() * v[2].getValue();
            Assert.assertEquals(ref, lin.getValue(), 1.0e-15 * FastMath.abs(ref));
            Assert.assertEquals(v[0].getValue(), lin.getPartialDerivative(1, 0, 0, 0), 1.0e-15 * FastMath.abs(v[0].getValue()));
            Assert.assertEquals(v[1].getValue(), lin.getPartialDerivative(0, 1, 0, 0), 1.0e-15 * FastMath.abs(v[1].getValue()));
            Assert.assertEquals(v[2].getValue(), lin.getPartialDerivative(0, 0, 1, 0), 1.0e-15 * FastMath.abs(v[2].getValue()));

            lin = u[0].linearCombination(u[0], v[0], u[1], v[1], u[2], v[2], u[3], v[3]);
            ref = u[0].getValue() * v[0].getValue() +
                  u[1].getValue() * v[1].getValue() +
                  u[2].getValue() * v[2].getValue() +
                  u[3].getValue() * v[3].getValue();
            Assert.assertEquals(ref, lin.getValue(), 1.0e-15 * FastMath.abs(ref));
            Assert.assertEquals(v[0].getValue(), lin.getPartialDerivative(1, 0, 0, 0), 1.0e-15 * FastMath.abs(v[0].getValue()));
            Assert.assertEquals(v[1].getValue(), lin.getPartialDerivative(0, 1, 0, 0), 1.0e-15 * FastMath.abs(v[1].getValue()));
            Assert.assertEquals(v[2].getValue(), lin.getPartialDerivative(0, 0, 1, 0), 1.0e-15 * FastMath.abs(v[2].getValue()));
            Assert.assertEquals(v[3].getValue(), lin.getPartialDerivative(0, 0, 0, 1), 1.0e-15 * FastMath.abs(v[3].getValue()));

        }
    }

    @Test
    public void testLinearCombination2DoubleDS() {
        // we compare accurate versus naive dot product implementations
        // on regular vectors (i.e. not extreme cases like in the previous test)
        Well1024a random = new Well1024a(0xc6af886975069f11l);

        DSFactory factory = new DSFactory(4, 1);
        for (int i = 0; i < 10000; ++i) {
            final double[] u = new double[4];
            final DerivativeStructure[] v = new DerivativeStructure[factory.getCompiler().getFreeParameters()];
            for (int j = 0; j < u.length; ++j) {
                u[j] = 1e17 * random.nextDouble();
                v[j] = factory.variable(j, 1e17 * random.nextDouble());
            }

            DerivativeStructure lin = v[0].linearCombination(u[0], v[0], u[1], v[1]);
            double ref = u[0] * v[0].getValue() +
                         u[1] * v[1].getValue();
            Assert.assertEquals(ref, lin.getValue(), 1.0e-15 * FastMath.abs(ref));
            Assert.assertEquals(u[0], lin.getPartialDerivative(1, 0, 0, 0), 1.0e-15 * FastMath.abs(v[0].getValue()));
            Assert.assertEquals(u[1], lin.getPartialDerivative(0, 1, 0, 0), 1.0e-15 * FastMath.abs(v[1].getValue()));

            lin = v[0].linearCombination(u[0], v[0], u[1], v[1], u[2], v[2]);
            ref = u[0] * v[0].getValue() +
                  u[1] * v[1].getValue() +
                  u[2] * v[2].getValue();
            Assert.assertEquals(ref, lin.getValue(), 1.0e-15 * FastMath.abs(ref));
            Assert.assertEquals(u[0], lin.getPartialDerivative(1, 0, 0, 0), 1.0e-15 * FastMath.abs(v[0].getValue()));
            Assert.assertEquals(u[1], lin.getPartialDerivative(0, 1, 0, 0), 1.0e-15 * FastMath.abs(v[1].getValue()));
            Assert.assertEquals(u[2], lin.getPartialDerivative(0, 0, 1, 0), 1.0e-15 * FastMath.abs(v[2].getValue()));

            lin = v[0].linearCombination(u[0], v[0], u[1], v[1], u[2], v[2], u[3], v[3]);
            ref = u[0] * v[0].getValue() +
                  u[1] * v[1].getValue() +
                  u[2] * v[2].getValue() +
                  u[3] * v[3].getValue();
            Assert.assertEquals(ref, lin.getValue(), 1.0e-15 * FastMath.abs(ref));
            Assert.assertEquals(u[0], lin.getPartialDerivative(1, 0, 0, 0), 1.0e-15 * FastMath.abs(v[0].getValue()));
            Assert.assertEquals(u[1], lin.getPartialDerivative(0, 1, 0, 0), 1.0e-15 * FastMath.abs(v[1].getValue()));
            Assert.assertEquals(u[2], lin.getPartialDerivative(0, 0, 1, 0), 1.0e-15 * FastMath.abs(v[2].getValue()));
            Assert.assertEquals(u[3], lin.getPartialDerivative(0, 0, 0, 1), 1.0e-15 * FastMath.abs(v[3].getValue()));

        }
    }

    @Test
    public void testSerialization() {
        DerivativeStructure a = new DSFactory(3, 2).variable(0, 1.3);
        DerivativeStructure b = (DerivativeStructure) UnitTestUtils.serializeAndRecover(a);
        Assert.assertEquals(a.getFreeParameters(), b.getFreeParameters());
        Assert.assertEquals(a.getOrder(), b.getOrder());
        checkEquals(a, b, 1.0e-15);
    }

    @Test
    public void testZero() {
        DerivativeStructure zero = new DSFactory(3, 2).variable(2, 17.0).getField().getZero();
        double[] a = zero.getAllDerivatives();
        Assert.assertEquals(10, a.length);
        for (int i = 0; i < a.length; ++i) {
            Assert.assertEquals(0.0, a[i], 1.0e-15);
        }
    }

    @Test
    public void testOne() {
        DerivativeStructure one = new DSFactory(3, 2).variable(2, 17.0).getField().getOne();
        double[] a = one.getAllDerivatives();
        Assert.assertEquals(10, a.length);
        for (int i = 0; i < a.length; ++i) {
            Assert.assertEquals(i == 0 ? 1.0 : 0.0, a[i], 1.0e-15);
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
            map.put(new DSFactory(parameters, order).constant(17.0).getField(), 0);
        }

        // despite we have created numerous factories,
        // there should be only one field for each pair parameters/order
        Assert.assertEquals(pairs.size(), map.size());
        @SuppressWarnings("unchecked")
        Field<DerivativeStructure> first = (Field<DerivativeStructure>) map.entrySet().iterator().next().getKey();
        Assert.assertTrue(first.equals(first));
        Assert.assertFalse(first.equals(Binary64Field.getInstance()));

    }

    @Test
    public void testRebaseConditions() {
        final DSFactory f32 = new DSFactory(3, 2);
        final DSFactory f22 = new DSFactory(2, 2);
        final DSFactory f31 = new DSFactory(3, 1);
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

    @Test
    public void testRebaseNoVariables() {
        final DerivativeStructure x = new DSFactory(0, 2).constant(1.0);
        Assert.assertSame(x, x.rebase());
    }

    @Test
    public void testRebaseValueMoreIntermediateThanBase() {
        doTestRebaseValue(createBaseVariables(new DSFactory(2, 4), 1.5, -2.0),
                          q -> new DerivativeStructure[] {
                              q[0].add(q[1].multiply(3)),
                              q[0].log(),
                              q[1].divide(q[0].sin())
                          },
                          new DSFactory(3, 4),
                          p -> p[0].add(p[1].divide(p[2])),
                          1.0e-15);
    }

    @Test
    public void testRebaseValueLessIntermediateThanBase() {
        doTestRebaseValue(createBaseVariables(new DSFactory(3, 4), 1.5, -2.0, 0.5),
                          q -> new DerivativeStructure[] {
                              q[0].add(q[1].multiply(3)),
                              q[0].add(q[1]).subtract(q[2])
                          },
                          new DSFactory(2, 4),
                          p -> p[0].multiply(p[1]),
                          1.0e-15);
    }

    @Test
    public void testRebaseValueEqualIntermediateAndBase() {
        doTestRebaseValue(createBaseVariables(new DSFactory(2, 4), 1.5, -2.0),
                          q -> new DerivativeStructure[] {
                              q[0].add(q[1].multiply(3)),
                              q[0].add(q[1])
                          },
                          new DSFactory(2, 4),
                          p -> p[0].multiply(p[1]),
                          1.0e-15);
    }

    private void doTestRebaseValue(final DerivativeStructure[] q,
                                   final CalculusFieldMultivariateVectorFunction<DerivativeStructure> qToP,
                                   final DSFactory factoryP,
                                   final CalculusFieldMultivariateFunction<DerivativeStructure> f,
                                   final double tol) {

        // intermediate variables as functions of base variables
        final DerivativeStructure[] pBase = qToP.value(q);
        
        // reference function
        final DerivativeStructure ref = f.value(pBase);

        // intermediate variables as independent variables
        final DerivativeStructure[] pIntermediate = creatIntermediateVariables(factoryP, pBase);

        // function of the intermediate variables
        final DerivativeStructure fI = f.value(pIntermediate);

        // function rebased to base variables
        final DerivativeStructure rebased = fI.rebase(pBase);

        Assert.assertEquals(q[0].getFreeParameters(),                   ref.getFreeParameters());
        Assert.assertEquals(q[0].getOrder(),                            ref.getOrder());
        Assert.assertEquals(factoryP.getCompiler().getFreeParameters(), fI.getFreeParameters());
        Assert.assertEquals(factoryP.getCompiler().getOrder(),          fI.getOrder());
        Assert.assertEquals(ref.getFreeParameters(),                    rebased.getFreeParameters());
        Assert.assertEquals(ref.getOrder(),                             rebased.getOrder());

        checkEquals(ref, rebased, tol);

        // compare with Taylor map based implementation
        checkEquals(composeWithTaylorMap(fI, pBase), rebased, tol);

    }

    @Test
    public void testOrdersSum() {
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 4; ++j) {
                DSCompiler compiler = DSCompiler.getCompiler(i, j);
                for (int k = 0; k < compiler.getSize(); ++k) {
                    Assert.assertEquals(IntStream.of(compiler.getPartialDerivativeOrders(k)).sum(),
                                        compiler.getPartialDerivativeOrdersSum(k));
                }
            }
        }
    }

    @Test
    public void testDivisionVersusAlternativeImplementation() {
        final DSFactory factory = new DSFactory(3, 10);
        final DerivativeStructure lhs = FastMath.cos(factory.variable(1, 2.5));
        final DerivativeStructure rhs = factory.variable(2, -4).multiply(factory.variable(0, -3.));
        compareDivisionToVersionViaPower(lhs, rhs, 1e-12);
    }

    @Test
    public void testReciprocalVersusAlternativeImplementation() {
        final DSFactory factory = new DSFactory(2, 15);
        final DerivativeStructure operand = factory.variable(0, 1.).
                add(factory.variable(1, 0.).multiply(2.));
        compareReciprocalToVersionViaPower(operand, 1e-15);
    }

    @Test
    public void testSqrtVersusRootN() {
        final DSFactory factory = new DSFactory(2, 8);
        DerivativeStructure ds = factory.variable(1, -2.).multiply(factory.variable(0, 1.));
        Assert.assertArrayEquals(ds.sqrt().getAllDerivatives(), ds.rootN(2).getAllDerivatives(), 1e-10);
    }

    @Test
    public void testRunTimeClass() {
        Field<DerivativeStructure> field = new DSFactory(3, 2).constant(0.0).getField();
        Assert.assertEquals(DerivativeStructure.class, field.getRuntimeClass());
    }

    private void checkF0F1(DerivativeStructure ds, double value, double...derivatives) {

        // check dimension
        Assert.assertEquals(derivatives.length, ds.getFreeParameters());

        // check value, directly and also as 0th order derivative
        Assert.assertEquals(value, ds.getValue(), 1.0e-15);
        Assert.assertEquals(value, ds.getPartialDerivative(new int[ds.getFreeParameters()]), 1.0e-15);

        // check first order derivatives
        for (int i = 0; i < derivatives.length; ++i) {
            int[] orders = new int[derivatives.length];
            orders[i] = 1;
            Assert.assertEquals(derivatives[i], ds.getPartialDerivative(orders), 1.0e-15);
        }

    }

    public static void checkEquals(DerivativeStructure ds1, DerivativeStructure ds2, double epsilon) {

        // check dimension
        Assert.assertEquals(ds1.getFreeParameters(), ds2.getFreeParameters());
        Assert.assertEquals(ds1.getOrder(), ds2.getOrder());

        int[] derivatives = new int[ds1.getFreeParameters()];
        int sum = 0;
        while (true) {

            if (sum <= ds1.getOrder()) {
                Assert.assertEquals(ds1.getPartialDerivative(derivatives),
                                    ds2.getPartialDerivative(derivatives),
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

    /** Compose the present derivatives on the right with a compatible so-called Taylor map i.e. an array of other
     * partial derivatives. The output has the same number of independent variables than the right-hand side. */
    private DerivativeStructure composeWithTaylorMap(final DerivativeStructure lhs, final DerivativeStructure[] rhs) {

        // turn right-hand side of composition into Taylor expansions without constant term
        final DSFactory rhsFactory = rhs[0].getFactory();
        final TaylorExpansion[] rhsAsExpansions = new TaylorExpansion[rhs.length];
        for (int k = 0; k < rhs.length; k++) {
            final DerivativeStructure copied = new DerivativeStructure(rhsFactory, rhs[k].getAllDerivatives());
            copied.setDerivativeComponent(0, 0.);
            rhsAsExpansions[k] = new TaylorExpansion(copied);
        }
        // turn left-hand side of composition into Taylor expansion
        final TaylorExpansion lhsAsExpansion = new TaylorExpansion(lhs);

        // initialize quantities
        TaylorExpansion te = new TaylorExpansion(rhsFactory.constant(lhs.getValue()));
        TaylorExpansion[][] powers = new TaylorExpansion[rhs.length][lhs.getOrder()];  // for lazy storage of powers

        // compose the Taylor expansions
        final double[] coefficients = lhsAsExpansion.coefficients;
        for (int j = 1; j < coefficients.length; j++) {
            if (coefficients[j] != 0.) {  // filter out null terms
                TaylorExpansion inter = new TaylorExpansion(rhsFactory.constant(coefficients[j]));
                final int[] orders = lhs.getFactory().getCompiler().getPartialDerivativeOrders(j);
                for (int i = 0; i < orders.length; i++) {
                    if (orders[i] != 0) {  // only consider non-trivial powers
                        if (powers[i][orders[i] - 1] == null) {
                            // this power has not been computed yet
                            final DerivativeStructure ds = new DerivativeStructure(rhsFactory, rhs[i].getAllDerivatives());
                            ds.setDerivativeComponent(0, 0.);
                            TaylorExpansion inter2 = new TaylorExpansion(ds);
                            for (int k = 1; k < orders[i]; k++) {
                                inter2 = inter2.multiply(rhsAsExpansions[i]);
                            }
                            powers[i][orders[i] - 1] = inter2;
                        }
                        inter = inter.multiply(powers[i][orders[i] - 1]);
                    }
                }
                te = te.add(inter);
            }
        }

        // convert into derivatives object
        return te.buildDsEquivalent();
    }

    /** Class to map partial derivatives to corresponding Taylor expansion. */
    private static class TaylorExpansion {

        /** Polynomial coefficients of the Taylor expansion in the local canonical basis. */
        final double[] coefficients;

        final double[] factorials;
        final DSFactory dsFactory;

       /** Constructor. */
        TaylorExpansion(final DerivativeStructure ds) {
            final double[] data = ds.getAllDerivatives();
            this.dsFactory = ds.getFactory();
            this.coefficients = new double[data.length];

            // compute relevant factorials (would be more efficient to compute products of factorials to map Taylor
            // expansions and partial derivatives)
            this.factorials = new double[ds.getOrder() + 1];
            Arrays.fill(this.factorials, 1.);
            for (int i = 2; i < this.factorials.length; i++) {
                this.factorials[i] = this.factorials[i - 1] * (double) (i);  // avoid limit of 20! in ArithmeticUtils
            }

            // transform partial derivatives into coefficients of Taylor expansion
            for (int j = 0; j < data.length; j++) {
                this.coefficients[j] = data[j];
                if (this.coefficients[j] != 0.) {
                    int[] orders = ds.getFactory().getCompiler().getPartialDerivativeOrders(j);
                    for (int order : orders) {
                        this.coefficients[j] /= this.factorials[order];
                    }
                }
            }
        }

        /** Builder for the corresponding {@link DerivativeStructure}. */
        public DerivativeStructure buildDsEquivalent() throws MathIllegalArgumentException {
            final DSCompiler dsc = this.dsFactory.getCompiler();
            final double[] data = new double[this.coefficients.length];
            for (int j = 0; j < data.length; j++) {
                data[j] = this.coefficients[j];
                if (data[j] != 0.) {
                    int[] orders = dsc.getPartialDerivativeOrders(j);
                    for (int order : orders) {
                        data[j] *= this.factorials[order];
                    }
                }
            }
            return new DerivativeStructure(this.dsFactory, data);
        }

        TaylorExpansion add(final TaylorExpansion te) {
            return new TaylorExpansion(this.buildDsEquivalent().add(te.buildDsEquivalent()));
        }

        TaylorExpansion multiply(final TaylorExpansion te) {
            return new TaylorExpansion(this.buildDsEquivalent().multiply(te.buildDsEquivalent()));
        }

    }

    private DerivativeStructure[] createBaseVariables(final DSFactory factory, double... q) {
        final DerivativeStructure[] qDS = new DerivativeStructure[q.length];
        for (int i = 0; i < q.length; ++i) {
            qDS[i] = factory.variable(i, q[i]);
        }
        return qDS;
    }

    private DerivativeStructure[] creatIntermediateVariables(final DSFactory factory, DerivativeStructure... pBase) {
        final DerivativeStructure[] pIntermediate = new DerivativeStructure[pBase.length];
        for (int i = 0; i < pBase.length; ++i) {
            pIntermediate[i] = factory.variable(i, pBase[i].getValue());
        }
        return pIntermediate;
    }

    private void compareDivisionToVersionViaPower(final DerivativeStructure lhs, final DerivativeStructure rhs,
                                                  final double tolerance) {
        DSCompiler compiler = lhs.getFactory().getCompiler();
        final double[] result = new double[compiler.getSize()];
        compiler.multiply(lhs.getAllDerivatives(), 0, rhs.reciprocal().getAllDerivatives(), 0, result, 0);
        final double[] result2 = result.clone();
        compiler.divide(lhs.getAllDerivatives(), 0, rhs.getAllDerivatives(), 0, result2, 0);
        Assert.assertArrayEquals(result, result2, tolerance);
    }

    private void compareReciprocalToVersionViaPower(final DerivativeStructure operand, final double tolerance) {
        DSCompiler compiler = operand.getFactory().getCompiler();
        final double[] result = new double[compiler.getSize()];
        compiler.pow(operand.getAllDerivatives(), 0, -1, result, 0);
        final double[] result2 = result.clone();
        compiler.reciprocal(operand.getAllDerivatives(), 0, result2, 0);
        Assert.assertArrayEquals(result, result2, tolerance);
    }

}
