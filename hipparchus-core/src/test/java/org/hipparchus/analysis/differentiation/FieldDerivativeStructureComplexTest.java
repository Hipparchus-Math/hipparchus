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

import java.util.Arrays;
import java.util.List;

import org.hipparchus.Field;
import org.hipparchus.UnitTestUtils;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexField;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Precision;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test for class {@link FieldDerivativeStructure} on {@link Complex}.
 */
public class FieldDerivativeStructureComplexTest extends FieldDerivativeStructureAbstractTest<Complex> {

    @Override
    protected Field<Complex> getField() {
        return ComplexField.getInstance();
    }

    @Override
    @Test
    public void testComposeField() {
        doTestComposeField(new double[] { 1.0e-100, 5.0e-14, 2.0e-13, 3.0e-13, 2.0e-13, 1.0e-100 });
    }

    @Override
    @Test
    public void testComposePrimitive() {
        doTestComposePrimitive(new double[] { 1.0e-100, 5.0e-14, 2.0e-13, 3.0e-13, 2.0e-13, 1.0e-100 });
    }

    @Override
    @Test
    public void testHypotNoOverflow() {
        doTestHypotNoOverflow(180);
    }

    @Override
    @Test
    public void testLinearCombinationReference() {
        doTestLinearCombinationReference(x -> build(x), 4.8e-16, 1.0);
    }

    @Override
    @Test
    public void testLinearCombination1DSDS() {
        doTestLinearCombination1DSDS(1.0e-15);
    }

    @Override
    @Test
    public void testLinearCombination1FieldDS() {
        doTestLinearCombination1FieldDS(1.0e-15);
    }

    @Override
    @Test
    public void testLinearCombination1DoubleDS() {
        doTestLinearCombination1DoubleDS(1.0e-15);
    }

    @Override
    @Test
    public void testAtan2() {
        double[] epsilon = new double[] { 9.0e-16, 3.0e-15, 2.9e-14, 1.0e-12, 8.0e-11 };
        for (int maxOrder = 0; maxOrder < 5; ++maxOrder) {
            final FDSFactory<Complex> factory = buildFactory(2, maxOrder);
            for (double x = -1.7; x < 2; x += 0.2) {
                FieldDerivativeStructure<Complex> dsX = factory.variable(0, x);
                for (double y = -1.7; y < 2; y += 0.2) {
                    FieldDerivativeStructure<Complex> dsY = factory.variable(1, y);
                    FieldDerivativeStructure<Complex> atan2 = FieldDerivativeStructure.atan2(dsY, dsX);
                    FieldDerivativeStructure<Complex> ref = dsY.divide(dsX).atan();
                    if (x < 0) {
                        ref = (y < 0) ? ref.subtract(FastMath.PI) : ref.add(FastMath.PI);
                    }
                    double fullTurns = MathUtils.normalizeAngle(atan2.getValue().getReal(), ref.getValue().getReal()) - atan2.getValue().getReal();
                    atan2 = atan2.add(fullTurns);
                    FieldDerivativeStructure<Complex> zero = atan2.subtract(ref);
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

    @Override
    @Test
    public void testAtan2SpecialCases() {
        Assert.assertTrue(build(+0.0).atan2(build(+0.0)).isNaN());
        Assert.assertTrue(build(-0.0).atan2(build(+0.0)).isNaN());
        Assert.assertTrue(build(+0.0).atan2(build(-0.0)).isNaN());
        Assert.assertTrue(build(-0.0).atan2(build(-0.0)).isNaN());
    }

    @Test
    public void testAtan2SpecialCasesDerivatives() {

        final FDSFactory<Complex> factory = buildFactory(2, 2);
        FieldDerivativeStructure<Complex> pp =
                FieldDerivativeStructure.atan2(factory.variable(1, buildScalar(+0.0)), factory.variable(1, buildScalar(+0.0)));
        Assert.assertTrue(pp.getValue().isNaN());

        FieldDerivativeStructure<Complex> pn =
                FieldDerivativeStructure.atan2(factory.variable(1, buildScalar(+0.0)), factory.variable(1, buildScalar(-0.0)));
        Assert.assertTrue(pn.getValue().isNaN());

        FieldDerivativeStructure<Complex> np =
                FieldDerivativeStructure.atan2(factory.variable(1, buildScalar(-0.0)), factory.variable(1, buildScalar(+0.0)));
        Assert.assertTrue(np.getValue().isNaN());

        FieldDerivativeStructure<Complex> nn =
                FieldDerivativeStructure.atan2(factory.variable(1, buildScalar(-0.0)), factory.variable(1, buildScalar(-0.0)));
        Assert.assertTrue(nn.getValue().isNaN());

    }

    @Override
    @Test
    public void testPowDoubleDS() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {

            final FDSFactory<Complex> factory = buildFactory(3, maxOrder);
            FieldDerivativeStructure<Complex> x = factory.variable(0, 0.1);
            FieldDerivativeStructure<Complex> y = factory.variable(1, 0.2);
            FieldDerivativeStructure<Complex> z = factory.variable(2, 0.3);
            List<FieldDerivativeStructure<Complex>> list = Arrays.asList(x, y, z,
                                                                           x.add(y).add(z),
                                                                           x.multiply(y).multiply(z));

            for (FieldDerivativeStructure<Complex> ds : list) {
                // the special case a = 0 is included here
                for (double a : new double[] { 0.0, 0.1, 1.0, 2.0, 5.0 }) {
                    FieldDerivativeStructure<Complex> reference = (a == 0) ?
                                                    x.getField().getZero() :
                                                    factory.constant(a).pow(ds);
                    FieldDerivativeStructure<Complex> result = FieldDerivativeStructure.pow(a, ds);
                    checkEquals(reference, result, 2.0e-14 * FastMath.abs(reference.getReal()));
                }

            }

            // negative base: -1^x can be evaluated for integers only, so value is sometimes OK, derivatives are always NaN
            FieldDerivativeStructure<Complex> negEvenInteger = FieldDerivativeStructure.pow(-2.0, factory.variable(0, 2.0));
            Assert.assertEquals(4.0, negEvenInteger.getReal(), 1.0e-15);
            Assert.assertTrue(Double.isNaN(negEvenInteger.getPartialDerivative(1, 0, 0).getReal()));
            FieldDerivativeStructure<Complex> negOddInteger = FieldDerivativeStructure.pow(-2.0, factory.variable(0, 3.0));
            Assert.assertEquals(-8.0, negOddInteger.getReal(), 1.0e-15);
            Assert.assertTrue(Double.isNaN(negOddInteger.getPartialDerivative(1, 0, 0).getReal()));
            FieldDerivativeStructure<Complex> negNonInteger = FieldDerivativeStructure.pow(-2.0, factory.variable(0, 2.001));
            Assert.assertEquals(4.0027537969708465469, negNonInteger.getValue().getReal(), 2.0e-15);
            Assert.assertEquals(0.012575063293019489803, negNonInteger.getValue().getImaginary(), 2.0e-15);
            Assert.assertTrue(Double.isNaN(negNonInteger.getPartialDerivative(1, 0, 0).getReal()));

            FieldDerivativeStructure<Complex> zeroNeg = FieldDerivativeStructure.pow(0.0, factory.variable(0, -1.0));
            Assert.assertTrue(Double.isNaN(zeroNeg.getReal()));
            Assert.assertTrue(Double.isNaN(zeroNeg.getPartialDerivative(1, 0, 0).getReal()));
            FieldDerivativeStructure<Complex> posNeg = FieldDerivativeStructure.pow(2.0, factory.variable(0, -2.0));
            Assert.assertEquals(1.0 / 4.0, posNeg.getReal(), 1.0e-15);
            Assert.assertEquals(FastMath.log(2.0) / 4.0, posNeg.getPartialDerivative(1, 0, 0).getReal(), 1.0e-15);

            // very special case: 0^0 where the power is a primitive
            FieldDerivativeStructure<Complex> zeroDsZeroDouble = factory.variable(0, 0.0).pow(0.0);
            boolean first = true;
            for (final Complex d : zeroDsZeroDouble.getAllDerivatives()) {
                if (first) {
                    Assert.assertEquals(1.0, d.getReal(), Precision.EPSILON);
                    first = false;
                } else {
                    Assert.assertEquals(0.0, d.getReal(), Precision.SAFE_MIN);
                }
            }
            FieldDerivativeStructure<Complex> zeroDsZeroInt = factory.variable(0, 0.0).pow(0);
            first = true;
            for (final Complex d : zeroDsZeroInt.getAllDerivatives()) {
                if (first) {
                    Assert.assertEquals(1.0, d.getReal(), Precision.EPSILON);
                    first = false;
                } else {
                    Assert.assertEquals(0.0, d.getReal(), Precision.SAFE_MIN);
                }
            }

            // 0^p with p smaller than 1.0
            FieldDerivativeStructure<Complex> u = factory.variable(1, -0.0).pow(0.25);
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

    @Override
    @Test
    public void testLog10() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            if (x <= 0) {
                Assert.assertTrue(Double.isNaN(FastMath.log10(x)));
                Assert.assertFalse(build(x).log10().getValue().isNaN());
            } else {
                checkRelative(FastMath.log10(x), build(x).log10());
            }
        }
    }

    @Override
    @Test
    public void testRootN() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            for (int n = 1; n < 5; ++n) {
                if (x < 0) {
                    // special case for Complex
                    final double doubleRoot = new Binary64(x).rootN(n).getReal();
                    if (n % 2 == 0) {
                        Assert.assertTrue(Double.isNaN(doubleRoot));
                    } else {
                        Assert.assertTrue(doubleRoot < 0);
                    }
                    Assert.assertEquals(FastMath.PI / n, build(x).rootN(n).getValue().getArgument(), 1.0e-15);
                } else {
                    checkRelative(FastMath.pow(x, 1.0 / n), build(x).rootN(n));
                }
            }
        }
    }

    @Override
    @Test
    public void testRootNSingularity() {
        doTestRootNSingularity(false);
    }

    @Override
    @Test
    public void testCbrt() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            if ( x < 0) {
                // special case for Complex
                Assert.assertTrue(FastMath.cbrt(x) < 0);
                Assert.assertEquals(FastMath.PI / 3, build(x).cbrt().getValue().getArgument(), 1.0e-15);
            } else {
                checkRelative(FastMath.cbrt(x), build(x).cbrt());
            }
        }
    }

    @Test
    public void testCbrtComplex() {
        Complex z = new Complex(15, 2);
        UnitTestUtils.assertEquals(z, z.multiply(z).multiply(z).cbrt(), 1.0e-14);
        Complex branchCutPlus = new Complex(-8.0, +0.0);
        Complex cbrtPlus = branchCutPlus.cbrt();
        UnitTestUtils.assertEquals(branchCutPlus, cbrtPlus.multiply(cbrtPlus).multiply(cbrtPlus), 1.0e-14);
        Assert.assertEquals(1.0, cbrtPlus.getReal(), 1.0e-15);
        Assert.assertEquals(FastMath.sqrt(3.0), cbrtPlus.getImaginary(), 1.0e-15);
        Complex branchCutMinus = new Complex(-8.0, -0.0);
        Complex cbrtMinus = branchCutMinus.cbrt();
        UnitTestUtils.assertEquals(branchCutMinus, cbrtMinus.multiply(cbrtMinus).multiply(cbrtMinus), 1.0e-14);
        Assert.assertEquals(1.0, cbrtMinus.getReal(), 1.0e-15);
        Assert.assertEquals(-FastMath.sqrt(3.0), cbrtMinus.getImaginary(), 1.0e-15);
    }

    @Override
    @Test
    public void testPowField() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            for (double y = 0.1; y < 4; y += 0.2) {
                if ( x < 0) {
                    // special case for Complex
                    Assert.assertTrue(Double.isNaN(FastMath.pow(x, y)));
                    Assert.assertFalse(build(x).pow(build(y)).isNaN());
                } else {
                    checkRelative(FastMath.pow(x, y), build(x).pow(build(y)));
                }
            }
        }
    }

    @Override
    @Test
    public void testPowDouble() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            for (double y = 0.1; y < 4; y += 0.2) {
                if ( x < 0) {
                    // special case for Complex
                    Assert.assertTrue(Double.isNaN(FastMath.pow(x, y)));
                    Assert.assertFalse(build(x).pow(y).isNaN());
                } else {
                    checkRelative(FastMath.pow(x, y), build(x).pow(y));
                }
            }
        }
    }

}
