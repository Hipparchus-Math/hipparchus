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
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test for class {@link UnivariateDerivative1}.
 */
public class UnivariateDerivative1Test extends UnivariateDerivativeAbstractTest<UnivariateDerivative1> {

    @Override
    protected UnivariateDerivative1 build(final double x) {
        return new UnivariateDerivative1(x, 1.0);
    }

    @Override
    protected int getMaxOrder() {
        return 1;
    }

    @Test
    public void testGetFirstDerivative() {
        UnivariateDerivative1 ud1 = new UnivariateDerivative1(-0.5, 2.5);
        Assert.assertEquals(-0.5, ud1.getReal(), 1.0e-15);
        Assert.assertEquals(-0.5, ud1.getValue(), 1.0e-15);
        Assert.assertEquals(+2.5, ud1.getFirstDerivative(), 1.0e-15);
    }

    @Test
    public void testConversion() {
        UnivariateDerivative1 udA = new UnivariateDerivative1(-0.5, 2.5);
        DerivativeStructure ds = udA.toDerivativeStructure();
        Assert.assertEquals(1, ds.getFreeParameters());
        Assert.assertEquals(1, ds.getOrder());
        Assert.assertEquals(-0.5, ds.getValue(), 1.0e-15);
        Assert.assertEquals(-0.5, ds.getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals( 2.5, ds.getPartialDerivative(1), 1.0e-15);
        UnivariateDerivative1 udB = new UnivariateDerivative1(ds);
        Assert.assertNotSame(udA, udB);
        Assert.assertEquals(udA, udB);
        try {
            new UnivariateDerivative1(new DSFactory(2, 2).variable(0, 1.0));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
        }
        try {
            new UnivariateDerivative1(new DSFactory(1, 2).variable(0, 1.0));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
        }
    }

    @Test
    public void testDoublePow() {
        Assert.assertSame(build(3).getField().getZero(), UnivariateDerivative1.pow(0.0, build(1.5)));
        UnivariateDerivative1 ud = UnivariateDerivative1.pow(2.0, build(1.5));
        DSFactory factory = new DSFactory(1, 1);
        DerivativeStructure ds = factory.constant(2.0).pow(factory.variable(0, 1.5));
        Assert.assertEquals(ds.getValue(), ud.getValue(), 1.0e-15);
        Assert.assertEquals(ds.getPartialDerivative(1), ud.getFirstDerivative(), 1.0e-15);
    }

    @Test
    public void testTaylor() {
        Assert.assertEquals(2.5, new UnivariateDerivative1(2, 1).taylor(0.5), 1.0e-15);
    }

    @Test
    public void testHashcode() {
        Assert.assertEquals(2108686789, new UnivariateDerivative1(2, 1).hashCode());
    }

    @Test
    public void testEquals() {
        UnivariateDerivative1 ud1 = new UnivariateDerivative1(12, -34);
        Assert.assertEquals(ud1, ud1);
        Assert.assertNotEquals(ud1, "");
        Assert.assertEquals(ud1, new UnivariateDerivative1(12, -34));
        Assert.assertNotEquals(ud1, new UnivariateDerivative1(21, -34));
        Assert.assertNotEquals(ud1, new UnivariateDerivative1(12, -43));
        Assert.assertNotEquals(ud1, new UnivariateDerivative1(21, -43));
    }


    @Test
    public void testComparableFirstTerm() {
        // GIVEN
        final UnivariateDerivative1 ud1a = new UnivariateDerivative1(12, -34);
        final UnivariateDerivative1 ud1b = new UnivariateDerivative1(2, 0);
        // WHEN
        final int actualComparison = ud1a.compareTo(ud1b);
        // THEN
        final int expectedComparison = 1;
        Assert.assertEquals(expectedComparison, actualComparison);
    }

    @Test
    public void testComparableSecondTerm() {
        // GIVEN
        final UnivariateDerivative1 ud1a = new UnivariateDerivative1(12, -34);
        final UnivariateDerivative1 ud1b = new UnivariateDerivative1(12, 0);
        // WHEN
        final int actualComparison = ud1a.compareTo(ud1b);
        // THEN
        final int expectedComparison = -1;
        Assert.assertEquals(expectedComparison, actualComparison);
    }

    @Test
    public void testRunTimeClass() {
        Field<UnivariateDerivative1> field = build(0.0).getField();
        Assert.assertEquals(UnivariateDerivative1.class, field.getRuntimeClass());
    }

}
