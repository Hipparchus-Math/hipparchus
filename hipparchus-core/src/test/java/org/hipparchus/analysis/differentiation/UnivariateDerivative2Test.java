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
 * Test for class {@link UnivariateDerivative2}.
 */
public class UnivariateDerivative2Test extends UnivariateDerivativeAbstractTest<UnivariateDerivative2> {

    @Override
    protected UnivariateDerivative2 build(final double x) {
        return new UnivariateDerivative2(x, 1.0, 0.0);
    }

    @Override
    protected int getMaxOrder() {
        return 2;
    }

    @Test
    public void testGetFirstAndSecondDerivative() {
        UnivariateDerivative2 ud1 = new UnivariateDerivative2(-0.5, 2.5, 4.5);
        Assert.assertEquals(-0.5, ud1.getReal(), 1.0e-15);
        Assert.assertEquals(-0.5, ud1.getValue(), 1.0e-15);
        Assert.assertEquals(+2.5, ud1.getFirstDerivative(), 1.0e-15);
        Assert.assertEquals(+4.5, ud1.getSecondDerivative(), 1.0e-15);
    }

    @Test
    public void testConversion() {
        UnivariateDerivative2 udA = new UnivariateDerivative2(-0.5, 2.5, 4.5);
        DerivativeStructure ds = udA.toDerivativeStructure();
        Assert.assertEquals(1, ds.getFreeParameters());
        Assert.assertEquals(2, ds.getOrder());
        Assert.assertEquals(-0.5, ds.getValue(), 1.0e-15);
        Assert.assertEquals(-0.5, ds.getPartialDerivative(0), 1.0e-15);
        Assert.assertEquals( 2.5, ds.getPartialDerivative(1), 1.0e-15);
        Assert.assertEquals( 4.5, ds.getPartialDerivative(2), 1.0e-15);
        UnivariateDerivative2 udB = new UnivariateDerivative2(ds);
        Assert.assertNotSame(udA, udB);
        Assert.assertEquals(udA, udB);
        try {
            new UnivariateDerivative2(new DSFactory(2, 2).variable(0, 1.0));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
        }
        try {
            new UnivariateDerivative2(new DSFactory(1, 1).variable(0, 1.0));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
        }
    }

    @Test
    public void testDoublePow() {
        Assert.assertSame(build(3).getField().getZero(), UnivariateDerivative2.pow(0.0, build(1.5)));
        UnivariateDerivative2 ud = UnivariateDerivative2.pow(2.0, build(1.5));
        DSFactory factory = new DSFactory(1, 2);
        DerivativeStructure ds = factory.constant(2.0).pow(factory.variable(0, 1.5));
        Assert.assertEquals(ds.getValue(), ud.getValue(), 1.0e-15);
        Assert.assertEquals(ds.getPartialDerivative(1), ud.getFirstDerivative(), 1.0e-15);
        Assert.assertEquals(ds.getPartialDerivative(2), ud.getSecondDerivative(), 1.0e-15);
    }

    @Test
    public void testTaylor() {
        Assert.assertEquals(-0.125, new UnivariateDerivative2(1, -3, 4).taylor(0.75), 1.0e-15);
    }

    @Test
    public void testHashcode() {
        Assert.assertEquals(-1025507011, new UnivariateDerivative2(2, 1, -1).hashCode());
    }

    @Test
    public void testEquals() {
        UnivariateDerivative2 ud2 = new UnivariateDerivative2(12, -34, 56);
        Assert.assertEquals(ud2, ud2);
        Assert.assertNotEquals(ud2, "");
        Assert.assertEquals(ud2, new UnivariateDerivative2(12, -34, 56));
        Assert.assertNotEquals(ud2, new UnivariateDerivative2(21, -34, 56));
        Assert.assertNotEquals(ud2, new UnivariateDerivative2(12, -43, 56));
        Assert.assertNotEquals(ud2, new UnivariateDerivative2(12, -34, 65));
        Assert.assertNotEquals(ud2, new UnivariateDerivative2(21, -43, 65));
    }

    @Test
    public void testComparableFirstTerm() {
        // GIVEN
        final UnivariateDerivative2 ud2a = new UnivariateDerivative2(12, -34, 25);
        final UnivariateDerivative2 ud2b = new UnivariateDerivative2(2, 0, 25);
        // WHEN
        final int actualComparison = ud2a.compareTo(ud2b);
        // THEN
        final int expectedComparison = 1;
        Assert.assertEquals(expectedComparison, actualComparison);
    }

    @Test
    public void testComparableSecondTerm() {
        // GIVEN
        final UnivariateDerivative2 ud2a = new UnivariateDerivative2(12, -34, 25);
        final UnivariateDerivative2 ud2b = new UnivariateDerivative2(12, 0, 25);
        // WHEN
        final int actualComparison = ud2a.compareTo(ud2b);
        // THEN
        final int expectedComparison = -1;
        Assert.assertEquals(expectedComparison, actualComparison);
    }

    @Test
    public void testComparableThirdTerm() {
        // GIVEN
        final UnivariateDerivative2 ud2a = new UnivariateDerivative2(12, -34, 25);
        final UnivariateDerivative2 ud2b = new UnivariateDerivative2(12, -34, 25);
        // WHEN
        final int actualComparison = ud2a.compareTo(ud2b);
        // THEN
        final int expectedComparison = 0;
        Assert.assertEquals(expectedComparison, actualComparison);
    }

    @Test
    public void testRunTimeClass() {
        Field<UnivariateDerivative2> field = build(0.0).getField();
        Assert.assertEquals(UnivariateDerivative2.class, field.getRuntimeClass());
    }

}
