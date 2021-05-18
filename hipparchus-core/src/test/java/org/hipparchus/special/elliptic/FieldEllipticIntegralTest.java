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
package org.hipparchus.special.elliptic;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.Decimal64Field;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class FieldEllipticIntegralTest {

    @Test
    public void testNoConvergence() {
        doTestNoConvergence(Decimal64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestNoConvergence(final Field<T> field) {
        try {
            new EllipticIntegral(Double.NaN).getBigK();
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedCoreFormats.CONVERGENCE_FAILED, mise.getSpecifier());
        }
    }

    @Test
    public void testComplementary() {
        doTestComplementary(Decimal64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestComplementary(final Field<T> field) {
        for (double k = 0.01; k < 1; k += 0.01) {
            FieldEllipticIntegral<T> ei1 = build(field, k);
            FieldEllipticIntegral<T> ei2 = build(field, FastMath.sqrt(1 - k * k));
            Assert.assertEquals(ei1.getBigK().getReal(), ei2.getBigKPrime().getReal(), FastMath.ulp(ei1.getBigK()).getReal());
        }
    }

    @Test
    public void testAbramowitzStegunExample3() {
        doTestAbramowitzStegunExample3(Decimal64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestAbramowitzStegunExample3(final Field<T> field) {
        final FieldEllipticIntegral<T> ei = build(field, FastMath.sqrt(80.0 / 81.0));
        Assert.assertEquals(80.0 / 81.0, ei.getK().multiply(ei.getK()).getReal(), 1.0e-15);
        Assert.assertEquals(3.591545001, ei.getBigK().getReal(), 2.0e-9);
    }

    private <T extends CalculusFieldElement<T>> FieldEllipticIntegral<T> build(final Field<T> field, final double k) {
        return new FieldEllipticIntegral<>(field.getZero().newInstance(k));
    }

}
