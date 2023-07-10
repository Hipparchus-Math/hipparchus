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
package org.hipparchus.util;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.MathIllegalStateException;
import org.junit.Assert;
import org.junit.Test;

/**
 * Tests for ContinuedFraction.
 */
public class FieldContinuedFractionTest {

    @Test
    public void testGoldenRatio() {
        FieldContinuedFraction cf = new FieldContinuedFraction() {

            @Override
            public <T extends CalculusFieldElement<T>> T getA(final int n, final T x) {
                return x.getField().getOne();
            }

            @Override
            public <T extends CalculusFieldElement<T>> T getB(final int n, final T x) {
                return x.getField().getOne();
            }

        };

        Binary64 gr = cf.evaluate(new Binary64(0.0), 10e-9);
        Assert.assertEquals(1.61803399, gr.getReal(), 10e-9);
    }

    @Test(expected = MathIllegalStateException.class)
    public void testNonConvergentContinuedFraction() {
        FieldContinuedFraction cf = new FieldContinuedFraction() {

            @Override
            public <T extends CalculusFieldElement<T>> T getA(final int n, final T x) {
                return x.getField().getOne();
            }

            @Override
            public <T extends CalculusFieldElement<T>> T getB(final int n, final T x) {
                return x.getField().getOne();
            }

        };

        cf.evaluate(new Binary64(0.0), 10e-9, 10);
    }

    @Test(expected = MathIllegalStateException.class)
    public void testInfinityDivergence() {
        FieldContinuedFraction cf = new FieldContinuedFraction() {

            @Override
            public <T extends CalculusFieldElement<T>> T getA(final int n, final T x) {
                return x.getField().getOne().divide(n);
            }

            @Override
            public <T extends CalculusFieldElement<T>> T getB(final int n, final T x) {
                return x.getField().getOne();
            }

        };

        cf.evaluate(new Binary64(1));
    }
}
