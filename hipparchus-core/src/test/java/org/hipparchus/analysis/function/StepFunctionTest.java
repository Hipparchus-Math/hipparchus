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

package org.hipparchus.analysis.function;

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Test for class {@link StepFunction}.
 */
class StepFunctionTest {
    private final double EPS = Math.ulp(1d);

    @Test
    void testPreconditions1() {
        assertThrows(NullArgumentException.class, () -> new StepFunction(null, new double[]{0, -1, -2}));
    }

    @Test
    void testPreconditions2() {
        assertThrows(NullArgumentException.class, () -> new StepFunction(new double[]{0, 1}, null));
    }

    @Test
    void testPreconditions3() {
        assertThrows(MathIllegalArgumentException.class, () -> new StepFunction(new double[]{0}, new double[]{}));
    }

    @Test
    void testPreconditions4() {
        assertThrows(MathIllegalArgumentException.class, () -> new StepFunction(new double[]{}, new double[]{0}));
    }

    @Test
    void testPreconditions5() {
        assertThrows(MathIllegalArgumentException.class, () -> new StepFunction(new double[]{0, 1}, new double[]{0, -1, -2}));
    }

    @Test
    void testPreconditions6() {
        assertThrows(MathIllegalArgumentException.class, () -> new StepFunction(new double[]{1, 0, 1}, new double[]{0, -1, -2}));
    }

    @Test
    void testSomeValues() {
        final double[] x = { -2, -0.5, 0, 1.9, 7.4, 21.3 };
        final double[] y = { 4, -1, -5.5, 0.4, 5.8, 51.2 };

        final UnivariateFunction f = new StepFunction(x, y);

        assertEquals(4, f.value(Double.NEGATIVE_INFINITY), EPS);
        assertEquals(4, f.value(-10), EPS);
        assertEquals(-1, f.value(-0.4), EPS);
        assertEquals(-5.5, f.value(0), EPS);
        assertEquals(0.4, f.value(2), EPS);
        assertEquals(5.8, f.value(10), EPS);
        assertEquals(51.2, f.value(30), EPS);
        assertEquals(51.2, f.value(Double.POSITIVE_INFINITY), EPS);
    }

    @Test
    void testEndpointBehavior() {
        final double[] x = {0, 1, 2, 3};
        final double[] xp = {-8, 1, 2, 3};
        final double[] y = {1, 2, 3, 4};
        final UnivariateFunction f = new StepFunction(x, y);
        final UnivariateFunction fp = new StepFunction(xp, y);
        assertEquals(f.value(-8), fp.value(-8), EPS);
        assertEquals(f.value(-10), fp.value(-10), EPS);
        assertEquals(f.value(0), fp.value(0), EPS);
        assertEquals(f.value(0.5), fp.value(0.5), EPS);
        for (int i = 0; i < x.length; i++) {
           assertEquals(y[i], f.value(x[i]), EPS);
           if (i > 0) {
               assertEquals(y[i - 1], f.value(x[i] - 0.5), EPS);
           } else {
               assertEquals(y[0], f.value(x[i] - 0.5), EPS);
           }
        }
    }

    @Test
    void testHeaviside() {
        final UnivariateFunction h = new StepFunction(new double[] {-1, 0},
                                                          new double[] {0, 1});

        assertEquals(0, h.value(Double.NEGATIVE_INFINITY), 0);
        assertEquals(0, h.value(-Double.MAX_VALUE), 0);
        assertEquals(0, h.value(-2), 0);
        assertEquals(0, h.value(-Double.MIN_VALUE), 0);
        assertEquals(1, h.value(0), 0);
        assertEquals(1, h.value(2), 0);
        assertEquals(1, h.value(Double.POSITIVE_INFINITY), 0);
    }
}
