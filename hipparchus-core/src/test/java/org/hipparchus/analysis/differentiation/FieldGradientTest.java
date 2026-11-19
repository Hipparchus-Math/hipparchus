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

import org.hipparchus.complex.Complex;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertEquals;

class FieldGradientTest {

    @Test
    void testStackVariable() {
        // GIVEN
        final FieldGradient<Complex> gradient = new FieldGradient<>(Complex.ONE, new Complex(2));
        // WHEN
        final FieldGradient<Complex> gradientWithMoreVariable = gradient.stackVariable();
        // THEN
        Assertions.assertEquals(gradient.getValue(), gradientWithMoreVariable.getValue());
        Assertions.assertEquals(gradient.getFreeParameters() + 1, gradientWithMoreVariable.getFreeParameters());
        Assertions.assertEquals(0., gradientWithMoreVariable.getGradient()[gradient.getFreeParameters()].getReal());
        Assertions.assertArrayEquals(gradient.getGradient(), Arrays.copyOfRange(gradientWithMoreVariable.getGradient(),
                0, gradient.getFreeParameters()));
    }

    @Test
    void testNorm() {
        final FieldGradient<Complex> gradient = new FieldGradient<>(Complex.ONE, new Complex(2));
        assertEquals(3., gradient.norm());
    }
}
