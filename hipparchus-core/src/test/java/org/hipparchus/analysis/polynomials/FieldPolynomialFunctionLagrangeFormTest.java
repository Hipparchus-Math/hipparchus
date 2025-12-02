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
package org.hipparchus.analysis.polynomials;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.MathArrays;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import static org.junit.jupiter.api.Assertions.*;

class FieldPolynomialFunctionLagrangeFormTest {

    @Test
    void testExceptionSize() {
        // GIVEN
        final Binary64Field field = Binary64Field.getInstance();
        final Binary64[] x = MathArrays.buildArray(field, 1);
        final Binary64[] y = x.clone();
        // WHEN & THEN
        assertThrows(MathIllegalArgumentException.class, () -> new FieldPolynomialFunctionLagrangeForm<>(x, y));
    }

    @Test
    void testExceptionOrder() {
        // GIVEN
        final Binary64Field field = Binary64Field.getInstance();
        final Binary64[] x = MathArrays.buildArray(field, 2);
        x[0] = field.getOne();
        x[1] = field.getZero();
        final Binary64[] y = x.clone();
        // WHEN & THEN
        assertThrows(MathIllegalArgumentException.class, () -> new FieldPolynomialFunctionLagrangeForm<>(x, y));
    }

    @Test
    void testGetter() {
        // GIVEN
        final Binary64Field field = Binary64Field.getInstance();
        Binary64[] x = { field.getOne(), field.getOne().twice() };
        Binary64[] y = { field.getOne().negate(), field.getOne().multiply(4.) };
        final FieldPolynomialFunctionLagrangeForm<Binary64> lagrangeForm = new FieldPolynomialFunctionLagrangeForm<>(x, y);
        // WHEN & THEN
        assertArrayEquals(x, lagrangeForm.getInterpolatingPoints());
        assertArrayEquals(y, lagrangeForm.getInterpolatingValues());
        final Binary64[] expected = MathArrays.buildArray(field, 2);
        expected[0] = field.getOne().multiply(-6.);
        expected[1] = field.getOne().multiply(5.);
        assertArrayEquals(expected, lagrangeForm.getCoefficients());
    }

    @Test
    void testGetCoefficients() {
        // GIVEN
        final Binary64Field field = Binary64Field.getInstance();
        Binary64[] x = { field.getOne(), field.getOne().twice() };
        Binary64[] y = { field.getOne().negate(), field.getOne().multiply(4.) };
        final FieldPolynomialFunctionLagrangeForm<Binary64> lagrangeForm = new FieldPolynomialFunctionLagrangeForm<>(x, y);
        // WHEN
        final Binary64[] coefficients = lagrangeForm.getCoefficients();
        // WHEN & THEN
        assertArrayEquals(lagrangeForm.getCoefficients(), coefficients);
    }

    @ParameterizedTest
    @ValueSource(ints = {2, 3, 4, 5, 6, 7, 8, 9, 10})
    void testValueDegree(final int n) {
        // GIVEN
        final Binary64Field field = Binary64Field.getInstance();
        final Binary64[] x = MathArrays.buildArray(field, n);
        final Binary64[] y = MathArrays.buildArray(field, n);
        for (int i = 0; i < n; i++) {
            x[i] = field.getOne().multiply(i);
            y[i] = field.getOne().multiply(2 * i + 1);
        }
        final FieldPolynomialFunctionLagrangeForm<Binary64> lagrangeForm = new FieldPolynomialFunctionLagrangeForm<>(x, y);
        final Binary64 five = field.getOne().multiply(5);
        // WHEN
        final Binary64 fieldValue = lagrangeForm.value(five);
        // WHEN & THEN
        final double[] xDouble = new double[n];
        final double[] yDouble = new double[n];
        for (int i = 0; i < n; i++) {
            xDouble[i] = x[i].getReal();
            yDouble[i] = y[i].getReal();
        }
        assertEquals(PolynomialFunctionLagrangeForm.evaluate(xDouble, yDouble, five.getReal()),
                fieldValue.getReal());
    }

    @ParameterizedTest
    @ValueSource(doubles = {2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.6, 11})
    void testValuePoint(final double z) {
        // GIVEN
        final Binary64Field field = Binary64Field.getInstance();
        final Binary64[] x = MathArrays.buildArray(field, 10);
        final Binary64[] y = MathArrays.buildArray(field, 10);
        for (int i = 0; i < x.length; i++) {
            x[i] = field.getOne().multiply(i);
            y[i] = field.getOne().multiply(2 * i + 1);
        }
        final FieldPolynomialFunctionLagrangeForm<Binary64> lagrangeForm = new FieldPolynomialFunctionLagrangeForm<>(x, y);
        // WHEN
        final Binary64 fieldValue = lagrangeForm.value(new Binary64(z));
        // WHEN & THEN
        final double[] xDouble = new double[x.length];
        final double[] yDouble = xDouble.clone();
        for (int i = 0; i < xDouble.length; i++) {
            xDouble[i] = x[i].getReal();
            yDouble[i] = y[i].getReal();
        }
        assertEquals(PolynomialFunctionLagrangeForm.evaluate(xDouble, yDouble, z), fieldValue.getReal());
    }
}
