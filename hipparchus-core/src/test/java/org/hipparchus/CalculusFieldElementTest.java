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

package org.hipparchus;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class CalculusFieldElementTest {

    @Test
    void testMultiplyInt() {
        // GIVEN
        final double value = 3.;
        final int factor = 2;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.multiply(factor);
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(value * factor);
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testAddDouble() {
        // GIVEN
        final double value1 = 1.;
        final double value2 = 2.;
        final TestCalculusFieldElement testElement1 = new TestCalculusFieldElement(value1);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement1.add(value2);
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(value1 + value2);
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testSubtractDouble() {
        // GIVEN
        final double value1 = 1.;
        final double value2 = 2.;
        final TestCalculusFieldElement testElement1 = new TestCalculusFieldElement(value1);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement1.subtract(value2);
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(value1 - value2);
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testMultiplyDouble() {
        // GIVEN
        final double value1 = 3.;
        final double value2 = 2.;
        final TestCalculusFieldElement testElement1 = new TestCalculusFieldElement(value1);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement1.multiply(value2);
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(value1 * value2);
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testDivideDouble() {
        // GIVEN
        final double value1 = 3.;
        final double value2 = 2.;
        final TestCalculusFieldElement testElement1 = new TestCalculusFieldElement(value1);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement1.divide(value2);
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(value1 / value2);
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testSubtract() {
        // GIVEN
        final double value1 = 1.;
        final double value2 = 2.;
        final TestCalculusFieldElement testElement1 = new TestCalculusFieldElement(value1);
        final TestCalculusFieldElement testElement2 = new TestCalculusFieldElement(value2);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement1.subtract(testElement2);
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(value1 - value2);
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testDivide() {
        // GIVEN
        final double value1 = 3.;
        final double value2 = 2.;
        final TestCalculusFieldElement testElement1 = new TestCalculusFieldElement(value1);
        final TestCalculusFieldElement testElement2 = new TestCalculusFieldElement(value2);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement1.divide(testElement2);
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(value1 / value2);
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testSquare() {
        // GIVEN
        final double value = 3.;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.square();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(value * value);
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testSqrt() {
        // GIVEN
        final double value = 3.;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.sqrt();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.sqrt(value));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testCbrt() {
        // GIVEN
        final double value = 3.;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.cbrt();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.cbrt(value));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testRootN() {
        // GIVEN
        final int n = 4;
        final double value = 3.;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.rootN(n);
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.pow(value, 1. / n));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testPowInt() {
        // GIVEN
        final int exponent = 4;
        final double value = 3.;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.pow(exponent);
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.pow(value, exponent));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testPowDouble() {
        // GIVEN
        final double exponent = 4.5;
        final double value = 3.;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.pow(exponent);
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.pow(value, exponent));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testSin() {
        // GIVEN
        final double value = 3.;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.sin();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.sin(value));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testCos() {
        // GIVEN
        final double value = 3.;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.cos();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.cos(value));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testTan() {
        // GIVEN
        final double value = 3.;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.tan();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.tan(value));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testSinh() {
        // GIVEN
        final double value = 3.;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.sinh();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.sinh(value));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testCosh() {
        // GIVEN
        final double value = 3.;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.cosh();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.cosh(value));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testTanh() {
        // GIVEN
        final double value = 3.;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.tanh();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.tanh(value));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testSign() {
        // GIVEN
        final double value = 3.6;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.sign();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.signum(value));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testUlp() {
        // GIVEN
        final double value = 3.6;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.ulp();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.ulp(value));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testFloor() {
        // GIVEN
        final double value = 3.6;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.floor();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.floor(value));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testCeil() {
        // GIVEN
        final double value = 3.6;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.ceil();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.ceil(value));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testRint() {
        // GIVEN
        final double value = 3.6;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.rint();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.rint(value));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testTwice() {
        // GIVEN
        final double value = 3.6;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.twice();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(value * 2);
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testHalf() {
        // GIVEN
        final double value = 3.6;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.half();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(value / 2);
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testToDegrees() {
        // GIVEN
        final double value = 3.;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.toDegrees();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.toDegrees(value));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testToRadians() {
        // GIVEN
        final double value = 3.;
        final TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement.toRadians();
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(FastMath.toRadians(value));
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testLinearCombinationDouble2() {
        // GIVEN
        final double value1 = 3.;
        final double value2 = 2.;
        final TestCalculusFieldElement testElement1 = new TestCalculusFieldElement(value1);
        final TestCalculusFieldElement testElement2 = new TestCalculusFieldElement(value2);
        final double coeff1 = -5.;
        final double coeff2 = 4.;
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement1.linearCombination(coeff1, testElement1,
                coeff2, testElement2);
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(
                coeff1 * value1 + coeff2 * value2);
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testLinearCombinationDouble3() {
        // GIVEN
        final double value1 = 3.;
        final double value2 = 2.;
        final double value3 = -1;
        final TestCalculusFieldElement testElement1 = new TestCalculusFieldElement(value1);
        final TestCalculusFieldElement testElement2 = new TestCalculusFieldElement(value2);
        final TestCalculusFieldElement testElement3 = new TestCalculusFieldElement(value3);
        final double coeff1 = -5.;
        final double coeff2 = 4.;
        final double coeff3 = 6.;
        // WHEN
        final TestCalculusFieldElement actualOperation = testElement1.linearCombination(coeff1, testElement1,
                coeff2, testElement2, coeff3, testElement3);
        // THEN
        final TestCalculusFieldElement expectedOperation = new TestCalculusFieldElement(
                coeff1 * value1 + coeff2 * value2 + coeff3 * value3);
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testGetExponent() {
        // GIVEN
        final double value = 3.5;
        final  TestCalculusFieldElement testElement = new TestCalculusFieldElement(value);
        // WHEN
        final int actualOperation = testElement.getExponent();
        // THEN
        final int expectedOperation = FastMath.getExponent(value);
        assertEquals(expectedOperation, actualOperation);
    }

    @Test
    void testIsFinite() {
        assertTrue(new TestCalculusFieldElement(1.).isFinite());
        assertFalse(new TestCalculusFieldElement(Double.NaN).isFinite());
        assertFalse(new TestCalculusFieldElement(Double.POSITIVE_INFINITY).isFinite());
    }

    @Test
    void testIsInfinite() {
        assertFalse(new TestCalculusFieldElement(1.).isInfinite());
        assertTrue(new TestCalculusFieldElement(Double.POSITIVE_INFINITY).isInfinite());
    }

    @Test
    void testIsNan() {
        assertFalse(new TestCalculusFieldElement(1.).isNaN());
        assertTrue(new TestCalculusFieldElement(Double.NaN).isNaN());
    }

    @Test
    void testNorm() {
        assertEquals(0., new TestCalculusFieldElement(0.).norm(), 0.0);
    }

    private static class TestCalculusFieldElement implements CalculusFieldElement<TestCalculusFieldElement> {

        private final double value;

        TestCalculusFieldElement (double value) {
            this.value = value;
        }

        @Override
        public TestCalculusFieldElement getAddendum() {
            return subtract(getReal());
        }

        @Override
        public TestCalculusFieldElement newInstance(double value) {
            return new TestCalculusFieldElement(value);
        }

        @Override
        public TestCalculusFieldElement scalb(int n) {
            return null;
        }

        @Override
        public TestCalculusFieldElement ulp() {
            return new TestCalculusFieldElement(FastMath.ulp(value));
        }

        @Override
        public TestCalculusFieldElement hypot(TestCalculusFieldElement y) throws MathIllegalArgumentException {
            return null;
        }

        @Override
        public TestCalculusFieldElement pow(TestCalculusFieldElement e) throws MathIllegalArgumentException {
            return new TestCalculusFieldElement(FastMath.pow(value, e.value));
        }

        @Override
        public TestCalculusFieldElement exp() {
            return null;
        }

        @Override
        public TestCalculusFieldElement expm1() {
            return null;
        }

        @Override
        public TestCalculusFieldElement log() {
            return null;
        }

        @Override
        public TestCalculusFieldElement log1p() {
            return null;
        }

        @Override
        public TestCalculusFieldElement log10() {
            return null;
        }

        @Override
        public TestCalculusFieldElement sin() {
            return new TestCalculusFieldElement(FastMath.sin(value));
        }

        @Override
        public TestCalculusFieldElement cos() {
            return new TestCalculusFieldElement(FastMath.cos(value));
        }

        @Override
        public TestCalculusFieldElement acos() {
            return null;
        }

        @Override
        public TestCalculusFieldElement asin() {
            return null;
        }

        @Override
        public TestCalculusFieldElement atan() {
            return null;
        }

        @Override
        public TestCalculusFieldElement atan2(TestCalculusFieldElement x) throws MathIllegalArgumentException {
            return null;
        }

        @Override
        public TestCalculusFieldElement sinh() {
            return new TestCalculusFieldElement(FastMath.sinh(value));
        }

        @Override
        public TestCalculusFieldElement cosh() {
            return new TestCalculusFieldElement(FastMath.cosh(value));
        }

        @Override
        public TestCalculusFieldElement acosh() {
            return null;
        }

        @Override
        public TestCalculusFieldElement asinh() {
            return null;
        }

        @Override
        public TestCalculusFieldElement atanh() {
            return null;
        }

        @Override
        public TestCalculusFieldElement linearCombination(TestCalculusFieldElement[] a, TestCalculusFieldElement[] b) throws MathIllegalArgumentException {
            return null;
        }

        @Override
        public TestCalculusFieldElement linearCombination(TestCalculusFieldElement a1, TestCalculusFieldElement b1, TestCalculusFieldElement a2, TestCalculusFieldElement b2) {
            return newInstance(a1.value * b1.value + a2.value * b2.value);
        }

        @Override
        public TestCalculusFieldElement linearCombination(TestCalculusFieldElement a1, TestCalculusFieldElement b1, TestCalculusFieldElement a2, TestCalculusFieldElement b2, TestCalculusFieldElement a3, TestCalculusFieldElement b3) {
            return newInstance(a1.value * b1.value + a2.value * b2.value + a3.value * b3.value);
        }

        @Override
        public TestCalculusFieldElement linearCombination(TestCalculusFieldElement a1, TestCalculusFieldElement b1, TestCalculusFieldElement a2, TestCalculusFieldElement b2, TestCalculusFieldElement a3, TestCalculusFieldElement b3, TestCalculusFieldElement a4, TestCalculusFieldElement b4) {
            return newInstance(a1.value * b1.value + a2.value * b2.value + a3.value * b3.value +
                    a4.value * b4.value);
        }

        @Override
        public TestCalculusFieldElement ceil() {
            return new TestCalculusFieldElement(FastMath.ceil(value));
        }

        @Override
        public TestCalculusFieldElement floor() {
            return new TestCalculusFieldElement(FastMath.floor(value));
        }

        @Override
        public TestCalculusFieldElement rint() {
            return new TestCalculusFieldElement(FastMath.rint(value));
        }

        @Override
        public TestCalculusFieldElement remainder(double a) {
            return null;
        }

        @Override
        public TestCalculusFieldElement remainder(TestCalculusFieldElement a) {
            return null;
        }

        @Override
        public TestCalculusFieldElement sign() {
            return new TestCalculusFieldElement(FastMath.signum(value));
        }

        @Override
        public TestCalculusFieldElement copySign(TestCalculusFieldElement sign) {
            return null;
        }

        @Override
        public TestCalculusFieldElement abs() {
            return new TestCalculusFieldElement(FastMath.abs(value));
        }

        @Override
        public double getReal() {
            return value;
        }

        @Override
        public TestCalculusFieldElement add(TestCalculusFieldElement a) throws NullArgumentException {
            return new TestCalculusFieldElement(value + a.value);
        }

        @Override
        public TestCalculusFieldElement negate() {
            return new TestCalculusFieldElement(-value);
        }

        @Override
        public TestCalculusFieldElement multiply(TestCalculusFieldElement a) throws NullArgumentException {
            return new TestCalculusFieldElement(value * a.value);
        }

        @Override
        public TestCalculusFieldElement reciprocal() throws MathRuntimeException {
            return new TestCalculusFieldElement(1. / value);
        }

        @Override
        public Field<TestCalculusFieldElement> getField() {
            return null;
        }

        @Override
        public boolean equals(Object obj) {
            if (obj instanceof TestCalculusFieldElement) {
                return Double.compare(value, ((TestCalculusFieldElement) obj).value) == 0;
            } else {
                return false;
            }
        }
    }

}
