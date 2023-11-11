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
package org.hipparchus.analysis.differentiation;

import org.hipparchus.Field;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class DerivativeTest {

    private static final double TOLERANCE = 1e-10;

    @Test
    public void testGetReal() {
        // GIVEN
        final double expectedOperation = 0.5;
        final TestDerivative testDerivative = new TestDerivative(expectedOperation);
        // WHEN
        final double actualOperation = testDerivative.getReal();
        // THEN
        Assert.assertEquals(expectedOperation, actualOperation, 0.);
    }

    @Test
    public void testLog10() {
        // GIVEN
        final double value = 0.5;
        final TestDerivative testDerivative = new TestDerivative(value);
        // WHEN
        final TestDerivative actualOperation = testDerivative.log10();
        // THEN
        final TestDerivative expectedOperation = new TestDerivative(FastMath.log10(value));
        Assert.assertEquals(expectedOperation.getValue(), actualOperation.getValue(), TOLERANCE);
    }

    @Test
    public void testPow() {
        // GIVEN
        final double value1 = 2.;
        final double value2 = 3.;
        final TestDerivative testDerivative1 = new TestDerivative(value1);
        final TestDerivative testDerivative2 = new TestDerivative(value2);
        // WHEN
        final TestDerivative actualOperation = testDerivative1.pow(testDerivative2);
        // THEN
        final TestDerivative expectedOperation = new TestDerivative(FastMath.pow(value1, value2));
        Assert.assertEquals(expectedOperation.getValue(), actualOperation.getValue(), TOLERANCE);
    }

    @Test
    public void testCosh() {
        // GIVEN
        final double value = 0.5;
        final TestDerivative testDerivative = new TestDerivative(value);
        // WHEN
        final TestDerivative actualOperation = testDerivative.cosh();
        // THEN
        final TestDerivative expectedOperation = new TestDerivative(FastMath.cosh(value));
        Assert.assertEquals(expectedOperation, actualOperation);
    }

    @Test
    public void testSinh() {
        // GIVEN
        final double value = 0.5;
        final TestDerivative testDerivative = new TestDerivative(value);
        // WHEN
        final TestDerivative actualOperation = testDerivative.sinh();
        // THEN
        final TestDerivative expectedOperation = new TestDerivative(FastMath.sinh(value));
        Assert.assertEquals(expectedOperation, actualOperation);
    }

    @Test
    public void testAcos() {
        // GIVEN
        final double value = 0.5;
        final TestDerivative testDerivative = new TestDerivative(value);
        // WHEN
        final TestDerivative actualOperation = testDerivative.acos();
        // THEN
        final TestDerivative expectedOperation = new TestDerivative(FastMath.acos(value));
        Assert.assertEquals(expectedOperation.getValue(), actualOperation.getValue(), TOLERANCE);
    }

    @Test
    public void testFloor() {
        // GIVEN
        final double value = 0.5;
        final TestDerivative testDerivative = new TestDerivative(value);
        // WHEN
        final TestDerivative actualOperation = testDerivative.floor();
        // THEN
        final TestDerivative expectedOperation = new TestDerivative(FastMath.floor(value));
        Assert.assertEquals(expectedOperation, actualOperation);
    }

    @Test
    public void testCeil() {
        // GIVEN
        final double value = 0.5;
        final TestDerivative testDerivative = new TestDerivative(value);
        // WHEN
        final TestDerivative actualOperation = testDerivative.ceil();
        // THEN
        final TestDerivative expectedOperation = new TestDerivative(FastMath.ceil(value));
        Assert.assertEquals(expectedOperation.getValue(), actualOperation.getValue(), TOLERANCE);
    }

    @Test
    public void testRint() {
        // GIVEN
        final double value = 0.5;
        final TestDerivative testDerivative = new TestDerivative(value);
        // WHEN
        final TestDerivative actualOperation = testDerivative.rint();
        // THEN
        final TestDerivative expectedOperation = new TestDerivative(FastMath.rint(value));
        Assert.assertEquals(expectedOperation, actualOperation);
    }

    @Test
    public void testSign() {
        // GIVEN
        final double value = 0.5;
        final TestDerivative testDerivative = new TestDerivative(value);
        // WHEN
        final TestDerivative actualOperation = testDerivative.sign();
        // THEN
        final TestDerivative expectedOperation = new TestDerivative(FastMath.signum(value));
        Assert.assertEquals(expectedOperation.getValue(), actualOperation.getValue(), TOLERANCE);
    }

    static class TestDerivative implements Derivative<TestDerivative> {

        private final double value;

        TestDerivative (final double value) {
            this.value = value;
        }

        @Override
        public TestDerivative newInstance(double value) {
            return new TestDerivative(value);
        }

        @Override
        public TestDerivative scalb(int n) {
            return null;
        }

        @Override
        public TestDerivative hypot(TestDerivative y) throws MathIllegalArgumentException {
            return null;
        }

        @Override
        public TestDerivative exp() {
            return new TestDerivative(FastMath.exp(value));
        }

        @Override
        public TestDerivative expm1() {
            return null;
        }

        @Override
        public TestDerivative log() {
            return new TestDerivative(FastMath.log(value));
        }

        @Override
        public TestDerivative log1p() {
            return null;
        }

        @Override
        public TestDerivative cos() {
            return null;
        }

        @Override
        public TestDerivative sin() {
            return null;
        }

        @Override
        public TestDerivative asin() {
            return new TestDerivative(FastMath.asin(value));
        }

        @Override
        public TestDerivative atan() {
            return null;
        }

        @Override
        public TestDerivative atan2(TestDerivative x) throws MathIllegalArgumentException {
            return null;
        }

        @Override
        public TestDerivative acosh() {
            return null;
        }

        @Override
        public TestDerivative asinh() {
            return null;
        }

        @Override
        public TestDerivative atanh() {
            return null;
        }

        @Override
        public TestDerivative linearCombination(TestDerivative[] a, TestDerivative[] b) throws MathIllegalArgumentException {
            return null;
        }

        @Override
        public TestDerivative linearCombination(TestDerivative a1, TestDerivative b1, TestDerivative a2, TestDerivative b2) {
            return null;
        }

        @Override
        public TestDerivative linearCombination(TestDerivative a1, TestDerivative b1, TestDerivative a2, TestDerivative b2,
                                                TestDerivative a3, TestDerivative b3) {
            return null;
        }

        @Override
        public TestDerivative linearCombination(TestDerivative a1, TestDerivative b1, TestDerivative a2, TestDerivative b2,
                                                TestDerivative a3, TestDerivative b3, TestDerivative a4, TestDerivative b4) {
            return null;
        }

        @Override
        public TestDerivative remainder(double a) {
            return null;
        }

        @Override
        public TestDerivative remainder(TestDerivative a) {
            return null;
        }

        @Override
        public TestDerivative copySign(TestDerivative sign) {
            return null;
        }

        @Override
        public TestDerivative abs() {
            return null;
        }

        @Override
        public TestDerivative add(TestDerivative a) throws NullArgumentException {
            return new TestDerivative(value + a.value);
        }

        @Override
        public TestDerivative negate() {
            return new TestDerivative(-value);
        }

        @Override
        public TestDerivative multiply(TestDerivative a) throws NullArgumentException {
            return new TestDerivative(value * a.value);
        }

        @Override
        public TestDerivative reciprocal() throws MathRuntimeException {
            return new TestDerivative(1. / value);
        }

        @Override
        public Field<TestDerivative> getField() {
            return null;
        }

        @Override
        public int getFreeParameters() {
            return 0;
        }

        @Override
        public int getOrder() {
            return 0;
        }

        @Override
        public double getValue() {
            return value;
        }

        @Override
        public double getPartialDerivative(int... orders) throws MathIllegalArgumentException {
            return 0;
        }

        @Override
        public TestDerivative compose(double... f) throws MathIllegalArgumentException {
            return null;
        }

        @Override
        public boolean equals(Object obj) {
            if (obj instanceof DerivativeTest.TestDerivative) {
                return Double.compare(value, ((DerivativeTest.TestDerivative) obj).value) == 0;
            } else {
                return false;
            }
        }
    }

}