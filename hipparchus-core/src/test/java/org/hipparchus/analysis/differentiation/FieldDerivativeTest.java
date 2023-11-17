package org.hipparchus.analysis.differentiation;

import org.hipparchus.Field;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.analysis.differentiation.DerivativeTest.TestDerivative;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class FieldDerivativeTest {

    private static final double TOLERANCE = 1e-10;

    @Test
    public void testGetReal() {
        // GIVEN
        final double expectedOperation = 0.5;
        final TestDerivative testDerivative = new TestDerivative(expectedOperation);
        final TestFieldDerivative testFieldDerivative = new TestFieldDerivative(testDerivative);
        // WHEN
        final double actualOperation = testFieldDerivative.getReal();
        // THEN
        Assert.assertEquals(expectedOperation, actualOperation, 0.);
    }

    @Test
    public void testPow() {
        // GIVEN
        final TestDerivative testDerivative = new TestDerivative(2.);
        final TestFieldDerivative testFieldDerivative = new TestFieldDerivative(testDerivative);
        final TestDerivative testDerivative2 = new TestDerivative(3.);
        final TestFieldDerivative testFieldDerivative2 = new TestFieldDerivative(testDerivative2);
        // WHEN
        final TestFieldDerivative actualOperation = testFieldDerivative.pow(testFieldDerivative2);
        // THEN
        final TestFieldDerivative expectedOperation = new TestFieldDerivative(new TestDerivative(FastMath
                .pow(testDerivative.getValue(), testDerivative2.getValue())));
        Assert.assertEquals(expectedOperation.getValue().getValue(), actualOperation.getValue().getValue(), TOLERANCE);
    }

    @Test
    public void testLog10() {
        // GIVEN
        final TestDerivative testDerivative = new TestDerivative(0.5);
        final TestFieldDerivative testFieldDerivative = new TestFieldDerivative(testDerivative);
        // WHEN
        final TestFieldDerivative actualOperation = testFieldDerivative.log10();
        // THEN
        final TestFieldDerivative expectedOperation = new TestFieldDerivative(new TestDerivative(FastMath
                .log10(testDerivative.getValue())));
        Assert.assertEquals(expectedOperation.getValue().getValue(), actualOperation.getValue().getValue(), TOLERANCE);
    }

    @Test
    public void testCosh() {
        // GIVEN
        final TestDerivative testDerivative = new TestDerivative(0.5);
        final TestFieldDerivative testFieldDerivative = new TestFieldDerivative(testDerivative);
        // WHEN
        final TestFieldDerivative actualOperation = testFieldDerivative.cosh();
        // THEN
        final TestFieldDerivative expectedOperation = new TestFieldDerivative(new TestDerivative(FastMath
                .cosh(testDerivative.getValue())));
        Assert.assertEquals(expectedOperation.getValue().getValue(), actualOperation.getValue().getValue(), TOLERANCE);
    }

    @Test
    public void testSinh() {
        // GIVEN
        final TestDerivative testDerivative = new TestDerivative(0.5);
        final TestFieldDerivative testFieldDerivative = new TestFieldDerivative(testDerivative);
        // WHEN
        final TestFieldDerivative actualOperation = testFieldDerivative.sinh();
        // THEN
        final TestFieldDerivative expectedOperation = new TestFieldDerivative(new TestDerivative(FastMath
                .sinh(testDerivative.getValue())));
        Assert.assertEquals(expectedOperation.getValue().getValue(), actualOperation.getValue().getValue(), TOLERANCE);
    }

    @Test
    public void testAcos() {
        // GIVEN
        final TestDerivative testDerivative = new TestDerivative(0.5);
        final TestFieldDerivative testFieldDerivative = new TestFieldDerivative(testDerivative);
        // WHEN
        final TestFieldDerivative actualOperation = testFieldDerivative.acos();
        // THEN
        final TestFieldDerivative expectedOperation = new TestFieldDerivative(new TestDerivative(FastMath
                .acos(testDerivative.getValue())));
        Assert.assertEquals(expectedOperation.getValue().getValue(), actualOperation.getValue().getValue(), TOLERANCE);
    }

    @Test
    public void testFloor() {
        // GIVEN
        final TestDerivative testDerivative = new TestDerivative(0.5);
        final TestFieldDerivative testFieldDerivative = new TestFieldDerivative(testDerivative);
        // WHEN
        final TestFieldDerivative actualOperation = testFieldDerivative.floor();
        // THEN
        final TestFieldDerivative expectedOperation = new TestFieldDerivative(new TestDerivative(FastMath
                .floor(testDerivative.getValue())));
        Assert.assertEquals(expectedOperation.getValue().getValue(), actualOperation.getValue().getValue(), TOLERANCE);
    }

    @Test
    public void testCeil() {
        // GIVEN
        final TestDerivative testDerivative = new TestDerivative(0.5);
        final TestFieldDerivative testFieldDerivative = new TestFieldDerivative(testDerivative);
        // WHEN
        final TestFieldDerivative actualOperation = testFieldDerivative.ceil();
        // THEN
        final TestFieldDerivative expectedOperation = new TestFieldDerivative(new TestDerivative(FastMath
                .ceil(testDerivative.getValue())));
        Assert.assertEquals(expectedOperation.getValue().getValue(), actualOperation.getValue().getValue(), TOLERANCE);
    }

    @Test
    public void testRint() {
        // GIVEN
        final TestDerivative testDerivative = new TestDerivative(0.5);
        final TestFieldDerivative testFieldDerivative = new TestFieldDerivative(testDerivative);
        // WHEN
        final TestFieldDerivative actualOperation = testFieldDerivative.rint();
        // THEN
        final TestFieldDerivative expectedOperation = new TestFieldDerivative(new TestDerivative(FastMath
                .rint(testDerivative.getValue())));
        Assert.assertEquals(expectedOperation.getValue().getValue(), actualOperation.getValue().getValue(), TOLERANCE);
    }

    @Test
    public void testSign() {
        // GIVEN
        final TestDerivative testDerivative = new TestDerivative(0.5);
        final TestFieldDerivative testFieldDerivative = new TestFieldDerivative(testDerivative);
        // WHEN
        final TestFieldDerivative actualOperation = testFieldDerivative.sign();
        // THEN
        final TestFieldDerivative expectedOperation = new TestFieldDerivative(new TestDerivative(FastMath
                .signum(testDerivative.getValue())));
        Assert.assertEquals(expectedOperation, actualOperation);
    }

    private static class TestFieldDerivative implements FieldDerivative<TestDerivative, TestFieldDerivative> {

        private final TestDerivative value;

        TestFieldDerivative (TestDerivative value) {
            this.value = value;
        }

        @Override
        public TestFieldDerivative newInstance(double value) {
            return newInstance(new TestDerivative(value));
        }

        @Override
        public TestFieldDerivative scalb(int n) {
            return null;
        }

        @Override
        public TestFieldDerivative hypot(TestFieldDerivative y) throws MathIllegalArgumentException {
            return null;
        }

        @Override
        public TestFieldDerivative exp() {
            return newInstance(value.exp());
        }

        @Override
        public TestFieldDerivative expm1() {
            return null;
        }

        @Override
        public TestFieldDerivative log() {
            return newInstance(value.log());
        }

        @Override
        public TestFieldDerivative log1p() {
            return null;
        }

        @Override
        public TestFieldDerivative cos() {
            return null;
        }

        @Override
        public TestFieldDerivative sin() {
            return null;
        }

        @Override
        public TestFieldDerivative asin() {
            return newInstance(value.asin());
        }

        @Override
        public TestFieldDerivative atan() {
            return null;
        }

        @Override
        public TestFieldDerivative atan2(TestFieldDerivative x) throws MathIllegalArgumentException {
            return null;
        }

        @Override
        public TestFieldDerivative acosh() {
            return null;
        }

        @Override
        public TestFieldDerivative asinh() {
            return null;
        }

        @Override
        public TestFieldDerivative atanh() {
            return null;
        }

        @Override
        public TestFieldDerivative linearCombination(TestFieldDerivative[] a, TestFieldDerivative[] b) throws MathIllegalArgumentException {
            return null;
        }

        @Override
        public TestFieldDerivative linearCombination(TestFieldDerivative a1, TestFieldDerivative b1, TestFieldDerivative a2, TestFieldDerivative b2) {
            return null;
        }

        @Override
        public TestFieldDerivative linearCombination(TestFieldDerivative a1, TestFieldDerivative b1, TestFieldDerivative a2, TestFieldDerivative b2, TestFieldDerivative a3, TestFieldDerivative b3) {
            return null;
        }

        @Override
        public TestFieldDerivative linearCombination(TestFieldDerivative a1, TestFieldDerivative b1, TestFieldDerivative a2, TestFieldDerivative b2, TestFieldDerivative a3, TestFieldDerivative b3, TestFieldDerivative a4, TestFieldDerivative b4) {
            return null;
        }

        @Override
        public TestFieldDerivative remainder(double a) {
            return null;
        }

        @Override
        public TestFieldDerivative remainder(TestFieldDerivative a) {
            return null;
        }

        @Override
        public TestFieldDerivative copySign(TestFieldDerivative sign) {
            return null;
        }

        @Override
        public TestFieldDerivative abs() {
            return null;
        }

        @Override
        public TestFieldDerivative add(TestFieldDerivative a) throws NullArgumentException {
            return new TestFieldDerivative(value.add(a.value));
        }

        @Override
        public TestFieldDerivative negate() {
            return new TestFieldDerivative(value.negate());
        }

        @Override
        public TestFieldDerivative multiply(TestFieldDerivative a) throws NullArgumentException {
            return new TestFieldDerivative(value.multiply(a.value));
        }

        @Override
        public TestFieldDerivative reciprocal() throws MathRuntimeException {
            return new TestFieldDerivative(value.reciprocal());
        }

        @Override
        public Field<TestFieldDerivative> getField() {
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
        public TestDerivative getValue() {
            return value;
        }

        @Override
        public TestDerivative getPartialDerivative(int... orders) throws MathIllegalArgumentException {
            return null;
        }

        @Override
        public TestFieldDerivative newInstance(TestDerivative value) {
            return new TestFieldDerivative(value);
        }

        @Override
        public boolean equals(Object obj) {
            if (obj instanceof TestFieldDerivative) {
                return value.equals(((TestFieldDerivative) obj).value);
            } else {
                return false;
            }
        }
    }

}