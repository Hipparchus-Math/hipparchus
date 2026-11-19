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

package org.hipparchus.complex;

import org.hipparchus.CalculusFieldElementAbstractTest;
import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Precision;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;


/**
 */
class ComplexTest extends CalculusFieldElementAbstractTest<Complex> {


    private final double  inf = Double.POSITIVE_INFINITY;
    private final double  neginf = Double.NEGATIVE_INFINITY;
    private final double  nan = Double.NaN;
    private final double  pi = FastMath.PI;
    private final Complex oneInf = new Complex(1, inf);
    private final Complex oneNegInf = new Complex(1, neginf);
    private final Complex infOne = new Complex(inf, 1);
    private final Complex infZero = new Complex(inf, 0);
    private final Complex infNaN = new Complex(inf, nan);
    private final Complex infNegInf = new Complex(inf, neginf);
    private final Complex infInf = new Complex(inf, inf);
    private final Complex negInfInf = new Complex(neginf, inf);
    private final Complex negInfZero = new Complex(neginf, 0);
    private final Complex negInfOne = new Complex(neginf, 1);
    private final Complex negInfNaN = new Complex(neginf, nan);
    private final Complex negInfNegInf = new Complex(neginf, neginf);
    private final Complex oneNaN = new Complex(1, nan);
    private final Complex zeroInf = new Complex(0, inf);
    private final Complex zeroNaN = new Complex(0, nan);
    private final Complex nanInf = new Complex(nan, inf);
    private final Complex nanNegInf = new Complex(nan, neginf);
    private final Complex nanZero = new Complex(nan, 0);

    @Override
    protected Complex build(final double x) {
        return new Complex(x, 0.0);
    }

    @Test
    void testConstructor() {
        Complex z = new Complex(3.0, 4.0);
        assertEquals(3.0, z.getReal(), 1.0e-5);
        assertEquals(4.0, z.getImaginary(), 1.0e-5);
        assertEquals(3.0, z.getRealPart(), 1.0e-5);
        assertEquals(4.0, z.getImaginaryPart(), 1.0e-5);
    }

    @Test
    void testConstructorNaN() {
        Complex z = new Complex(3.0, Double.NaN);
        assertTrue(z.isNaN());

        z = new Complex(nan, 4.0);
        assertTrue(z.isNaN());

        z = new Complex(3.0, 4.0);
        assertFalse(z.isNaN());
    }

    @Test
    void testNorm() {
        Complex z = new Complex(3.0, 4.0);
        assertEquals(5.0, z.norm(), 1.0e-5);
    }

    @Test
    void testNormNaN() {
        assertTrue(Double.isNaN(Complex.NaN.norm()));
        Complex z = new Complex(inf, nan);
        assertTrue(Double.isNaN(z.norm()));
    }

    @Test
    void testNormInfinite() {
        Complex z = Complex.NaN.newInstance(inf);
        assertEquals(inf, z.norm(), 0);
        z = new Complex(0, neginf);
        assertEquals(inf, z.norm(), 0);
        z = new Complex(inf, neginf);
        assertEquals(inf, z.norm(), 0);
    }

    @Test
    void testAdd() {
        Complex x = new Complex(3.0, 4.0);
        Complex y = new Complex(5.0, 6.0);
        Complex z = x.add(y);
        assertEquals(8.0, z.getReal(), 1.0e-5);
        assertEquals(10.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    void testAddNaN() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.add(Complex.NaN);
        assertSame(Complex.NaN, z);
        z = new Complex(1, nan);
        Complex w = x.add(z);
        assertSame(Complex.NaN, w);
        assertSame(Complex.NaN, Complex.NaN.add(Double.NaN));
    }

    @Test
    void testAddInf() {
        Complex x = new Complex(1, 1);
        Complex z = new Complex(inf, 0);
        Complex w = x.add(z);
        assertEquals(1, w.getImaginary(), 0);
        assertEquals(inf, w.getReal(), 0);

        x = new Complex(neginf, 0);
        assertTrue(Double.isNaN(x.add(z).getReal()));
    }

    @Test
    void testScalarAdd() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = 2.0;
        Complex yComplex = new Complex(yDouble);
        assertEquals(x.add(yComplex), x.add(yDouble));
    }

    @Test
    void testScalarAddNaN() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = Double.NaN;
        Complex yComplex = new Complex(yDouble);
        assertEquals(x.add(yComplex), x.add(yDouble));
    }

    @Test
    void testScalarAddInf() {
        Complex x = new Complex(1, 1);
        double yDouble = Double.POSITIVE_INFINITY;

        Complex yComplex = new Complex(yDouble);
        assertEquals(x.add(yComplex), x.add(yDouble));

        x = new Complex(neginf, 0);
        assertEquals(x.add(yComplex), x.add(yDouble));
    }

    @Test
    void testConjugate() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.conjugate();
        assertEquals(3.0, z.getReal(), 1.0e-5);
        assertEquals(-4.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    void testConjugateNaN() {
        Complex z = Complex.NaN.conjugate();
        assertTrue(z.isNaN());
    }

    @Test
    void testConjugateInfiinite() {
        Complex z = new Complex(0, inf);
        assertEquals(neginf, z.conjugate().getImaginary(), 0);
        z = new Complex(0, neginf);
        assertEquals(inf, z.conjugate().getImaginary(), 0);
    }

    @Test
    void testDivide() {
        Complex x = new Complex(3.0, 4.0);
        Complex y = new Complex(5.0, 6.0);
        Complex z = x.divide(y);
        assertEquals(39.0 / 61.0, z.getReal(), 1.0e-5);
        assertEquals(2.0 / 61.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    void testDivideReal() {
        Complex x = new Complex(2d, 3d);
        Complex y = new Complex(2d, 0d);
        assertEquals(new Complex(1d, 1.5), x.divide(y));

    }

    @Test
    void testDivideImaginary() {
        Complex x = new Complex(2d, 3d);
        Complex y = new Complex(0d, 2d);
        assertEquals(new Complex(1.5d, -1d), x.divide(y));
    }

    @Test
    void testDivideInf() {
        Complex x = new Complex(3, 4);
        Complex w = new Complex(neginf, inf);
        assertEquals(Complex.ZERO, x.divide(w));

        Complex z = w.divide(x);
        assertTrue(Double.isNaN(z.getReal()));
        assertEquals(inf, z.getImaginary(), 0);

        w = new Complex(inf, inf);
        z = w.divide(x);
        assertTrue(Double.isNaN(z.getImaginary()));
        assertEquals(inf, z.getReal(), 0);

        w = new Complex(1, inf);
        z = w.divide(w);
        assertTrue(Double.isNaN(z.getReal()));
        assertTrue(Double.isNaN(z.getImaginary()));
    }

    @Test
    void testDivideZero() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.divide(Complex.ZERO);
        // Assertions.assertEquals(z, Complex.INF); // See MATH-657
        assertEquals(Complex.NaN, z);
    }

    @Test
    void testDivideZeroZero() {
        Complex x = new Complex(0.0, 0.0);
        Complex z = x.divide(Complex.ZERO);
        assertEquals(Complex.NaN, z);
    }

    @Test
    void testDivideNaN() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.divide(Complex.NaN);
        assertTrue(z.isNaN());
        assertTrue(Complex.NaN.divide(Complex.NaN).isNaN());
    }

    @Test
    void testDivideNaNInf() {
       Complex z = oneInf.divide(Complex.ONE);
       assertTrue(Double.isNaN(z.getReal()));
       assertEquals(inf, z.getImaginary(), 0);

       z = negInfNegInf.divide(oneNaN);
       assertTrue(Double.isNaN(z.getReal()));
       assertTrue(Double.isNaN(z.getImaginary()));

       z = negInfInf.divide(Complex.ONE);
       assertTrue(Double.isNaN(z.getReal()));
       assertTrue(Double.isNaN(z.getImaginary()));
    }

    @Test
    void testScalarDivide() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = 2.0;
        Complex yComplex = new Complex(yDouble);
        assertEquals(x.divide(yComplex), x.divide(yDouble));
    }

    @Test
    void testScalarDivideNaN() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = Double.NaN;
        Complex yComplex = new Complex(yDouble);
        assertEquals(x.divide(yComplex), x.divide(yDouble));
        assertTrue(Complex.NaN.divide(Double.NaN).isNaN());
    }

    @Test
    void testScalarDivideInf() {
        Complex x = new Complex(1,1);
        double yDouble = Double.POSITIVE_INFINITY;
        Complex yComplex = new Complex(yDouble);
        UnitTestUtils.customAssertEquals(x.divide(yComplex), x.divide(yDouble), 0);

        yDouble = Double.NEGATIVE_INFINITY;
        yComplex = new Complex(yDouble);
        UnitTestUtils.customAssertEquals(x.divide(yComplex), x.divide(yDouble), 0);

        x = new Complex(1, Double.NEGATIVE_INFINITY);
        UnitTestUtils.customAssertEquals(x.divide(yComplex), x.divide(yDouble), 0);
    }

    @Test
    void testScalarDivideZero() {
        Complex x = new Complex(1,1);
        UnitTestUtils.customAssertEquals(x.divide(Complex.ZERO), x.divide(0), 0);
    }

    @Test
    void testReciprocal() {
        Complex z = new Complex(5.0, 6.0);
        Complex act = z.reciprocal();
        double expRe = 5.0 / 61.0;
        double expIm = -6.0 / 61.0;
        assertEquals(expRe, act.getReal(), FastMath.ulp(expRe));
        assertEquals(expIm, act.getImaginary(), FastMath.ulp(expIm));
    }

    @Test
    void testReciprocalReal() {
        Complex z = new Complex(-2.0, 0.0);
        assertTrue(Complex.equals(new Complex(-0.5, 0.0), z.reciprocal()));
    }

    @Test
    void testReciprocalImaginary() {
        Complex z = new Complex(0.0, -2.0);
        assertEquals(new Complex(0.0, 0.5), z.reciprocal());
    }

    @Test
    void testReciprocalInf() {
        Complex z = new Complex(neginf, inf);
        assertEquals(Complex.ZERO, z.reciprocal());

        z = new Complex(1, inf).reciprocal();
        assertEquals(Complex.ZERO, z);
    }

    @Test
    void testReciprocalZero() {
        assertEquals(Complex.INF, Complex.ZERO.reciprocal());
    }

    @Test
    void testReciprocalNaN() {
        assertTrue(Complex.NaN.reciprocal().isNaN());
    }

    @Test
    void testMultiply() {
        Complex x = new Complex(3.0, 4.0);
        Complex y = new Complex(5.0, 6.0);
        Complex z = x.multiply(y);
        assertEquals(-9.0, z.getReal(), 1.0e-5);
        assertEquals(38.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    void testMultiplyNaN() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.multiply(Complex.NaN);
        assertSame(Complex.NaN, z);
        z = Complex.NaN.multiply(5);
        assertSame(Complex.NaN, z);
    }

    @Test
    void testMultiplyInfInf() {
        // Assertions.assertTrue(infInf.multiply(infInf).isNaN()); // MATH-620
        assertTrue(infInf.multiply(infInf).isInfinite());
    }

    @Test
    void testMultiplyNaNInf() {
        Complex z = new Complex(1,1);
        Complex w = z.multiply(infOne);
        assertEquals(w.getReal(), inf, 0);
        assertEquals(w.getImaginary(), inf, 0);

        // [MATH-164]
        assertEquals(Complex.INF, new Complex( 1, 0).multiply(infInf));
        assertEquals(Complex.INF, new Complex(-1, 0).multiply(infInf));
        assertEquals(Complex.INF, new Complex( 1, 0).multiply(negInfZero));

        w = oneInf.multiply(oneNegInf);
        assertEquals(w.getReal(), inf, 0);
        assertEquals(w.getImaginary(), inf, 0);

        w = negInfNegInf.multiply(oneNaN);
        assertTrue(Double.isNaN(w.getReal()));
        assertTrue(Double.isNaN(w.getImaginary()));

        z = new Complex(1, neginf);
        assertSame(Complex.INF, z.square());
    }

    @Test
    void testScalarMultiply() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = 2.0;
        Complex yComplex = new Complex(yDouble);
        assertEquals(x.multiply(yComplex), x.multiply(yDouble));
        int zInt = -5;
        Complex zComplex = new Complex(zInt);
        assertEquals(x.multiply(zComplex), x.multiply(zInt));
    }

    @Test
    void testScalarMultiplyNaN() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = Double.NaN;
        Complex yComplex = new Complex(yDouble);
        assertEquals(x.multiply(yComplex), x.multiply(yDouble));
        assertTrue(new Complex(Double.NaN, 0).multiply(5).isNaN());
        assertTrue(new Complex(0, Double.NaN).multiply(5).isNaN());
        assertTrue(Complex.NaN.multiply(5).isNaN());
        assertTrue(new Complex(Double.NaN, 0).multiply(5.0).isNaN());
        assertTrue(new Complex(0, Double.NaN).multiply(5.0).isNaN());
        assertTrue(Complex.NaN.multiply(5.0).isNaN());
        assertTrue(Complex.ONE.multiply(Double.NaN).isNaN());
    }

    @Test
    void testScalarMultiplyInf() {
        Complex x = new Complex(1, 1);
        double yDouble = Double.POSITIVE_INFINITY;
        Complex yComplex = new Complex(yDouble);
        assertEquals(x.multiply(yComplex), x.multiply(yDouble));

        yDouble = Double.NEGATIVE_INFINITY;
        yComplex = new Complex(yDouble);
        assertEquals(x.multiply(yComplex), x.multiply(yDouble));

        assertTrue(new Complex(Double.POSITIVE_INFINITY, 0).multiply(5).isInfinite());
        assertTrue(new Complex(0, Double.POSITIVE_INFINITY).multiply(5).isInfinite());
        assertTrue(Complex.INF.multiply(5).isInfinite());
        assertTrue(new Complex(Double.POSITIVE_INFINITY, 0).multiply(5.0).isInfinite());
        assertTrue(new Complex(0, Double.POSITIVE_INFINITY).multiply(5.0).isInfinite());
        assertTrue(Complex.INF.multiply(5.0).isInfinite());
        assertTrue(Complex.ONE.multiply(Double.POSITIVE_INFINITY).isInfinite());
    }

    @Test
    void testNegate() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.negate();
        assertEquals(-3.0, z.getReal(), 1.0e-5);
        assertEquals(-4.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    void testNegateNaN() {
        Complex z = Complex.NaN.negate();
        assertTrue(z.isNaN());
    }

    @Test
    void testSubtract() {
        Complex x = new Complex(3.0, 4.0);
        Complex y = new Complex(5.0, 6.0);
        Complex z = x.subtract(y);
        assertEquals(-2.0, z.getReal(), 1.0e-5);
        assertEquals(-2.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    void testSubtractNaN() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.subtract(Complex.NaN);
        assertSame(Complex.NaN, z);
        z = new Complex(1, nan);
        Complex w = x.subtract(z);
        assertSame(Complex.NaN, w);
        assertSame(Complex.NaN, Complex.NaN.subtract(Complex.NaN));
        assertSame(Complex.NaN, Complex.NaN.subtract(Double.NaN));
    }

    @Test
    void testSubtractInf() {
        Complex x = new Complex(1, 1);
        Complex z = new Complex(neginf, 0);
        Complex w = x.subtract(z);
        assertEquals(1, w.getImaginary(), 0);
        assertEquals(inf, w.getReal(), 0);

        x = new Complex(neginf, 0);
        assertTrue(Double.isNaN(x.subtract(z).getReal()));
    }

    @Test
    void testScalarSubtract() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = 2.0;
        Complex yComplex = new Complex(yDouble);
        assertEquals(x.subtract(yComplex), x.subtract(yDouble));
    }

    @Test
    void testScalarSubtractNaN() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = Double.NaN;
        Complex yComplex = new Complex(yDouble);
        assertEquals(x.subtract(yComplex), x.subtract(yDouble));
    }

    @Test
    void testScalarSubtractInf() {
        Complex x = new Complex(1, 1);
        double yDouble = Double.POSITIVE_INFINITY;
        Complex yComplex = new Complex(yDouble);
        assertEquals(x.subtract(yComplex), x.subtract(yDouble));

        x = new Complex(neginf, 0);
        assertEquals(x.subtract(yComplex), x.subtract(yDouble));
    }

    @Test
    void testEqualsNull() {
        Complex x = new Complex(3.0, 4.0);
        assertNotEquals(null, x);
    }

    @Test
    void testFloatingPointEqualsPrecondition1() {
        assertThrows(NullPointerException.class, () -> Complex.equals(new Complex(3.0, 4.0), null, 3));
    }

    @Test
    void testFloatingPointEqualsPrecondition2() {
        assertThrows(NullPointerException.class, () -> Complex.equals(null, new Complex(3.0, 4.0), 3));
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    void testEqualsClass() {
        Complex x = new Complex(3.0, 4.0);
        assertNotEquals(this, x);
    }

    @Test
    void testEqualsSame() {
        Complex x = new Complex(3.0, 4.0);
        assertEquals(x, x);
    }

    @Test
    void testFloatingPointEquals() {
        double re = -3.21;
        double im = 456789e10;

        final Complex x = new Complex(re, im);
        Complex y = new Complex(re, im);

        assertEquals(x, y);
        assertTrue(Complex.equals(x, y));

        final int maxUlps = 5;
        for (int i = 0; i < maxUlps; i++) {
            re = FastMath.nextUp(re);
            im = FastMath.nextUp(im);
        }
        y = new Complex(re, im);
        assertTrue(Complex.equals(x, y, maxUlps));

        re = FastMath.nextUp(re);
        im = FastMath.nextUp(im);
        y = new Complex(re, im);
        assertFalse(Complex.equals(x, y, maxUlps));
    }

    @Test
    void testFloatingPointEqualsNaN() {
        Complex c = new Complex(Double.NaN, 1);
        assertFalse(Complex.equals(c, c));

        c = new Complex(1, Double.NaN);
        assertFalse(Complex.equals(c, c));
    }

    @Test
    void testFloatingPointEqualsWithAllowedDelta() {
        final double re = 153.0000;
        final double im = 152.9375;
        final double tol1 = 0.0625;
        final Complex x = new Complex(re, im);
        final Complex y = new Complex(re + tol1, im + tol1);
        assertTrue(Complex.equals(x, y, tol1));

        final double tol2 = 0.0624;
        assertFalse(Complex.equals(x, y, tol2));
    }

    @Test
    void testFloatingPointEqualsWithAllowedDeltaNaN() {
        final Complex x = new Complex(0, Double.NaN);
        final Complex y = new Complex(Double.NaN, 0);
        assertFalse(Complex.equals(x, Complex.ZERO, 0.1));
        assertFalse(Complex.equals(x, x, 0.1));
        assertFalse(Complex.equals(x, y, 0.1));
    }

    @Test
    void testFloatingPointEqualsWithRelativeTolerance() {
        final double tol = 1e-4;
        final double re = 1;
        final double im = 1e10;

        final double f = 1 + tol;
        final Complex x = new Complex(re, im);
        final Complex y = new Complex(re * f, im * f);
        assertTrue(Complex.equalsWithRelativeTolerance(x, y, tol));
    }

    @Test
    void testFloatingPointEqualsWithRelativeToleranceNaN() {
        final Complex x = new Complex(0, Double.NaN);
        final Complex y = new Complex(Double.NaN, 0);
        assertFalse(Complex.equalsWithRelativeTolerance(x, Complex.ZERO, 0.1));
        assertFalse(Complex.equalsWithRelativeTolerance(x, x, 0.1));
        assertFalse(Complex.equalsWithRelativeTolerance(x, y, 0.1));
    }

    @Test
    void testEqualsTrue() {
        Complex x = new Complex(3.0, 4.0);
        Complex y = new Complex(3.0, 4.0);
        assertEquals(x, y);
    }

    @Test
    void testEqualsRealDifference() {
        Complex x = new Complex(0.0, 0.0);
        Complex y = new Complex(0.0 + Double.MIN_VALUE, 0.0);
        assertNotEquals(x, y);
    }

    @Test
    void testEqualsImaginaryDifference() {
        Complex x = new Complex(0.0, 0.0);
        Complex y = new Complex(0.0, 0.0 + Double.MIN_VALUE);
        assertNotEquals(x, y);
    }

    @Test
    void testEqualsNaN() {
        Complex realNaN = new Complex(Double.NaN, 0.0);
        Complex imaginaryNaN = new Complex(0.0, Double.NaN);
        Complex complexNaN = Complex.NaN;
        assertEquals(realNaN, imaginaryNaN);
        assertEquals(imaginaryNaN, complexNaN);
        assertEquals(realNaN, complexNaN);
    }

    @Test
    void testHashCode() {
        Complex x = new Complex(0.0, 0.0);
        Complex y = new Complex(0.0, 0.0 + Double.MIN_VALUE);
        assertNotEquals(x.hashCode(), y.hashCode());
        y = new Complex(0.0 + Double.MIN_VALUE, 0.0);
        assertNotEquals(x.hashCode(), y.hashCode());
        Complex realNaN = new Complex(Double.NaN, 0.0);
        Complex imaginaryNaN = new Complex(0.0, Double.NaN);
        assertEquals(realNaN.hashCode(), imaginaryNaN.hashCode());
        assertEquals(imaginaryNaN.hashCode(), Complex.NaN.hashCode());

        // MATH-1118
        // "equals" and "hashCode" must be compatible: if two objects have
        // different hash codes, "equals" must return false.
        final String msg = "'equals' not compatible with 'hashCode'";

        x = new Complex(0.0, 0.0);
        y = new Complex(0.0, -0.0);
        assertTrue(x.hashCode() != y.hashCode());
        assertNotEquals(x, y, msg);

        x = new Complex(0.0, 0.0);
        y = new Complex(-0.0, 0.0);
        assertTrue(x.hashCode() != y.hashCode());
        assertNotEquals(x, y, msg);
    }

    @Test
    void testToDegreesComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(FastMath.toDegrees(z.getReal()), FastMath.toDegrees(z.getImaginary()));
        UnitTestUtils.customAssertEquals(expected, z.toDegrees(), 1.0e-15);
    }

    @Test
    void testToRadiansComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(FastMath.toRadians(z.getReal()), FastMath.toRadians(z.getImaginary()));
        UnitTestUtils.customAssertEquals(expected, z.toRadians(), 1.0e-15);
    }

    @Test
    void testAcosComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(0.936812, -2.30551);
        UnitTestUtils.customAssertEquals(expected, z.acos(), 1.0e-5);
        UnitTestUtils.customAssertEquals(new Complex(FastMath.acos(0), 0),
                                         Complex.ZERO.acos(), 1.0e-12);
    }

    @Test
    void testAcosNaN() {
        assertTrue(Complex.NaN.acos().isNaN());
    }

    @Test
    void testAcosInf() {
        UnitTestUtils.customAssertSame(Complex.NaN, oneInf.acos());
        UnitTestUtils.customAssertSame(Complex.NaN, oneNegInf.acos());
        UnitTestUtils.customAssertSame(Complex.NaN, infOne.acos());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfOne.acos());
        UnitTestUtils.customAssertSame(Complex.NaN, infInf.acos());
        UnitTestUtils.customAssertSame(Complex.NaN, infNegInf.acos());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfInf.acos());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfNegInf.acos());
    }

    @Test
    void testAcosBranchCuts() {
        UnitTestUtils.customAssertEquals(new Complex(3.141592653589793238462, -0.76103968373182660633),
                                         FastMath.acos(new Complex(-1.3038404810405297, +0.0)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(new Complex(3.141592653589793238462, +0.76103968373182660633),
                                         FastMath.acos(new Complex(-1.3038404810405297, -0.0)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(new Complex(0.0, -0.76103968373182660633),
                                         FastMath.acos(new Complex(1.3038404810405297, +0.0)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(new Complex(0.0, +0.76103968373182660633),
                                         FastMath.acos(new Complex(1.3038404810405297, -0.0)),
                                         1.0e-14);
    }

    @Test
    void testAsinComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(0.633984, 2.30551);
        UnitTestUtils.customAssertEquals(expected, z.asin(), 1.0e-5);
    }

    @Test
    void testAsinNaN() {
        assertTrue(Complex.NaN.asin().isNaN());
    }

    @Test
    void testAsinInf() {
        UnitTestUtils.customAssertSame(Complex.NaN, oneInf.asin());
        UnitTestUtils.customAssertSame(Complex.NaN, oneNegInf.asin());
        UnitTestUtils.customAssertSame(Complex.NaN, infOne.asin());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfOne.asin());
        UnitTestUtils.customAssertSame(Complex.NaN, infInf.asin());
        UnitTestUtils.customAssertSame(Complex.NaN, infNegInf.asin());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfInf.asin());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfNegInf.asin());
    }

    @Test
    void testAsinBranchCuts() {
        UnitTestUtils.customAssertEquals(new Complex(-1.57079632679489661923, +0.76103968373182660633),
                                         FastMath.asin(new Complex(-1.3038404810405297, +0.0)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(new Complex(-1.57079632679489661923, -0.76103968373182660633),
                                         FastMath.asin(new Complex(-1.3038404810405297, -0.0)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(new Complex(1.57079632679489661923, +0.76103968373182660633),
                                         FastMath.asin(new Complex(1.3038404810405297, +0.0)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(new Complex(1.57079632679489661923, -0.76103968373182660633),
                                         FastMath.asin(new Complex(1.3038404810405297, -0.0)),
                                         1.0e-14);
    }

    @Test
    void testAtanComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(1.44831, 0.158997);
        UnitTestUtils.customAssertEquals(expected, z.atan(), 1.0e-5);
    }

    @Test
    void testAtanInf() {
        UnitTestUtils.customAssertSame(Complex.NaN, oneInf.atan());
        UnitTestUtils.customAssertSame(Complex.NaN, oneNegInf.atan());
        UnitTestUtils.customAssertSame(Complex.NaN, infOne.atan());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfOne.atan());
        UnitTestUtils.customAssertSame(Complex.NaN, infInf.atan());
        UnitTestUtils.customAssertSame(Complex.NaN, infNegInf.atan());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfInf.atan());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfNegInf.atan());
    }

    @Test
    void testAtanI() {
        assertTrue(Complex.I.atan().isNaN());
    }

    @Test
    void testAtanNaN() {
        assertTrue(Complex.NaN.atan().isNaN());
    }

    @Test
    void testAtanBranchCuts() {
        UnitTestUtils.customAssertEquals(new Complex(+1.5707963267948966192, +1.0986122886681096913),
                                         FastMath.atan(new Complex(+0.0, 1.25)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(new Complex(-1.5707963267948966192, +1.0986122886681096913),
                                         FastMath.atan(new Complex(-0.0, 1.25)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(new Complex(+1.5707963267948966192, -1.0986122886681096913),
                                         FastMath.atan(new Complex(+0.0, -1.25)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(new Complex(-1.5707963267948966192, -1.0986122886681096913),
                                         FastMath.atan(new Complex(-0.0, -1.25)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(new Complex(0.0, +0.25541281188299534160),
                                         FastMath.atan(new Complex(+0.0, 0.25)),
                                         1.0e-14);
        assertTrue(FastMath.copySign(1.0, FastMath.atan(new Complex(+0.0, 0.25)).getReal()) > 0.0);
        UnitTestUtils.customAssertEquals(new Complex(0.0, +0.25541281188299534160),
                                         FastMath.atan(new Complex(-0.0, 0.25)),
                                         1.0e-14);
        assertTrue(FastMath.copySign(1.0, FastMath.atan(new Complex(-0.0, 0.25)).getReal()) < 0.0);
        UnitTestUtils.customAssertEquals(new Complex(0.0, -0.25541281188299534160),
                                         FastMath.atan(new Complex(+0.0, -0.25)),
                                         1.0e-14);
        assertTrue(FastMath.copySign(1.0, FastMath.atan(new Complex(+0.0, -0.25)).getReal()) > 0.0);
        UnitTestUtils.customAssertEquals(new Complex(0.0, -0.25541281188299534160),
                                         FastMath.atan(new Complex(-0.0, -0.25)),
                                         1.0e-14);
        assertTrue(FastMath.copySign(1.0, FastMath.atan(new Complex(-0.0, -0.25)).getReal()) < 0.0);
    }

    @Test
    public void testAtanReal() {
        final Complex zP = new Complex(0.8734729023516287, 0.0);
        final Complex aP = new Complex(0.717964439926383,  0.0);
        Assertions.assertEquals(aP, zP.atan());
        Assertions.assertEquals(1.0, FastMath.copySign(1.0, zP.atan().getImaginary()), 1.0e-15);
        final Complex zM = new Complex(0.8734729023516287, -0.0);
        final Complex aM = new Complex(0.717964439926383,  -0.0);
        Assertions.assertEquals(aM, zM.atan());
        Assertions.assertEquals(-1.0, FastMath.copySign(1.0, zM.atan().getImaginary()), 1.0e-15);
    }

    @Test
    @Override
    public void testAtan2() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                final Complex z = build(x).atan2(build(y));
                final double  r = FastMath.atan2(x, y);
                checkRelative(r, new Complex(MathUtils.normalizeAngle(z.getReal(), r), z.getImaginary()));
            }
        }
    }

    @Test
    void testAtan2Complex() {
        for (double r1 : Arrays.asList(-3, 3)) {
            for (double i1 : Arrays.asList(-2, 0, 2)) {
                final Complex c1 = new Complex(r1, i1);
                for (double r2 : Arrays.asList(-1, 1)) {
                    for (double i2 : Arrays.asList(-5, 0, 5)) {
                        final Complex c2 = new Complex(r2, i2);
                        UnitTestUtils.customAssertEquals(c1.divide(c2), c1.atan2(c2).tan(), 1.0e-14);
                        final Complex atan   = c1.divide(c2).atan();
                        final Complex atan2  = c1.atan2(c2);
                        final double  deltaR = FastMath.abs(atan.getReal() - atan2.getReal()) / FastMath.PI;
                        assertTrue(FastMath.abs(deltaR - FastMath.rint(deltaR)) < 1.0e-14);
                        assertEquals(atan.getImaginary(), atan2.getImaginary(), 1.0e-14);
                    }
                }
            }
        }
    }

    @Test
    void testAtan2Real() {
        for (double r1 : Arrays.asList(-3, 3)) {
            final Complex c1 = new Complex(r1, 0);
            for (double r2 : Arrays.asList(-1, 1)) {
                final Complex c2 = new Complex(r2, 0);
                assertEquals(FastMath.atan2(r1, r2),
                                    MathUtils.normalizeAngle(c1.atan2(c2).getReal(), 0.0),
                                    1.0e-14);
            }
        }
    }

    @Override
    @Test
    public void testAtan2SpecialCases() {
        assertTrue(build(+0.0).atan2(build(+0.0)).isNaN());
        assertTrue(build(-0.0).atan2(build(+0.0)).isNaN());
        assertTrue(build(+0.0).atan2(build(-0.0)).isNaN());
        assertTrue(build(-0.0).atan2(build(-0.0)).isNaN());
    }

    @Test
    void testCosComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(-27.03495, -3.851153);
        UnitTestUtils.customAssertEquals(expected, z.cos(), 1.0e-5);
    }

    @Test
    void testCosNaN() {
        assertTrue(Complex.NaN.cos().isNaN());
    }

    @Test
    void testCosInf() {
        UnitTestUtils.customAssertSame(infNegInf, oneInf.cos());
        UnitTestUtils.customAssertSame(infInf, oneNegInf.cos());
        UnitTestUtils.customAssertSame(Complex.NaN, infOne.cos());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfOne.cos());
        UnitTestUtils.customAssertSame(Complex.NaN, infInf.cos());
        UnitTestUtils.customAssertSame(Complex.NaN, infNegInf.cos());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfInf.cos());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfNegInf.cos());
    }

    @Test
    void testCoshComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(-6.58066, -7.58155);
        UnitTestUtils.customAssertEquals(expected, z.cosh(), 1.0e-5);
    }

    @Test
    void testCoshNaN() {
        assertTrue(Complex.NaN.cosh().isNaN());
    }

    @Test
    void testCoshInf() {
        UnitTestUtils.customAssertSame(Complex.NaN, oneInf.cosh());
        UnitTestUtils.customAssertSame(Complex.NaN, oneNegInf.cosh());
        UnitTestUtils.customAssertSame(infInf, infOne.cosh());
        UnitTestUtils.customAssertSame(infNegInf, negInfOne.cosh());
        UnitTestUtils.customAssertSame(Complex.NaN, infInf.cosh());
        UnitTestUtils.customAssertSame(Complex.NaN, infNegInf.cosh());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfInf.cosh());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfNegInf.cosh());
    }

    @Test
    void testExpComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(-13.12878, -15.20078);
        UnitTestUtils.customAssertEquals(expected, z.exp(), 1.0e-5);
        UnitTestUtils.customAssertEquals(Complex.ONE,
                                         Complex.ZERO.exp(), 10e-12);
        Complex iPi = Complex.I.multiply(new Complex(pi,0));
        UnitTestUtils.customAssertEquals(Complex.ONE.negate(),
                                         iPi.exp(), 10e-12);
    }

    @Test
    void testExpNaN() {
        assertTrue(Complex.NaN.exp().isNaN());
    }

    @Test
    void testExpInf1() {
        UnitTestUtils.customAssertSame(Complex.NaN, oneInf.exp());
    }

    @Test
    void testExpInf2() {
        UnitTestUtils.customAssertSame(Complex.NaN, oneNegInf.exp());
    }

    @Test
    void testExpInf3() {
        UnitTestUtils.customAssertSame(infInf, infOne.exp());
    }

    @Test
    void testExpInf4() {
        final Complex exp = negInfOne.exp();
        UnitTestUtils.customAssertSame(Complex.ZERO, exp);
    }

    @Test
    void testExpInf5() {
        UnitTestUtils.customAssertSame(Complex.NaN, infInf.exp());
    }

    @Test
    void testExpInf6() {
        UnitTestUtils.customAssertSame(Complex.NaN, infNegInf.exp());
    }

    @Test
    void testExpInf7() {
        UnitTestUtils.customAssertSame(Complex.NaN, negInfInf.exp());
    }

    @Test
    void testExpInf8() {
        UnitTestUtils.customAssertSame(Complex.NaN, negInfNegInf.exp());
    }

    @Test
    void testExpM1() {
        UnitTestUtils.customAssertSame(Complex.NaN, Complex.NaN.expm1());
        final double testValue = FastMath.scalb(1.0, -30);
        assertEquals(FastMath.expm1(testValue), new Complex(testValue, 0).expm1().getReal(), 1.0e-30);
        assertTrue(FastMath.expm1(testValue) - new Complex(testValue).exp().subtract(1.0).getReal() > 4.0e-19);
        assertEquals(0.0, new Complex(0, testValue).expm1().getReal(), 1.0e-30);
        assertEquals(0.0, new Complex(0, testValue).expm1().getImaginary(), 1.0e-30);
    }

    @Test
    void testLogComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(1.60944, 0.927295);
        UnitTestUtils.customAssertEquals(expected, z.log(), 1.0e-5);
    }

    @Test
    void testLogNaN() {
        assertTrue(Complex.NaN.log().isNaN());
    }

    @Test
    void testLogInf() {
        UnitTestUtils.customAssertEquals(new Complex(inf, pi / 2),
                                         oneInf.log(), 10e-12);
        UnitTestUtils.customAssertEquals(new Complex(inf, -pi / 2),
                                         oneNegInf.log(), 10e-12);
        UnitTestUtils.customAssertEquals(infZero, infOne.log(), 10e-12);
        UnitTestUtils.customAssertEquals(new Complex(inf, pi),
                                         negInfOne.log(), 10e-12);
        UnitTestUtils.customAssertEquals(new Complex(inf, pi / 4),
                                         infInf.log(), 10e-12);
        UnitTestUtils.customAssertEquals(new Complex(inf, -pi / 4),
                                         infNegInf.log(), 10e-12);
        UnitTestUtils.customAssertEquals(new Complex(inf, 3d * pi / 4),
                                         negInfInf.log(), 10e-12);
        UnitTestUtils.customAssertEquals(new Complex(inf, - 3d * pi / 4),
                                         negInfNegInf.log(), 10e-12);
    }

    @Test
    void testLogBranchCut() {
        UnitTestUtils.customAssertEquals(new Complex(0.6931471805599453094, +3.1415926535897932384),
                                         FastMath.log(new Complex(-2.0, +0.0)),
                                         10e-12);
        UnitTestUtils.customAssertEquals(new Complex(0.6931471805599453094, -3.1415926535897932384),
                                         FastMath.log(new Complex(-2.0, -0.0)),
                                         10e-12);
    }

    @Test
    void testLogZero() {
        UnitTestUtils.customAssertSame(negInfZero, Complex.ZERO.log());
    }

    @Test
    void testLog1P() {
        Complex z = new Complex(2, 4);
        Complex expected = new Complex(1.60944, 0.927295);
        UnitTestUtils.customAssertEquals(expected, z.log1p(), 1.0e-5);
    }

    @Test
    void testLog10Complex() {
        UnitTestUtils.customAssertEquals(new Complex(2.0, 0.0), new Complex(100, 0).log10(), 1.0e-15);
        UnitTestUtils.customAssertEquals(new Complex(2.0, 0.5 * FastMath.PI / FastMath.log(10)), new Complex(0, 100).log10(), 1.0e-15);
    }

    @Test
    @Override
    public void testLog10() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            if (x < 0) {
                // special case for Complex
                assertTrue(Double.isNaN(FastMath.log10(x)));
                assertFalse(build(x).log10().isNaN());
            } else {
                checkRelative(FastMath.log10(x), build(x).log10());
            }
        }
    }

    @Test
    @Override
    public void testPowField() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            for (double y = 0.1; y < 4; y += 0.2) {
                if ( x < 0) {
                    // special case for Complex
                    assertTrue(Double.isNaN(FastMath.pow(x, y)));
                    assertFalse(build(x).pow(build(y)).isNaN());
                } else {
                    checkRelative(FastMath.pow(x, y), build(x).pow(build(y)));
                }
            }
        }
    }

    @Test
    @Override
    public void testPowDouble() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            for (double y = 0.1; y < 4; y += 0.2) {
                if ( x < 0) {
                    // special case for Complex
                    assertTrue(Double.isNaN(FastMath.pow(x, y)));
                    assertFalse(build(x).pow(y).isNaN());
                } else {
                    checkRelative(FastMath.pow(x, y), build(x).pow(y));
                }
            }
        }
    }

    @Test
    void testPow() {
        Complex x = new Complex(3, 4);
        Complex y = new Complex(5, 6);
        Complex expected = new Complex(-1.860893, 11.83677);
        UnitTestUtils.customAssertEquals(expected, x.pow(y), 1.0e-5);
        UnitTestUtils.customAssertEquals(new Complex(-46, 9).divide(2197), new Complex(2, -3).pow(new Complex(-3, 0)), 1.0e-15);
        UnitTestUtils.customAssertEquals(new Complex(-1, 0).divide(8), new Complex(-2, 0).pow(new Complex(-3, 0)), 1.0e-15);
        UnitTestUtils.customAssertEquals(new Complex(0, 2),
                                         new Complex(-4, 0).pow(new Complex(0.5, 0)),
                                         1.0e-15);
    }

    @Test
    void testPowNaNBase() {
        Complex x = new Complex(3, 4);
        assertTrue(Complex.NaN.pow(x).isNaN());
    }

    @Test
    void testPowNaNExponent() {
        Complex x = new Complex(3, 4);
        assertTrue(x.pow(Complex.NaN).isNaN());
    }

    @Test
    void testPowInf() {
        UnitTestUtils.customAssertSame(Complex.NaN, Complex.ONE.pow(oneInf));
        UnitTestUtils.customAssertSame(Complex.NaN, Complex.ONE.pow(oneNegInf));
        UnitTestUtils.customAssertSame(Complex.NaN, Complex.ONE.pow(infOne));
        UnitTestUtils.customAssertSame(Complex.NaN, Complex.ONE.pow(infInf));
        UnitTestUtils.customAssertSame(Complex.NaN, Complex.ONE.pow(infNegInf));
        UnitTestUtils.customAssertSame(Complex.NaN, Complex.ONE.pow(negInfInf));
        UnitTestUtils.customAssertSame(Complex.NaN, Complex.ONE.pow(negInfNegInf));
        UnitTestUtils.customAssertSame(Complex.INF, infOne.pow(Complex.ONE));
        UnitTestUtils.customAssertSame(Complex.INF, negInfOne.pow(Complex.ONE));
        UnitTestUtils.customAssertSame(Complex.INF, infInf.pow(Complex.ONE));
        UnitTestUtils.customAssertSame(Complex.INF, infNegInf.pow(Complex.ONE));
        UnitTestUtils.customAssertSame(Complex.INF, negInfInf.pow(Complex.ONE));
        UnitTestUtils.customAssertSame(Complex.INF, negInfNegInf.pow(Complex.ONE));
        UnitTestUtils.customAssertSame(Complex.NaN, negInfNegInf.pow(infNegInf));
        UnitTestUtils.customAssertSame(Complex.NaN, negInfNegInf.pow(negInfNegInf));
        UnitTestUtils.customAssertSame(Complex.NaN, negInfNegInf.pow(infInf));
        UnitTestUtils.customAssertSame(Complex.NaN, infInf.pow(infNegInf));
        UnitTestUtils.customAssertSame(Complex.NaN, infInf.pow(negInfNegInf));
        UnitTestUtils.customAssertSame(Complex.NaN, infInf.pow(infInf));
        UnitTestUtils.customAssertSame(Complex.NaN, infNegInf.pow(infNegInf));
        UnitTestUtils.customAssertSame(Complex.NaN, infNegInf.pow(negInfNegInf));
        UnitTestUtils.customAssertSame(Complex.NaN, infNegInf.pow(infInf));
    }

    @Test
    void testPowZero() {
        UnitTestUtils.customAssertEquals(Complex.ZERO,
                                         Complex.ZERO.pow(Complex.ONE), 1.0e-12);
        UnitTestUtils.customAssertSame(Complex.ONE,
                                       Complex.ZERO.pow(Complex.ZERO));
        UnitTestUtils.customAssertSame(Complex.NaN,
                                       Complex.ZERO.pow(Complex.I));
        UnitTestUtils.customAssertEquals(Complex.ONE,
                                         Complex.ONE.pow(Complex.ZERO), 10e-12);
        UnitTestUtils.customAssertEquals(Complex.ONE,
                                         Complex.I.pow(Complex.ZERO), 10e-12);
        UnitTestUtils.customAssertEquals(Complex.ONE,
                                         new Complex(-1, 3).pow(Complex.ZERO), 10e-12);
    }

    @Test
    void testZeroPow() {
        UnitTestUtils.customAssertEquals(Complex.ZERO, Complex.ZERO.pow(2.0), 1.0e-5);
    }

    @Test
    void testScalarPow() {
        Complex x = new Complex(3, 4);
        double yDouble = 5.0;
        Complex yComplex = new Complex(yDouble);
        assertEquals(x.pow(yComplex), x.pow(yDouble));
        UnitTestUtils.customAssertEquals(Complex.ONE.negate(), Complex.ONE.negate().pow(0.5).pow(2), 1.0e-15);
        UnitTestUtils.customAssertEquals(new Complex(2, 0), new Complex(4, 0).pow(0.5), 1.0e-15);
        UnitTestUtils.customAssertEquals(new Complex(2, 0), new Complex(4, 0).pow(new Complex(0.5, 0)), 1.0e-15);
    }

    @Test
    void testScalarPowNaNBase() {
        Complex x = Complex.NaN;
        double yDouble = 5.0;
        Complex yComplex = new Complex(yDouble);
        assertEquals(x.pow(yComplex), x.pow(yDouble));
    }

    @Test
    void testScalarPowNaNExponent() {
        Complex x = new Complex(3, 4);
        double yDouble = Double.NaN;
        Complex yComplex = new Complex(yDouble);
        assertEquals(x.pow(yComplex), x.pow(yDouble));
    }

    @Test
    void testScalarPowInf() {
       UnitTestUtils.customAssertSame(Complex.NaN, Complex.ONE.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.customAssertSame(Complex.NaN, Complex.ONE.pow(Double.NEGATIVE_INFINITY));
       UnitTestUtils.customAssertSame(Complex.INF, infOne.pow(1.0));
       UnitTestUtils.customAssertSame(Complex.INF, negInfOne.pow(1.0));
       UnitTestUtils.customAssertSame(Complex.INF, infInf.pow(1.0));
       UnitTestUtils.customAssertSame(Complex.INF, infNegInf.pow(1.0));
       UnitTestUtils.customAssertSame(Complex.INF, negInfInf.pow(10));
       UnitTestUtils.customAssertSame(Complex.INF, negInfNegInf.pow(1.0));
       UnitTestUtils.customAssertSame(Complex.NaN, negInfNegInf.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.customAssertSame(Complex.NaN, negInfNegInf.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.customAssertSame(Complex.NaN, infInf.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.customAssertSame(Complex.NaN, infInf.pow(Double.NEGATIVE_INFINITY));
       UnitTestUtils.customAssertSame(Complex.NaN, infNegInf.pow(Double.NEGATIVE_INFINITY));
       UnitTestUtils.customAssertSame(Complex.NaN, infNegInf.pow(Double.POSITIVE_INFINITY));
   }

    @Test
    void testScalarPowZero() {
       UnitTestUtils.customAssertEquals(Complex.ZERO, Complex.ZERO.pow(1.0), 1.0e-12);
       UnitTestUtils.customAssertSame(Complex.ONE, Complex.ZERO.pow(0.0));
       UnitTestUtils.customAssertEquals(Complex.ONE, Complex.ONE.pow(0.0), 10e-12);
       UnitTestUtils.customAssertEquals(Complex.ONE, Complex.I.pow(0.0), 10e-12);
       UnitTestUtils.customAssertEquals(Complex.ONE, new Complex(-1, 3).pow(0.0), 10e-12);
   }

    @Test
    void testpowNull() {
        assertThrows(NullArgumentException.class, () -> Complex.ONE.pow(null));
    }

    @Test
    void testSinComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(3.853738, -27.01681);
        UnitTestUtils.customAssertEquals(expected, z.sin(), 1.0e-5);
    }

    @Test
    void testSinInf() {
        UnitTestUtils.customAssertSame(infInf, oneInf.sin());
        UnitTestUtils.customAssertSame(infNegInf, oneNegInf.sin());
        UnitTestUtils.customAssertSame(Complex.NaN, infOne.sin());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfOne.sin());
        UnitTestUtils.customAssertSame(Complex.NaN, infInf.sin());
        UnitTestUtils.customAssertSame(Complex.NaN, infNegInf.sin());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfInf.sin());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfNegInf.sin());
    }

    @Test
    void testSinNaN() {
        assertTrue(Complex.NaN.sin().isNaN());
    }

    @Test
    void testSinhComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(-6.54812, -7.61923);
        UnitTestUtils.customAssertEquals(expected, z.sinh(), 1.0e-5);
    }

    @Test
    void testSinhNaN() {
        assertTrue(Complex.NaN.sinh().isNaN());
    }

    @Test
    void testSinhInf() {
        UnitTestUtils.customAssertSame(Complex.NaN, oneInf.sinh());
        UnitTestUtils.customAssertSame(Complex.NaN, oneNegInf.sinh());
        UnitTestUtils.customAssertSame(infInf, infOne.sinh());
        UnitTestUtils.customAssertSame(negInfInf, negInfOne.sinh());
        UnitTestUtils.customAssertSame(Complex.NaN, infInf.sinh());
        UnitTestUtils.customAssertSame(Complex.NaN, infNegInf.sinh());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfInf.sinh());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfNegInf.sinh());
    }

    @Test
    void testAsinhComplex() {
        for (double x = -2; x <= 2; x += 0.125) {
            for (double y = -2; y <= 2; y += 0.125) {
                final Complex z = new Complex(x, y);
                UnitTestUtils.customAssertEquals(z, z.asinh().sinh(), 1.0e-14);
            }
        }
    }

    @Test
    void testAsinhBranchCuts() {
        UnitTestUtils.customAssertEquals(new Complex(FastMath.log(2 + FastMath.sqrt(3)), 0.5 * FastMath.PI),
                                         new Complex(+0.0, 2.0).asinh(),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(new Complex(-FastMath.log(2 + FastMath.sqrt(3)), 0.5 * FastMath.PI),
                                         new Complex(-0.0, 2.0).asinh(),
                                         1.0e-14);
    }

    @Test
    void testAcoshComplex() {
        for (double x = -2; x <= 2; x += 0.125) {
            for (double y = -2; y <= 2; y += 0.125) {
                final Complex z = new Complex(x, y);
                UnitTestUtils.customAssertEquals(z, z.acosh().cosh(), 1.0e-14);
            }
        }
    }

    @Test
    void testAcoshBranchCuts() {
        UnitTestUtils.customAssertEquals(new Complex(FastMath.log(2 + FastMath.sqrt(3)), +FastMath.PI),
                                         new Complex(-2.0, +0.0).acosh(),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(new Complex(FastMath.log(2 + FastMath.sqrt(3)), -FastMath.PI),
                                         new Complex(-2.0, -0.0).acosh(),
                                         1.0e-14);
    }

    @Test
    void testAtanhComplex() {
        for (double x = -2; x <= 2; x += 0.125) {
            for (double y = -2; y <= 2; y += 0.125) {
                final Complex z = new Complex(x, y);
                if (FastMath.abs(x) == 1.0 && y == 0.0) {
                    assertTrue(z.atanh().isInfinite());
                } else {
                    UnitTestUtils.customAssertEquals(z, z.atanh().tanh(), 1.0e-14);
                }
            }
        }
    }

    @Test
    void testAtanhBranchCuts() {
        UnitTestUtils.customAssertEquals(new Complex(-0.5 * FastMath.log(3), +0.5 * FastMath.PI),
                                         new Complex(-2.0, +0.0).atanh(),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(new Complex(-0.5 * FastMath.log(3), -0.5 * FastMath.PI),
                                         new Complex(-2.0, -0.0).atanh(),
                                         1.0e-14);
    }

    @Test
    void testSqrtRealPositive() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(2, 1);
        UnitTestUtils.customAssertEquals(expected, z.sqrt(), 1.0e-5);
    }

    @Test
    void testSqrtRealZero() {
        Complex z = new Complex(0.0, 4);
        Complex expected = new Complex(1.41421, 1.41421);
        UnitTestUtils.customAssertEquals(expected, z.sqrt(), 1.0e-5);
    }

    @Test
    void testSqrtZero() {
        UnitTestUtils.customAssertEquals(Complex.ZERO, Complex.ZERO.sqrt(), 1.0e-15);
    }

    @Test
    void testSqrtRealNegative() {
        Complex z = new Complex(-3.0, 4);
        Complex expected = new Complex(1, 2);
        UnitTestUtils.customAssertEquals(expected, z.sqrt(), 1.0e-5);
    }

    @Test
    void testSqrtImaginaryZero() {
        Complex z = new Complex(-3.0, 0.0);
        Complex expected = new Complex(0.0, 1.73205);
        UnitTestUtils.customAssertEquals(expected, z.sqrt(), 1.0e-5);
    }

    @Test
    void testSqrtImaginaryNegative() {
        Complex z = new Complex(-3.0, -4.0);
        Complex expected = new Complex(1.0, -2.0);
        UnitTestUtils.customAssertEquals(expected, z.sqrt(), 1.0e-5);
    }

    @Test
    void testSqrtPolar() {
        double r = 1;
        for (int i = 0; i < 5; i++) {
            r += i;
            double theta = 0;
            for (int j =0; j < 11; j++) {
                theta += pi /12;
                Complex z = ComplexUtils.polar2Complex(r, theta);
                Complex sqrtz = ComplexUtils.polar2Complex(FastMath.sqrt(r), theta / 2);
                UnitTestUtils.customAssertEquals(sqrtz, z.sqrt(), 10e-12);
            }
        }
    }

    @Test
    void testSqrtNaN() {
        assertTrue(Complex.NaN.sqrt().isNaN());
    }

    @Test
    void testSqrtInf() {
        UnitTestUtils.customAssertSame(infNaN, oneInf.sqrt());
        UnitTestUtils.customAssertSame(infNaN, oneNegInf.sqrt());
        UnitTestUtils.customAssertSame(infZero, infOne.sqrt());
        UnitTestUtils.customAssertSame(zeroInf, negInfOne.sqrt());
        UnitTestUtils.customAssertSame(infNaN, infInf.sqrt());
        UnitTestUtils.customAssertSame(infNaN, infNegInf.sqrt());
        UnitTestUtils.customAssertSame(nanInf, negInfInf.sqrt());
        UnitTestUtils.customAssertSame(nanNegInf, negInfNegInf.sqrt());
    }

    @Test
    void testSqrt1z() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(4.08033, -2.94094);
        UnitTestUtils.customAssertEquals(expected, z.sqrt1z(), 1.0e-5);
    }

    @Test
    void testSqrt1zNaN() {
        assertTrue(Complex.NaN.sqrt1z().isNaN());
    }

    @Test
    @Override
    public void testCbrt() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            if ( x < 0) {
                // special case for Complex
                assertTrue(FastMath.cbrt(x) < 0);
                assertEquals(FastMath.PI / 3, build(x).cbrt().getArgument(), 1.0e-15);
            } else {
                checkRelative(FastMath.cbrt(x), build(x).cbrt());
            }
        }
    }

    @Test
    void testCbrtComplex() {
        Complex z = new Complex(15, 2);
        UnitTestUtils.customAssertEquals(z, z.square().multiply(z).cbrt(), 1.0e-14);
        Complex branchCutPlus = new Complex(-8.0, +0.0);
        Complex cbrtPlus = branchCutPlus.cbrt();
        UnitTestUtils.customAssertEquals(branchCutPlus, cbrtPlus.multiply(cbrtPlus).multiply(cbrtPlus), 1.0e-14);
        assertEquals(1.0, cbrtPlus.getReal(), 1.0e-15);
        assertEquals(FastMath.sqrt(3.0), cbrtPlus.getImaginary(), 1.0e-15);
        Complex branchCutMinus = new Complex(-8.0, -0.0);
        Complex cbrtMinus = branchCutMinus.cbrt();
        UnitTestUtils.customAssertEquals(branchCutMinus, cbrtMinus.multiply(cbrtMinus).multiply(cbrtMinus), 1.0e-14);
        assertEquals(1.0, cbrtMinus.getReal(), 1.0e-15);
        assertEquals(-FastMath.sqrt(3.0), cbrtMinus.getImaginary(), 1.0e-15);
    }

    @Test
    @Override
    public void testRootN() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            for (int n = 1; n < 5; ++n) {
                if (x < 0) {
                    // special case for Complex
                    final double doubleRoot = new Binary64(x).rootN(n).getReal();
                    if (n % 2 == 0) {
                        assertTrue(Double.isNaN(doubleRoot));
                    } else {
                        assertTrue(doubleRoot < 0);
                    }
                    assertEquals(FastMath.PI / n, build(x).rootN(n).getArgument(), 1.0e-15);
                } else {
                    checkRelative(FastMath.pow(x, 1.0 / n), build(x).rootN(n));
                }
            }
        }
    }

    @Test
    void testRootNComplex() {
        Complex z = new Complex(15, 2);
        UnitTestUtils.customAssertEquals(z, z.square().multiply(z).rootN(3), 1.0e-14);
        Complex branchCutPlus = new Complex(-8.0, +0.0);
        Complex cbrtPlus = branchCutPlus.rootN(3);
        UnitTestUtils.customAssertEquals(branchCutPlus, cbrtPlus.multiply(cbrtPlus).multiply(cbrtPlus), 1.0e-14);
        assertEquals(1.0, cbrtPlus.getReal(), 1.0e-15);
        assertEquals(FastMath.sqrt(3.0), cbrtPlus.getImaginary(), 1.0e-15);
        Complex branchCutMinus = new Complex(-8.0, -0.0);
        Complex cbrtMinus = branchCutMinus.rootN(3);
        UnitTestUtils.customAssertEquals(branchCutMinus, cbrtMinus.multiply(cbrtMinus).multiply(cbrtMinus), 1.0e-14);
        assertEquals(1.0, cbrtMinus.getReal(), 1.0e-15);
        assertEquals(-FastMath.sqrt(3.0), cbrtMinus.getImaginary(), 1.0e-15);
    }

    @Test
    void testTanComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(-0.000187346, 0.999356);
        UnitTestUtils.customAssertEquals(expected, z.tan(), 1.0e-5);
        /* Check that no overflow occurs (MATH-722) */
        Complex actual = new Complex(3.0, 1E10).tan();
        expected = new Complex(0, 1);
        UnitTestUtils.customAssertEquals(expected, actual, 1.0e-5);
        actual = new Complex(3.0, -1E10).tan();
        expected = new Complex(0, -1);
        UnitTestUtils.customAssertEquals(expected, actual, 1.0e-5);
    }

    @Test
    void testTanNaN() {
        assertTrue(Complex.NaN.tan().isNaN());
    }

    @Test
    void testTanInf() {
        UnitTestUtils.customAssertSame(Complex.valueOf(0.0, 1.0), oneInf.tan());
        UnitTestUtils.customAssertSame(Complex.valueOf(0.0, -1.0), oneNegInf.tan());
        UnitTestUtils.customAssertSame(Complex.NaN, infOne.tan());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfOne.tan());
        UnitTestUtils.customAssertSame(Complex.NaN, infInf.tan());
        UnitTestUtils.customAssertSame(Complex.NaN, infNegInf.tan());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfInf.tan());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfNegInf.tan());
    }

    @Test
    void testTanCritical() {
        UnitTestUtils.customAssertSame(infNaN, new Complex(pi/2, 0).tan());
        UnitTestUtils.customAssertSame(negInfNaN, new Complex(-pi/2, 0).tan());
    }

    @Test
    void testTanhComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(1.00071, 0.00490826);
        UnitTestUtils.customAssertEquals(expected, z.tanh(), 1.0e-5);
        /* Check that no overflow occurs (MATH-722) */
        Complex actual = new Complex(1E10, 3.0).tanh();
        expected = new Complex(1, 0);
        UnitTestUtils.customAssertEquals(expected, actual, 1.0e-5);
        actual = new Complex(-1E10, 3.0).tanh();
        expected = new Complex(-1, 0);
        UnitTestUtils.customAssertEquals(expected, actual, 1.0e-5);
    }

    @Test
    void testTanhNaN() {
        assertTrue(Complex.NaN.tanh().isNaN());
    }

    @Test
    void testTanhInf() {
        UnitTestUtils.customAssertSame(Complex.NaN, oneInf.tanh());
        UnitTestUtils.customAssertSame(Complex.NaN, oneNegInf.tanh());
        UnitTestUtils.customAssertSame(Complex.valueOf(1.0, 0.0), infOne.tanh());
        UnitTestUtils.customAssertSame(Complex.valueOf(-1.0, 0.0), negInfOne.tanh());
        UnitTestUtils.customAssertSame(Complex.NaN, infInf.tanh());
        UnitTestUtils.customAssertSame(Complex.NaN, infNegInf.tanh());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfInf.tanh());
        UnitTestUtils.customAssertSame(Complex.NaN, negInfNegInf.tanh());
    }

    @Test
    void testTanhCritical() {
        UnitTestUtils.customAssertSame(nanInf, new Complex(0, pi/2).tanh());
    }

    /** test issue MATH-221 */
    @Test
    void testMath221() {
        assertTrue(Complex.equals(new Complex(0,-1),
                                         new Complex(0,1).multiply(new Complex(-1,0))));
    }

    /**
     * Test: computing <b>third roots</b> of z.
     * <pre>
     * <code>
     * <b>z = -2 + 2 * i</b>
     *   ⇒ z_0 =  1      +          i
     *   ⇒ z_1 = -1.3660 + 0.3660 * i
     *   ⇒ z_2 =  0.3660 - 1.3660 * i
     * </code>
     * </pre>
     */
    @Test
    void testNthRoot_normal_thirdRoot() {
        // The complex number we want to compute all third-roots for.
        Complex z = new Complex(-2,2);
        // The List holding all third roots
        List<Complex> thirdRootsOfZ = z.nthRoot(3);
        // Returned Collection must not be empty!
        assertEquals(3, thirdRootsOfZ.size());
        // test z_0
        assertEquals(1.0,                  thirdRootsOfZ.get(0).getReal(),      1.0e-5);
        assertEquals(1.0,                  thirdRootsOfZ.get(0).getImaginary(), 1.0e-5);
        // test z_1
        assertEquals(-1.3660254037844386,  thirdRootsOfZ.get(1).getReal(),      1.0e-5);
        assertEquals(0.36602540378443843,  thirdRootsOfZ.get(1).getImaginary(), 1.0e-5);
        // test z_2
        assertEquals(0.366025403784439,    thirdRootsOfZ.get(2).getReal(),      1.0e-5);
        assertEquals(-1.3660254037844384,  thirdRootsOfZ.get(2).getImaginary(), 1.0e-5);
    }


    /**
     * Test: computing <b>fourth roots</b> of z.
     * <pre>
     * <code>
     * <b>z = 5 - 2 * i</b>
     *   ⇒ z_0 =  1.5164 - 0.1446 * i
     *   ⇒ z_1 =  0.1446 + 1.5164 * i
     *   ⇒ z_2 = -1.5164 + 0.1446 * i
     *   ⇒ z_3 = -1.5164 - 0.1446 * i
     * </code>
     * </pre>
     */
    @Test
    void testNthRoot_normal_fourthRoot() {
        // The complex number we want to compute all third-roots for.
        Complex z = new Complex(5,-2);
        // The List holding all fourth roots
        List<Complex> fourthRootsOfZ = z.nthRoot(4);
        // Returned Collection must not be empty!
        assertEquals(4, fourthRootsOfZ.size());
        // test z_0
        assertEquals(1.5164629308487783,     fourthRootsOfZ.get(0).getReal(),      1.0e-5);
        assertEquals(-0.14469266210702247,   fourthRootsOfZ.get(0).getImaginary(), 1.0e-5);
        // test z_1
        assertEquals(0.14469266210702256,    fourthRootsOfZ.get(1).getReal(),      1.0e-5);
        assertEquals(1.5164629308487783,     fourthRootsOfZ.get(1).getImaginary(), 1.0e-5);
        // test z_2
        assertEquals(-1.5164629308487783,    fourthRootsOfZ.get(2).getReal(),      1.0e-5);
        assertEquals(0.14469266210702267,    fourthRootsOfZ.get(2).getImaginary(), 1.0e-5);
        // test z_3
        assertEquals(-0.14469266210702275,   fourthRootsOfZ.get(3).getReal(),      1.0e-5);
        assertEquals(-1.5164629308487783,    fourthRootsOfZ.get(3).getImaginary(), 1.0e-5);
    }

    /**
     * Test: computing <b>third roots</b> of z.
     * <pre>
     * <code>
     * <b>z = 8</b>
     *   ⇒ z_0 =  2
     *   ⇒ z_1 = -1 + 1.73205 * i
     *   ⇒ z_2 = -1 - 1.73205 * i
     * </code>
     * </pre>
     */
    @Test
    void testNthRoot_cornercase_thirdRoot_imaginaryPartEmpty() {
        // The number 8 has three third roots. One we all already know is the number 2.
        // But there are two more complex roots.
        Complex z = new Complex(8,0);
        // The List holding all third roots
        List<Complex> thirdRootsOfZ = z.nthRoot(3);
        // Returned Collection must not be empty!
        assertEquals(3, thirdRootsOfZ.size());
        // test z_0
        assertEquals(2.0,                thirdRootsOfZ.get(0).getReal(),      1.0e-5);
        assertEquals(0.0,                thirdRootsOfZ.get(0).getImaginary(), 1.0e-5);
        // test z_1
        assertEquals(-1.0,               thirdRootsOfZ.get(1).getReal(),      1.0e-5);
        assertEquals(1.7320508075688774, thirdRootsOfZ.get(1).getImaginary(), 1.0e-5);
        // test z_2
        assertEquals(-1.0,               thirdRootsOfZ.get(2).getReal(),      1.0e-5);
        assertEquals(-1.732050807568877, thirdRootsOfZ.get(2).getImaginary(), 1.0e-5);
    }


    /**
     * Test: computing <b>third roots</b> of z with real part 0.
     * <pre>
     * <code>
     * <b>z = 2 * i</b>
     *   ⇒ z_0 =  1.0911 + 0.6299 * i
     *   ⇒ z_1 = -1.0911 + 0.6299 * i
     *   ⇒ z_2 = -2.3144 - 1.2599 * i
     * </code>
     * </pre>
     */
    @Test
    void testNthRoot_cornercase_thirdRoot_realPartZero() {
        // complex number with only imaginary part
        Complex z = new Complex(0,2);
        // The List holding all third roots
        List<Complex> thirdRootsOfZ = z.nthRoot(3);
        // Returned Collection must not be empty!
        assertEquals(3, thirdRootsOfZ.size());
        // test z_0
        assertEquals(1.0911236359717216,      thirdRootsOfZ.get(0).getReal(),      1.0e-5);
        assertEquals(0.6299605249474365,      thirdRootsOfZ.get(0).getImaginary(), 1.0e-5);
        // test z_1
        assertEquals(-1.0911236359717216,     thirdRootsOfZ.get(1).getReal(),      1.0e-5);
        assertEquals(0.6299605249474365,      thirdRootsOfZ.get(1).getImaginary(), 1.0e-5);
        // test z_2
        assertEquals(-2.3144374213981936E-16, thirdRootsOfZ.get(2).getReal(),      1.0e-5);
        assertEquals(-1.2599210498948732,     thirdRootsOfZ.get(2).getImaginary(), 1.0e-5);
    }

    /**
     * Test cornercases with NaN and Infinity.
     */
    @Test
    void testNthRoot_cornercase_NAN_Inf() {
        // NaN + finite -> NaN
        List<Complex> roots = oneNaN.nthRoot(3);
        assertEquals(1,roots.size());
        assertEquals(Complex.NaN, roots.get(0));

        roots = nanZero.nthRoot(3);
        assertEquals(1,roots.size());
        assertEquals(Complex.NaN, roots.get(0));

        // NaN + infinite -> NaN
        roots = nanInf.nthRoot(3);
        assertEquals(1,roots.size());
        assertEquals(Complex.NaN, roots.get(0));

        // finite + infinite -> Inf
        roots = oneInf.nthRoot(3);
        assertEquals(1,roots.size());
        assertEquals(Complex.INF, roots.get(0));

        // infinite + infinite -> Inf
        roots = negInfInf.nthRoot(3);
        assertEquals(1,roots.size());
        assertEquals(Complex.INF, roots.get(0));
    }

    @Test
    void testNthRootError() {
        try {
            Complex.ONE.nthRoot(-1);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.CANNOT_COMPUTE_NTH_ROOT_FOR_NEGATIVE_N,
                                miae.getSpecifier());
        }
    }

    @Test
    void testIsMathematicalInteger() {
        doTestIsMathmeaticalInteger(-0.0, true);
        doTestIsMathmeaticalInteger(+0.0, true);
        doTestIsMathmeaticalInteger(-0.5, false);
        doTestIsMathmeaticalInteger(+0.5, false);
        doTestIsMathmeaticalInteger(Double.NaN, false);
        doTestIsMathmeaticalInteger(Double.POSITIVE_INFINITY, false);
        doTestIsMathmeaticalInteger(Double.NEGATIVE_INFINITY, false);
        doTestIsMathmeaticalInteger(Double.MIN_NORMAL, false);
        doTestIsMathmeaticalInteger(Double.MIN_VALUE, false);
    }

    private void doTestIsMathmeaticalInteger(double imaginary, boolean expectedForInteger) {
        assertFalse(new Complex(Double.NaN, imaginary).isMathematicalInteger());
        assertFalse(new Complex(Double.POSITIVE_INFINITY, imaginary).isMathematicalInteger());
        assertFalse(new Complex(Double.NEGATIVE_INFINITY, imaginary).isMathematicalInteger());
        assertFalse(new Complex(Double.MIN_NORMAL, imaginary).isMathematicalInteger());
        assertFalse(new Complex(Double.MIN_VALUE, imaginary).isMathematicalInteger());

        assertEquals(expectedForInteger, new Complex(-0.0, imaginary).isMathematicalInteger());
        assertEquals(expectedForInteger, new Complex(+0.0, imaginary).isMathematicalInteger());

        for (int i = -1000; i < 1000; ++i) {
            final double d = i;
            assertEquals(expectedForInteger, new Complex(d, imaginary).isMathematicalInteger());
            assertFalse(new Complex(FastMath.nextAfter(d, Double.POSITIVE_INFINITY), imaginary).isMathematicalInteger());
            assertFalse(new Complex(FastMath.nextAfter(d, Double.NEGATIVE_INFINITY), imaginary).isMathematicalInteger());
        }

        double minNoFractional = 0x1L << 52;
        assertEquals(expectedForInteger, new Complex(minNoFractional, imaginary).isMathematicalInteger());
        assertFalse(new Complex(minNoFractional - 0.5, imaginary).isMathematicalInteger());
        assertEquals(expectedForInteger, new Complex(minNoFractional + 0.5, imaginary).isMathematicalInteger());

    }

    /**
     * Test standard values
     */
    @Test
    void testGetArgument() {
        Complex z = new Complex(1, 0);
        assertEquals(0.0, z.getArgument(), 1.0e-12);

        z = new Complex(1, 1);
        assertEquals(FastMath.PI/4, z.getArgument(), 1.0e-12);

        z = new Complex(0, 1);
        assertEquals(FastMath.PI/2, z.getArgument(), 1.0e-12);

        z = new Complex(-1, 1);
        assertEquals(3 * FastMath.PI/4, z.getArgument(), 1.0e-12);

        z = new Complex(-1, 0);
        assertEquals(FastMath.PI, z.getArgument(), 1.0e-12);

        z = new Complex(-1, -1);
        assertEquals(-3 * FastMath.PI/4, z.getArgument(), 1.0e-12);

        z = new Complex(0, -1);
        assertEquals(-FastMath.PI/2, z.getArgument(), 1.0e-12);

        z = new Complex(1, -1);
        assertEquals(-FastMath.PI/4, z.getArgument(), 1.0e-12);

    }

    /**
     * Verify atan2-style handling of infinite parts
     */
    @Test
    void testGetArgumentInf() {
        assertEquals(FastMath.PI/4, infInf.getArgument(), 1.0e-12);
        assertEquals(FastMath.PI/2, oneInf.getArgument(), 1.0e-12);
        assertEquals(0.0, infOne.getArgument(), 1.0e-12);
        assertEquals(FastMath.PI/2, zeroInf.getArgument(), 1.0e-12);
        assertEquals(0.0, infZero.getArgument(), 1.0e-12);
        assertEquals(FastMath.PI, negInfOne.getArgument(), 1.0e-12);
        assertEquals(-3.0*FastMath.PI/4, negInfNegInf.getArgument(), 1.0e-12);
        assertEquals(-FastMath.PI/2, oneNegInf.getArgument(), 1.0e-12);
    }

    /**
     * Verify that either part NaN results in NaN
     */
    @Test
    void testGetArgumentNaN() {
        assertTrue(Double.isNaN(nanZero.getArgument()));
        assertTrue(Double.isNaN(zeroNaN.getArgument()));
        assertTrue(Double.isNaN(Complex.NaN.getArgument()));
    }

    @Test
    void testValueOf() {
        assertEquals(2.0, Complex.valueOf(2.0).getReal(), 1.0e-15);
        assertEquals(0.0, Complex.valueOf(2.0).getImaginary(), 1.0e-15);
        assertTrue(Complex.valueOf(Double.NaN).isNaN());
        assertTrue(Complex.valueOf(Double.POSITIVE_INFINITY).isInfinite());
        assertEquals( 2.0, Complex.valueOf(2.0, -1.0).getReal(), 1.0e-15);
        assertEquals(-1.0, Complex.valueOf(2.0, -1.0).getImaginary(), 1.0e-15);
        assertTrue(Complex.valueOf(Double.NaN, 0.0).isNaN());
        assertTrue(Complex.valueOf(Double.POSITIVE_INFINITY, 0.0).isInfinite());
        assertTrue(Complex.valueOf(Double.NaN, -1.0).isNaN());
        assertTrue(Complex.valueOf(Double.POSITIVE_INFINITY, -1.0).isInfinite());
        assertTrue(Complex.valueOf(0.0, Double.NaN).isNaN());
        assertTrue(Complex.valueOf(0.0, Double.POSITIVE_INFINITY).isInfinite());
        assertTrue(Complex.valueOf(-1.0, Double.NaN).isNaN());
        assertTrue(Complex.valueOf(-1.0, Double.POSITIVE_INFINITY).isInfinite());
    }

    @Test
    void testField() {
        assertEquals(ComplexField.getInstance(), Complex.ZERO.getField());
    }

    @Test
    void testToString() {
        assertEquals("(1.0, -2.0)", new Complex(1, -2).toString());
    }

    @Test
    void testScalbComplex() {
        assertEquals(0.125,  new Complex(2.0, 1.0).scalb(-4).getReal(), 1.0e-15);
        assertEquals(0.0625, new Complex(2.0, 1.0).scalb(-4).getImaginary(), 1.0e-15);
    }

    @Test
    void testHypotComplex() {
        assertEquals(5.8269600298808519855, new Complex(3, 4).hypot(new Complex(5, 6)).getReal(), 1.0e-15);
        assertEquals(7.2078750814528590485, new Complex(3, 4).hypot(new Complex(5, 6)).getImaginary(), 1.0e-15);
    }

    @Test
    void testCeilComplex() {
        for (double x = -3.9; x < 3.9; x += 0.05) {
            for (double y = -3.9; y < 3.9; y += 0.05) {
                final Complex z = new Complex(x, y);
                assertEquals(FastMath.ceil(x), z.ceil().getReal(), 1.0e-15);
                assertEquals(FastMath.ceil(y), z.ceil().getImaginary(), 1.0e-15);
            }
        }
    }

    @Test
    void testFloorComplex() {
        for (double x = -3.9; x < 3.9; x += 0.05) {
            for (double y = -3.9; y < 3.9; y += 0.05) {
                final Complex z = new Complex(x, y);
                assertEquals(FastMath.floor(x), z.floor().getReal(), 1.0e-15);
                assertEquals(FastMath.floor(y), z.floor().getImaginary(), 1.0e-15);
            }
        }
    }

    @Test
    void testRintComplex() {
        for (double x = -3.9; x < 3.9; x += 0.05) {
            for (double y = -3.9; y < 3.9; y += 0.05) {
                final Complex z = new Complex(x, y);
                assertEquals(FastMath.rint(x), z.rint().getReal(), 1.0e-15);
                assertEquals(FastMath.rint(y), z.rint().getImaginary(), 1.0e-15);
            }
        }
    }

    @Test
    void testRemainderComplexComplex() {
        for (double x1 = -3.9; x1 < 3.9; x1 += 0.125) {
            for (double y1 = -3.9; y1 < 3.9; y1 += 0.125) {
                final Complex z1 = new Complex(x1, y1);
                for (double x2 = -3.92; x2 < 3.9; x2 += 0.125) {
                    for (double y2 = -3.92; y2 < 3.9; y2 += 0.125) {
                        final Complex z2 = new Complex(x2, y2);
                        final Complex r  = z1.remainder(z2);
                        final Complex q  = z1.subtract(r).divide(z2);
                        assertTrue(r.norm() <= z2.norm());
                        assertEquals(FastMath.rint(q.getReal()), q.getReal(), 2.0e-14);
                        assertEquals(FastMath.rint(q.getImaginary()), q.getImaginary(), 2.0e-14);
                    }
                }
            }
        }
    }

    @Test
    void testRemainderComplexDouble() {
        for (double x1 = -3.9; x1 < 3.9; x1 += 0.125) {
            for (double y1 = -3.9; y1 < 3.9; y1 += 0.125) {
                final Complex z1 = new Complex(x1, y1);
                for (double a = -3.92; a < 3.9; a += 0.125) {
                        final Complex r  = z1.remainder(a);
                        final Complex q  = z1.subtract(r).divide(a);
                        assertTrue(r.norm() <= FastMath.abs(a));
                        assertEquals(FastMath.rint(q.getReal()), q.getReal(), 2.0e-14);
                        assertEquals(FastMath.rint(q.getImaginary()), q.getImaginary(), 2.0e-14);
                }
            }
        }
    }

    @Test
    void testRemainderAxKr() {
        checkRemainder(new Complex(14, -5), new Complex(3, 4), new Complex(-1.0,  0.0));
        checkRemainder(new Complex(26, 120), new Complex(37, 226), new Complex(-11.0, -106.0));
        checkRemainder(new Complex(9.4, 6), new Complex(1.0, 1.0), new Complex(-0.6, 0.0));
        checkRemainder(new Complex(-5.89, 0.33), new Complex(2.4, -0.123), new Complex(-1.09, 0.084));
    }

    private void checkRemainder(final Complex c1, final Complex c2, final Complex expectedRemainder) {

        final Complex remainder = c1.remainder(c2);
        assertEquals(expectedRemainder.getReal(),      remainder.getReal(),      1.0e-15);
        assertEquals(expectedRemainder.getImaginary(), remainder.getImaginary(), 1.0e-15);

        final Complex crossCheck = c1.subtract(remainder).divide(c2);
        assertTrue(Precision.isMathematicalInteger(crossCheck.getReal()));
        assertTrue(Precision.isMathematicalInteger(crossCheck.getImaginary()));

    }

    @Test
    void testCopySignFieldComplex() {
        for (double x1 = -3.9; x1 < 3.9; x1 += 0.08) {
            for (double y1 = -3.9; y1 < 3.9; y1 += 0.08) {
                final Complex z1 = new Complex(x1, y1);
                for (double x2 = -3.9; x2 < 3.9; x2 += 0.08) {
                    for (double y2 = -3.9; y2 < 3.9; y2 += 0.08) {
                        final Complex z2 = new Complex(x2, y2);
                        assertEquals(FastMath.copySign(x1, x2), z1.copySign(z2).getReal(), 1.0e-15);
                        assertEquals(FastMath.copySign(y1, y2), z1.copySign(z2).getImaginary(), 1.0e-15);
                    }
                }
            }
        }
    }

    @Test
    void testCopySignDoubleComplex() {
        for (double x1 = -3.9; x1 < 3.9; x1 += 0.05) {
            for (double y1 = -3.9; y1 < 3.9; y1 += 0.05) {
                final Complex z1 = new Complex(x1, y1);
                for (double r = -3.9; r < 3.9; r += 0.05) {
                    assertEquals(FastMath.copySign(x1, r), z1.copySign(r).getReal(), 1.0e-15);
                    assertEquals(FastMath.copySign(y1, r), z1.copySign(r).getImaginary(), 1.0e-15);
                }
            }
        }
    }

    @Test
    void testSignComplex() {
        for (double x = -3.9; x < 3.9; x += 0.05) {
            for (double y = -3.9; y < 3.9; y += 0.05) {
                final Complex z = new Complex(x, y);
                assertEquals(1.0, z.sign().norm(), 1.0e-15);
                assertEquals(FastMath.copySign(1, FastMath.signum(x)), FastMath.copySign(1, z.sign().getRealPart()), 1.0e-15);
                assertEquals(FastMath.copySign(1, FastMath.signum(y)), FastMath.copySign(1, z.sign().getImaginaryPart()), 1.0e-15);
            }
        }
        assertTrue(Complex.NaN.sign().isNaN());
        for (int sR : Arrays.asList(-1, +1)) {
            for (int sI : Arrays.asList(-1, +1)) {
                Complex z = new Complex(FastMath.copySign(0, sR), FastMath.copySign(0, sI));
                assertTrue(z.isZero());
                Complex zSign = z.sign();
                assertTrue(zSign.isZero());
                assertEquals(sR, FastMath.copySign(1, zSign.getRealPart()), 1.0e-15);
                assertEquals(sI, FastMath.copySign(1, zSign.getImaginaryPart()), 1.0e-15);
            }
        }
    }

    @Test
    void testLinearCombination1() {
        final Complex[] a = new Complex[] {
            new Complex(-1321008684645961.0 / 268435456.0,
                        +5774608829631843.0 / 268435456.0),
            new Complex(-7645843051051357.0 / 8589934592.0,
                        0.0)
        };
        final Complex[] b = new Complex[] {
            new Complex(-5712344449280879.0 / 2097152.0,
                        -4550117129121957.0 / 2097152.0),
            new Complex(8846951984510141.0 / 131072.0,
                        0.0)
        };

        final Complex abSumInline = Complex.ZERO.linearCombination(a[0], b[0],
                                                                  a[1], b[1]);
        final Complex abSumArray = Complex.ZERO.linearCombination(a, b);

        UnitTestUtils.customAssertEquals(abSumInline, abSumArray, 0);
        UnitTestUtils.customAssertEquals(-1.8551294182586248737720779899, abSumInline.getReal(), 1.0e-15);

        final Complex naive = a[0].multiply(b[0]).add(a[1].multiply(b[1]));
        assertTrue(naive.subtract(abSumInline).norm() > 1.5);

    }

    @Test
    void testSignedZeroEquality() {

        assertFalse(new Complex(-0.0, 1.0).isZero());
        assertFalse(new Complex(+0.0, 1.0).isZero());
        assertFalse(new Complex( 1.0, -0.0).isZero());
        assertFalse(new Complex( 1.0, +0.0).isZero());

        assertTrue(new Complex(-0.0, -0.0).isZero());
        assertTrue(new Complex(-0.0, +0.0).isZero());
        assertTrue(new Complex(+0.0, -0.0).isZero());
        assertTrue(new Complex(+0.0, +0.0).isZero());

        assertNotEquals(Complex.ZERO, new Complex(-0.0, -0.0));
        assertNotEquals(Complex.ZERO, new Complex(-0.0, +0.0));
        assertNotEquals(Complex.ZERO, new Complex(+0.0, -0.0));
        assertEquals(Complex.ZERO, new Complex(+0.0, +0.0));

    }

    @Test
    void testSerial() {
        Complex z = new Complex(3.0, 4.0);
        assertEquals(z, UnitTestUtils.serializeAndRecover(z));
        Complex ncmplx = (Complex)UnitTestUtils.serializeAndRecover(oneNaN);
        assertEquals(nanZero, ncmplx);
        assertTrue(ncmplx.isNaN());
        Complex infcmplx = (Complex)UnitTestUtils.serializeAndRecover(infInf);
        assertEquals(infInf, infcmplx);
        assertTrue(infcmplx.isInfinite());
        TestComplex tz = new TestComplex(3.0, 4.0);
        assertEquals(tz, UnitTestUtils.serializeAndRecover(tz));
        TestComplex ntcmplx = (TestComplex)UnitTestUtils.serializeAndRecover(new TestComplex(oneNaN));
        assertEquals(nanZero, ntcmplx);
        assertTrue(ntcmplx.isNaN());
        TestComplex inftcmplx = (TestComplex)UnitTestUtils.serializeAndRecover(new TestComplex(infInf));
        assertEquals(infInf, inftcmplx);
        assertTrue(inftcmplx.isInfinite());
    }

    @Test
    public void testIsSmall() {
        Assertions.assertTrue(new Complex(1.0e-13, 1.0e-13).isSmall(new Complex(1.0, 1.0), 1.0e-12));
        Assertions.assertFalse(new Complex(1.0e-11, 1.0e-13).isSmall(new Complex(1.0, 1.0), 1.0e-12));
        Assertions.assertFalse(new Complex(1.0e-13, 1.0e-11).isSmall(new Complex(1.0, 1.0), 1.0e-12));
    }

    /**
     * Class to test extending Complex
     */
    public static class TestComplex extends Complex {

        /**
         * Serialization identifier.
         */
        private static final long serialVersionUID = 3268726724160389237L;

        public TestComplex(double real, double imaginary) {
            super(real, imaginary);
        }

        public TestComplex(Complex other) {
            this(other.getReal(), other.getImaginary());
        }

        @Override
        protected TestComplex createComplex(double real, double imaginary) {
            return new TestComplex(real, imaginary);
        }

    }
}
