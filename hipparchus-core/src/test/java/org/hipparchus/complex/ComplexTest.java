/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus.complex;

import java.util.List;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;


/**
 */
public class ComplexTest {


    private double inf = Double.POSITIVE_INFINITY;
    private double neginf = Double.NEGATIVE_INFINITY;
    private double nan = Double.NaN;
    private double pi = FastMath.PI;
    private Complex oneInf = new Complex(1, inf);
    private Complex oneNegInf = new Complex(1, neginf);
    private Complex infOne = new Complex(inf, 1);
    private Complex infZero = new Complex(inf, 0);
    private Complex infNaN = new Complex(inf, nan);
    private Complex infNegInf = new Complex(inf, neginf);
    private Complex infInf = new Complex(inf, inf);
    private Complex negInfInf = new Complex(neginf, inf);
    private Complex negInfZero = new Complex(neginf, 0);
    private Complex negInfOne = new Complex(neginf, 1);
    private Complex negInfNaN = new Complex(neginf, nan);
    private Complex negInfNegInf = new Complex(neginf, neginf);
    private Complex oneNaN = new Complex(1, nan);
    private Complex zeroInf = new Complex(0, inf);
    private Complex zeroNaN = new Complex(0, nan);
    private Complex nanInf = new Complex(nan, inf);
    private Complex nanNegInf = new Complex(nan, neginf);
    private Complex nanZero = new Complex(nan, 0);

    @Test
    public void testConstructor() {
        Complex z = new Complex(3.0, 4.0);
        Assert.assertEquals(3.0, z.getReal(), 1.0e-5);
        Assert.assertEquals(4.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    public void testConstructorNaN() {
        Complex z = new Complex(3.0, Double.NaN);
        Assert.assertTrue(z.isNaN());

        z = new Complex(nan, 4.0);
        Assert.assertTrue(z.isNaN());

        z = new Complex(3.0, 4.0);
        Assert.assertFalse(z.isNaN());
    }

    @Test
    public void testAbs() {
        Complex z = new Complex(3.0, 4.0);
        Assert.assertEquals(5.0, z.abs(), 1.0e-5);
    }

    @Test
    public void testAbsNaN() {
        Assert.assertTrue(Double.isNaN(Complex.NaN.abs()));
        Complex z = new Complex(inf, nan);
        Assert.assertTrue(Double.isNaN(z.abs()));
    }

    @Test
    public void testAbsInfinite() {
        Complex z = new Complex(inf, 0);
        Assert.assertEquals(inf, z.abs(), 0);
        z = new Complex(0, neginf);
        Assert.assertEquals(inf, z.abs(), 0);
        z = new Complex(inf, neginf);
        Assert.assertEquals(inf, z.abs(), 0);
    }

    @Test
    public void testAdd() {
        Complex x = new Complex(3.0, 4.0);
        Complex y = new Complex(5.0, 6.0);
        Complex z = x.add(y);
        Assert.assertEquals(8.0, z.getReal(), 1.0e-5);
        Assert.assertEquals(10.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    public void testAddNaN() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.add(Complex.NaN);
        Assert.assertSame(Complex.NaN, z);
        z = new Complex(1, nan);
        Complex w = x.add(z);
        Assert.assertSame(Complex.NaN, w);
    }

    @Test
    public void testAddInf() {
        Complex x = new Complex(1, 1);
        Complex z = new Complex(inf, 0);
        Complex w = x.add(z);
        Assert.assertEquals(w.getImaginary(), 1, 0);
        Assert.assertEquals(inf, w.getReal(), 0);

        x = new Complex(neginf, 0);
        Assert.assertTrue(Double.isNaN(x.add(z).getReal()));
    }


    @Test
    public void testScalarAdd() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = 2.0;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.add(yComplex), x.add(yDouble));
    }

    @Test
    public void testScalarAddNaN() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = Double.NaN;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.add(yComplex), x.add(yDouble));
    }

    @Test
    public void testScalarAddInf() {
        Complex x = new Complex(1, 1);
        double yDouble = Double.POSITIVE_INFINITY;

        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.add(yComplex), x.add(yDouble));

        x = new Complex(neginf, 0);
        Assert.assertEquals(x.add(yComplex), x.add(yDouble));
    }

    @Test
    public void testConjugate() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.conjugate();
        Assert.assertEquals(3.0, z.getReal(), 1.0e-5);
        Assert.assertEquals(-4.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    public void testConjugateNaN() {
        Complex z = Complex.NaN.conjugate();
        Assert.assertTrue(z.isNaN());
    }

    @Test
    public void testConjugateInfiinite() {
        Complex z = new Complex(0, inf);
        Assert.assertEquals(neginf, z.conjugate().getImaginary(), 0);
        z = new Complex(0, neginf);
        Assert.assertEquals(inf, z.conjugate().getImaginary(), 0);
    }

    @Test
    public void testDivide() {
        Complex x = new Complex(3.0, 4.0);
        Complex y = new Complex(5.0, 6.0);
        Complex z = x.divide(y);
        Assert.assertEquals(39.0 / 61.0, z.getReal(), 1.0e-5);
        Assert.assertEquals(2.0 / 61.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    public void testDivideReal() {
        Complex x = new Complex(2d, 3d);
        Complex y = new Complex(2d, 0d);
        Assert.assertEquals(new Complex(1d, 1.5), x.divide(y));

    }

    @Test
    public void testDivideImaginary() {
        Complex x = new Complex(2d, 3d);
        Complex y = new Complex(0d, 2d);
        Assert.assertEquals(new Complex(1.5d, -1d), x.divide(y));
    }

    @Test
    public void testDivideInf() {
        Complex x = new Complex(3, 4);
        Complex w = new Complex(neginf, inf);
        Assert.assertTrue(x.divide(w).equals(Complex.ZERO));

        Complex z = w.divide(x);
        Assert.assertTrue(Double.isNaN(z.getReal()));
        Assert.assertEquals(inf, z.getImaginary(), 0);

        w = new Complex(inf, inf);
        z = w.divide(x);
        Assert.assertTrue(Double.isNaN(z.getImaginary()));
        Assert.assertEquals(inf, z.getReal(), 0);

        w = new Complex(1, inf);
        z = w.divide(w);
        Assert.assertTrue(Double.isNaN(z.getReal()));
        Assert.assertTrue(Double.isNaN(z.getImaginary()));
    }

    @Test
    public void testDivideZero() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.divide(Complex.ZERO);
        // Assert.assertEquals(z, Complex.INF); // See MATH-657
        Assert.assertEquals(z, Complex.NaN);
    }

    @Test
    public void testDivideZeroZero() {
        Complex x = new Complex(0.0, 0.0);
        Complex z = x.divide(Complex.ZERO);
        Assert.assertEquals(z, Complex.NaN);
    }

    @Test
    public void testDivideNaN() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.divide(Complex.NaN);
        Assert.assertTrue(z.isNaN());
    }

    @Test
    public void testDivideNaNInf() {
       Complex z = oneInf.divide(Complex.ONE);
       Assert.assertTrue(Double.isNaN(z.getReal()));
       Assert.assertEquals(inf, z.getImaginary(), 0);

       z = negInfNegInf.divide(oneNaN);
       Assert.assertTrue(Double.isNaN(z.getReal()));
       Assert.assertTrue(Double.isNaN(z.getImaginary()));

       z = negInfInf.divide(Complex.ONE);
       Assert.assertTrue(Double.isNaN(z.getReal()));
       Assert.assertTrue(Double.isNaN(z.getImaginary()));
    }

    @Test
    public void testScalarDivide() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = 2.0;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.divide(yComplex), x.divide(yDouble));
    }

    @Test
    public void testScalarDivideNaN() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = Double.NaN;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.divide(yComplex), x.divide(yDouble));
    }

    @Test
    public void testScalarDivideInf() {
        Complex x = new Complex(1,1);
        double yDouble = Double.POSITIVE_INFINITY;
        Complex yComplex = new Complex(yDouble);
        UnitTestUtils.assertEquals(x.divide(yComplex), x.divide(yDouble), 0);

        yDouble = Double.NEGATIVE_INFINITY;
        yComplex = new Complex(yDouble);
        UnitTestUtils.assertEquals(x.divide(yComplex), x.divide(yDouble), 0);

        x = new Complex(1, Double.NEGATIVE_INFINITY);
        UnitTestUtils.assertEquals(x.divide(yComplex), x.divide(yDouble), 0);
    }

    @Test
    public void testScalarDivideZero() {
        Complex x = new Complex(1,1);
        UnitTestUtils.assertEquals(x.divide(Complex.ZERO), x.divide(0), 0);
    }

    @Test
    public void testReciprocal() {
        Complex z = new Complex(5.0, 6.0);
        Complex act = z.reciprocal();
        double expRe = 5.0 / 61.0;
        double expIm = -6.0 / 61.0;
        Assert.assertEquals(expRe, act.getReal(), FastMath.ulp(expRe));
        Assert.assertEquals(expIm, act.getImaginary(), FastMath.ulp(expIm));
    }

    @Test
    public void testReciprocalReal() {
        Complex z = new Complex(-2.0, 0.0);
        Assert.assertTrue(Complex.equals(new Complex(-0.5, 0.0), z.reciprocal()));
    }

    @Test
    public void testReciprocalImaginary() {
        Complex z = new Complex(0.0, -2.0);
        Assert.assertEquals(new Complex(0.0, 0.5), z.reciprocal());
    }

    @Test
    public void testReciprocalInf() {
        Complex z = new Complex(neginf, inf);
        Assert.assertTrue(z.reciprocal().equals(Complex.ZERO));

        z = new Complex(1, inf).reciprocal();
        Assert.assertEquals(z, Complex.ZERO);
    }

    @Test
    public void testReciprocalZero() {
        Assert.assertEquals(Complex.ZERO.reciprocal(), Complex.INF);
    }

    @Test
    public void testReciprocalNaN() {
        Assert.assertTrue(Complex.NaN.reciprocal().isNaN());
    }

    @Test
    public void testMultiply() {
        Complex x = new Complex(3.0, 4.0);
        Complex y = new Complex(5.0, 6.0);
        Complex z = x.multiply(y);
        Assert.assertEquals(-9.0, z.getReal(), 1.0e-5);
        Assert.assertEquals(38.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    public void testMultiplyNaN() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.multiply(Complex.NaN);
        Assert.assertSame(Complex.NaN, z);
        z = Complex.NaN.multiply(5);
        Assert.assertSame(Complex.NaN, z);
    }

    @Test
    public void testMultiplyInfInf() {
        // Assert.assertTrue(infInf.multiply(infInf).isNaN()); // MATH-620
        Assert.assertTrue(infInf.multiply(infInf).isInfinite());
    }

    @Test
    public void testMultiplyNaNInf() {
        Complex z = new Complex(1,1);
        Complex w = z.multiply(infOne);
        Assert.assertEquals(w.getReal(), inf, 0);
        Assert.assertEquals(w.getImaginary(), inf, 0);

        // [MATH-164]
        Assert.assertTrue(new Complex( 1,0).multiply(infInf).equals(Complex.INF));
        Assert.assertTrue(new Complex(-1,0).multiply(infInf).equals(Complex.INF));
        Assert.assertTrue(new Complex( 1,0).multiply(negInfZero).equals(Complex.INF));

        w = oneInf.multiply(oneNegInf);
        Assert.assertEquals(w.getReal(), inf, 0);
        Assert.assertEquals(w.getImaginary(), inf, 0);

        w = negInfNegInf.multiply(oneNaN);
        Assert.assertTrue(Double.isNaN(w.getReal()));
        Assert.assertTrue(Double.isNaN(w.getImaginary()));

        z = new Complex(1, neginf);
        Assert.assertSame(Complex.INF, z.multiply(z));
    }

    @Test
    public void testScalarMultiply() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = 2.0;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.multiply(yComplex), x.multiply(yDouble));
        int zInt = -5;
        Complex zComplex = new Complex(zInt);
        Assert.assertEquals(x.multiply(zComplex), x.multiply(zInt));
    }

    @Test
    public void testScalarMultiplyNaN() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = Double.NaN;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.multiply(yComplex), x.multiply(yDouble));
    }

    @Test
    public void testScalarMultiplyInf() {
        Complex x = new Complex(1, 1);
        double yDouble = Double.POSITIVE_INFINITY;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.multiply(yComplex), x.multiply(yDouble));

        yDouble = Double.NEGATIVE_INFINITY;
        yComplex = new Complex(yDouble);
        Assert.assertEquals(x.multiply(yComplex), x.multiply(yDouble));
    }

    @Test
    public void testNegate() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.negate();
        Assert.assertEquals(-3.0, z.getReal(), 1.0e-5);
        Assert.assertEquals(-4.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    public void testNegateNaN() {
        Complex z = Complex.NaN.negate();
        Assert.assertTrue(z.isNaN());
    }

    @Test
    public void testSubtract() {
        Complex x = new Complex(3.0, 4.0);
        Complex y = new Complex(5.0, 6.0);
        Complex z = x.subtract(y);
        Assert.assertEquals(-2.0, z.getReal(), 1.0e-5);
        Assert.assertEquals(-2.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    public void testSubtractNaN() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.subtract(Complex.NaN);
        Assert.assertSame(Complex.NaN, z);
        z = new Complex(1, nan);
        Complex w = x.subtract(z);
        Assert.assertSame(Complex.NaN, w);
    }

    @Test
    public void testSubtractInf() {
        Complex x = new Complex(1, 1);
        Complex z = new Complex(neginf, 0);
        Complex w = x.subtract(z);
        Assert.assertEquals(w.getImaginary(), 1, 0);
        Assert.assertEquals(inf, w.getReal(), 0);

        x = new Complex(neginf, 0);
        Assert.assertTrue(Double.isNaN(x.subtract(z).getReal()));
    }

    @Test
    public void testScalarSubtract() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = 2.0;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.subtract(yComplex), x.subtract(yDouble));
    }

    @Test
    public void testScalarSubtractNaN() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = Double.NaN;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.subtract(yComplex), x.subtract(yDouble));
    }

    @Test
    public void testScalarSubtractInf() {
        Complex x = new Complex(1, 1);
        double yDouble = Double.POSITIVE_INFINITY;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.subtract(yComplex), x.subtract(yDouble));

        x = new Complex(neginf, 0);
        Assert.assertEquals(x.subtract(yComplex), x.subtract(yDouble));
    }


    @Test
    public void testEqualsNull() {
        Complex x = new Complex(3.0, 4.0);
        Assert.assertFalse(x.equals(null));
    }

    @Test(expected=NullPointerException.class)
    public void testFloatingPointEqualsPrecondition1() {
        Complex.equals(new Complex(3.0, 4.0), null, 3);
    }
    @Test(expected=NullPointerException.class)
    public void testFloatingPointEqualsPrecondition2() {
        Complex.equals(null, new Complex(3.0, 4.0), 3);
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    public void testEqualsClass() {
        Complex x = new Complex(3.0, 4.0);
        Assert.assertFalse(x.equals(this));
    }

    @Test
    public void testEqualsSame() {
        Complex x = new Complex(3.0, 4.0);
        Assert.assertTrue(x.equals(x));
    }

    @Test
    public void testFloatingPointEquals() {
        double re = -3.21;
        double im = 456789e10;

        final Complex x = new Complex(re, im);
        Complex y = new Complex(re, im);

        Assert.assertTrue(x.equals(y));
        Assert.assertTrue(Complex.equals(x, y));

        final int maxUlps = 5;
        for (int i = 0; i < maxUlps; i++) {
            re = FastMath.nextUp(re);
            im = FastMath.nextUp(im);
        }
        y = new Complex(re, im);
        Assert.assertTrue(Complex.equals(x, y, maxUlps));

        re = FastMath.nextUp(re);
        im = FastMath.nextUp(im);
        y = new Complex(re, im);
        Assert.assertFalse(Complex.equals(x, y, maxUlps));
    }

    @Test
    public void testFloatingPointEqualsNaN() {
        Complex c = new Complex(Double.NaN, 1);
        Assert.assertFalse(Complex.equals(c, c));

        c = new Complex(1, Double.NaN);
        Assert.assertFalse(Complex.equals(c, c));
    }

    @Test
    public void testFloatingPointEqualsWithAllowedDelta() {
        final double re = 153.0000;
        final double im = 152.9375;
        final double tol1 = 0.0625;
        final Complex x = new Complex(re, im);
        final Complex y = new Complex(re + tol1, im + tol1);
        Assert.assertTrue(Complex.equals(x, y, tol1));

        final double tol2 = 0.0624;
        Assert.assertFalse(Complex.equals(x, y, tol2));
    }

    @Test
    public void testFloatingPointEqualsWithAllowedDeltaNaN() {
        final Complex x = new Complex(0, Double.NaN);
        final Complex y = new Complex(Double.NaN, 0);
        Assert.assertFalse(Complex.equals(x, Complex.ZERO, 0.1));
        Assert.assertFalse(Complex.equals(x, x, 0.1));
        Assert.assertFalse(Complex.equals(x, y, 0.1));
    }

    @Test
    public void testFloatingPointEqualsWithRelativeTolerance() {
        final double tol = 1e-4;
        final double re = 1;
        final double im = 1e10;

        final double f = 1 + tol;
        final Complex x = new Complex(re, im);
        final Complex y = new Complex(re * f, im * f);
        Assert.assertTrue(Complex.equalsWithRelativeTolerance(x, y, tol));
    }

    @Test
    public void testFloatingPointEqualsWithRelativeToleranceNaN() {
        final Complex x = new Complex(0, Double.NaN);
        final Complex y = new Complex(Double.NaN, 0);
        Assert.assertFalse(Complex.equalsWithRelativeTolerance(x, Complex.ZERO, 0.1));
        Assert.assertFalse(Complex.equalsWithRelativeTolerance(x, x, 0.1));
        Assert.assertFalse(Complex.equalsWithRelativeTolerance(x, y, 0.1));
    }

    @Test
    public void testEqualsTrue() {
        Complex x = new Complex(3.0, 4.0);
        Complex y = new Complex(3.0, 4.0);
        Assert.assertTrue(x.equals(y));
    }

    @Test
    public void testEqualsRealDifference() {
        Complex x = new Complex(0.0, 0.0);
        Complex y = new Complex(0.0 + Double.MIN_VALUE, 0.0);
        Assert.assertFalse(x.equals(y));
    }

    @Test
    public void testEqualsImaginaryDifference() {
        Complex x = new Complex(0.0, 0.0);
        Complex y = new Complex(0.0, 0.0 + Double.MIN_VALUE);
        Assert.assertFalse(x.equals(y));
    }

    @Test
    public void testEqualsNaN() {
        Complex realNaN = new Complex(Double.NaN, 0.0);
        Complex imaginaryNaN = new Complex(0.0, Double.NaN);
        Complex complexNaN = Complex.NaN;
        Assert.assertTrue(realNaN.equals(imaginaryNaN));
        Assert.assertTrue(imaginaryNaN.equals(complexNaN));
        Assert.assertTrue(realNaN.equals(complexNaN));
    }

    @Test
    public void testHashCode() {
        Complex x = new Complex(0.0, 0.0);
        Complex y = new Complex(0.0, 0.0 + Double.MIN_VALUE);
        Assert.assertFalse(x.hashCode()==y.hashCode());
        y = new Complex(0.0 + Double.MIN_VALUE, 0.0);
        Assert.assertFalse(x.hashCode()==y.hashCode());
        Complex realNaN = new Complex(Double.NaN, 0.0);
        Complex imaginaryNaN = new Complex(0.0, Double.NaN);
        Assert.assertEquals(realNaN.hashCode(), imaginaryNaN.hashCode());
        Assert.assertEquals(imaginaryNaN.hashCode(), Complex.NaN.hashCode());

        // MATH-1118
        // "equals" and "hashCode" must be compatible: if two objects have
        // different hash codes, "equals" must return false.
        final String msg = "'equals' not compatible with 'hashCode'";

        x = new Complex(0.0, 0.0);
        y = new Complex(0.0, -0.0);
        Assert.assertTrue(x.hashCode() != y.hashCode());
        Assert.assertFalse(msg, x.equals(y));

        x = new Complex(0.0, 0.0);
        y = new Complex(-0.0, 0.0);
        Assert.assertTrue(x.hashCode() != y.hashCode());
        Assert.assertFalse(msg, x.equals(y));
    }

    @Test
    public void testAcos() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(0.936812, -2.30551);
        UnitTestUtils.assertEquals(expected, z.acos(), 1.0e-5);
        UnitTestUtils.assertEquals(new Complex(FastMath.acos(0), 0),
                Complex.ZERO.acos(), 1.0e-12);
    }

    @Test
    public void testAcosInf() {
        UnitTestUtils.assertSame(Complex.NaN, oneInf.acos());
        UnitTestUtils.assertSame(Complex.NaN, oneNegInf.acos());
        UnitTestUtils.assertSame(Complex.NaN, infOne.acos());
        UnitTestUtils.assertSame(Complex.NaN, negInfOne.acos());
        UnitTestUtils.assertSame(Complex.NaN, infInf.acos());
        UnitTestUtils.assertSame(Complex.NaN, infNegInf.acos());
        UnitTestUtils.assertSame(Complex.NaN, negInfInf.acos());
        UnitTestUtils.assertSame(Complex.NaN, negInfNegInf.acos());
    }

    @Test
    public void testAcosNaN() {
        Assert.assertTrue(Complex.NaN.acos().isNaN());
    }

    @Test
    public void testAsin() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(0.633984, 2.30551);
        UnitTestUtils.assertEquals(expected, z.asin(), 1.0e-5);
    }

    @Test
    public void testAsinNaN() {
        Assert.assertTrue(Complex.NaN.asin().isNaN());
    }

    @Test
    public void testAsinInf() {
        UnitTestUtils.assertSame(Complex.NaN, oneInf.asin());
        UnitTestUtils.assertSame(Complex.NaN, oneNegInf.asin());
        UnitTestUtils.assertSame(Complex.NaN, infOne.asin());
        UnitTestUtils.assertSame(Complex.NaN, negInfOne.asin());
        UnitTestUtils.assertSame(Complex.NaN, infInf.asin());
        UnitTestUtils.assertSame(Complex.NaN, infNegInf.asin());
        UnitTestUtils.assertSame(Complex.NaN, negInfInf.asin());
        UnitTestUtils.assertSame(Complex.NaN, negInfNegInf.asin());
    }


    @Test
    public void testAtan() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(1.44831, 0.158997);
        UnitTestUtils.assertEquals(expected, z.atan(), 1.0e-5);
    }

    @Test
    public void testAtanInf() {
        UnitTestUtils.assertSame(Complex.NaN, oneInf.atan());
        UnitTestUtils.assertSame(Complex.NaN, oneNegInf.atan());
        UnitTestUtils.assertSame(Complex.NaN, infOne.atan());
        UnitTestUtils.assertSame(Complex.NaN, negInfOne.atan());
        UnitTestUtils.assertSame(Complex.NaN, infInf.atan());
        UnitTestUtils.assertSame(Complex.NaN, infNegInf.atan());
        UnitTestUtils.assertSame(Complex.NaN, negInfInf.atan());
        UnitTestUtils.assertSame(Complex.NaN, negInfNegInf.atan());
    }

    @Test
    public void testAtanI() {
        Assert.assertTrue(Complex.I.atan().isNaN());
    }

    @Test
    public void testAtanNaN() {
        Assert.assertTrue(Complex.NaN.atan().isNaN());
    }

    @Test
    public void testCos() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(-27.03495, -3.851153);
        UnitTestUtils.assertEquals(expected, z.cos(), 1.0e-5);
    }

    @Test
    public void testCosNaN() {
        Assert.assertTrue(Complex.NaN.cos().isNaN());
    }

    @Test
    public void testCosInf() {
        UnitTestUtils.assertSame(infNegInf, oneInf.cos());
        UnitTestUtils.assertSame(infInf, oneNegInf.cos());
        UnitTestUtils.assertSame(Complex.NaN, infOne.cos());
        UnitTestUtils.assertSame(Complex.NaN, negInfOne.cos());
        UnitTestUtils.assertSame(Complex.NaN, infInf.cos());
        UnitTestUtils.assertSame(Complex.NaN, infNegInf.cos());
        UnitTestUtils.assertSame(Complex.NaN, negInfInf.cos());
        UnitTestUtils.assertSame(Complex.NaN, negInfNegInf.cos());
    }

    @Test
    public void testCosh() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(-6.58066, -7.58155);
        UnitTestUtils.assertEquals(expected, z.cosh(), 1.0e-5);
    }

    @Test
    public void testCoshNaN() {
        Assert.assertTrue(Complex.NaN.cosh().isNaN());
    }

    @Test
    public void testCoshInf() {
        UnitTestUtils.assertSame(Complex.NaN, oneInf.cosh());
        UnitTestUtils.assertSame(Complex.NaN, oneNegInf.cosh());
        UnitTestUtils.assertSame(infInf, infOne.cosh());
        UnitTestUtils.assertSame(infNegInf, negInfOne.cosh());
        UnitTestUtils.assertSame(Complex.NaN, infInf.cosh());
        UnitTestUtils.assertSame(Complex.NaN, infNegInf.cosh());
        UnitTestUtils.assertSame(Complex.NaN, negInfInf.cosh());
        UnitTestUtils.assertSame(Complex.NaN, negInfNegInf.cosh());
    }

    @Test
    public void testExp() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(-13.12878, -15.20078);
        UnitTestUtils.assertEquals(expected, z.exp(), 1.0e-5);
        UnitTestUtils.assertEquals(Complex.ONE,
                Complex.ZERO.exp(), 10e-12);
        Complex iPi = Complex.I.multiply(new Complex(pi,0));
        UnitTestUtils.assertEquals(Complex.ONE.negate(),
                iPi.exp(), 10e-12);
    }

    @Test
    public void testExpNaN() {
        Assert.assertTrue(Complex.NaN.exp().isNaN());
    }

    @Test
    public void testExpInf1() {
        UnitTestUtils.assertSame(Complex.NaN, oneInf.exp());
    }

    @Test
    public void testExpInf2() {
        UnitTestUtils.assertSame(Complex.NaN, oneNegInf.exp());
    }

    @Test
    public void testExpInf3() {
        UnitTestUtils.assertSame(infInf, infOne.exp());
    }

    @Test
    public void testExpInf4() {
        final Complex exp = negInfOne.exp();
        UnitTestUtils.assertSame(Complex.ZERO, exp);
    }

    @Test
    public void testExpInf5() {
        UnitTestUtils.assertSame(Complex.NaN, infInf.exp());
    }

    @Test
    public void testExpInf6() {
        UnitTestUtils.assertSame(Complex.NaN, infNegInf.exp());
    }

    @Test
    public void testExpInf7() {
        UnitTestUtils.assertSame(Complex.NaN, negInfInf.exp());
    }

    @Test
    public void testExpInf8() {
        UnitTestUtils.assertSame(Complex.NaN, negInfNegInf.exp());
    }

    @Test
    public void testLog() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(1.60944, 0.927295);
        UnitTestUtils.assertEquals(expected, z.log(), 1.0e-5);
    }

    @Test
    public void testLogNaN() {
        Assert.assertTrue(Complex.NaN.log().isNaN());
    }

    @Test
    public void testLogInf() {
        UnitTestUtils.assertEquals(new Complex(inf, pi / 2),
                oneInf.log(), 10e-12);
        UnitTestUtils.assertEquals(new Complex(inf, -pi / 2),
                oneNegInf.log(), 10e-12);
        UnitTestUtils.assertEquals(infZero, infOne.log(), 10e-12);
        UnitTestUtils.assertEquals(new Complex(inf, pi),
                negInfOne.log(), 10e-12);
        UnitTestUtils.assertEquals(new Complex(inf, pi / 4),
                infInf.log(), 10e-12);
        UnitTestUtils.assertEquals(new Complex(inf, -pi / 4),
                infNegInf.log(), 10e-12);
        UnitTestUtils.assertEquals(new Complex(inf, 3d * pi / 4),
                negInfInf.log(), 10e-12);
        UnitTestUtils.assertEquals(new Complex(inf, - 3d * pi / 4),
                negInfNegInf.log(), 10e-12);
    }

    @Test
    public void testLogZero() {
        UnitTestUtils.assertSame(negInfZero, Complex.ZERO.log());
    }

    @Test
    public void testPow() {
        Complex x = new Complex(3, 4);
        Complex y = new Complex(5, 6);
        Complex expected = new Complex(-1.860893, 11.83677);
        UnitTestUtils.assertEquals(expected, x.pow(y), 1.0e-5);
    }

    @Test
    public void testPowNaNBase() {
        Complex x = new Complex(3, 4);
        Assert.assertTrue(Complex.NaN.pow(x).isNaN());
    }

    @Test
    public void testPowNaNExponent() {
        Complex x = new Complex(3, 4);
        Assert.assertTrue(x.pow(Complex.NaN).isNaN());
    }

   @Test
   public void testPowInf() {
       UnitTestUtils.assertSame(Complex.NaN,Complex.ONE.pow(oneInf));
       UnitTestUtils.assertSame(Complex.NaN,Complex.ONE.pow(oneNegInf));
       UnitTestUtils.assertSame(Complex.NaN,Complex.ONE.pow(infOne));
       UnitTestUtils.assertSame(Complex.NaN,Complex.ONE.pow(infInf));
       UnitTestUtils.assertSame(Complex.NaN,Complex.ONE.pow(infNegInf));
       UnitTestUtils.assertSame(Complex.NaN,Complex.ONE.pow(negInfInf));
       UnitTestUtils.assertSame(Complex.NaN,Complex.ONE.pow(negInfNegInf));
       UnitTestUtils.assertSame(Complex.NaN,infOne.pow(Complex.ONE));
       UnitTestUtils.assertSame(Complex.NaN,negInfOne.pow(Complex.ONE));
       UnitTestUtils.assertSame(Complex.NaN,infInf.pow(Complex.ONE));
       UnitTestUtils.assertSame(Complex.NaN,infNegInf.pow(Complex.ONE));
       UnitTestUtils.assertSame(Complex.NaN,negInfInf.pow(Complex.ONE));
       UnitTestUtils.assertSame(Complex.NaN,negInfNegInf.pow(Complex.ONE));
       UnitTestUtils.assertSame(Complex.NaN,negInfNegInf.pow(infNegInf));
       UnitTestUtils.assertSame(Complex.NaN,negInfNegInf.pow(negInfNegInf));
       UnitTestUtils.assertSame(Complex.NaN,negInfNegInf.pow(infInf));
       UnitTestUtils.assertSame(Complex.NaN,infInf.pow(infNegInf));
       UnitTestUtils.assertSame(Complex.NaN,infInf.pow(negInfNegInf));
       UnitTestUtils.assertSame(Complex.NaN,infInf.pow(infInf));
       UnitTestUtils.assertSame(Complex.NaN,infNegInf.pow(infNegInf));
       UnitTestUtils.assertSame(Complex.NaN,infNegInf.pow(negInfNegInf));
       UnitTestUtils.assertSame(Complex.NaN,infNegInf.pow(infInf));
   }

   @Test
   public void testPowZero() {
       UnitTestUtils.assertSame(Complex.NaN,
               Complex.ZERO.pow(Complex.ONE));
       UnitTestUtils.assertSame(Complex.NaN,
               Complex.ZERO.pow(Complex.ZERO));
       UnitTestUtils.assertSame(Complex.NaN,
               Complex.ZERO.pow(Complex.I));
       UnitTestUtils.assertEquals(Complex.ONE,
               Complex.ONE.pow(Complex.ZERO), 10e-12);
       UnitTestUtils.assertEquals(Complex.ONE,
               Complex.I.pow(Complex.ZERO), 10e-12);
       UnitTestUtils.assertEquals(Complex.ONE,
               new Complex(-1, 3).pow(Complex.ZERO), 10e-12);
   }

    @Test
    public void testScalarPow() {
        Complex x = new Complex(3, 4);
        double yDouble = 5.0;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.pow(yComplex), x.pow(yDouble));
    }

    @Test
    public void testScalarPowNaNBase() {
        Complex x = Complex.NaN;
        double yDouble = 5.0;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.pow(yComplex), x.pow(yDouble));
    }

    @Test
    public void testScalarPowNaNExponent() {
        Complex x = new Complex(3, 4);
        double yDouble = Double.NaN;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.pow(yComplex), x.pow(yDouble));
    }

   @Test
   public void testScalarPowInf() {
       UnitTestUtils.assertSame(Complex.NaN,Complex.ONE.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.assertSame(Complex.NaN,Complex.ONE.pow(Double.NEGATIVE_INFINITY));
       UnitTestUtils.assertSame(Complex.NaN,infOne.pow(1.0));
       UnitTestUtils.assertSame(Complex.NaN,negInfOne.pow(1.0));
       UnitTestUtils.assertSame(Complex.NaN,infInf.pow(1.0));
       UnitTestUtils.assertSame(Complex.NaN,infNegInf.pow(1.0));
       UnitTestUtils.assertSame(Complex.NaN,negInfInf.pow(10));
       UnitTestUtils.assertSame(Complex.NaN,negInfNegInf.pow(1.0));
       UnitTestUtils.assertSame(Complex.NaN,negInfNegInf.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.assertSame(Complex.NaN,negInfNegInf.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.assertSame(Complex.NaN,infInf.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.assertSame(Complex.NaN,infInf.pow(Double.NEGATIVE_INFINITY));
       UnitTestUtils.assertSame(Complex.NaN,infNegInf.pow(Double.NEGATIVE_INFINITY));
       UnitTestUtils.assertSame(Complex.NaN,infNegInf.pow(Double.POSITIVE_INFINITY));
   }

   @Test
   public void testScalarPowZero() {
       UnitTestUtils.assertSame(Complex.NaN, Complex.ZERO.pow(1.0));
       UnitTestUtils.assertSame(Complex.NaN, Complex.ZERO.pow(0.0));
       UnitTestUtils.assertEquals(Complex.ONE, Complex.ONE.pow(0.0), 10e-12);
       UnitTestUtils.assertEquals(Complex.ONE, Complex.I.pow(0.0), 10e-12);
       UnitTestUtils.assertEquals(Complex.ONE, new Complex(-1, 3).pow(0.0), 10e-12);
   }

    @Test(expected=NullArgumentException.class)
    public void testpowNull() {
        Complex.ONE.pow(null);
    }

    @Test
    public void testSin() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(3.853738, -27.01681);
        UnitTestUtils.assertEquals(expected, z.sin(), 1.0e-5);
    }

    @Test
    public void testSinInf() {
        UnitTestUtils.assertSame(infInf, oneInf.sin());
        UnitTestUtils.assertSame(infNegInf, oneNegInf.sin());
        UnitTestUtils.assertSame(Complex.NaN, infOne.sin());
        UnitTestUtils.assertSame(Complex.NaN, negInfOne.sin());
        UnitTestUtils.assertSame(Complex.NaN, infInf.sin());
        UnitTestUtils.assertSame(Complex.NaN, infNegInf.sin());
        UnitTestUtils.assertSame(Complex.NaN, negInfInf.sin());
        UnitTestUtils.assertSame(Complex.NaN, negInfNegInf.sin());
    }

    @Test
    public void testSinNaN() {
        Assert.assertTrue(Complex.NaN.sin().isNaN());
    }

    @Test
    public void testSinh() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(-6.54812, -7.61923);
        UnitTestUtils.assertEquals(expected, z.sinh(), 1.0e-5);
    }

    @Test
    public void testSinhNaN() {
        Assert.assertTrue(Complex.NaN.sinh().isNaN());
    }

    @Test
    public void testSinhInf() {
        UnitTestUtils.assertSame(Complex.NaN, oneInf.sinh());
        UnitTestUtils.assertSame(Complex.NaN, oneNegInf.sinh());
        UnitTestUtils.assertSame(infInf, infOne.sinh());
        UnitTestUtils.assertSame(negInfInf, negInfOne.sinh());
        UnitTestUtils.assertSame(Complex.NaN, infInf.sinh());
        UnitTestUtils.assertSame(Complex.NaN, infNegInf.sinh());
        UnitTestUtils.assertSame(Complex.NaN, negInfInf.sinh());
        UnitTestUtils.assertSame(Complex.NaN, negInfNegInf.sinh());
    }

    @Test
    public void testSqrtRealPositive() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(2, 1);
        UnitTestUtils.assertEquals(expected, z.sqrt(), 1.0e-5);
    }

    @Test
    public void testSqrtRealZero() {
        Complex z = new Complex(0.0, 4);
        Complex expected = new Complex(1.41421, 1.41421);
        UnitTestUtils.assertEquals(expected, z.sqrt(), 1.0e-5);
    }

    @Test
    public void testSqrtRealNegative() {
        Complex z = new Complex(-3.0, 4);
        Complex expected = new Complex(1, 2);
        UnitTestUtils.assertEquals(expected, z.sqrt(), 1.0e-5);
    }

    @Test
    public void testSqrtImaginaryZero() {
        Complex z = new Complex(-3.0, 0.0);
        Complex expected = new Complex(0.0, 1.73205);
        UnitTestUtils.assertEquals(expected, z.sqrt(), 1.0e-5);
    }

    @Test
    public void testSqrtImaginaryNegative() {
        Complex z = new Complex(-3.0, -4.0);
        Complex expected = new Complex(1.0, -2.0);
        UnitTestUtils.assertEquals(expected, z.sqrt(), 1.0e-5);
    }

    @Test
    public void testSqrtPolar() {
        double r = 1;
        for (int i = 0; i < 5; i++) {
            r += i;
            double theta = 0;
            for (int j =0; j < 11; j++) {
                theta += pi /12;
                Complex z = ComplexUtils.polar2Complex(r, theta);
                Complex sqrtz = ComplexUtils.polar2Complex(FastMath.sqrt(r), theta / 2);
                UnitTestUtils.assertEquals(sqrtz, z.sqrt(), 10e-12);
            }
        }
    }

    @Test
    public void testSqrtNaN() {
        Assert.assertTrue(Complex.NaN.sqrt().isNaN());
    }

    @Test
    public void testSqrtInf() {
        UnitTestUtils.assertSame(infNaN, oneInf.sqrt());
        UnitTestUtils.assertSame(infNaN, oneNegInf.sqrt());
        UnitTestUtils.assertSame(infZero, infOne.sqrt());
        UnitTestUtils.assertSame(zeroInf, negInfOne.sqrt());
        UnitTestUtils.assertSame(infNaN, infInf.sqrt());
        UnitTestUtils.assertSame(infNaN, infNegInf.sqrt());
        UnitTestUtils.assertSame(nanInf, negInfInf.sqrt());
        UnitTestUtils.assertSame(nanNegInf, negInfNegInf.sqrt());
    }

    @Test
    public void testSqrt1z() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(4.08033, -2.94094);
        UnitTestUtils.assertEquals(expected, z.sqrt1z(), 1.0e-5);
    }

    @Test
    public void testSqrt1zNaN() {
        Assert.assertTrue(Complex.NaN.sqrt1z().isNaN());
    }

    @Test
    public void testTan() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(-0.000187346, 0.999356);
        UnitTestUtils.assertEquals(expected, z.tan(), 1.0e-5);
        /* Check that no overflow occurs (MATH-722) */
        Complex actual = new Complex(3.0, 1E10).tan();
        expected = new Complex(0, 1);
        UnitTestUtils.assertEquals(expected, actual, 1.0e-5);
        actual = new Complex(3.0, -1E10).tan();
        expected = new Complex(0, -1);
        UnitTestUtils.assertEquals(expected, actual, 1.0e-5);
    }

    @Test
    public void testTanNaN() {
        Assert.assertTrue(Complex.NaN.tan().isNaN());
    }

    @Test
    public void testTanInf() {
        UnitTestUtils.assertSame(Complex.valueOf(0.0, 1.0), oneInf.tan());
        UnitTestUtils.assertSame(Complex.valueOf(0.0, -1.0), oneNegInf.tan());
        UnitTestUtils.assertSame(Complex.NaN, infOne.tan());
        UnitTestUtils.assertSame(Complex.NaN, negInfOne.tan());
        UnitTestUtils.assertSame(Complex.NaN, infInf.tan());
        UnitTestUtils.assertSame(Complex.NaN, infNegInf.tan());
        UnitTestUtils.assertSame(Complex.NaN, negInfInf.tan());
        UnitTestUtils.assertSame(Complex.NaN, negInfNegInf.tan());
    }

   @Test
   public void testTanCritical() {
        UnitTestUtils.assertSame(infNaN, new Complex(pi/2, 0).tan());
        UnitTestUtils.assertSame(negInfNaN, new Complex(-pi/2, 0).tan());
    }

    @Test
    public void testTanh() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(1.00071, 0.00490826);
        UnitTestUtils.assertEquals(expected, z.tanh(), 1.0e-5);
        /* Check that no overflow occurs (MATH-722) */
        Complex actual = new Complex(1E10, 3.0).tanh();
        expected = new Complex(1, 0);
        UnitTestUtils.assertEquals(expected, actual, 1.0e-5);
        actual = new Complex(-1E10, 3.0).tanh();
        expected = new Complex(-1, 0);
        UnitTestUtils.assertEquals(expected, actual, 1.0e-5);
    }

    @Test
    public void testTanhNaN() {
        Assert.assertTrue(Complex.NaN.tanh().isNaN());
    }

    @Test
    public void testTanhInf() {
        UnitTestUtils.assertSame(Complex.NaN, oneInf.tanh());
        UnitTestUtils.assertSame(Complex.NaN, oneNegInf.tanh());
        UnitTestUtils.assertSame(Complex.valueOf(1.0, 0.0), infOne.tanh());
        UnitTestUtils.assertSame(Complex.valueOf(-1.0, 0.0), negInfOne.tanh());
        UnitTestUtils.assertSame(Complex.NaN, infInf.tanh());
        UnitTestUtils.assertSame(Complex.NaN, infNegInf.tanh());
        UnitTestUtils.assertSame(Complex.NaN, negInfInf.tanh());
        UnitTestUtils.assertSame(Complex.NaN, negInfNegInf.tanh());
    }

    @Test
    public void testTanhCritical() {
        UnitTestUtils.assertSame(nanInf, new Complex(0, pi/2).tanh());
    }

    /** test issue MATH-221 */
    @Test
    public void testMath221() {
        Assert.assertTrue(Complex.equals(new Complex(0,-1),
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
    public void testNthRoot_normal_thirdRoot() {
        // The complex number we want to compute all third-roots for.
        Complex z = new Complex(-2,2);
        // The List holding all third roots
        Complex[] thirdRootsOfZ = z.nthRoot(3).toArray(new Complex[0]);
        // Returned Collection must not be empty!
        Assert.assertEquals(3, thirdRootsOfZ.length);
        // test z_0
        Assert.assertEquals(1.0,                  thirdRootsOfZ[0].getReal(),      1.0e-5);
        Assert.assertEquals(1.0,                  thirdRootsOfZ[0].getImaginary(), 1.0e-5);
        // test z_1
        Assert.assertEquals(-1.3660254037844386,  thirdRootsOfZ[1].getReal(),      1.0e-5);
        Assert.assertEquals(0.36602540378443843,  thirdRootsOfZ[1].getImaginary(), 1.0e-5);
        // test z_2
        Assert.assertEquals(0.366025403784439,    thirdRootsOfZ[2].getReal(),      1.0e-5);
        Assert.assertEquals(-1.3660254037844384,  thirdRootsOfZ[2].getImaginary(), 1.0e-5);
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
    public void testNthRoot_normal_fourthRoot() {
        // The complex number we want to compute all third-roots for.
        Complex z = new Complex(5,-2);
        // The List holding all fourth roots
        Complex[] fourthRootsOfZ = z.nthRoot(4).toArray(new Complex[0]);
        // Returned Collection must not be empty!
        Assert.assertEquals(4, fourthRootsOfZ.length);
        // test z_0
        Assert.assertEquals(1.5164629308487783,     fourthRootsOfZ[0].getReal(),      1.0e-5);
        Assert.assertEquals(-0.14469266210702247,   fourthRootsOfZ[0].getImaginary(), 1.0e-5);
        // test z_1
        Assert.assertEquals(0.14469266210702256,    fourthRootsOfZ[1].getReal(),      1.0e-5);
        Assert.assertEquals(1.5164629308487783,     fourthRootsOfZ[1].getImaginary(), 1.0e-5);
        // test z_2
        Assert.assertEquals(-1.5164629308487783,    fourthRootsOfZ[2].getReal(),      1.0e-5);
        Assert.assertEquals(0.14469266210702267,    fourthRootsOfZ[2].getImaginary(), 1.0e-5);
        // test z_3
        Assert.assertEquals(-0.14469266210702275,   fourthRootsOfZ[3].getReal(),      1.0e-5);
        Assert.assertEquals(-1.5164629308487783,    fourthRootsOfZ[3].getImaginary(), 1.0e-5);
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
    public void testNthRoot_cornercase_thirdRoot_imaginaryPartEmpty() {
        // The number 8 has three third roots. One we all already know is the number 2.
        // But there are two more complex roots.
        Complex z = new Complex(8,0);
        // The List holding all third roots
        Complex[] thirdRootsOfZ = z.nthRoot(3).toArray(new Complex[0]);
        // Returned Collection must not be empty!
        Assert.assertEquals(3, thirdRootsOfZ.length);
        // test z_0
        Assert.assertEquals(2.0,                thirdRootsOfZ[0].getReal(),      1.0e-5);
        Assert.assertEquals(0.0,                thirdRootsOfZ[0].getImaginary(), 1.0e-5);
        // test z_1
        Assert.assertEquals(-1.0,               thirdRootsOfZ[1].getReal(),      1.0e-5);
        Assert.assertEquals(1.7320508075688774, thirdRootsOfZ[1].getImaginary(), 1.0e-5);
        // test z_2
        Assert.assertEquals(-1.0,               thirdRootsOfZ[2].getReal(),      1.0e-5);
        Assert.assertEquals(-1.732050807568877, thirdRootsOfZ[2].getImaginary(), 1.0e-5);
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
    public void testNthRoot_cornercase_thirdRoot_realPartZero() {
        // complex number with only imaginary part
        Complex z = new Complex(0,2);
        // The List holding all third roots
        Complex[] thirdRootsOfZ = z.nthRoot(3).toArray(new Complex[0]);
        // Returned Collection must not be empty!
        Assert.assertEquals(3, thirdRootsOfZ.length);
        // test z_0
        Assert.assertEquals(1.0911236359717216,      thirdRootsOfZ[0].getReal(),      1.0e-5);
        Assert.assertEquals(0.6299605249474365,      thirdRootsOfZ[0].getImaginary(), 1.0e-5);
        // test z_1
        Assert.assertEquals(-1.0911236359717216,     thirdRootsOfZ[1].getReal(),      1.0e-5);
        Assert.assertEquals(0.6299605249474365,      thirdRootsOfZ[1].getImaginary(), 1.0e-5);
        // test z_2
        Assert.assertEquals(-2.3144374213981936E-16, thirdRootsOfZ[2].getReal(),      1.0e-5);
        Assert.assertEquals(-1.2599210498948732,     thirdRootsOfZ[2].getImaginary(), 1.0e-5);
    }

    /**
     * Test cornercases with NaN and Infinity.
     */
    @Test
    public void testNthRoot_cornercase_NAN_Inf() {
        // NaN + finite -> NaN
        List<Complex> roots = oneNaN.nthRoot(3);
        Assert.assertEquals(1,roots.size());
        Assert.assertEquals(Complex.NaN, roots.get(0));

        roots = nanZero.nthRoot(3);
        Assert.assertEquals(1,roots.size());
        Assert.assertEquals(Complex.NaN, roots.get(0));

        // NaN + infinite -> NaN
        roots = nanInf.nthRoot(3);
        Assert.assertEquals(1,roots.size());
        Assert.assertEquals(Complex.NaN, roots.get(0));

        // finite + infinite -> Inf
        roots = oneInf.nthRoot(3);
        Assert.assertEquals(1,roots.size());
        Assert.assertEquals(Complex.INF, roots.get(0));

        // infinite + infinite -> Inf
        roots = negInfInf.nthRoot(3);
        Assert.assertEquals(1,roots.size());
        Assert.assertEquals(Complex.INF, roots.get(0));
    }

    /**
     * Test standard values
     */
    @Test
    public void testGetArgument() {
        Complex z = new Complex(1, 0);
        Assert.assertEquals(0.0, z.getArgument(), 1.0e-12);

        z = new Complex(1, 1);
        Assert.assertEquals(FastMath.PI/4, z.getArgument(), 1.0e-12);

        z = new Complex(0, 1);
        Assert.assertEquals(FastMath.PI/2, z.getArgument(), 1.0e-12);

        z = new Complex(-1, 1);
        Assert.assertEquals(3 * FastMath.PI/4, z.getArgument(), 1.0e-12);

        z = new Complex(-1, 0);
        Assert.assertEquals(FastMath.PI, z.getArgument(), 1.0e-12);

        z = new Complex(-1, -1);
        Assert.assertEquals(-3 * FastMath.PI/4, z.getArgument(), 1.0e-12);

        z = new Complex(0, -1);
        Assert.assertEquals(-FastMath.PI/2, z.getArgument(), 1.0e-12);

        z = new Complex(1, -1);
        Assert.assertEquals(-FastMath.PI/4, z.getArgument(), 1.0e-12);

    }

    /**
     * Verify atan2-style handling of infinite parts
     */
    @Test
    public void testGetArgumentInf() {
        Assert.assertEquals(FastMath.PI/4, infInf.getArgument(), 1.0e-12);
        Assert.assertEquals(FastMath.PI/2, oneInf.getArgument(), 1.0e-12);
        Assert.assertEquals(0.0, infOne.getArgument(), 1.0e-12);
        Assert.assertEquals(FastMath.PI/2, zeroInf.getArgument(), 1.0e-12);
        Assert.assertEquals(0.0, infZero.getArgument(), 1.0e-12);
        Assert.assertEquals(FastMath.PI, negInfOne.getArgument(), 1.0e-12);
        Assert.assertEquals(-3.0*FastMath.PI/4, negInfNegInf.getArgument(), 1.0e-12);
        Assert.assertEquals(-FastMath.PI/2, oneNegInf.getArgument(), 1.0e-12);
    }

    /**
     * Verify that either part NaN results in NaN
     */
    @Test
    public void testGetArgumentNaN() {
        Assert.assertTrue(Double.isNaN(nanZero.getArgument()));
        Assert.assertTrue(Double.isNaN(zeroNaN.getArgument()));
        Assert.assertTrue(Double.isNaN(Complex.NaN.getArgument()));
    }

    @Test
    public void testSerial() {
        Complex z = new Complex(3.0, 4.0);
        Assert.assertEquals(z, UnitTestUtils.serializeAndRecover(z));
        Complex ncmplx = (Complex)UnitTestUtils.serializeAndRecover(oneNaN);
        Assert.assertEquals(nanZero, ncmplx);
        Assert.assertTrue(ncmplx.isNaN());
        Complex infcmplx = (Complex)UnitTestUtils.serializeAndRecover(infInf);
        Assert.assertEquals(infInf, infcmplx);
        Assert.assertTrue(infcmplx.isInfinite());
        TestComplex tz = new TestComplex(3.0, 4.0);
        Assert.assertEquals(tz, UnitTestUtils.serializeAndRecover(tz));
        TestComplex ntcmplx = (TestComplex)UnitTestUtils.serializeAndRecover(new TestComplex(oneNaN));
        Assert.assertEquals(nanZero, ntcmplx);
        Assert.assertTrue(ntcmplx.isNaN());
        TestComplex inftcmplx = (TestComplex)UnitTestUtils.serializeAndRecover(new TestComplex(infInf));
        Assert.assertEquals(infInf, inftcmplx);
        Assert.assertTrue(inftcmplx.isInfinite());
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

        public TestComplex(Complex other){
            this(other.getReal(), other.getImaginary());
        }

        @Override
        protected TestComplex createComplex(double real, double imaginary){
            return new TestComplex(real, imaginary);
        }

    }
}
