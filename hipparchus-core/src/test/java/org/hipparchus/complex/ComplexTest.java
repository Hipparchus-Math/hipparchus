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

import java.util.Arrays;
import java.util.List;

import org.hipparchus.CalculusFieldElementAbstractTest;
import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Precision;
import org.junit.Assert;
import org.junit.Test;


/**
 */
public class ComplexTest extends CalculusFieldElementAbstractTest<Complex> {


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

    @Override
    protected Complex build(final double x) {
        return new Complex(x, 0.0);
    }

    @Test
    public void testConstructor() {
        Complex z = new Complex(3.0, 4.0);
        Assert.assertEquals(3.0, z.getReal(), 1.0e-5);
        Assert.assertEquals(4.0, z.getImaginary(), 1.0e-5);
        Assert.assertEquals(3.0, z.getRealPart(), 1.0e-5);
        Assert.assertEquals(4.0, z.getImaginaryPart(), 1.0e-5);
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
    public void testNorm() {
        Complex z = new Complex(3.0, 4.0);
        Assert.assertEquals(5.0, z.norm(), 1.0e-5);
    }

    @Test
    public void testNormNaN() {
        Assert.assertTrue(Double.isNaN(Complex.NaN.norm()));
        Complex z = new Complex(inf, nan);
        Assert.assertTrue(Double.isNaN(z.norm()));
    }

    @Test
    public void testNormInfinite() {
        Complex z = Complex.NaN.newInstance(inf);
        Assert.assertEquals(inf, z.norm(), 0);
        z = new Complex(0, neginf);
        Assert.assertEquals(inf, z.norm(), 0);
        z = new Complex(inf, neginf);
        Assert.assertEquals(inf, z.norm(), 0);
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
        Assert.assertSame(Complex.NaN, Complex.NaN.add(Double.NaN));
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
        Assert.assertTrue(Complex.NaN.divide(Complex.NaN).isNaN());
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
        Assert.assertTrue(Complex.NaN.divide(Double.NaN).isNaN());
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
        Assert.assertTrue(new Complex(Double.NaN, 0).multiply(5).isNaN());
        Assert.assertTrue(new Complex(0, Double.NaN).multiply(5).isNaN());
        Assert.assertTrue(Complex.NaN.multiply(5).isNaN());
        Assert.assertTrue(new Complex(Double.NaN, 0).multiply(5.0).isNaN());
        Assert.assertTrue(new Complex(0, Double.NaN).multiply(5.0).isNaN());
        Assert.assertTrue(Complex.NaN.multiply(5.0).isNaN());
        Assert.assertTrue(Complex.ONE.multiply(Double.NaN).isNaN());
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

        Assert.assertTrue(new Complex(Double.POSITIVE_INFINITY, 0).multiply(5).isInfinite());
        Assert.assertTrue(new Complex(0, Double.POSITIVE_INFINITY).multiply(5).isInfinite());
        Assert.assertTrue(Complex.INF.multiply(5).isInfinite());
        Assert.assertTrue(new Complex(Double.POSITIVE_INFINITY, 0).multiply(5.0).isInfinite());
        Assert.assertTrue(new Complex(0, Double.POSITIVE_INFINITY).multiply(5.0).isInfinite());
        Assert.assertTrue(Complex.INF.multiply(5.0).isInfinite());
        Assert.assertTrue(Complex.ONE.multiply(Double.POSITIVE_INFINITY).isInfinite());
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
        Assert.assertSame(Complex.NaN, Complex.NaN.subtract(Complex.NaN));
        Assert.assertSame(Complex.NaN, Complex.NaN.subtract(Double.NaN));
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
    public void testToDegreesComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(FastMath.toDegrees(z.getReal()), FastMath.toDegrees(z.getImaginary()));
        UnitTestUtils.assertEquals(expected, z.toDegrees(), 1.0e-15);
    }

    @Test
    public void testToRadiansComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(FastMath.toRadians(z.getReal()), FastMath.toRadians(z.getImaginary()));
        UnitTestUtils.assertEquals(expected, z.toRadians(), 1.0e-15);
    }

    @Test
    public void testAcosComplex() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(0.936812, -2.30551);
        UnitTestUtils.assertEquals(expected, z.acos(), 1.0e-5);
        UnitTestUtils.assertEquals(new Complex(FastMath.acos(0), 0),
                Complex.ZERO.acos(), 1.0e-12);
    }

    @Test
    public void testAcosNaN() {
        Assert.assertTrue(Complex.NaN.acos().isNaN());
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
    public void testAcosBranchCuts() {
        UnitTestUtils.assertEquals(new Complex(3.141592653589793238462, -0.76103968373182660633),
                                   FastMath.acos(new Complex(-1.3038404810405297, +0.0)),
                                   1.0e-14);
        UnitTestUtils.assertEquals(new Complex(3.141592653589793238462, +0.76103968373182660633),
                                   FastMath.acos(new Complex(-1.3038404810405297, -0.0)),
                                   1.0e-14);
        UnitTestUtils.assertEquals(new Complex(0.0, -0.76103968373182660633),
                                   FastMath.acos(new Complex(1.3038404810405297, +0.0)),
                                   1.0e-14);
        UnitTestUtils.assertEquals(new Complex(0.0, +0.76103968373182660633),
                                   FastMath.acos(new Complex(1.3038404810405297, -0.0)),
                                   1.0e-14);
    }

    @Test
    public void testAsinComplex() {
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
    public void testAsinBranchCuts() {
        UnitTestUtils.assertEquals(new Complex(-1.57079632679489661923, +0.76103968373182660633),
                                   FastMath.asin(new Complex(-1.3038404810405297, +0.0)),
                                   1.0e-14);
        UnitTestUtils.assertEquals(new Complex(-1.57079632679489661923, -0.76103968373182660633),
                                   FastMath.asin(new Complex(-1.3038404810405297, -0.0)),
                                   1.0e-14);
        UnitTestUtils.assertEquals(new Complex(1.57079632679489661923, +0.76103968373182660633),
                                   FastMath.asin(new Complex(1.3038404810405297, +0.0)),
                                   1.0e-14);
        UnitTestUtils.assertEquals(new Complex(1.57079632679489661923, -0.76103968373182660633),
                                   FastMath.asin(new Complex(1.3038404810405297, -0.0)),
                                   1.0e-14);
    }

    @Test
    public void testAtanComplex() {
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
    public void testAtanBranchCuts() {
        UnitTestUtils.assertEquals(new Complex(+1.5707963267948966192, +1.0986122886681096913),
                                   FastMath.atan(new Complex(+0.0, 1.25)),
                                   1.0e-14);
        UnitTestUtils.assertEquals(new Complex(-1.5707963267948966192, +1.0986122886681096913),
                                   FastMath.atan(new Complex(-0.0, 1.25)),
                                   1.0e-14);
        UnitTestUtils.assertEquals(new Complex(+1.5707963267948966192, -1.0986122886681096913),
                                   FastMath.atan(new Complex(+0.0, -1.25)),
                                   1.0e-14);
        UnitTestUtils.assertEquals(new Complex(-1.5707963267948966192, -1.0986122886681096913),
                                   FastMath.atan(new Complex(-0.0, -1.25)),
                                   1.0e-14);
        UnitTestUtils.assertEquals(new Complex(0.0, +0.25541281188299534160),
                                   FastMath.atan(new Complex(+0.0, 0.25)),
                                   1.0e-14);
        Assert.assertTrue(FastMath.copySign(1.0, FastMath.atan(new Complex(+0.0, 0.25)).getReal()) > 0.0);
        UnitTestUtils.assertEquals(new Complex(0.0, +0.25541281188299534160),
                                   FastMath.atan(new Complex(-0.0, 0.25)),
                                   1.0e-14);
        Assert.assertTrue(FastMath.copySign(1.0, FastMath.atan(new Complex(-0.0, 0.25)).getReal()) < 0.0);
        UnitTestUtils.assertEquals(new Complex(0.0, -0.25541281188299534160),
                                   FastMath.atan(new Complex(+0.0, -0.25)),
                                   1.0e-14);
        Assert.assertTrue(FastMath.copySign(1.0, FastMath.atan(new Complex(+0.0, -0.25)).getReal()) > 0.0);
        UnitTestUtils.assertEquals(new Complex(0.0, -0.25541281188299534160),
                                   FastMath.atan(new Complex(-0.0, -0.25)),
                                   1.0e-14);
        Assert.assertTrue(FastMath.copySign(1.0, FastMath.atan(new Complex(-0.0, -0.25)).getReal()) < 0.0);
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
    public void testAtan2Complex() {
        for (double r1 : Arrays.asList(-3, 3)) {
            for (double i1 : Arrays.asList(-2, 0, 2)) {
                final Complex c1 = new Complex(r1, i1);
                for (double r2 : Arrays.asList(-1, 1)) {
                    for (double i2 : Arrays.asList(-5, 0, 5)) {
                        final Complex c2 = new Complex(r2, i2);
                        UnitTestUtils.assertEquals(c1.divide(c2), c1.atan2(c2).tan(), 1.0e-14);
                        final Complex atan   = c1.divide(c2).atan();
                        final Complex atan2  = c1.atan2(c2);
                        final double  deltaR = FastMath.abs(atan.getReal() - atan2.getReal()) / FastMath.PI;
                        Assert.assertTrue(FastMath.abs(deltaR - FastMath.rint(deltaR)) < 1.0e-14);
                        Assert.assertEquals(atan.getImaginary(), atan2.getImaginary(), 1.0e-14);
                    }
                }
            }
        }
    }

    @Test
    public void testAtan2Real() {
        for (double r1 : Arrays.asList(-3, 3)) {
            final Complex c1 = new Complex(r1, 0);
            for (double r2 : Arrays.asList(-1, 1)) {
                final Complex c2 = new Complex(r2, 0);
                Assert.assertEquals(FastMath.atan2(r1, r2),
                                    MathUtils.normalizeAngle(c1.atan2(c2).getReal(), 0.0),
                                    1.0e-14);
            }
        }
    }

    @Override
    @Test
    public void testAtan2SpecialCases() {
        Assert.assertTrue(build(+0.0).atan2(build(+0.0)).isNaN());
        Assert.assertTrue(build(-0.0).atan2(build(+0.0)).isNaN());
        Assert.assertTrue(build(+0.0).atan2(build(-0.0)).isNaN());
        Assert.assertTrue(build(-0.0).atan2(build(-0.0)).isNaN());
    }

    @Test
    public void testCosComplex() {
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
    public void testCoshComplex() {
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
    public void testExpComplex() {
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
    public void testExpM1() {
        UnitTestUtils.assertSame(Complex.NaN, Complex.NaN.expm1());
        final double testValue = FastMath.scalb(1.0, -30);
        Assert.assertEquals(FastMath.expm1(testValue), new Complex(testValue, 0).expm1().getReal(), 1.0e-30);
        Assert.assertTrue(FastMath.expm1(testValue) - new Complex(testValue).exp().subtract(1.0).getReal() > 4.0e-19);
        Assert.assertEquals(0.0, new Complex(0, testValue).expm1().getReal(), 1.0e-30);
        Assert.assertEquals(0.0, new Complex(0, testValue).expm1().getImaginary(), 1.0e-30);
    }

    @Test
    public void testLogComplex() {
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
    public void testLogBranchCut() {
        UnitTestUtils.assertEquals(new Complex(0.6931471805599453094, +3.1415926535897932384),
                                   FastMath.log(new Complex(-2.0, +0.0)),
                                   10e-12);
        UnitTestUtils.assertEquals(new Complex(0.6931471805599453094, -3.1415926535897932384),
                                   FastMath.log(new Complex(-2.0, -0.0)),
                                   10e-12);
    }

    @Test
    public void testLogZero() {
        UnitTestUtils.assertSame(negInfZero, Complex.ZERO.log());
    }

    @Test
    public void testLog1P() {
        Complex z = new Complex(2, 4);
        Complex expected = new Complex(1.60944, 0.927295);
        UnitTestUtils.assertEquals(expected, z.log1p(), 1.0e-5);
    }

    @Test
    public void testLog10Complex() {
        UnitTestUtils.assertEquals(new Complex(2.0, 0.0), new Complex(100, 0).log10(), 1.0e-15);
        UnitTestUtils.assertEquals(new Complex(2.0, 0.5 * FastMath.PI / FastMath.log(10)), new Complex(0, 100).log10(), 1.0e-15);
    }

    @Test
    @Override
    public void testLog10() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            if (x < 0) {
                // special case for Complex
                Assert.assertTrue(Double.isNaN(FastMath.log10(x)));
                Assert.assertFalse(build(x).log10().isNaN());
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
                    Assert.assertTrue(Double.isNaN(FastMath.pow(x, y)));
                    Assert.assertFalse(build(x).pow(build(y)).isNaN());
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
                    Assert.assertTrue(Double.isNaN(FastMath.pow(x, y)));
                    Assert.assertFalse(build(x).pow(y).isNaN());
                } else {
                    checkRelative(FastMath.pow(x, y), build(x).pow(y));
                }
            }
        }
    }

    @Test
    public void testPow() {
        Complex x = new Complex(3, 4);
        Complex y = new Complex(5, 6);
        Complex expected = new Complex(-1.860893, 11.83677);
        UnitTestUtils.assertEquals(expected, x.pow(y), 1.0e-5);
        UnitTestUtils.assertEquals(new Complex(-46, 9).divide(2197), new Complex(2, -3).pow(new Complex(-3, 0)), 1.0e-15);
        UnitTestUtils.assertEquals(new Complex(-1, 0).divide(8), new Complex(-2, 0).pow(new Complex(-3, 0)), 1.0e-15);
        UnitTestUtils.assertEquals(new Complex(0, 2),
                                   new Complex(-4, 0).pow(new Complex(0.5, 0)),
                                   1.0e-15);
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
        UnitTestUtils.assertSame(Complex.NaN, Complex.ONE.pow(oneInf));
        UnitTestUtils.assertSame(Complex.NaN, Complex.ONE.pow(oneNegInf));
        UnitTestUtils.assertSame(Complex.NaN, Complex.ONE.pow(infOne));
        UnitTestUtils.assertSame(Complex.NaN, Complex.ONE.pow(infInf));
        UnitTestUtils.assertSame(Complex.NaN, Complex.ONE.pow(infNegInf));
        UnitTestUtils.assertSame(Complex.NaN, Complex.ONE.pow(negInfInf));
        UnitTestUtils.assertSame(Complex.NaN, Complex.ONE.pow(negInfNegInf));
        UnitTestUtils.assertSame(Complex.INF, infOne.pow(Complex.ONE));
        UnitTestUtils.assertSame(Complex.INF, negInfOne.pow(Complex.ONE));
        UnitTestUtils.assertSame(Complex.INF, infInf.pow(Complex.ONE));
        UnitTestUtils.assertSame(Complex.INF, infNegInf.pow(Complex.ONE));
        UnitTestUtils.assertSame(Complex.INF, negInfInf.pow(Complex.ONE));
        UnitTestUtils.assertSame(Complex.INF, negInfNegInf.pow(Complex.ONE));
        UnitTestUtils.assertSame(Complex.NaN, negInfNegInf.pow(infNegInf));
        UnitTestUtils.assertSame(Complex.NaN, negInfNegInf.pow(negInfNegInf));
        UnitTestUtils.assertSame(Complex.NaN, negInfNegInf.pow(infInf));
        UnitTestUtils.assertSame(Complex.NaN, infInf.pow(infNegInf));
        UnitTestUtils.assertSame(Complex.NaN, infInf.pow(negInfNegInf));
        UnitTestUtils.assertSame(Complex.NaN, infInf.pow(infInf));
        UnitTestUtils.assertSame(Complex.NaN, infNegInf.pow(infNegInf));
        UnitTestUtils.assertSame(Complex.NaN, infNegInf.pow(negInfNegInf));
        UnitTestUtils.assertSame(Complex.NaN, infNegInf.pow(infInf));
    }

    @Test
    public void testPowZero() {
        UnitTestUtils.assertEquals(Complex.ZERO,
                                 Complex.ZERO.pow(Complex.ONE), 1.0e-12);
        UnitTestUtils.assertSame(Complex.ONE,
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
    public void testZeroPow() {
        UnitTestUtils.assertEquals(Complex.ZERO, Complex.ZERO.pow(2.0), 1.0e-5);
    }

    @Test
    public void testScalarPow() {
        Complex x = new Complex(3, 4);
        double yDouble = 5.0;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.pow(yComplex), x.pow(yDouble));
        UnitTestUtils.assertEquals(Complex.ONE.negate(), Complex.ONE.negate().pow(0.5).pow(2), 1.0e-15);
        UnitTestUtils.assertEquals(new Complex(2, 0), new Complex(4, 0).pow(0.5), 1.0e-15);
        UnitTestUtils.assertEquals(new Complex(2, 0), new Complex(4, 0).pow(new Complex(0.5, 0)), 1.0e-15);
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
       UnitTestUtils.assertSame(Complex.NaN, Complex.ONE.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.assertSame(Complex.NaN, Complex.ONE.pow(Double.NEGATIVE_INFINITY));
       UnitTestUtils.assertSame(Complex.INF, infOne.pow(1.0));
       UnitTestUtils.assertSame(Complex.INF, negInfOne.pow(1.0));
       UnitTestUtils.assertSame(Complex.INF, infInf.pow(1.0));
       UnitTestUtils.assertSame(Complex.INF, infNegInf.pow(1.0));
       UnitTestUtils.assertSame(Complex.INF, negInfInf.pow(10));
       UnitTestUtils.assertSame(Complex.INF, negInfNegInf.pow(1.0));
       UnitTestUtils.assertSame(Complex.NaN, negInfNegInf.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.assertSame(Complex.NaN, negInfNegInf.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.assertSame(Complex.NaN, infInf.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.assertSame(Complex.NaN, infInf.pow(Double.NEGATIVE_INFINITY));
       UnitTestUtils.assertSame(Complex.NaN, infNegInf.pow(Double.NEGATIVE_INFINITY));
       UnitTestUtils.assertSame(Complex.NaN, infNegInf.pow(Double.POSITIVE_INFINITY));
   }

   @Test
   public void testScalarPowZero() {
       UnitTestUtils.assertEquals(Complex.ZERO, Complex.ZERO.pow(1.0), 1.0e-12);
       UnitTestUtils.assertSame(Complex.ONE, Complex.ZERO.pow(0.0));
       UnitTestUtils.assertEquals(Complex.ONE, Complex.ONE.pow(0.0), 10e-12);
       UnitTestUtils.assertEquals(Complex.ONE, Complex.I.pow(0.0), 10e-12);
       UnitTestUtils.assertEquals(Complex.ONE, new Complex(-1, 3).pow(0.0), 10e-12);
   }

    @Test(expected=NullArgumentException.class)
    public void testpowNull() {
        Complex.ONE.pow(null);
    }

    @Test
    public void testSinComplex() {
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
    public void testSinhComplex() {
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
    public void testAsinhComplex() {
        for (double x = -2; x <= 2; x += 0.125) {
            for (double y = -2; y <= 2; y += 0.125) {
                final Complex z = new Complex(x, y);
                UnitTestUtils.assertEquals(z, z.asinh().sinh(), 1.0e-14);
            }
        }
    }

    @Test
    public void testAsinhBranchCuts() {
        UnitTestUtils.assertEquals(new Complex(FastMath.log(2 + FastMath.sqrt(3)), 0.5 * FastMath.PI),
                                   new Complex(+0.0, 2.0).asinh(),
                                   1.0e-14);
        UnitTestUtils.assertEquals(new Complex(-FastMath.log(2 + FastMath.sqrt(3)), 0.5 * FastMath.PI),
                                   new Complex(-0.0, 2.0).asinh(),
                                   1.0e-14);
    }

    @Test
    public void testAcoshComplex() {
        for (double x = -2; x <= 2; x += 0.125) {
            for (double y = -2; y <= 2; y += 0.125) {
                final Complex z = new Complex(x, y);
                UnitTestUtils.assertEquals(z, z.acosh().cosh(), 1.0e-14);
            }
        }
    }

    @Test
    public void testAcoshBranchCuts() {
        UnitTestUtils.assertEquals(new Complex(FastMath.log(2 + FastMath.sqrt(3)), +FastMath.PI),
                                   new Complex(-2.0, +0.0).acosh(),
                                   1.0e-14);
        UnitTestUtils.assertEquals(new Complex(FastMath.log(2 + FastMath.sqrt(3)), -FastMath.PI),
                                   new Complex(-2.0, -0.0).acosh(),
                                   1.0e-14);
    }

    @Test
    public void testAtanhComplex() {
        for (double x = -2; x <= 2; x += 0.125) {
            for (double y = -2; y <= 2; y += 0.125) {
                final Complex z = new Complex(x, y);
                if (FastMath.abs(x) == 1.0 && y == 0.0) {
                    Assert.assertTrue(z.atanh().isInfinite());
                } else {
                    UnitTestUtils.assertEquals(z, z.atanh().tanh(), 1.0e-14);
                }
            }
        }
    }

    @Test
    public void testAtanhBranchCuts() {
        UnitTestUtils.assertEquals(new Complex(-0.5 * FastMath.log(3), +0.5 * FastMath.PI),
                                   new Complex(-2.0, +0.0).atanh(),
                                   1.0e-14);
        UnitTestUtils.assertEquals(new Complex(-0.5 * FastMath.log(3), -0.5 * FastMath.PI),
                                   new Complex(-2.0, -0.0).atanh(),
                                   1.0e-14);
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
    public void testSqrtZero() {
        UnitTestUtils.assertEquals(Complex.ZERO, Complex.ZERO.sqrt(), 1.0e-15);
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
    @Override
    public void testCbrt() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            if ( x < 0) {
                // special case for Complex
                Assert.assertTrue(FastMath.cbrt(x) < 0);
                Assert.assertEquals(FastMath.PI / 3, build(x).cbrt().getArgument(), 1.0e-15);
            } else {
                checkRelative(FastMath.cbrt(x), build(x).cbrt());
            }
        }
    }

    @Test
    public void testCbrtComplex() {
        Complex z = new Complex(15, 2);
        UnitTestUtils.assertEquals(z, z.multiply(z).multiply(z).cbrt(), 1.0e-14);
        Complex branchCutPlus = new Complex(-8.0, +0.0);
        Complex cbrtPlus = branchCutPlus.cbrt();
        UnitTestUtils.assertEquals(branchCutPlus, cbrtPlus.multiply(cbrtPlus).multiply(cbrtPlus), 1.0e-14);
        Assert.assertEquals(1.0, cbrtPlus.getReal(), 1.0e-15);
        Assert.assertEquals(FastMath.sqrt(3.0), cbrtPlus.getImaginary(), 1.0e-15);
        Complex branchCutMinus = new Complex(-8.0, -0.0);
        Complex cbrtMinus = branchCutMinus.cbrt();
        UnitTestUtils.assertEquals(branchCutMinus, cbrtMinus.multiply(cbrtMinus).multiply(cbrtMinus), 1.0e-14);
        Assert.assertEquals(1.0, cbrtMinus.getReal(), 1.0e-15);
        Assert.assertEquals(-FastMath.sqrt(3.0), cbrtMinus.getImaginary(), 1.0e-15);
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
                        Assert.assertTrue(Double.isNaN(doubleRoot));
                    } else {
                        Assert.assertTrue(doubleRoot < 0);
                    }
                    Assert.assertEquals(FastMath.PI / n, build(x).rootN(n).getArgument(), 1.0e-15);
                } else {
                    checkRelative(FastMath.pow(x, 1.0 / n), build(x).rootN(n));
                }
            }
        }
    }

    @Test
    public void testRootNComplex() {
        Complex z = new Complex(15, 2);
        UnitTestUtils.assertEquals(z, z.multiply(z).multiply(z).rootN(3), 1.0e-14);
        Complex branchCutPlus = new Complex(-8.0, +0.0);
        Complex cbrtPlus = branchCutPlus.rootN(3);
        UnitTestUtils.assertEquals(branchCutPlus, cbrtPlus.multiply(cbrtPlus).multiply(cbrtPlus), 1.0e-14);
        Assert.assertEquals(1.0, cbrtPlus.getReal(), 1.0e-15);
        Assert.assertEquals(FastMath.sqrt(3.0), cbrtPlus.getImaginary(), 1.0e-15);
        Complex branchCutMinus = new Complex(-8.0, -0.0);
        Complex cbrtMinus = branchCutMinus.rootN(3);
        UnitTestUtils.assertEquals(branchCutMinus, cbrtMinus.multiply(cbrtMinus).multiply(cbrtMinus), 1.0e-14);
        Assert.assertEquals(1.0, cbrtMinus.getReal(), 1.0e-15);
        Assert.assertEquals(-FastMath.sqrt(3.0), cbrtMinus.getImaginary(), 1.0e-15);
    }

    @Test
    public void testTanComplex() {
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
    public void testTanhComplex() {
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
        List<Complex> thirdRootsOfZ = z.nthRoot(3);
        // Returned Collection must not be empty!
        Assert.assertEquals(3, thirdRootsOfZ.size());
        // test z_0
        Assert.assertEquals(1.0,                  thirdRootsOfZ.get(0).getReal(),      1.0e-5);
        Assert.assertEquals(1.0,                  thirdRootsOfZ.get(0).getImaginary(), 1.0e-5);
        // test z_1
        Assert.assertEquals(-1.3660254037844386,  thirdRootsOfZ.get(1).getReal(),      1.0e-5);
        Assert.assertEquals(0.36602540378443843,  thirdRootsOfZ.get(1).getImaginary(), 1.0e-5);
        // test z_2
        Assert.assertEquals(0.366025403784439,    thirdRootsOfZ.get(2).getReal(),      1.0e-5);
        Assert.assertEquals(-1.3660254037844384,  thirdRootsOfZ.get(2).getImaginary(), 1.0e-5);
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
        List<Complex> fourthRootsOfZ = z.nthRoot(4);
        // Returned Collection must not be empty!
        Assert.assertEquals(4, fourthRootsOfZ.size());
        // test z_0
        Assert.assertEquals(1.5164629308487783,     fourthRootsOfZ.get(0).getReal(),      1.0e-5);
        Assert.assertEquals(-0.14469266210702247,   fourthRootsOfZ.get(0).getImaginary(), 1.0e-5);
        // test z_1
        Assert.assertEquals(0.14469266210702256,    fourthRootsOfZ.get(1).getReal(),      1.0e-5);
        Assert.assertEquals(1.5164629308487783,     fourthRootsOfZ.get(1).getImaginary(), 1.0e-5);
        // test z_2
        Assert.assertEquals(-1.5164629308487783,    fourthRootsOfZ.get(2).getReal(),      1.0e-5);
        Assert.assertEquals(0.14469266210702267,    fourthRootsOfZ.get(2).getImaginary(), 1.0e-5);
        // test z_3
        Assert.assertEquals(-0.14469266210702275,   fourthRootsOfZ.get(3).getReal(),      1.0e-5);
        Assert.assertEquals(-1.5164629308487783,    fourthRootsOfZ.get(3).getImaginary(), 1.0e-5);
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
        List<Complex> thirdRootsOfZ = z.nthRoot(3);
        // Returned Collection must not be empty!
        Assert.assertEquals(3, thirdRootsOfZ.size());
        // test z_0
        Assert.assertEquals(2.0,                thirdRootsOfZ.get(0).getReal(),      1.0e-5);
        Assert.assertEquals(0.0,                thirdRootsOfZ.get(0).getImaginary(), 1.0e-5);
        // test z_1
        Assert.assertEquals(-1.0,               thirdRootsOfZ.get(1).getReal(),      1.0e-5);
        Assert.assertEquals(1.7320508075688774, thirdRootsOfZ.get(1).getImaginary(), 1.0e-5);
        // test z_2
        Assert.assertEquals(-1.0,               thirdRootsOfZ.get(2).getReal(),      1.0e-5);
        Assert.assertEquals(-1.732050807568877, thirdRootsOfZ.get(2).getImaginary(), 1.0e-5);
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
        List<Complex> thirdRootsOfZ = z.nthRoot(3);
        // Returned Collection must not be empty!
        Assert.assertEquals(3, thirdRootsOfZ.size());
        // test z_0
        Assert.assertEquals(1.0911236359717216,      thirdRootsOfZ.get(0).getReal(),      1.0e-5);
        Assert.assertEquals(0.6299605249474365,      thirdRootsOfZ.get(0).getImaginary(), 1.0e-5);
        // test z_1
        Assert.assertEquals(-1.0911236359717216,     thirdRootsOfZ.get(1).getReal(),      1.0e-5);
        Assert.assertEquals(0.6299605249474365,      thirdRootsOfZ.get(1).getImaginary(), 1.0e-5);
        // test z_2
        Assert.assertEquals(-2.3144374213981936E-16, thirdRootsOfZ.get(2).getReal(),      1.0e-5);
        Assert.assertEquals(-1.2599210498948732,     thirdRootsOfZ.get(2).getImaginary(), 1.0e-5);
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

    @Test
    public void testNthRootError() {
        try {
            Complex.ONE.nthRoot(-1);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.CANNOT_COMPUTE_NTH_ROOT_FOR_NEGATIVE_N,
                                miae.getSpecifier());
        }
    }

    @Test
    public void testIsMathematicalInteger() {
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
        Assert.assertFalse(new Complex(Double.NaN, imaginary).isMathematicalInteger());
        Assert.assertFalse(new Complex(Double.POSITIVE_INFINITY, imaginary).isMathematicalInteger());
        Assert.assertFalse(new Complex(Double.NEGATIVE_INFINITY, imaginary).isMathematicalInteger());
        Assert.assertFalse(new Complex(Double.MIN_NORMAL, imaginary).isMathematicalInteger());
        Assert.assertFalse(new Complex(Double.MIN_VALUE, imaginary).isMathematicalInteger());

        Assert.assertEquals(expectedForInteger, new Complex(-0.0, imaginary).isMathematicalInteger());
        Assert.assertEquals(expectedForInteger, new Complex(+0.0, imaginary).isMathematicalInteger());

        for (int i = -1000; i < 1000; ++i) {
            final double d = i;
            Assert.assertEquals(expectedForInteger, new Complex(d, imaginary).isMathematicalInteger());
            Assert.assertFalse(new Complex(FastMath.nextAfter(d, Double.POSITIVE_INFINITY), imaginary).isMathematicalInteger());
            Assert.assertFalse(new Complex(FastMath.nextAfter(d, Double.NEGATIVE_INFINITY), imaginary).isMathematicalInteger());
        }

        double minNoFractional = 0x1l << 52;
        Assert.assertEquals(expectedForInteger, new Complex(minNoFractional, imaginary).isMathematicalInteger());
        Assert.assertFalse(new Complex(minNoFractional - 0.5, imaginary).isMathematicalInteger());
        Assert.assertEquals(expectedForInteger, new Complex(minNoFractional + 0.5, imaginary).isMathematicalInteger());

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
    public void testValueOf() {
        Assert.assertEquals(2.0, Complex.valueOf(2.0).getReal(), 1.0e-15);
        Assert.assertEquals(0.0, Complex.valueOf(2.0).getImaginary(), 1.0e-15);
        Assert.assertTrue(Complex.valueOf(Double.NaN).isNaN());
        Assert.assertTrue(Complex.valueOf(Double.POSITIVE_INFINITY).isInfinite());
        Assert.assertEquals( 2.0, Complex.valueOf(2.0, -1.0).getReal(), 1.0e-15);
        Assert.assertEquals(-1.0, Complex.valueOf(2.0, -1.0).getImaginary(), 1.0e-15);
        Assert.assertTrue(Complex.valueOf(Double.NaN, 0.0).isNaN());
        Assert.assertTrue(Complex.valueOf(Double.POSITIVE_INFINITY, 0.0).isInfinite());
        Assert.assertTrue(Complex.valueOf(Double.NaN, -1.0).isNaN());
        Assert.assertTrue(Complex.valueOf(Double.POSITIVE_INFINITY, -1.0).isInfinite());
        Assert.assertTrue(Complex.valueOf(0.0, Double.NaN).isNaN());
        Assert.assertTrue(Complex.valueOf(0.0, Double.POSITIVE_INFINITY).isInfinite());
        Assert.assertTrue(Complex.valueOf(-1.0, Double.NaN).isNaN());
        Assert.assertTrue(Complex.valueOf(-1.0, Double.POSITIVE_INFINITY).isInfinite());
    }

    @Test
    public void testField() {
        Assert.assertEquals(ComplexField.getInstance(), Complex.ZERO.getField());
    }

    @Test
    public void testToString() {
        Assert.assertEquals("(1.0, -2.0)", new Complex(1, -2).toString());
    }

    @Test
    public void testScalbComplex() {
        Assert.assertEquals(0.125,  new Complex(2.0, 1.0).scalb(-4).getReal(), 1.0e-15);
        Assert.assertEquals(0.0625, new Complex(2.0, 1.0).scalb(-4).getImaginary(), 1.0e-15);
    }

    @Test
    public void testHypotComplex() {
        Assert.assertEquals(5.8269600298808519855, new Complex(3, 4).hypot(new Complex(5, 6)).getReal(), 1.0e-15);
        Assert.assertEquals(7.2078750814528590485, new Complex(3, 4).hypot(new Complex(5, 6)).getImaginary(), 1.0e-15);
    }

    @Test
    public void testCeilComplex() {
        for (double x = -3.9; x < 3.9; x += 0.05) {
            for (double y = -3.9; y < 3.9; y += 0.05) {
                final Complex z = new Complex(x, y);
                Assert.assertEquals(FastMath.ceil(x), z.ceil().getReal(), 1.0e-15);
                Assert.assertEquals(FastMath.ceil(y), z.ceil().getImaginary(), 1.0e-15);
            }
        }
    }

    @Test
    public void testFloorComplex() {
        for (double x = -3.9; x < 3.9; x += 0.05) {
            for (double y = -3.9; y < 3.9; y += 0.05) {
                final Complex z = new Complex(x, y);
                Assert.assertEquals(FastMath.floor(x), z.floor().getReal(), 1.0e-15);
                Assert.assertEquals(FastMath.floor(y), z.floor().getImaginary(), 1.0e-15);
            }
        }
    }

    @Test
    public void testRintComplex() {
        for (double x = -3.9; x < 3.9; x += 0.05) {
            for (double y = -3.9; y < 3.9; y += 0.05) {
                final Complex z = new Complex(x, y);
                Assert.assertEquals(FastMath.rint(x), z.rint().getReal(), 1.0e-15);
                Assert.assertEquals(FastMath.rint(y), z.rint().getImaginary(), 1.0e-15);
            }
        }
    }

    @Test
    public void testRemainderComplexComplex() {
        for (double x1 = -3.9; x1 < 3.9; x1 += 0.125) {
            for (double y1 = -3.9; y1 < 3.9; y1 += 0.125) {
                final Complex z1 = new Complex(x1, y1);
                for (double x2 = -3.92; x2 < 3.9; x2 += 0.125) {
                    for (double y2 = -3.92; y2 < 3.9; y2 += 0.125) {
                        final Complex z2 = new Complex(x2, y2);
                        final Complex r  = z1.remainder(z2);
                        final Complex q  = z1.subtract(r).divide(z2);
                        Assert.assertTrue(r.norm() <= z2.norm());
                        Assert.assertEquals(FastMath.rint(q.getReal()), q.getReal(), 2.0e-14);
                        Assert.assertEquals(FastMath.rint(q.getImaginary()), q.getImaginary(), 2.0e-14);
                    }
                }
            }
        }
    }

    @Test
    public void testRemainderComplexDouble() {
        for (double x1 = -3.9; x1 < 3.9; x1 += 0.125) {
            for (double y1 = -3.9; y1 < 3.9; y1 += 0.125) {
                final Complex z1 = new Complex(x1, y1);
                for (double a = -3.92; a < 3.9; a += 0.125) {
                        final Complex r  = z1.remainder(a);
                        final Complex q  = z1.subtract(r).divide(a);
                        Assert.assertTrue(r.norm() <= FastMath.abs(a));
                        Assert.assertEquals(FastMath.rint(q.getReal()), q.getReal(), 2.0e-14);
                        Assert.assertEquals(FastMath.rint(q.getImaginary()), q.getImaginary(), 2.0e-14);
                }
            }
        }
    }

    @Test
    public void testRemainderAxKr() {
        checkRemainder(new Complex(14, -5), new Complex(3, 4), new Complex(-1.0,  0.0));
        checkRemainder(new Complex(26, 120), new Complex(37, 226), new Complex(-11.0, -106.0));
        checkRemainder(new Complex(9.4, 6), new Complex(1.0, 1.0), new Complex(-0.6, 0.0));
        checkRemainder(new Complex(-5.89, 0.33), new Complex(2.4, -0.123), new Complex(-1.09, 0.084));
    }

    private void checkRemainder(final Complex c1, final Complex c2, final Complex expectedRemainder) {

        final Complex remainder = c1.remainder(c2);
        Assert.assertEquals(expectedRemainder.getReal(),      remainder.getReal(),      1.0e-15);
        Assert.assertEquals(expectedRemainder.getImaginary(), remainder.getImaginary(), 1.0e-15);

        final Complex crossCheck = c1.subtract(remainder).divide(c2);
        Assert.assertTrue(Precision.isMathematicalInteger(crossCheck.getReal()));
        Assert.assertTrue(Precision.isMathematicalInteger(crossCheck.getImaginary()));

    }

    @Test
    public void testCopySignFieldComplex() {
        for (double x1 = -3.9; x1 < 3.9; x1 += 0.08) {
            for (double y1 = -3.9; y1 < 3.9; y1 += 0.08) {
                final Complex z1 = new Complex(x1, y1);
                for (double x2 = -3.9; x2 < 3.9; x2 += 0.08) {
                    for (double y2 = -3.9; y2 < 3.9; y2 += 0.08) {
                        final Complex z2 = new Complex(x2, y2);
                        Assert.assertEquals(FastMath.copySign(x1, x2), z1.copySign(z2).getReal(), 1.0e-15);
                        Assert.assertEquals(FastMath.copySign(y1, y2), z1.copySign(z2).getImaginary(), 1.0e-15);
                    }
                }
            }
        }
    }

    @Test
    public void testCopySignDoubleComplex() {
        for (double x1 = -3.9; x1 < 3.9; x1 += 0.05) {
            for (double y1 = -3.9; y1 < 3.9; y1 += 0.05) {
                final Complex z1 = new Complex(x1, y1);
                for (double r = -3.9; r < 3.9; r += 0.05) {
                    Assert.assertEquals(FastMath.copySign(x1, r), z1.copySign(r).getReal(), 1.0e-15);
                    Assert.assertEquals(FastMath.copySign(y1, r), z1.copySign(r).getImaginary(), 1.0e-15);
                }
            }
        }
    }

    @Test
    public void testSignComplex() {
        for (double x = -3.9; x < 3.9; x += 0.05) {
            for (double y = -3.9; y < 3.9; y += 0.05) {
                final Complex z = new Complex(x, y);
                Assert.assertEquals(1.0, z.sign().norm(), 1.0e-15);
                Assert.assertEquals(FastMath.copySign(1, FastMath.signum(x)), FastMath.copySign(1, z.sign().getRealPart()), 1.0e-15);
                Assert.assertEquals(FastMath.copySign(1, FastMath.signum(y)), FastMath.copySign(1, z.sign().getImaginaryPart()), 1.0e-15);
            }
        }
        Assert.assertTrue(Complex.NaN.sign().isNaN());
        for (int sR : Arrays.asList(-1, +1)) {
            for (int sI : Arrays.asList(-1, +1)) {
                Complex z = new Complex(FastMath.copySign(0, sR), FastMath.copySign(0, sI));
                Assert.assertTrue(z.isZero());
                Complex zSign = z.sign();
                Assert.assertTrue(zSign.isZero());
                Assert.assertEquals(sR, FastMath.copySign(1, zSign.getRealPart()), 1.0e-15);
                Assert.assertEquals(sI, FastMath.copySign(1, zSign.getImaginaryPart()), 1.0e-15);
            }
        }
    }

    @Test
    public void testLinearCombination1() {
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

        UnitTestUtils.assertEquals(abSumInline, abSumArray, 0);
        UnitTestUtils.assertEquals(-1.8551294182586248737720779899, abSumInline.getReal(), 1.0e-15);

        final Complex naive = a[0].multiply(b[0]).add(a[1].multiply(b[1]));
        Assert.assertTrue(naive.subtract(abSumInline).norm() > 1.5);

    }

    @Test
    public void testSignedZeroEquality() {

        Assert.assertFalse(new Complex(-0.0, 1.0).isZero());
        Assert.assertFalse(new Complex(+0.0, 1.0).isZero());
        Assert.assertFalse(new Complex( 1.0, -0.0).isZero());
        Assert.assertFalse(new Complex( 1.0, +0.0).isZero());

        Assert.assertTrue(new Complex(-0.0, -0.0).isZero());
        Assert.assertTrue(new Complex(-0.0, +0.0).isZero());
        Assert.assertTrue(new Complex(+0.0, -0.0).isZero());
        Assert.assertTrue(new Complex(+0.0, +0.0).isZero());

        Assert.assertFalse(new Complex(-0.0, -0.0).equals(Complex.ZERO));
        Assert.assertFalse(new Complex(-0.0, +0.0).equals(Complex.ZERO));
        Assert.assertFalse(new Complex(+0.0, -0.0).equals(Complex.ZERO));
        Assert.assertTrue(new Complex(+0.0, +0.0).equals(Complex.ZERO));

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

        public TestComplex(Complex other) {
            this(other.getReal(), other.getImaginary());
        }

        @Override
        protected TestComplex createComplex(double real, double imaginary) {
            return new TestComplex(real, imaginary);
        }

    }
}
