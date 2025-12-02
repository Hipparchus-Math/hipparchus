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
package org.hipparchus.complex;

import org.hipparchus.CalculusFieldElementAbstractTest;
import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Precision;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;


public class FieldComplexTest extends CalculusFieldElementAbstractTest<FieldComplex<Binary64>> {


    private FieldComplex<Binary64> oneInf       = build(1,                        Double.POSITIVE_INFINITY);
    private FieldComplex<Binary64> oneNegInf    = build(1,                        Double.NEGATIVE_INFINITY);
    private FieldComplex<Binary64> infOne       = build(Double.POSITIVE_INFINITY, 1);
    private FieldComplex<Binary64> infZero      = build(Double.POSITIVE_INFINITY, 0);
    private FieldComplex<Binary64> infNaN       = build(Double.POSITIVE_INFINITY, Double.NaN);
    private FieldComplex<Binary64> infNegInf    = build(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY);
    private FieldComplex<Binary64> infInf       = build(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
    private FieldComplex<Binary64> negInfInf    = build(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
    private FieldComplex<Binary64> negInfZero   = build(Double.NEGATIVE_INFINITY, 0);
    private FieldComplex<Binary64> negInfOne    = build(Double.NEGATIVE_INFINITY, 1);
    private FieldComplex<Binary64> negInfNaN    = build(Double.NEGATIVE_INFINITY, Double.NaN);
    private FieldComplex<Binary64> negInfNegInf = build(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
    private FieldComplex<Binary64> oneNaN       = build(1,                        Double.NaN);
    private FieldComplex<Binary64> zeroInf      = build(0,                        Double.POSITIVE_INFINITY);
    private FieldComplex<Binary64> zeroNaN      = build(0,                        Double.NaN);
    private FieldComplex<Binary64> nanInf       = build(Double.NaN,               Double.POSITIVE_INFINITY);
    private FieldComplex<Binary64> nanNegInf    = build(Double.NaN,               Double.NEGATIVE_INFINITY);
    private FieldComplex<Binary64> nanZero      = build(Double.NaN);

    @Override
    protected FieldComplex<Binary64> build(final double x) {
        return build(x, 0.0);
    }

    private FieldComplex<Binary64> build(final double real, double imaginary) {
        return new FieldComplex<>(new Binary64(real), new Binary64(imaginary));
    }

    @Test
    void testConstructor() {
        FieldComplex<Binary64> z = build(3.0, 4.0);
        assertEquals(3.0, z.getRealPart().getReal(), 1.0e-5);
        assertEquals(4.0, z.getImaginary().getReal(), 1.0e-5);
        assertEquals(3.0, z.getRealPart().getReal(), 1.0e-5);
        assertEquals(4.0, z.getImaginaryPart().getReal(), 1.0e-5);
    }

    @Test
    void testConstructorNaN() {
        FieldComplex<Binary64> z = build(3.0, Double.NaN);
        assertTrue(z.isNaN());

        z = build(Double.NaN, 4.0);
        assertTrue(z.isNaN());

        z = build(3.0, 4.0);
        assertFalse(z.isNaN());
    }

    @Test
    void testNorm() {
        FieldComplex<Binary64> z = build(3.0, 4.0);
        assertEquals(5.0, z.norm(), 1.0e-5);
    }

    @Test
    void testNormNaN() {
        assertTrue(Double.isNaN(FieldComplex.getNaN(Binary64Field.getInstance()).norm()));
        FieldComplex<Binary64> z = build(Double.POSITIVE_INFINITY, Double.NaN);
        assertTrue(Double.isNaN(z.norm()));
    }

    @Test
    void testNormInfinite() {
        FieldComplex<Binary64> z = FieldComplex.getNaN(Binary64Field.getInstance()).newInstance(Double.POSITIVE_INFINITY);
        assertEquals(Double.POSITIVE_INFINITY, z.norm(), 0);
        z = build(0, Double.NEGATIVE_INFINITY);
        assertEquals(Double.POSITIVE_INFINITY, z.norm(), 0);
        z = build(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY);
        assertEquals(Double.POSITIVE_INFINITY, z.norm(), 0);
    }

    @Test
    void testAdd() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        FieldComplex<Binary64> y = build(5.0, 6.0);
        FieldComplex<Binary64> z = x.add(y);
        assertEquals(8.0, z.getRealPart().getReal(), 1.0e-5);
        assertEquals(10.0, z.getImaginary().getReal(), 1.0e-5);
    }

    @Test
    void testAddT() {
        FieldComplex<Binary64> z = build(3.0, 4.0).add(new Binary64(5.0));
        assertEquals(8.0, z.getRealPart().getReal(), 1.0e-5);
        assertEquals(4.0, z.getImaginary().getReal(), 1.0e-5);
    }

    @Test
    void testAddNaN() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        FieldComplex<Binary64> z = x.add(FieldComplex.getNaN(Binary64Field.getInstance()));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), z);
        z = build(1, Double.NaN);
        FieldComplex<Binary64> w = x.add(z);
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), w);
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), FieldComplex.getNaN(Binary64Field.getInstance()).add(Double.NaN));
    }

    @Test
    void testAddInf() {
        FieldComplex<Binary64> x = build(1, 1);
        FieldComplex<Binary64> z = build(Double.POSITIVE_INFINITY, 0);
        FieldComplex<Binary64> w = x.add(z);
        assertEquals(1, w.getImaginary().getReal(), 0);
        assertEquals(Double.POSITIVE_INFINITY, w.getRealPart().getReal(), 0);

        x = build(Double.NEGATIVE_INFINITY, 0);
        assertTrue(Double.isNaN(x.add(z).getReal()));
    }

    @Test
    void testScalarAdd() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        double yDouble = 2.0;
        FieldComplex<Binary64> yComplex = build(yDouble);
        assertEquals(x.add(yComplex), x.add(yDouble));
        assertEquals(x.add(yComplex), x.add(new Binary64(yDouble)));
    }

    @Test
    void testScalarAddNaN() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        double yDouble = Double.NaN;
        FieldComplex<Binary64> yComplex = build(yDouble);
        assertEquals(x.add(yComplex), x.add(yDouble));
        assertEquals(x.add(yComplex), x.add(new Binary64(yDouble)));
        assertTrue(build(Double.NaN).add(0).isNaN());
        assertTrue(build(Double.NaN).add(Binary64.ZERO).isNaN());
    }

    @Test
    void testScalarAddInf() {
        FieldComplex<Binary64> x = build(1, 1);
        double yDouble = Double.POSITIVE_INFINITY;

        FieldComplex<Binary64> yComplex = build(yDouble);
        assertEquals(x.add(yComplex), x.add(yDouble));
        assertEquals(x.add(yComplex), x.add(new Binary64(yDouble)));

        x = build(Double.NEGATIVE_INFINITY, 0);
        assertEquals(x.add(yComplex), x.add(yDouble));
        assertEquals(x.add(yComplex), x.add(new Binary64(yDouble)));

    }

    @Test
    void testConjugate() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        FieldComplex<Binary64> z = x.conjugate();
        assertEquals(3.0, z.getRealPart().getReal(), 1.0e-5);
        assertEquals(-4.0, z.getImaginaryPart().getReal(), 1.0e-5);
    }

    @Test
    void testConjugateNaN() {
        FieldComplex<Binary64> z = FieldComplex.getNaN(Binary64Field.getInstance()).conjugate();
        assertTrue(z.isNaN());
    }

    @Test
    void testConjugateInfiinite() {
        FieldComplex<Binary64> z = build(0, Double.POSITIVE_INFINITY);
        assertEquals(Double.NEGATIVE_INFINITY, z.conjugate().getImaginary().getReal(), 0);
        z = build(0, Double.NEGATIVE_INFINITY);
        assertEquals(Double.POSITIVE_INFINITY, z.conjugate().getImaginary().getReal(), 0);
    }

    @Test
    void testDivide() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        FieldComplex<Binary64> y = build(5.0, 6.0);
        FieldComplex<Binary64> z = x.divide(y);
        assertEquals(39.0 / 61.0, z.getRealPart().getReal(), 1.0e-5);
        assertEquals(2.0 / 61.0, z.getImaginaryPart().getReal(), 1.0e-5);
    }

    @Test
    void testDivideReal() {
        FieldComplex<Binary64> x = build(2d, 3d);
        FieldComplex<Binary64> y = build(2d, 0d);
        assertEquals(build(1d, 1.5), x.divide(y));

    }

    @Test
    void testDivideImaginary() {
        FieldComplex<Binary64> x = build(2d, 3d);
        FieldComplex<Binary64> y = build(0d, 2d);
        assertEquals(build(1.5d, -1d), x.divide(y));
    }

    @Test
    void testDivideInf() {
        FieldComplex<Binary64> x = build(3, 4);
        FieldComplex<Binary64> w = build(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        assertEquals(x.divide(w), FieldComplex.getZero(Binary64Field.getInstance()));

        FieldComplex<Binary64> z = w.divide(x);
        assertTrue(Double.isNaN(z.getReal()));
        assertEquals(Double.POSITIVE_INFINITY, z.getImaginary().getReal(), 0);

        w = build(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
        z = w.divide(x);
        assertTrue(Double.isNaN(z.getImaginary().getReal()));
        assertEquals(Double.POSITIVE_INFINITY, z.getRealPart().getReal(), 0);

        w = build(1, Double.POSITIVE_INFINITY);
        z = w.divide(w);
        assertTrue(Double.isNaN(z.getReal()));
        assertTrue(Double.isNaN(z.getImaginary().getReal()));
    }

    @Test
    void testDivideZero() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        FieldComplex<Binary64> z = x.divide(FieldComplex.getZero(Binary64Field.getInstance()));
        // Assertions.assertEquals(z, FieldComplex.getInf(Binary64Field.getInstance())); // See MATH-657
        assertEquals(z, FieldComplex.getNaN(Binary64Field.getInstance()));
        assertTrue(build(3.0).divide(Binary64.ZERO).isNaN());
    }

    @Test
    void testDivideZeroZero() {
        FieldComplex<Binary64> x = build(0.0, 0.0);
        FieldComplex<Binary64> z = x.divide(FieldComplex.getZero(Binary64Field.getInstance()));
        assertEquals(z, FieldComplex.getNaN(Binary64Field.getInstance()));
    }

    @Test
    void testDivideNaN() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        FieldComplex<Binary64> z = x.divide(FieldComplex.getNaN(Binary64Field.getInstance()));
        assertTrue(z.isNaN());
        assertTrue(x.divide(Binary64.NAN).isNaN());
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).divide(Binary64.ONE).isNaN());
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).divide(FieldComplex.getNaN(Binary64Field.getInstance())).isNaN());
    }

    @Test
    void testDivideNaNInf() {
       FieldComplex<Binary64> z = oneInf.divide(FieldComplex.getOne(Binary64Field.getInstance()));
       assertTrue(Double.isNaN(z.getReal()));
       assertEquals(Double.POSITIVE_INFINITY, z.getImaginary().getReal(), 0);

       z = negInfNegInf.divide(oneNaN);
       assertTrue(Double.isNaN(z.getReal()));
       assertTrue(Double.isNaN(z.getImaginary().getReal()));

       z = negInfInf.divide(FieldComplex.getOne(Binary64Field.getInstance()));
       assertTrue(Double.isNaN(z.getReal()));
       assertTrue(Double.isNaN(z.getImaginary().getReal()));
    }

    @Test
    void testScalarDivide() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        double yDouble = 2.0;
        FieldComplex<Binary64> yComplex = build(yDouble);
        assertEquals(x.divide(yComplex), x.divide(yDouble));
        assertEquals(x.divide(yComplex), x.divide(new Binary64(yDouble)));
    }

    @Test
    void testScalarDivideNaN() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        double yDouble = Double.NaN;
        FieldComplex<Binary64> yComplex = build(yDouble);
        assertEquals(x.divide(yComplex), x.divide(yDouble));
        assertEquals(x.divide(yComplex), x.divide(new Binary64(yDouble)));
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).divide(Double.NaN).isNaN());
    }

    @Test
    void testScalarDivideInf() {
        FieldComplex<Binary64> x = build(1,1);
        double yDouble = Double.POSITIVE_INFINITY;
        FieldComplex<Binary64> yComplex = build(yDouble);
        UnitTestUtils.customAssertEquals(x.divide(yComplex), x.divide(yDouble), 0);
        UnitTestUtils.customAssertEquals(x.divide(yComplex), x.divide(new Binary64(yDouble)), 0);

        yDouble = Double.NEGATIVE_INFINITY;
        yComplex = build(yDouble);
        UnitTestUtils.customAssertEquals(x.divide(yComplex), x.divide(yDouble), 0);
        UnitTestUtils.customAssertEquals(x.divide(yComplex), x.divide(new Binary64(yDouble)), 0);

        x = build(1, Double.NEGATIVE_INFINITY);
        UnitTestUtils.customAssertEquals(x.divide(yComplex), x.divide(yDouble), 0);
        UnitTestUtils.customAssertEquals(x.divide(yComplex), x.divide(new Binary64(yDouble)), 0);

    }

    @Test
    void testScalarDivideZero() {
        FieldComplex<Binary64> x = build(1,1);
        UnitTestUtils.customAssertEquals(x.divide(FieldComplex.getZero(Binary64Field.getInstance())), x.divide(0), 0);
        UnitTestUtils.customAssertEquals(x.divide(FieldComplex.getZero(Binary64Field.getInstance())), x.divide(new Binary64(0)), 0);
    }

    @Test
    void testReciprocal() {
        FieldComplex<Binary64> z = build(5.0, 6.0);
        FieldComplex<Binary64> act = z.reciprocal();
        double expRe = 5.0 / 61.0;
        double expIm = -6.0 / 61.0;
        assertEquals(expRe, act.getRealPart().getReal(), FastMath.ulp(expRe));
        assertEquals(expIm, act.getImaginaryPart().getReal(), FastMath.ulp(expIm));
    }

    @Test
    void testReciprocalReal() {
        FieldComplex<Binary64> z = build(-2.0, 0.0);
        assertTrue(FieldComplex.equals(build(-0.5, 0.0), z.reciprocal()));
    }

    @Test
    void testReciprocalImaginary() {
        FieldComplex<Binary64> z = build(0.0, -2.0);
        assertEquals(build(0.0, 0.5), z.reciprocal());
    }

    @Test
    void testReciprocalInf() {
        FieldComplex<Binary64> z = build(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        assertEquals(z.reciprocal(), FieldComplex.getZero(Binary64Field.getInstance()));

        z = build(1, Double.POSITIVE_INFINITY).reciprocal();
        assertEquals(z, FieldComplex.getZero(Binary64Field.getInstance()));
    }

    @Test
    void testReciprocalZero() {
        assertEquals(FieldComplex.getZero(Binary64Field.getInstance()).reciprocal(), FieldComplex.getInf(Binary64Field.getInstance()));
    }

    @Test
    void testReciprocalNaN() {
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).reciprocal().isNaN());
    }

    @Test
    void testMultiply() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        FieldComplex<Binary64> y = build(5.0, 6.0);
        FieldComplex<Binary64> z = x.multiply(y);
        assertEquals(-9.0, z.getRealPart().getReal(), 1.0e-5);
        assertEquals(38.0, z.getImaginaryPart().getReal(), 1.0e-5);
    }

    @Test
    void testMultiplyT() {
        FieldComplex<Binary64> z = build(3.0, 4.0).multiply(new Binary64(5.0));
        assertEquals(15.0, z.getRealPart().getReal(), 1.0e-5);
        assertEquals(20.0, z.getImaginaryPart().getReal(), 1.0e-5);
    }

    @Test
    void testMultiplyNaN() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        FieldComplex<Binary64> z = x.multiply(FieldComplex.getNaN(Binary64Field.getInstance()));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), z);
        z = FieldComplex.getNaN(Binary64Field.getInstance()).multiply(5);
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), z);
    }

    @Test
    void testMultiplyInfInf() {
        // Assertions.assertTrue(infInf.multiply(infInf).isNaN()); // MATH-620
        assertTrue(infInf.multiply(infInf).isInfinite());
    }

    @Test
    void testMultiplyNaNInf() {
        FieldComplex<Binary64> z = build(1,1);
        FieldComplex<Binary64> w = z.multiply(infOne);
        assertEquals(Double.POSITIVE_INFINITY, w.getRealPart().getReal(), 0);
        assertEquals(Double.POSITIVE_INFINITY, w.getImaginaryPart().getReal(), 0);

        // [MATH-164]
        assertEquals(build(1, 0).multiply(infInf), FieldComplex.getInf(Binary64Field.getInstance()));
        assertEquals(build(-1, 0).multiply(infInf), FieldComplex.getInf(Binary64Field.getInstance()));
        assertEquals(build(1, 0).multiply(negInfZero), FieldComplex.getInf(Binary64Field.getInstance()));

        w = oneInf.multiply(oneNegInf);
        assertEquals(Double.POSITIVE_INFINITY, w.getRealPart().getReal(), 0);
        assertEquals(Double.POSITIVE_INFINITY, w.getImaginaryPart().getReal(), 0);

        w = negInfNegInf.multiply(oneNaN);
        assertTrue(w.getRealPart().isNaN());
        assertTrue(w.getImaginaryPart().isNaN());

        z = build(1, Double.NEGATIVE_INFINITY);
        UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), z.square());
    }

    @Test
    void testScalarMultiply() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        double yDouble = 2.0;
        FieldComplex<Binary64> yComplex = build(yDouble);
        assertEquals(x.multiply(yComplex), x.multiply(yDouble));
        assertEquals(x.multiply(yComplex), x.multiply(new Binary64(yDouble)));
        int zInt = -5;
        FieldComplex<Binary64> zComplex = build(zInt);
        assertEquals(x.multiply(zComplex), x.multiply(zInt));
    }

    @Test
    void testScalarMultiplyNaN() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        double yDouble = Double.NaN;
        FieldComplex<Binary64> yComplex = build(yDouble);
        assertEquals(x.multiply(yComplex), x.multiply(yDouble));
        assertEquals(x.multiply(yComplex), x.multiply(new Binary64(yDouble)));
        assertTrue(build(Double.NaN, 0).multiply(5).isNaN());
        assertTrue(build(Double.NaN, 0).multiply(new Binary64(5)).isNaN());
        assertTrue(build(0, Double.NaN).multiply(5).isNaN());
        assertTrue(build(0, Double.NaN).multiply(new Binary64(5)).isNaN());
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).multiply(5).isNaN());
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).multiply(new Binary64(5)).isNaN());
        assertTrue(build(Double.NaN, 0).multiply(5.0).isNaN());
        assertTrue(build(Double.NaN, 0).multiply(new Binary64(5)).isNaN());
        assertTrue(build(0, Double.NaN).multiply(5.0).isNaN());
        assertTrue(build(0, Double.NaN).multiply(new Binary64(5)).isNaN());
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).multiply(5.0).isNaN());
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).multiply(new Binary64(5.0)).isNaN());
        assertTrue(FieldComplex.getOne(Binary64Field.getInstance()).multiply(Double.NaN).isNaN());
        assertTrue(FieldComplex.getOne(Binary64Field.getInstance()).multiply(new Binary64(Double.NaN)).isNaN());
    }

    @Test
    void testScalarMultiplyInf() {
        FieldComplex<Binary64> x = build(1, 1);
        double yDouble = Double.POSITIVE_INFINITY;
        FieldComplex<Binary64> yComplex = build(yDouble);
        assertEquals(x.multiply(yComplex), x.multiply(yDouble));
        assertEquals(x.multiply(yComplex), x.multiply(new Binary64(yDouble)));

        yDouble = Double.NEGATIVE_INFINITY;
        yComplex = build(yDouble);
        assertEquals(x.multiply(yComplex), x.multiply(yDouble));
        assertEquals(x.multiply(yComplex), x.multiply(new Binary64(yDouble)));

        assertTrue(build(Double.POSITIVE_INFINITY, 0).multiply(5).isInfinite());
        assertTrue(build(Double.POSITIVE_INFINITY, 0).multiply(new Binary64(5)).isInfinite());
        assertTrue(build(0, Double.POSITIVE_INFINITY).multiply(5).isInfinite());
        assertTrue(build(0, Double.POSITIVE_INFINITY).multiply(new Binary64(5)).isInfinite());
        assertTrue(FieldComplex.getInf(Binary64Field.getInstance()).multiply(5).isInfinite());
        assertTrue(FieldComplex.getInf(Binary64Field.getInstance()).multiply(new Binary64(5)).isInfinite());
        assertTrue(build(Double.POSITIVE_INFINITY, 0).multiply(5.0).isInfinite());
        assertTrue(build(Double.POSITIVE_INFINITY, 0).multiply(new Binary64(5.0)).isInfinite());
        assertTrue(build(0, Double.POSITIVE_INFINITY).multiply(5.0).isInfinite());
        assertTrue(build(0, Double.POSITIVE_INFINITY).multiply(new Binary64(5.0)).isInfinite());
        assertTrue(FieldComplex.getInf(Binary64Field.getInstance()).multiply(5.0).isInfinite());
        assertTrue(FieldComplex.getInf(Binary64Field.getInstance()).multiply(new Binary64(5.0)).isInfinite());
        assertTrue(FieldComplex.getOne(Binary64Field.getInstance()).multiply(Double.POSITIVE_INFINITY).isInfinite());
        assertTrue(FieldComplex.getOne(Binary64Field.getInstance()).multiply(new Binary64(Double.POSITIVE_INFINITY)).isInfinite());
    }

    @Test
    void testNegate() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        FieldComplex<Binary64> z = x.negate();
        assertEquals(-3.0, z.getRealPart().getReal(), 1.0e-5);
        assertEquals(-4.0, z.getImaginaryPart().getReal(), 1.0e-5);
    }

    @Test
    void testNegateNaN() {
        FieldComplex<Binary64> z = FieldComplex.getNaN(Binary64Field.getInstance()).negate();
        assertTrue(z.isNaN());
    }

    @Test
    void testSubtract() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        FieldComplex<Binary64> y = build(5.0, 6.0);
        FieldComplex<Binary64> z = x.subtract(y);
        assertEquals(-2.0, z.getRealPart().getReal(), 1.0e-5);
        assertEquals(-2.0, z.getImaginaryPart().getReal(), 1.0e-5);
    }

    @Test
    void testSubtractT() {
        FieldComplex<Binary64> z = build(3.0, 4.0).subtract(new Binary64(5.0));
        assertEquals(-2.0, z.getRealPart().getReal(), 1.0e-5);
        assertEquals( 4.0, z.getImaginary().getReal(), 1.0e-5);
    }

    @Test
    void testSubtractNaN() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        FieldComplex<Binary64> z = x.subtract(FieldComplex.getNaN(Binary64Field.getInstance()));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), z);
        z = build(1, Double.NaN);
        FieldComplex<Binary64> w = x.subtract(z);
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), w);
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), FieldComplex.getNaN(Binary64Field.getInstance()).subtract(FieldComplex.getNaN(Binary64Field.getInstance())));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), FieldComplex.getNaN(Binary64Field.getInstance()).subtract(Double.NaN));
    }

    @Test
    void testSubtractInf() {
        FieldComplex<Binary64> x = build(1, 1);
        FieldComplex<Binary64> z = build(Double.NEGATIVE_INFINITY, 0);
        FieldComplex<Binary64> w = x.subtract(z);
        assertEquals(1, w.getImaginaryPart().getReal(), 0);
        assertEquals(Double.POSITIVE_INFINITY, w.getRealPart().getReal(), 0);

        x = build(Double.NEGATIVE_INFINITY, 0);
        assertTrue(Double.isNaN(x.subtract(z).getReal()));
    }

    @Test
    void testScalarSubtract() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        double yDouble = 2.0;
        FieldComplex<Binary64> yComplex = build(yDouble);
        assertEquals(x.subtract(yComplex), x.subtract(yDouble));
        assertEquals(x.subtract(yComplex), x.subtract(new Binary64(yDouble)));
    }

    @Test
    void testScalarSubtractNaN() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        double yDouble = Double.NaN;
        FieldComplex<Binary64> yComplex = build(yDouble);
        assertEquals(x.subtract(yComplex), x.subtract(yDouble));
        assertEquals(x.subtract(yComplex), x.subtract(new Binary64(yDouble)));
        assertTrue(build(Double.NaN).subtract(0).isNaN());
        assertTrue(build(Double.NaN).subtract(Binary64.ZERO).isNaN());
    }

    @Test
    void testScalarSubtractInf() {
        FieldComplex<Binary64> x = build(1, 1);
        double yDouble = Double.POSITIVE_INFINITY;
        FieldComplex<Binary64> yComplex = build(yDouble);
        assertEquals(x.subtract(yComplex), x.subtract(yDouble));
        assertEquals(x.subtract(yComplex), x.subtract(new Binary64(yDouble)));

        x = build(Double.NEGATIVE_INFINITY, 0);
        assertEquals(x.subtract(yComplex), x.subtract(yDouble));
        assertEquals(x.subtract(yComplex), x.subtract(new Binary64(yDouble)));
    }

    @Test
    void testEqualsNull() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        assertNotEquals(null, x);
    }

    @Test
    void testFloatingPointEqualsPrecondition1() {
        assertThrows(NullPointerException.class, () -> {
            FieldComplex.equals(build(3.0, 4.0), null, 3);
        });
    }

    @Test
    void testFloatingPointEqualsPrecondition2() {
        assertThrows(NullPointerException.class, () -> {
            FieldComplex.equals(null, build(3.0, 4.0), 3);
        });
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    void testEqualsClass() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        assertNotEquals(x, this);
    }

    @Test
    void testEqualsSame() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        assertEquals(x, x);
    }

    @Test
    void testFloatingPointEquals() {
        double re = -3.21;
        double im = 456789e10;

        final FieldComplex<Binary64> x = build(re, im);
        FieldComplex<Binary64> y = build(re, im);

        assertEquals(x, y);
        assertTrue(FieldComplex.equals(x, y));

        final int maxUlps = 5;
        for (int i = 0; i < maxUlps; i++) {
            re = FastMath.nextUp(re);
            im = FastMath.nextUp(im);
        }
        y = build(re, im);
        assertTrue(FieldComplex.equals(x, y, maxUlps));

        re = FastMath.nextUp(re);
        im = FastMath.nextUp(im);
        y = build(re, im);
        assertFalse(FieldComplex.equals(x, y, maxUlps));
    }

    @Test
    void testFloatingPointEqualsNaN() {
        FieldComplex<Binary64> c = build(Double.NaN, 1);
        assertFalse(FieldComplex.equals(c, c));

        c = build(1, Double.NaN);
        assertFalse(FieldComplex.equals(c, c));
    }

    @Test
    void testFloatingPointEqualsWithAllowedDelta() {
        final double re = 153.0000;
        final double im = 152.9375;
        final double tol1 = 0.0625;
        final FieldComplex<Binary64> x = build(re, im);
        final FieldComplex<Binary64> y = build(re + tol1, im + tol1);
        assertTrue(FieldComplex.equals(x, y, tol1));

        final double tol2 = 0.0624;
        assertFalse(FieldComplex.equals(x, y, tol2));
    }

    @Test
    void testFloatingPointEqualsWithAllowedDeltaNaN() {
        final FieldComplex<Binary64> x = build(0, Double.NaN);
        final FieldComplex<Binary64> y = build(Double.NaN, 0);
        assertFalse(FieldComplex.equals(x, FieldComplex.getZero(Binary64Field.getInstance()), 0.1));
        assertFalse(FieldComplex.equals(x, x, 0.1));
        assertFalse(FieldComplex.equals(x, y, 0.1));
    }

    @Test
    void testFloatingPointEqualsWithRelativeTolerance() {
        final double tol = 1e-4;
        final double re = 1;
        final double im = 1e10;

        final double f = 1 + tol;
        final FieldComplex<Binary64> x = build(re, im);
        final FieldComplex<Binary64> y = build(re * f, im * f);
        assertTrue(FieldComplex.equalsWithRelativeTolerance(x, y, tol));
    }

    @Test
    void testFloatingPointEqualsWithRelativeToleranceNaN() {
        final FieldComplex<Binary64> x = build(0, Double.NaN);
        final FieldComplex<Binary64> y = build(Double.NaN, 0);
        assertFalse(FieldComplex.equalsWithRelativeTolerance(x, FieldComplex.getZero(Binary64Field.getInstance()), 0.1));
        assertFalse(FieldComplex.equalsWithRelativeTolerance(x, x, 0.1));
        assertFalse(FieldComplex.equalsWithRelativeTolerance(x, y, 0.1));
    }

    @Test
    void testEqualsTrue() {
        FieldComplex<Binary64> x = build(3.0, 4.0);
        FieldComplex<Binary64> y = build(3.0, 4.0);
        assertEquals(x, y);
    }

    @Test
    void testEqualsRealDifference() {
        FieldComplex<Binary64> x = build(0.0, 0.0);
        FieldComplex<Binary64> y = build(0.0 + Double.MIN_VALUE, 0.0);
        assertNotEquals(x, y);
    }

    @Test
    void testEqualsImaginaryDifference() {
        FieldComplex<Binary64> x = build(0.0, 0.0);
        FieldComplex<Binary64> y = build(0.0, 0.0 + Double.MIN_VALUE);
        assertNotEquals(x, y);
    }

    @Test
    void testEqualsNaN() {
        FieldComplex<Binary64> realNaN = build(Double.NaN, 0.0);
        FieldComplex<Binary64> imaginaryNaN = build(0.0, Double.NaN);
        FieldComplex<Binary64> complexNaN = FieldComplex.getNaN(Binary64Field.getInstance());
        assertEquals(realNaN, imaginaryNaN);
        assertEquals(imaginaryNaN, complexNaN);
        assertEquals(realNaN, complexNaN);
    }

    @Test
    void testHashCode() {
        FieldComplex<Binary64> x = build(0.0, 0.0);
        FieldComplex<Binary64> y = build(0.0, 0.0 + Double.MIN_VALUE);
        assertNotEquals(x.hashCode(), y.hashCode());
        y = build(0.0 + Double.MIN_VALUE, 0.0);
        assertNotEquals(x.hashCode(), y.hashCode());
        FieldComplex<Binary64> realNaN = build(Double.NaN, 0.0);
        FieldComplex<Binary64> imaginaryNaN = build(0.0, Double.NaN);
        assertEquals(realNaN.hashCode(), imaginaryNaN.hashCode());
        assertEquals(imaginaryNaN.hashCode(), FieldComplex.getNaN(Binary64Field.getInstance()).hashCode());

        // MATH-1118
        // "equals" and "hashCode" must be compatible: if two objects have
        // different hash codes, "equals" must return false.
        final String msg = "'equals' not compatible with 'hashCode'";

        x = build(0.0, 0.0);
        y = build(0.0, -0.0);
        assertTrue(x.hashCode() != y.hashCode());
        assertNotEquals(x, y, msg);

        x = build(0.0, 0.0);
        y = build(-0.0, 0.0);
        assertTrue(x.hashCode() != y.hashCode());
        assertNotEquals(x, y, msg);
    }

    @Test
    void testToDegreesComplex() {
        FieldComplex<Binary64> z = build(3, 4);
        FieldComplex<Binary64> expected = build(FastMath.toDegrees(z.getRealPart().getReal()), FastMath.toDegrees(z.getImaginaryPart().getReal()));
        UnitTestUtils.customAssertEquals(expected, z.toDegrees(), 1.0e-15);
    }

    @Test
    void testToRadiansComplex() {
        FieldComplex<Binary64> z = build(3, 4);
        FieldComplex<Binary64> expected = build(FastMath.toRadians(z.getRealPart().getReal()), FastMath.toRadians(z.getImaginaryPart().getReal()));
        UnitTestUtils.customAssertEquals(expected, z.toRadians(), 1.0e-15);
    }

    @Test
    void testAcosComplex() {
        FieldComplex<Binary64> z = build(3, 4);
        FieldComplex<Binary64> expected = build(0.936812, -2.30551);
        UnitTestUtils.customAssertEquals(expected, z.acos(), 1.0e-5);
        UnitTestUtils.customAssertEquals(build(FastMath.acos(0), 0),
                                         FieldComplex.getZero(Binary64Field.getInstance()).acos(), 1.0e-12);
    }

    @Test
    void testAcosNaN() {
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).acos().isNaN());
    }

    @Test
    void testAcosInf() {
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), oneInf.acos());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), oneNegInf.acos());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infOne.acos());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfOne.acos());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.acos());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.acos());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfInf.acos());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.acos());
    }

    @Test
    void testAcosBranchCuts() {
        UnitTestUtils.customAssertEquals(build(3.141592653589793238462, -0.76103968373182660633),
                                         FastMath.acos(build(-1.3038404810405297, +0.0)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(build(3.141592653589793238462, +0.76103968373182660633),
                                         FastMath.acos(build(-1.3038404810405297, -0.0)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(build(0.0, -0.76103968373182660633),
                                         FastMath.acos(build(1.3038404810405297, +0.0)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(build(0.0, +0.76103968373182660633),
                                         FastMath.acos(build(1.3038404810405297, -0.0)),
                                         1.0e-14);
    }

    @Test
    void testAsinComplex() {
        FieldComplex<Binary64> z = build(3, 4);
        FieldComplex<Binary64> expected = build(0.633984, 2.30551);
        UnitTestUtils.customAssertEquals(expected, z.asin(), 1.0e-5);
    }

    @Test
    void testAsinNaN() {
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).asin().isNaN());
    }

    @Test
    void testAsinInf() {
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), oneInf.asin());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), oneNegInf.asin());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infOne.asin());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfOne.asin());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.asin());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.asin());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfInf.asin());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.asin());
    }

    @Test
    void testAsinBranchCuts() {
        UnitTestUtils.customAssertEquals(build(-1.57079632679489661923, +0.76103968373182660633),
                                         FastMath.asin(build(-1.3038404810405297, +0.0)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(build(-1.57079632679489661923, -0.76103968373182660633),
                                         FastMath.asin(build(-1.3038404810405297, -0.0)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(build(1.57079632679489661923, +0.76103968373182660633),
                                         FastMath.asin(build(1.3038404810405297, +0.0)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(build(1.57079632679489661923, -0.76103968373182660633),
                                         FastMath.asin(build(1.3038404810405297, -0.0)),
                                         1.0e-14);
    }

    @Test
    void testAtanComplex() {
        FieldComplex<Binary64> z = build(3, 4);
        FieldComplex<Binary64> expected = build(1.44831, 0.158997);
        UnitTestUtils.customAssertEquals(expected, z.atan(), 1.0e-5);
    }

    @Test
    void testAtanInf() {
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), oneInf.atan());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), oneNegInf.atan());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infOne.atan());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfOne.atan());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.atan());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.atan());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfInf.atan());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.atan());
    }

    @Test
    void testAtanI() {
        assertTrue(FieldComplex.getI(Binary64Field.getInstance()).atan().isNaN());
    }

    @Test
    void testAtanNaN() {
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).atan().isNaN());
    }

    @Test
    void testAtanBranchCuts() {
        UnitTestUtils.customAssertEquals(build(+1.5707963267948966192, +1.0986122886681096913),
                                         FastMath.atan(build(+0.0, 1.25)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(build(-1.5707963267948966192, +1.0986122886681096913),
                                         FastMath.atan(build(-0.0, 1.25)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(build(+1.5707963267948966192, -1.0986122886681096913),
                                         FastMath.atan(build(+0.0, -1.25)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(build(-1.5707963267948966192, -1.0986122886681096913),
                                         FastMath.atan(build(-0.0, -1.25)),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(build(0.0, +0.25541281188299534160),
                                         FastMath.atan(build(+0.0, 0.25)),
                                         1.0e-14);
        assertTrue(FastMath.copySign(1.0, FastMath.atan(build(+0.0, 0.25)).getReal()) > 0.0);
        UnitTestUtils.customAssertEquals(build(0.0, +0.25541281188299534160),
                                         FastMath.atan(build(-0.0, 0.25)),
                                         1.0e-14);
        assertTrue(FastMath.copySign(1.0, FastMath.atan(build(-0.0, 0.25)).getReal()) < 0.0);
        UnitTestUtils.customAssertEquals(build(0.0, -0.25541281188299534160),
                                         FastMath.atan(build(+0.0, -0.25)),
                                         1.0e-14);
        assertTrue(FastMath.copySign(1.0, FastMath.atan(build(+0.0, -0.25)).getReal()) > 0.0);
        UnitTestUtils.customAssertEquals(build(0.0, -0.25541281188299534160),
                                         FastMath.atan(build(-0.0, -0.25)),
                                         1.0e-14);
        assertTrue(FastMath.copySign(1.0, FastMath.atan(build(-0.0, -0.25)).getReal()) < 0.0);
    }

    @Test
    public void testAtanReal() {
        final FieldComplex<Binary64> zP = build(0.8734729023516287, 0.0);
        final FieldComplex<Binary64> aP = build(0.717964439926383,  0.0);
        Assertions.assertEquals(aP, zP.atan());
        Assertions.assertEquals(1.0, FastMath.copySign(new Binary64(1.0), zP.atan().getImaginary()).getReal(), 1.0e-15);
        final FieldComplex<Binary64> zM = build(0.8734729023516287, -0.0);
        final FieldComplex<Binary64> aM = build(0.717964439926383,  -0.0);
        Assertions.assertEquals(aM, zM.atan());
        Assertions.assertEquals(-1.0, FastMath.copySign(new Binary64(1.0), zM.atan().getImaginary()).getReal(), 1.0e-15);
    }

    @Test
    @Override
    public void testAtan2() {
        for (double x = -3; x < 3; x += 0.2) {
            for (double y = -3; y < 3; y += 0.2) {
                final FieldComplex<Binary64> z = build(x).atan2(build(y));
                final double  r = FastMath.atan2(x, y);
                checkRelative(r, build(MathUtils.normalizeAngle(z.getRealPart().getReal(), r), z.getImaginaryPart().getReal()));
            }
        }
    }

    @Test
    void testAtan2Complex() {
        for (double r1 : Arrays.asList(-3, 3)) {
            for (double i1 : Arrays.asList(-2, 0, 2)) {
                final FieldComplex<Binary64> c1 = build(r1, i1);
                for (double r2 : Arrays.asList(-1, 1)) {
                    for (double i2 : Arrays.asList(-5, 0, 5)) {
                        final FieldComplex<Binary64> c2 = build(r2, i2);
                        UnitTestUtils.customAssertEquals(c1.divide(c2), c1.atan2(c2).tan(), 1.0e-14);
                        final FieldComplex<Binary64> atan   = c1.divide(c2).atan();
                        final FieldComplex<Binary64> atan2  = c1.atan2(c2);
                        final double  deltaR = FastMath.abs(atan.getReal() - atan2.getReal()) / FastMath.PI;
                        assertTrue(FastMath.abs(deltaR - FastMath.rint(deltaR)) < 1.0e-14);
                        assertEquals(atan.getImaginaryPart().getReal(), atan2.getImaginaryPart().getReal(), 1.0e-14);
                    }
                }
            }
        }
    }

    @Test
    void testAtan2Real() {
        for (double r1 : Arrays.asList(-3, 3)) {
            final FieldComplex<Binary64> c1 = build(r1, 0);
            for (double r2 : Arrays.asList(-1, 1)) {
                final FieldComplex<Binary64> c2 = build(r2, 0);
                assertEquals(FastMath.atan2(r1, r2),
                                    MathUtils.normalizeAngle(c1.atan2(c2).getRealPart().getReal(), 0.0),
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
        FieldComplex<Binary64> z = build(3, 4);
        FieldComplex<Binary64> expected = build(-27.03495, -3.851153);
        UnitTestUtils.customAssertEquals(expected, z.cos(), 1.0e-5);
    }

    @Test
    void testCosNaN() {
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).cos().isNaN());
    }

    @Test
    void testCosInf() {
        UnitTestUtils.customAssertSame(infNegInf, oneInf.cos());
        UnitTestUtils.customAssertSame(infInf, oneNegInf.cos());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infOne.cos());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfOne.cos());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.cos());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.cos());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfInf.cos());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.cos());
    }

    @Test
    void testCoshComplex() {
        FieldComplex<Binary64> z = build(3, 4);
        FieldComplex<Binary64> expected = build(-6.58066, -7.58155);
        UnitTestUtils.customAssertEquals(expected, z.cosh(), 1.0e-5);
    }

    @Test
    void testCoshNaN() {
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).cosh().isNaN());
    }

    @Test
    void testCoshInf() {
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), oneInf.cosh());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), oneNegInf.cosh());
        UnitTestUtils.customAssertSame(infInf, infOne.cosh());
        UnitTestUtils.customAssertSame(infNegInf, negInfOne.cosh());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.cosh());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.cosh());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfInf.cosh());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.cosh());
    }

    @Test
    void testExpComplex() {
        FieldComplex<Binary64> z = build(3, 4);
        FieldComplex<Binary64> expected = build(-13.12878, -15.20078);
        UnitTestUtils.customAssertEquals(expected, z.exp(), 1.0e-5);
        UnitTestUtils.customAssertEquals(FieldComplex.getOne(Binary64Field.getInstance()),
                                         FieldComplex.getZero(Binary64Field.getInstance()).exp(), 10e-12);
        FieldComplex<Binary64> iPi = FieldComplex.getI(Binary64Field.getInstance()).multiply(build(FastMath.PI, 0));
        UnitTestUtils.customAssertEquals(FieldComplex.getOne(Binary64Field.getInstance()).negate(),
                                         iPi.exp(), 10e-12);
    }

    @Test
    void testExpNaN() {
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).exp().isNaN());
    }

    @Test
    void testExpInf1() {
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), oneInf.exp());
    }

    @Test
    void testExpInf2() {
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), oneNegInf.exp());
    }

    @Test
    void testExpInf3() {
        UnitTestUtils.customAssertSame(infInf, infOne.exp());
    }

    @Test
    void testExpInf4() {
        final FieldComplex<Binary64> exp = negInfOne.exp();
        UnitTestUtils.customAssertSame(FieldComplex.getZero(Binary64Field.getInstance()), exp);
    }

    @Test
    void testExpInf5() {
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.exp());
    }

    @Test
    void testExpInf6() {
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.exp());
    }

    @Test
    void testExpInf7() {
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfInf.exp());
    }

    @Test
    void testExpInf8() {
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.exp());
    }

    @Test
    void testExpM1() {
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), FieldComplex.getNaN(Binary64Field.getInstance()).expm1());
        final double testValue = FastMath.scalb(1.0, -30);
        assertEquals(FastMath.expm1(testValue), build(testValue, 0).expm1().getRealPart().getReal(), 1.0e-30);
        assertTrue(FastMath.expm1(testValue) - build(testValue).exp().subtract(1.0).getReal() > 4.0e-19);
        assertEquals(0.0, build(0, testValue).expm1().getRealPart().getReal(), 1.0e-30);
        assertEquals(0.0, build(0, testValue).expm1().getImaginaryPart().getReal(), 1.0e-30);
    }

    @Test
    void testLogComplex() {
        FieldComplex<Binary64> z = build(3, 4);
        FieldComplex<Binary64> expected = build(1.60944, 0.927295);
        UnitTestUtils.customAssertEquals(expected, z.log(), 1.0e-5);
    }

    @Test
    void testLogNaN() {
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).log().isNaN());
    }

    @Test
    void testLogInf() {
        UnitTestUtils.customAssertEquals(build(Double.POSITIVE_INFINITY, FastMath.PI / 2),
                                         oneInf.log(), 10e-12);
        UnitTestUtils.customAssertEquals(build(Double.POSITIVE_INFINITY, -FastMath.PI / 2),
                                         oneNegInf.log(), 10e-12);
        UnitTestUtils.customAssertEquals(infZero, infOne.log(), 10e-12);
        UnitTestUtils.customAssertEquals(build(Double.POSITIVE_INFINITY, FastMath.PI),
                                         negInfOne.log(), 10e-12);
        UnitTestUtils.customAssertEquals(build(Double.POSITIVE_INFINITY, FastMath.PI / 4),
                                         infInf.log(), 10e-12);
        UnitTestUtils.customAssertEquals(build(Double.POSITIVE_INFINITY, -FastMath.PI / 4),
                                         infNegInf.log(), 10e-12);
        UnitTestUtils.customAssertEquals(build(Double.POSITIVE_INFINITY, 3d * FastMath.PI / 4),
                                         negInfInf.log(), 10e-12);
        UnitTestUtils.customAssertEquals(build(Double.POSITIVE_INFINITY, - 3d * FastMath.PI / 4),
                                         negInfNegInf.log(), 10e-12);
    }

    @Test
    void testLogZero() {
        UnitTestUtils.customAssertSame(negInfZero, FieldComplex.getZero(Binary64Field.getInstance()).log());
    }

    @Test
    void testLog1P() {
        FieldComplex<Binary64> z = build(2, 4);
        FieldComplex<Binary64> expected = build(1.60944, 0.927295);
        UnitTestUtils.customAssertEquals(expected, z.log1p(), 1.0e-5);
    }

    @Test
    void testLog10Complex() {
        UnitTestUtils.customAssertEquals(build(2.0, 0.0), build(100, 0).log10(), 1.0e-15);
        UnitTestUtils.customAssertEquals(build(2.0, 0.5 * FastMath.PI / FastMath.log(10)), build(0, 100).log10(), 1.0e-15);
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

    public void testPowT() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            for (double y = 0.1; y < 4; y += 0.2) {
                if ( x < 0) {
                    // special case for Complex
                    assertTrue(Double.isNaN(FastMath.pow(x, y)));
                    assertFalse(build(x).pow(new Binary64(y)).isNaN());
                } else {
                    checkRelative(FastMath.pow(x, y), build(x).pow(new Binary64(y)));
                }
            }
        }
    }

    @Test
    void testPow() {
        FieldComplex<Binary64> x = build(3, 4);
        FieldComplex<Binary64> y = build(5, 6);
        FieldComplex<Binary64> expected = build(-1.860893, 11.83677);
        UnitTestUtils.customAssertEquals(expected, x.pow(y), 1.0e-5);
        UnitTestUtils.customAssertEquals(build(-46, 9).divide(2197), build(2, -3).pow(build(-3, 0)), 1.0e-15);
        UnitTestUtils.customAssertEquals(build(-1, 0).divide(8), build(-2, 0).pow(build(-3, 0)), 1.0e-15);
        UnitTestUtils.customAssertEquals(build(0, 2),
                                         build(-4, 0).pow(build(0.5, 0)),
                                         1.0e-15);
    }

    @Test
    void testPowNaNBase() {
        FieldComplex<Binary64> x = build(3, 4);
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).pow(x).isNaN());
    }

    @Test
    void testPowNaNExponent() {
        FieldComplex<Binary64> x = build(3, 4);
        assertTrue(x.pow(FieldComplex.getNaN(Binary64Field.getInstance())).isNaN());
    }

    @Test
    void testPowInf() {
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), FieldComplex.getOne(Binary64Field.getInstance()).pow(oneInf));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), FieldComplex.getOne(Binary64Field.getInstance()).pow(oneNegInf));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), FieldComplex.getOne(Binary64Field.getInstance()).pow(infOne));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), FieldComplex.getOne(Binary64Field.getInstance()).pow(infInf));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), FieldComplex.getOne(Binary64Field.getInstance()).pow(infNegInf));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), FieldComplex.getOne(Binary64Field.getInstance()).pow(negInfInf));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), FieldComplex.getOne(Binary64Field.getInstance()).pow(negInfNegInf));
        UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), infOne.pow(FieldComplex.getOne(Binary64Field.getInstance())));
        UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), negInfOne.pow(FieldComplex.getOne(Binary64Field.getInstance())));
        UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), infInf.pow(FieldComplex.getOne(Binary64Field.getInstance())));
        UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), infNegInf.pow(FieldComplex.getOne(Binary64Field.getInstance())));
        UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), negInfInf.pow(FieldComplex.getOne(Binary64Field.getInstance())));
        UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), negInfNegInf.pow(FieldComplex.getOne(Binary64Field.getInstance())));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.pow(infNegInf));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.pow(negInfNegInf));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.pow(infInf));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.pow(infNegInf));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.pow(negInfNegInf));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.pow(infInf));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.pow(infNegInf));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.pow(negInfNegInf));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.pow(infInf));
    }

    @Test
    void testPowZero() {
        UnitTestUtils.customAssertEquals(FieldComplex.getZero(Binary64Field.getInstance()),
                                         FieldComplex.getZero(Binary64Field.getInstance()).pow(FieldComplex.getOne(Binary64Field.getInstance())), 1.0e-12);
        UnitTestUtils.customAssertSame(FieldComplex.getOne(Binary64Field.getInstance()),
                                       FieldComplex.getZero(Binary64Field.getInstance()).pow(FieldComplex.getZero(Binary64Field.getInstance())));
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()),
                                       FieldComplex.getZero(Binary64Field.getInstance()).pow(FieldComplex.getI(Binary64Field.getInstance())));
        UnitTestUtils.customAssertEquals(FieldComplex.getOne(Binary64Field.getInstance()),
                                         FieldComplex.getOne(Binary64Field.getInstance()).pow(FieldComplex.getZero(Binary64Field.getInstance())), 10e-12);
        UnitTestUtils.customAssertEquals(FieldComplex.getOne(Binary64Field.getInstance()),
                                         FieldComplex.getI(Binary64Field.getInstance()).pow(FieldComplex.getZero(Binary64Field.getInstance())), 10e-12);
        UnitTestUtils.customAssertEquals(FieldComplex.getOne(Binary64Field.getInstance()),
                                         build(-1, 3).pow(FieldComplex.getZero(Binary64Field.getInstance())), 10e-12);
    }

    @Test
    void testZeroPow() {
        UnitTestUtils.customAssertEquals(FieldComplex.getZero(Binary64Field.getInstance()), FieldComplex.getZero(Binary64Field.getInstance()).pow(2.0), 1.0e-5);
    }

    @Test
    void testScalarPow() {
        FieldComplex<Binary64> x = build(3, 4);
        double yDouble = 5.0;
        FieldComplex<Binary64> yComplex = build(yDouble);
        assertEquals(x.pow(yComplex), x.pow(yDouble));
        assertEquals(x.pow(yComplex), x.pow(new Binary64(yDouble)));
        UnitTestUtils.customAssertEquals(FieldComplex.getOne(Binary64Field.getInstance()).negate(),
                                         FieldComplex.getOne(Binary64Field.getInstance()).negate().pow(0.5).pow(2),
                                         1.0e-15);
        UnitTestUtils.customAssertEquals(FieldComplex.getOne(Binary64Field.getInstance()).negate(),
                                         FieldComplex.getOne(Binary64Field.getInstance()).negate().pow(new Binary64(0.5)).pow(new Binary64(2)),
                                         1.0e-15);
        UnitTestUtils.customAssertEquals(build(2, 0), build(4, 0).pow(0.5), 1.0e-15);
        UnitTestUtils.customAssertEquals(build(2, 0), build(4, 0).pow(new Binary64(0.5)), 1.0e-15);
        UnitTestUtils.customAssertEquals(build(2, 0), build(4, 0).pow(build(0.5, 0)), 1.0e-15);
    }

    @Test
    void testScalarPowNaNBase() {
        FieldComplex<Binary64> x = FieldComplex.getNaN(Binary64Field.getInstance());
        double yDouble = 5.0;
        FieldComplex<Binary64> yComplex = build(yDouble);
        assertEquals(x.pow(yComplex), x.pow(yDouble));
        assertEquals(x.pow(yComplex), x.pow(new Binary64(yDouble)));
    }

    @Test
    void testScalarPowNaNExponent() {
        FieldComplex<Binary64> x = build(3, 4);
        double yDouble = Double.NaN;
        FieldComplex<Binary64> yComplex = build(yDouble);
        assertEquals(x.pow(yComplex), x.pow(yDouble));
        assertEquals(x.pow(yComplex), x.pow(new Binary64(yDouble)));
    }

    @Test
    void testScalarPowInf() {
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), FieldComplex.getOne(Binary64Field.getInstance()).pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), FieldComplex.getOne(Binary64Field.getInstance()).pow(new Binary64(Double.POSITIVE_INFINITY)));
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), FieldComplex.getOne(Binary64Field.getInstance()).pow(Double.NEGATIVE_INFINITY));
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), FieldComplex.getOne(Binary64Field.getInstance()).pow(new Binary64(Double.NEGATIVE_INFINITY)));
       UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), infOne.pow(1.0));
       UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), infOne.pow(new Binary64(1.0)));
       UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), negInfOne.pow(1.0));
       UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), negInfOne.pow(new Binary64(1.0)));
       UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), infInf.pow(1.0));
       UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), infInf.pow(new Binary64(1.0)));
       UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), infNegInf.pow(1.0));
       UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), infNegInf.pow(new Binary64(1.0)));
       UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), negInfInf.pow(10));
       UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), negInfInf.pow(new Binary64(10)));
       UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), negInfNegInf.pow(1.0));
       UnitTestUtils.customAssertSame(FieldComplex.getInf(Binary64Field.getInstance()), negInfNegInf.pow(new Binary64(1.0)));
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.pow(new Binary64(Double.POSITIVE_INFINITY)));
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.pow(new Binary64(Double.POSITIVE_INFINITY)));
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.pow(new Binary64(Double.POSITIVE_INFINITY)));
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.pow(Double.NEGATIVE_INFINITY));
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.pow(new Binary64(Double.NEGATIVE_INFINITY)));
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.pow(Double.NEGATIVE_INFINITY));
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.pow(new Binary64(Double.NEGATIVE_INFINITY)));
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.pow(Double.POSITIVE_INFINITY));
       UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.pow(new Binary64(Double.POSITIVE_INFINITY)));
   }

    @Test
    void testScalarPowZero() {
       UnitTestUtils.customAssertEquals(FieldComplex.getZero(Binary64Field.getInstance()), FieldComplex.getZero(Binary64Field.getInstance()).pow(1.0), 1.0e-12);
       UnitTestUtils.customAssertEquals(FieldComplex.getZero(Binary64Field.getInstance()), FieldComplex.getZero(Binary64Field.getInstance()).pow(new Binary64(1.0)), 1.0e-12);
       UnitTestUtils.customAssertSame(FieldComplex.getOne(Binary64Field.getInstance()), FieldComplex.getZero(Binary64Field.getInstance()).pow(0.0));
       UnitTestUtils.customAssertSame(FieldComplex.getOne(Binary64Field.getInstance()), FieldComplex.getZero(Binary64Field.getInstance()).pow(new Binary64(0.0)));
       UnitTestUtils.customAssertEquals(FieldComplex.getOne(Binary64Field.getInstance()), FieldComplex.getOne(Binary64Field.getInstance()).pow(0.0), 10e-12);
       UnitTestUtils.customAssertEquals(FieldComplex.getOne(Binary64Field.getInstance()), FieldComplex.getOne(Binary64Field.getInstance()).pow(new Binary64(0.0)), 10e-12);
       UnitTestUtils.customAssertEquals(FieldComplex.getOne(Binary64Field.getInstance()), FieldComplex.getI(Binary64Field.getInstance()).pow(0.0), 10e-12);
       UnitTestUtils.customAssertEquals(FieldComplex.getOne(Binary64Field.getInstance()), FieldComplex.getI(Binary64Field.getInstance()).pow(new Binary64(0.0)), 10e-12);
       UnitTestUtils.customAssertEquals(FieldComplex.getOne(Binary64Field.getInstance()), build(-1, 3).pow(0.0), 10e-12);
       UnitTestUtils.customAssertEquals(FieldComplex.getOne(Binary64Field.getInstance()), build(-1, 3).pow(new Binary64(0.0)), 10e-12);
   }

    @Test
    void testpowNull() {
        assertThrows(NullArgumentException.class, () -> {
            FieldComplex.getOne(Binary64Field.getInstance()).pow((FieldComplex<Binary64>) null);
        });
    }

    @Test
    void testSinComplex() {
        FieldComplex<Binary64> z = build(3, 4);
        FieldComplex<Binary64> expected = build(3.853738, -27.01681);
        UnitTestUtils.customAssertEquals(expected, z.sin(), 1.0e-5);
    }

    @Test
    void testSinInf() {
        UnitTestUtils.customAssertSame(infInf, oneInf.sin());
        UnitTestUtils.customAssertSame(infNegInf, oneNegInf.sin());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infOne.sin());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfOne.sin());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.sin());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.sin());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfInf.sin());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.sin());
    }

    @Test
    void testSinNaN() {
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).sin().isNaN());
    }

    @Test
    void testSinhComplex() {
        FieldComplex<Binary64> z = build(3, 4);
        FieldComplex<Binary64> expected = build(-6.54812, -7.61923);
        UnitTestUtils.customAssertEquals(expected, z.sinh(), 1.0e-5);
    }

    @Test
    void testSinhNaN() {
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).sinh().isNaN());
    }

    @Test
    void testSinhInf() {
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), oneInf.sinh());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), oneNegInf.sinh());
        UnitTestUtils.customAssertSame(infInf, infOne.sinh());
        UnitTestUtils.customAssertSame(negInfInf, negInfOne.sinh());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.sinh());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.sinh());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfInf.sinh());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.sinh());
    }

    @Test
    void testAsinhComplex() {
        for (double x = -2; x <= 2; x += 0.125) {
            for (double y = -2; y <= 2; y += 0.125) {
                final FieldComplex<Binary64> z = build(x, y);
                UnitTestUtils.customAssertEquals(z, z.asinh().sinh(), 1.0e-14);
            }
        }
    }

    @Test
    void testAsinhBranchCuts() {
        UnitTestUtils.customAssertEquals(build(FastMath.log(2 + FastMath.sqrt(3)), 0.5 * FastMath.PI),
                                         build(+0.0, 2.0).asinh(),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(build(-FastMath.log(2 + FastMath.sqrt(3)), 0.5 * FastMath.PI),
                                         build(-0.0, 2.0).asinh(),
                                         1.0e-14);
    }

    @Test
    void testAcoshComplex() {
        for (double x = -2; x <= 2; x += 0.125) {
            for (double y = -2; y <= 2; y += 0.125) {
                final FieldComplex<Binary64> z = build(x, y);
                UnitTestUtils.customAssertEquals(z, z.acosh().cosh(), 1.0e-14);
            }
        }
    }

    @Test
    void testAcoshBranchCuts() {
        UnitTestUtils.customAssertEquals(build(FastMath.log(2 + FastMath.sqrt(3)), +FastMath.PI),
                                         build(-2.0, +0.0).acosh(),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(build(FastMath.log(2 + FastMath.sqrt(3)), -FastMath.PI),
                                         build(-2.0, -0.0).acosh(),
                                         1.0e-14);
    }

    @Test
    void testAtanhComplex() {
        for (double x = -2; x <= 2; x += 0.125) {
            for (double y = -2; y <= 2; y += 0.125) {
                final FieldComplex<Binary64> z = build(x, y);
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
        UnitTestUtils.customAssertEquals(build(-0.5 * FastMath.log(3), +0.5 * FastMath.PI),
                                         build(-2.0, +0.0).atanh(),
                                         1.0e-14);
        UnitTestUtils.customAssertEquals(build(-0.5 * FastMath.log(3), -0.5 * FastMath.PI),
                                         build(-2.0, -0.0).atanh(),
                                         1.0e-14);
    }

    @Test
    void testSqrtRealPositive() {
        FieldComplex<Binary64> z = build(3, 4);
        FieldComplex<Binary64> expected = build(2, 1);
        UnitTestUtils.customAssertEquals(expected, z.sqrt(), 1.0e-5);
    }

    @Test
    void testSqrtRealZero() {
        FieldComplex<Binary64> z = build(0.0, 4);
        FieldComplex<Binary64> expected = build(1.41421, 1.41421);
        UnitTestUtils.customAssertEquals(expected, z.sqrt(), 1.0e-5);
    }

    @Test
    void testSqrtZero() {
        UnitTestUtils.customAssertEquals(FieldComplex.getZero(Binary64Field.getInstance()), FieldComplex.getZero(Binary64Field.getInstance()).sqrt(), 1.0e-15);
    }

    @Test
    void testSqrtRealNegative() {
        FieldComplex<Binary64> z = build(-3.0, 4);
        FieldComplex<Binary64> expected = build(1, 2);
        UnitTestUtils.customAssertEquals(expected, z.sqrt(), 1.0e-5);
    }

    @Test
    void testSqrtImaginaryZero() {
        FieldComplex<Binary64> z = build(-3.0, 0.0);
        FieldComplex<Binary64> expected = build(0.0, 1.73205);
        UnitTestUtils.customAssertEquals(expected, z.sqrt(), 1.0e-5);
    }

    @Test
    void testSqrtImaginaryNegative() {
        FieldComplex<Binary64> z = build(-3.0, -4.0);
        FieldComplex<Binary64> expected = build(1.0, -2.0);
        UnitTestUtils.customAssertEquals(expected, z.sqrt(), 1.0e-5);
    }

    @Test
    void testSqrtPolar() {
        Binary64 r = Binary64.ONE;
        for (int i = 0; i < 5; i++) {
            r = r.add(i);
            Binary64 theta = Binary64.ZERO;
            for (int j =0; j < 11; j++) {
                theta = theta.add(FastMath.PI / 12);
                FieldComplex<Binary64> z = ComplexUtils.polar2Complex(r, theta);
                FieldComplex<Binary64> sqrtz = ComplexUtils.polar2Complex(FastMath.sqrt(r), theta.half());
                UnitTestUtils.customAssertEquals(sqrtz, z.sqrt(), 10e-12);
            }
        }
    }

    @Test
    void testSqrtNaN() {
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).sqrt().isNaN());
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
        FieldComplex<Binary64> z = build(3, 4);
        FieldComplex<Binary64> expected = build(4.08033, -2.94094);
        UnitTestUtils.customAssertEquals(expected, z.sqrt1z(), 1.0e-5);
    }

    @Test
    void testSqrt1zNaN() {
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).sqrt1z().isNaN());
    }

    @Test
    @Override
    public void testCbrt() {
        for (double x = -0.9; x < 0.9; x += 0.05) {
            if ( x < 0) {
                // special case for Complex
                assertTrue(FastMath.cbrt(x) < 0);
                assertEquals(FastMath.PI / 3, build(x).cbrt().getArgument().getReal(), 1.0e-15);
            } else {
                checkRelative(FastMath.cbrt(x), build(x).cbrt());
            }
        }
    }

    @Test
    void testCbrtComplex() {
        FieldComplex<Binary64> z = build(15, 2);
        UnitTestUtils.customAssertEquals(z, z.square().multiply(z).cbrt(), 1.0e-14);
        FieldComplex<Binary64> branchCutPlus = build(-8.0, +0.0);
        FieldComplex<Binary64> cbrtPlus = branchCutPlus.cbrt();
        UnitTestUtils.customAssertEquals(branchCutPlus, cbrtPlus.multiply(cbrtPlus).multiply(cbrtPlus), 1.0e-14);
        assertEquals(1.0, cbrtPlus.getRealPart().getReal(), 1.0e-15);
        assertEquals(FastMath.sqrt(3.0), cbrtPlus.getImaginaryPart().getReal(), 1.0e-15);
        FieldComplex<Binary64> branchCutMinus = build(-8.0, -0.0);
        FieldComplex<Binary64> cbrtMinus = branchCutMinus.cbrt();
        UnitTestUtils.customAssertEquals(branchCutMinus, cbrtMinus.multiply(cbrtMinus).multiply(cbrtMinus), 1.0e-14);
        assertEquals(1.0, cbrtMinus.getRealPart().getReal(), 1.0e-15);
        assertEquals(-FastMath.sqrt(3.0), cbrtMinus.getImaginaryPart().getReal(), 1.0e-15);
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
                    assertEquals(FastMath.PI / n, build(x).rootN(n).getArgument().getReal(), 1.0e-15);
                } else {
                    checkRelative(FastMath.pow(x, 1.0 / n), build(x).rootN(n));
                }
            }
        }
    }

    @Test
    void testRootNComplex() {
        FieldComplex<Binary64> z = build(15, 2);
        UnitTestUtils.customAssertEquals(z, z.square().multiply(z).rootN(3), 1.0e-14);
        FieldComplex<Binary64> branchCutPlus = build(-8.0, +0.0);
        FieldComplex<Binary64> cbrtPlus = branchCutPlus.rootN(3);
        UnitTestUtils.customAssertEquals(branchCutPlus, cbrtPlus.multiply(cbrtPlus).multiply(cbrtPlus), 1.0e-14);
        assertEquals(1.0, cbrtPlus.getRealPart().getReal(), 1.0e-15);
        assertEquals(FastMath.sqrt(3.0), cbrtPlus.getImaginaryPart().getReal(), 1.0e-15);
        FieldComplex<Binary64> branchCutMinus = build(-8.0, -0.0);
        FieldComplex<Binary64> cbrtMinus = branchCutMinus.rootN(3);
        UnitTestUtils.customAssertEquals(branchCutMinus, cbrtMinus.multiply(cbrtMinus).multiply(cbrtMinus), 1.0e-14);
        assertEquals(1.0, cbrtMinus.getRealPart().getReal(), 1.0e-15);
        assertEquals(-FastMath.sqrt(3.0), cbrtMinus.getImaginaryPart().getReal(), 1.0e-15);
    }

    @Test
    void testTanComplex() {
        FieldComplex<Binary64> z = build(3, 4);
        FieldComplex<Binary64> expected = build(-0.000187346, 0.999356);
        UnitTestUtils.customAssertEquals(expected, z.tan(), 1.0e-5);
        /* Check that no overflow occurs (MATH-722) */
        FieldComplex<Binary64> actual = build(3.0, 1E10).tan();
        expected = build(0, 1);
        UnitTestUtils.customAssertEquals(expected, actual, 1.0e-5);
        actual = build(3.0, -1E10).tan();
        expected = build(0, -1);
        UnitTestUtils.customAssertEquals(expected, actual, 1.0e-5);
    }

    @Test
    void testTanNaN() {
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).tan().isNaN());
    }

    @Test
    void testTanInf() {
        UnitTestUtils.customAssertSame(FieldComplex.valueOf(new Binary64(0.0), new Binary64(1.0)), oneInf.tan());
        UnitTestUtils.customAssertSame(FieldComplex.valueOf(new Binary64(0.0), new Binary64(-1.0)), oneNegInf.tan());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infOne.tan());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfOne.tan());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.tan());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.tan());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfInf.tan());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.tan());
    }

    @Test
    void testTanCritical() {
        UnitTestUtils.customAssertSame(infNaN, build(MathUtils.SEMI_PI, 0).tan());
        UnitTestUtils.customAssertSame(negInfNaN, build(-MathUtils.SEMI_PI, 0).tan());
    }

    @Test
    void testTanhComplex() {
        FieldComplex<Binary64> z = build(3, 4);
        FieldComplex<Binary64> expected = build(1.00071, 0.00490826);
        UnitTestUtils.customAssertEquals(expected, z.tanh(), 1.0e-5);
        /* Check that no overflow occurs (MATH-722) */
        FieldComplex<Binary64> actual = build(1E10, 3.0).tanh();
        expected = build(1, 0);
        UnitTestUtils.customAssertEquals(expected, actual, 1.0e-5);
        actual = build(-1E10, 3.0).tanh();
        expected = build(-1, 0);
        UnitTestUtils.customAssertEquals(expected, actual, 1.0e-5);
    }

    @Test
    void testTanhNaN() {
        assertTrue(FieldComplex.getNaN(Binary64Field.getInstance()).tanh().isNaN());
    }

    @Test
    void testTanhInf() {
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), oneInf.tanh());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), oneNegInf.tanh());
        UnitTestUtils.customAssertSame(FieldComplex.valueOf(new Binary64(1.0), new Binary64(0.0)), infOne.tanh());
        UnitTestUtils.customAssertSame(FieldComplex.valueOf(new Binary64(-1.0), new Binary64(0.0)), negInfOne.tanh());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infInf.tanh());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), infNegInf.tanh());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfInf.tanh());
        UnitTestUtils.customAssertSame(FieldComplex.getNaN(Binary64Field.getInstance()), negInfNegInf.tanh());
    }

    @Test
    void testTanhCritical() {
        UnitTestUtils.customAssertSame(nanInf, build(0, MathUtils.SEMI_PI).tanh());
    }

    /** test issue MATH-221 */
    @Test
    void testMath221() {
        assertTrue(FieldComplex.equals(build(0,-1), build(0,1).multiply(build(-1,0))));
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
        FieldComplex<Binary64> z = build(-2,2);
        // The List holding all third roots
        List<FieldComplex<Binary64>> thirdRootsOfZ = z.nthRoot(3);
        // Returned Collection must not be empty!
        assertEquals(3, thirdRootsOfZ.size());
        // test z_0
        assertEquals(1.0,                  thirdRootsOfZ.get(0).getRealPart().getReal(),      1.0e-5);
        assertEquals(1.0,                  thirdRootsOfZ.get(0).getImaginaryPart().getReal(), 1.0e-5);
        // test z_1
        assertEquals(-1.3660254037844386,  thirdRootsOfZ.get(1).getRealPart().getReal(),      1.0e-5);
        assertEquals(0.36602540378443843,  thirdRootsOfZ.get(1).getImaginaryPart().getReal(), 1.0e-5);
        // test z_2
        assertEquals(0.366025403784439,    thirdRootsOfZ.get(2).getRealPart().getReal(),      1.0e-5);
        assertEquals(-1.3660254037844384,  thirdRootsOfZ.get(2).getImaginaryPart().getReal(), 1.0e-5);
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
        FieldComplex<Binary64> z = build(5,-2);
        // The List holding all fourth roots
        List<FieldComplex<Binary64>> fourthRootsOfZ = z.nthRoot(4);
        // Returned Collection must not be empty!
        assertEquals(4, fourthRootsOfZ.size());
        // test z_0
        assertEquals(1.5164629308487783,     fourthRootsOfZ.get(0).getRealPart().getReal(),      1.0e-5);
        assertEquals(-0.14469266210702247,   fourthRootsOfZ.get(0).getImaginaryPart().getReal(), 1.0e-5);
        // test z_1
        assertEquals(0.14469266210702256,    fourthRootsOfZ.get(1).getRealPart().getReal(),      1.0e-5);
        assertEquals(1.5164629308487783,     fourthRootsOfZ.get(1).getImaginaryPart().getReal(), 1.0e-5);
        // test z_2
        assertEquals(-1.5164629308487783,    fourthRootsOfZ.get(2).getRealPart().getReal(),      1.0e-5);
        assertEquals(0.14469266210702267,    fourthRootsOfZ.get(2).getImaginaryPart().getReal(), 1.0e-5);
        // test z_3
        assertEquals(-0.14469266210702275,   fourthRootsOfZ.get(3).getRealPart().getReal(),      1.0e-5);
        assertEquals(-1.5164629308487783,    fourthRootsOfZ.get(3).getImaginaryPart().getReal(), 1.0e-5);
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
        FieldComplex<Binary64> z = build(8,0);
        // The List holding all third roots
        List<FieldComplex<Binary64>> thirdRootsOfZ = z.nthRoot(3);
        // Returned Collection must not be empty!
        assertEquals(3, thirdRootsOfZ.size());
        // test z_0
        assertEquals(2.0,                thirdRootsOfZ.get(0).getRealPart().getReal(),      1.0e-5);
        assertEquals(0.0,                thirdRootsOfZ.get(0).getImaginaryPart().getReal(), 1.0e-5);
        // test z_1
        assertEquals(-1.0,               thirdRootsOfZ.get(1).getRealPart().getReal(),      1.0e-5);
        assertEquals(1.7320508075688774, thirdRootsOfZ.get(1).getImaginaryPart().getReal(), 1.0e-5);
        // test z_2
        assertEquals(-1.0,               thirdRootsOfZ.get(2).getRealPart().getReal(),      1.0e-5);
        assertEquals(-1.732050807568877, thirdRootsOfZ.get(2).getImaginaryPart().getReal(), 1.0e-5);
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
        FieldComplex<Binary64> z = build(0,2);
        // The List holding all third roots
        List<FieldComplex<Binary64>> thirdRootsOfZ = z.nthRoot(3);
        // Returned Collection must not be empty!
        assertEquals(3, thirdRootsOfZ.size());
        // test z_0
        assertEquals(1.0911236359717216,      thirdRootsOfZ.get(0).getRealPart().getReal(),      1.0e-5);
        assertEquals(0.6299605249474365,      thirdRootsOfZ.get(0).getImaginaryPart().getReal(), 1.0e-5);
        // test z_1
        assertEquals(-1.0911236359717216,     thirdRootsOfZ.get(1).getRealPart().getReal(),      1.0e-5);
        assertEquals(0.6299605249474365,      thirdRootsOfZ.get(1).getImaginaryPart().getReal(), 1.0e-5);
        // test z_2
        assertEquals(-2.3144374213981936E-16, thirdRootsOfZ.get(2).getRealPart().getReal(),      1.0e-5);
        assertEquals(-1.2599210498948732,     thirdRootsOfZ.get(2).getImaginaryPart().getReal(), 1.0e-5);
    }

    /**
     * Test cornercases with NaN and Infinity.
     */
    @Test
    void testNthRoot_cornercase_NAN_Inf() {
        // NaN + finite -> NaN
        List<FieldComplex<Binary64>> roots = oneNaN.nthRoot(3);
        assertEquals(1,roots.size());
        assertEquals(FieldComplex.getNaN(Binary64Field.getInstance()), roots.get(0));

        roots = nanZero.nthRoot(3);
        assertEquals(1,roots.size());
        assertEquals(FieldComplex.getNaN(Binary64Field.getInstance()), roots.get(0));

        // NaN + infinite -> NaN
        roots = nanInf.nthRoot(3);
        assertEquals(1,roots.size());
        assertEquals(FieldComplex.getNaN(Binary64Field.getInstance()), roots.get(0));

        // finite + infinite -> Inf
        roots = oneInf.nthRoot(3);
        assertEquals(1,roots.size());
        assertEquals(FieldComplex.getInf(Binary64Field.getInstance()), roots.get(0));

        // infinite + infinite -> Inf
        roots = negInfInf.nthRoot(3);
        assertEquals(1,roots.size());
        assertEquals(FieldComplex.getInf(Binary64Field.getInstance()), roots.get(0));
    }

    @Test
    void testNthRootError() {
        try {
            FieldComplex.getOne(Binary64Field.getInstance()).nthRoot(-1);
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
        assertFalse(build(Double.NaN, imaginary).isMathematicalInteger());
        assertFalse(build(Double.POSITIVE_INFINITY, imaginary).isMathematicalInteger());
        assertFalse(build(Double.NEGATIVE_INFINITY, imaginary).isMathematicalInteger());
        assertFalse(build(Double.MIN_NORMAL, imaginary).isMathematicalInteger());
        assertFalse(build(Double.MIN_VALUE, imaginary).isMathematicalInteger());

        assertEquals(expectedForInteger, build(-0.0, imaginary).isMathematicalInteger());
        assertEquals(expectedForInteger, build(+0.0, imaginary).isMathematicalInteger());

        for (int i = -1000; i < 1000; ++i) {
            final double d = i;
            assertEquals(expectedForInteger, build(d, imaginary).isMathematicalInteger());
            assertFalse(build(FastMath.nextAfter(d, Double.POSITIVE_INFINITY), imaginary).isMathematicalInteger());
            assertFalse(build(FastMath.nextAfter(d, Double.NEGATIVE_INFINITY), imaginary).isMathematicalInteger());
        }

        double minNoFractional = 0x1l << 52;
        assertEquals(expectedForInteger, build(minNoFractional, imaginary).isMathematicalInteger());
        assertFalse(build(minNoFractional - 0.5, imaginary).isMathematicalInteger());
        assertEquals(expectedForInteger, build(minNoFractional + 0.5, imaginary).isMathematicalInteger());

    }

    /**
     * Test standard values
     */
    @Test
    void testGetArgument() {
        FieldComplex<Binary64> z = build(1, 0);
        assertEquals(0.0, z.getArgument().getReal(), 1.0e-12);

        z = build(1, 1);
        assertEquals(FastMath.PI/4, z.getArgument().getReal(), 1.0e-12);

        z = build(0, 1);
        assertEquals(FastMath.PI/2, z.getArgument().getReal(), 1.0e-12);

        z = build(-1, 1);
        assertEquals(3 * FastMath.PI/4, z.getArgument().getReal(), 1.0e-12);

        z = build(-1, 0);
        assertEquals(FastMath.PI, z.getArgument().getReal(), 1.0e-12);

        z = build(-1, -1);
        assertEquals(-3 * FastMath.PI/4, z.getArgument().getReal(), 1.0e-12);

        z = build(0, -1);
        assertEquals(-FastMath.PI/2, z.getArgument().getReal(), 1.0e-12);

        z = build(1, -1);
        assertEquals(-FastMath.PI/4, z.getArgument().getReal(), 1.0e-12);

    }

    /**
     * Verify atan2-style handling of infinite parts
     */
    @Test
    void testGetArgumentInf() {
        assertEquals(FastMath.PI/4, infInf.getArgument().getReal(), 1.0e-12);
        assertEquals(FastMath.PI/2, oneInf.getArgument().getReal(), 1.0e-12);
        assertEquals(0.0, infOne.getArgument().getReal(), 1.0e-12);
        assertEquals(FastMath.PI/2, zeroInf.getArgument().getReal(), 1.0e-12);
        assertEquals(0.0, infZero.getArgument().getReal(), 1.0e-12);
        assertEquals(FastMath.PI, negInfOne.getArgument().getReal(), 1.0e-12);
        assertEquals(-3.0*FastMath.PI/4, negInfNegInf.getArgument().getReal(), 1.0e-12);
        assertEquals(-FastMath.PI/2, oneNegInf.getArgument().getReal(), 1.0e-12);
    }

    /**
     * Verify that either part NaN results in NaN
     */
    @Test
    void testGetArgumentNaN() {
        assertTrue(Double.isNaN(nanZero.getArgument().getReal()));
        assertTrue(Double.isNaN(zeroNaN.getArgument().getReal()));
        assertTrue(Double.isNaN(FieldComplex.getNaN(Binary64Field.getInstance()).getArgument().getReal()));
    }

    @Test
    void testValueOf() {
        assertEquals(2.0, FieldComplex.valueOf(new Binary64(2.0)).getRealPart().getReal(), 1.0e-15);
        assertEquals(0.0, FieldComplex.valueOf(new Binary64(2.0)).getImaginaryPart().getReal(), 1.0e-15);
        assertTrue(FieldComplex.valueOf(new Binary64(Double.NaN)).isNaN());
        assertTrue(FieldComplex.valueOf(new Binary64(Double.POSITIVE_INFINITY)).isInfinite());
        assertEquals( 2.0, FieldComplex.valueOf(new Binary64(2.0), new Binary64(-1.0)).getRealPart().getReal(), 1.0e-15);
        assertEquals(-1.0, FieldComplex.valueOf(new Binary64(2.0), new Binary64(-1.0)).getImaginaryPart().getReal(), 1.0e-15);
        assertTrue(FieldComplex.valueOf(new Binary64(Double.NaN), new Binary64(0.0)).isNaN());
        assertTrue(FieldComplex.valueOf(new Binary64(Double.POSITIVE_INFINITY), new Binary64(0.0)).isInfinite());
        assertTrue(FieldComplex.valueOf(new Binary64(Double.NaN), new Binary64(-1.0)).isNaN());
        assertTrue(FieldComplex.valueOf(new Binary64(Double.POSITIVE_INFINITY), new Binary64(-1.0)).isInfinite());
        assertTrue(FieldComplex.valueOf(new Binary64(0.0), new Binary64(Double.NaN)).isNaN());
        assertTrue(FieldComplex.valueOf(new Binary64(0.0), new Binary64(Double.POSITIVE_INFINITY)).isInfinite());
        assertTrue(FieldComplex.valueOf(new Binary64(-1.0), new Binary64(Double.NaN)).isNaN());
        assertTrue(FieldComplex.valueOf(new Binary64(-1.0), new Binary64(Double.POSITIVE_INFINITY)).isInfinite());
    }

    @Test
    void testField() {
        assertEquals(FieldComplexField.getField(Binary64Field.getInstance()),
                            FieldComplex.getZero(Binary64Field.getInstance()).getField());
    }

    @Test
    void testToString() {
        assertEquals("(1.0, -2.0)", build(1, -2).toString());
    }

    @Test
    void testScalbComplex() {
        assertEquals(0.125,  build(2.0, 1.0).scalb(-4).getRealPart().getReal(), 1.0e-15);
        assertEquals(0.0625, build(2.0, 1.0).scalb(-4).getImaginaryPart().getReal(), 1.0e-15);
    }

    @Test
    void testHypotComplex() {
        assertEquals(5.8269600298808519855, build(3, 4).hypot(build(5, 6)).getRealPart().getReal(), 1.0e-15);
        assertEquals(7.2078750814528590485, build(3, 4).hypot(build(5, 6)).getImaginaryPart().getReal(), 1.0e-15);
    }

    @Test
    void testCeilComplex() {
        for (double x = -3.9; x < 3.9; x += 0.05) {
            for (double y = -3.9; y < 3.9; y += 0.05) {
                final FieldComplex<Binary64> z = build(x, y);
                assertEquals(FastMath.ceil(x), z.ceil().getRealPart().getReal(), 1.0e-15);
                assertEquals(FastMath.ceil(y), z.ceil().getImaginaryPart().getReal(), 1.0e-15);
            }
        }
    }

    @Test
    void testFloorComplex() {
        for (double x = -3.9; x < 3.9; x += 0.05) {
            for (double y = -3.9; y < 3.9; y += 0.05) {
                final FieldComplex<Binary64> z = build(x, y);
                assertEquals(FastMath.floor(x), z.floor().getRealPart().getReal(), 1.0e-15);
                assertEquals(FastMath.floor(y), z.floor().getImaginaryPart().getReal(), 1.0e-15);
            }
        }
    }

    @Test
    void testRintComplex() {
        for (double x = -3.9; x < 3.9; x += 0.05) {
            for (double y = -3.9; y < 3.9; y += 0.05) {
                final FieldComplex<Binary64> z = build(x, y);
                assertEquals(FastMath.rint(x), z.rint().getRealPart().getReal(), 1.0e-15);
                assertEquals(FastMath.rint(y), z.rint().getImaginaryPart().getReal(), 1.0e-15);
            }
        }
    }

    @Test
    void testRemainderComplexComplex() {
        for (double x1 = -3.9; x1 < 3.9; x1 += 0.125) {
            for (double y1 = -3.9; y1 < 3.9; y1 += 0.125) {
                final FieldComplex<Binary64> z1 = build(x1, y1);
                for (double x2 = -3.92; x2 < 3.9; x2 += 0.125) {
                    for (double y2 = -3.92; y2 < 3.9; y2 += 0.125) {
                        final FieldComplex<Binary64> z2 = build(x2, y2);
                        final FieldComplex<Binary64> r  = z1.remainder(z2);
                        final FieldComplex<Binary64> q  = z1.subtract(r).divide(z2);
                        assertTrue(r.norm() <= z2.norm());
                        assertEquals(FastMath.rint(q.getRealPart().getReal()), q.getRealPart().getReal(), 2.0e-14);
                        assertEquals(FastMath.rint(q.getImaginaryPart().getReal()), q.getImaginaryPart().getReal(), 2.0e-14);
                    }
                }
            }
        }
    }

    @Test
    void testRemainderComplexDouble() {
        for (double x1 = -3.9; x1 < 3.9; x1 += 0.125) {
            for (double y1 = -3.9; y1 < 3.9; y1 += 0.125) {
                final FieldComplex<Binary64> z1 = build(x1, y1);
                for (double a = -3.92; a < 3.9; a += 0.125) {
                        final FieldComplex<Binary64> r  = z1.remainder(a);
                        final FieldComplex<Binary64> q  = z1.subtract(r).divide(a);
                        assertTrue(r.norm() <= FastMath.abs(a));
                        assertEquals(FastMath.rint(q.getRealPart().getReal()), q.getRealPart().getReal(), 2.0e-14);
                        assertEquals(FastMath.rint(q.getImaginaryPart().getReal()), q.getImaginaryPart().getReal(), 2.0e-14);
                }
            }
        }
    }

    @Test
    void testRemainderAxKr() {
        checkRemainder(build(14, -5), build(3, 4), build(-1.0,  0.0));
        checkRemainder(build(26, 120), build(37, 226), build(-11.0, -106.0));
        checkRemainder(build(9.4, 6), build(1.0, 1.0), build(-0.6, 0.0));
        checkRemainder(build(-5.89, 0.33), build(2.4, -0.123), build(-1.09, 0.084));
    }

    private void checkRemainder(final FieldComplex<Binary64> c1, final FieldComplex<Binary64> c2, final FieldComplex<Binary64> expectedRemainder) {

        final FieldComplex<Binary64> remainder = c1.remainder(c2);
        assertEquals(expectedRemainder.getRealPart().getReal(),      remainder.getRealPart().getReal(),      1.0e-15);
        assertEquals(expectedRemainder.getImaginaryPart().getReal(), remainder.getImaginaryPart().getReal(), 1.0e-15);

        final FieldComplex<Binary64> crossCheck = c1.subtract(remainder).divide(c2);
        assertTrue(Precision.isMathematicalInteger(crossCheck.getRealPart().getReal()));
        assertTrue(Precision.isMathematicalInteger(crossCheck.getImaginaryPart().getReal()));

    }

    @Test
    void testCopySignFieldComplex() {
        for (double x1 = -3.9; x1 < 3.9; x1 += 0.08) {
            for (double y1 = -3.9; y1 < 3.9; y1 += 0.08) {
                final FieldComplex<Binary64> z1 = build(x1, y1);
                for (double x2 = -3.9; x2 < 3.9; x2 += 0.08) {
                    for (double y2 = -3.9; y2 < 3.9; y2 += 0.08) {
                        final FieldComplex<Binary64> z2 = build(x2, y2);
                        assertEquals(FastMath.copySign(x1, x2), z1.copySign(z2).getRealPart().getReal(), 1.0e-15);
                        assertEquals(FastMath.copySign(y1, y2), z1.copySign(z2).getImaginaryPart().getReal(), 1.0e-15);
                    }
                }
            }
        }
    }

    @Test
    void testCopySignDoubleComplex() {
        for (double x1 = -3.9; x1 < 3.9; x1 += 0.05) {
            for (double y1 = -3.9; y1 < 3.9; y1 += 0.05) {
                final FieldComplex<Binary64> z1 = build(x1, y1);
                for (double r = -3.9; r < 3.9; r += 0.05) {
                    assertEquals(FastMath.copySign(x1, r), z1.copySign(r).getRealPart().getReal(), 1.0e-15);
                    assertEquals(FastMath.copySign(y1, r), z1.copySign(r).getImaginaryPart().getReal(), 1.0e-15);
                }
            }
        }
    }

    @Test
    void testSignumComplex() {
        for (double x = -3.9; x < 3.9; x += 0.05) {
            for (double y = -3.9; y < 3.9; y += 0.05) {
                final FieldComplex<Binary64> z = build(x, y);
                assertEquals(1.0, z.sign().norm(), 1.0e-15);
                assertEquals(FastMath.copySign(1, FastMath.signum(x)), FastMath.copySign(Binary64.ONE, z.sign().getRealPart()).getReal(), 1.0e-15);
                assertEquals(FastMath.copySign(1, FastMath.signum(y)), FastMath.copySign(Binary64.ONE, z.sign().getImaginaryPart()).getReal(), 1.0e-15);
            }
        }
        assertTrue(Complex.NaN.sign().isNaN());
        for (int sR : Arrays.asList(-1, +1)) {
            for (int sI : Arrays.asList(-1, +1)) {
                FieldComplex<Binary64> z = build(FastMath.copySign(0, sR), FastMath.copySign(0, sI));
                assertTrue(z.isZero());
                FieldComplex<Binary64> zSign = z.sign();
                assertTrue(zSign.isZero());
                assertEquals(sR, FastMath.copySign(Binary64.ONE, zSign.getRealPart()).getReal(), 1.0e-15);
                assertEquals(sI, FastMath.copySign(Binary64.ONE, zSign.getImaginaryPart()).getReal(), 1.0e-15);
            }
        }
    }

    @Test
    void testLinearCombination1() {
        final FieldComplex<Binary64>[] a = MathArrays.buildArray(build(0.0).getField(), 2);
        a[0] = build(-1321008684645961.0 / 268435456.0, +5774608829631843.0 / 268435456.0);
        a[1] = build(-7645843051051357.0 / 8589934592.0, 0.0);
        final FieldComplex<Binary64>[] b = MathArrays.buildArray(build(0.0).getField(), 2);
        b[0] = build(-5712344449280879.0 / 2097152.0, -4550117129121957.0 / 2097152.0);
        b[1] = build(8846951984510141.0 / 131072.0, 0.0);

        final FieldComplex<Binary64> abSumInline = FieldComplex.getZero(Binary64Field.getInstance()).linearCombination(a[0], b[0],
                                                                                                                         a[1], b[1]);
        final FieldComplex<Binary64> abSumArray = FieldComplex.getZero(Binary64Field.getInstance()).linearCombination(a, b);

        UnitTestUtils.customAssertEquals(abSumInline, abSumArray, 0);
        UnitTestUtils.customAssertEquals(-1.8551294182586248737720779899, abSumInline.getRealPart().getReal(), 1.0e-15);

        final FieldComplex<Binary64> naive = a[0].multiply(b[0]).add(a[1].multiply(b[1]));
        assertTrue(naive.subtract(abSumInline).norm() > 1.5);

    }

    @Test
    void testSignedZeroEquality() {

        assertFalse(build(-0.0, 1.0).isZero());
        assertFalse(build(+0.0, 1.0).isZero());
        assertFalse(build( 1.0, -0.0).isZero());
        assertFalse(build( 1.0, +0.0).isZero());

        assertTrue(build(-0.0, -0.0).isZero());
        assertTrue(build(-0.0, +0.0).isZero());
        assertTrue(build(+0.0, -0.0).isZero());
        assertTrue(build(+0.0, +0.0).isZero());

        assertNotEquals(build(-0.0, -0.0), FieldComplex.getZero(Binary64Field.getInstance()));
        assertNotEquals(build(-0.0, +0.0), FieldComplex.getZero(Binary64Field.getInstance()));
        assertNotEquals(build(+0.0, -0.0), FieldComplex.getZero(Binary64Field.getInstance()));
        assertEquals(build(+0.0, +0.0), FieldComplex.getZero(Binary64Field.getInstance()));

    }

    /**
     * Class to test extending Complex
     */
    public static class TestComplex extends FieldComplex<Binary64> {

        public TestComplex(Binary64 real, Binary64 imaginary) {
            super(real, imaginary);
        }

        public TestComplex(FieldComplex<Binary64> other) {
            this(other.getRealPart(), other.getImaginaryPart());
        }

        @Override
        protected TestComplex createComplex(Binary64 real, Binary64 imaginary) {
            return new TestComplex(real, imaginary);
        }

    }
}
