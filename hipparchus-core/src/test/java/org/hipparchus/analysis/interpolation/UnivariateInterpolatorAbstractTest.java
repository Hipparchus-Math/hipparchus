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
package org.hipparchus.analysis.interpolation;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.UnitTestUtils;
import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.polynomials.FieldPolynomialFunction;
import org.hipparchus.analysis.polynomials.FieldPolynomialSplineFunction;
import org.hipparchus.analysis.polynomials.PolynomialFunction;
import org.hipparchus.analysis.polynomials.PolynomialSplineFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.Binary64;
import org.junit.Assert;
import org.junit.Test;

/**
 * Base test for interpolators.
 */
public abstract class UnivariateInterpolatorAbstractTest {

    /** error tolerance for spline interpolator value at knot points */
    protected double knotTolerance = 1E-12;

    /** error tolerance for interpolating polynomial coefficients */
    protected double coefficientTolerance = 1E-6;

    /** error tolerance for interpolated values */
    protected double interpolationTolerance = 1E-12;

    protected abstract UnivariateInterpolator buildDoubleInterpolator();

    protected abstract FieldUnivariateInterpolator buildFieldInterpolator();

    @Test
    public void testInterpolateLinearDegenerateTwoSegment()
        {
        double x[] = { 0.0, 0.5, 1.0 };
        double y[] = { 0.0, 0.5, 1.0 };
        UnivariateInterpolator i = buildDoubleInterpolator();
        UnivariateFunction f = i.interpolate(x, y);
        verifyInterpolation(f, x, y);

        // Verify coefficients using analytical values
        PolynomialFunction polynomials[] = ((PolynomialSplineFunction) f).getPolynomials();
        double target[] = {y[0], 1d};
        UnitTestUtils.assertEquals(polynomials[0].getCoefficients(), target, coefficientTolerance);
        target = new double[]{y[1], 1d};
        UnitTestUtils.assertEquals(polynomials[1].getCoefficients(), target, coefficientTolerance);

        // Check interpolation
        Assert.assertEquals(0.0,f.value(0.0), interpolationTolerance);
        Assert.assertEquals(0.4,f.value(0.4), interpolationTolerance);
        Assert.assertEquals(1.0,f.value(1.0), interpolationTolerance);
    }

    @Test
    public void testInterpolateLinearDegenerateTwoSegmentD64()
        {
        Binary64 x[] = buildD64(0.0, 0.5, 1.0);
        Binary64 y[] = buildD64(0.0, 0.5, 1.0);
        FieldUnivariateInterpolator i = buildFieldInterpolator();
        CalculusFieldUnivariateFunction<Binary64> f = i.interpolate(x, y);
        verifyInterpolation(f, x, y);

        // Verify coefficients using analytical values
        FieldPolynomialFunction<Binary64> polynomials[] = ((FieldPolynomialSplineFunction<Binary64>) f).getPolynomials();
        checkCoeffs(coefficientTolerance, polynomials[0], y[0].getReal(), 1.0);
        checkCoeffs(coefficientTolerance, polynomials[1], y[1].getReal(), 1.0);

        // Check interpolation
        Assert.assertEquals(0.0, f.value(new Binary64(0.0)).getReal(), interpolationTolerance);
        Assert.assertEquals(0.4, f.value(new Binary64(0.4)).getReal(), interpolationTolerance);
        Assert.assertEquals(1.0, f.value(new Binary64(1.0)).getReal(), interpolationTolerance);
    }

    @Test
    public void testInterpolateLinearDegenerateThreeSegment()
        {
        double x[] = { 0.0, 0.5, 1.0, 1.5 };
        double y[] = { 0.0, 0.5, 1.0, 1.5 };
        UnivariateInterpolator i = buildDoubleInterpolator();
        UnivariateFunction f = i.interpolate(x, y);
        verifyInterpolation(f, x, y);

        // Verify coefficients using analytical values
        PolynomialFunction polynomials[] = ((PolynomialSplineFunction) f).getPolynomials();
        double target[] = {y[0], 1d};
        UnitTestUtils.assertEquals(polynomials[0].getCoefficients(), target, coefficientTolerance);
        target = new double[]{y[1], 1d};
        UnitTestUtils.assertEquals(polynomials[1].getCoefficients(), target, coefficientTolerance);
        target = new double[]{y[2], 1d};
        UnitTestUtils.assertEquals(polynomials[2].getCoefficients(), target, coefficientTolerance);

        // Check interpolation
        Assert.assertEquals(0,f.value(0), interpolationTolerance);
        Assert.assertEquals(1.4,f.value(1.4), interpolationTolerance);
        Assert.assertEquals(1.5,f.value(1.5), interpolationTolerance);
    }

    @Test
    public void testInterpolateLinearDegenerateThreeSegmentD64()
        {
        Binary64 x[] = buildD64(0.0, 0.5, 1.0, 1.5);
        Binary64 y[] = buildD64(0.0, 0.5, 1.0, 1.5);
        FieldUnivariateInterpolator i = buildFieldInterpolator();
        CalculusFieldUnivariateFunction<Binary64> f = i.interpolate(x, y);
        verifyInterpolation(f, x, y);

        // Verify coefficients using analytical values
        FieldPolynomialFunction<Binary64> polynomials[] = ((FieldPolynomialSplineFunction<Binary64>) f).getPolynomials();
        checkCoeffs(coefficientTolerance, polynomials[0], y[0].getReal(), 1.0);
        checkCoeffs(coefficientTolerance, polynomials[1], y[1].getReal(), 1.0);
        checkCoeffs(coefficientTolerance, polynomials[2], y[2].getReal(), 1.0);

        // Check interpolation
        Assert.assertEquals(0,   f.value(new Binary64(0)).getReal(),   interpolationTolerance);
        Assert.assertEquals(1.4, f.value(new Binary64(1.4)).getReal(), interpolationTolerance);
        Assert.assertEquals(1.5, f.value(new Binary64(1.5)).getReal(), interpolationTolerance);
    }

    @Test
    public void testIllegalArguments() {
        UnivariateInterpolator i = buildDoubleInterpolator();
        try
        {
            double yval[] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
            i.interpolate( null, yval );
            Assert.fail( "Failed to detect x null pointer" );
        }
        catch ( NullArgumentException iae )
        {
            // Expected.
        }

        try
        {
            double xval[] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
            i.interpolate( xval, null );
            Assert.fail( "Failed to detect y null pointer" );
        }
        catch ( NullArgumentException iae )
        {
            // Expected.
        }

        // Data set arrays of different size.
        try {
            double xval[] = { 0.0, 1.0 };
            double yval[] = { 0.0, 1.0, 2.0 };
            i.interpolate(xval, yval);
            Assert.fail("Failed to detect data set array with different sizes.");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }
        // X values not sorted.
        try {
            double xval[] = { 0.0, 1.0, 0.5 };
            double yval[] = { 0.0, 1.0, 2.0 };
            i.interpolate(xval, yval);
            Assert.fail("Failed to detect unsorted arguments.");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }
        // Not enough data to interpolate.
        try {
            double xval[] = { 0.0 };
            double yval[] = { 0.0 };
            i.interpolate(xval, yval);
            Assert.fail("Failed to detect unsorted arguments.");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }
    }

    @Test
    public void testIllegalArgumentsD64() {
        FieldUnivariateInterpolator i = buildFieldInterpolator();
        try
        {
            Binary64 yval[] = buildD64(0.0, 1.0, 2.0, 3.0, 4.0);
            i.interpolate( null, yval );
            Assert.fail( "Failed to detect x null pointer" );
        }
        catch ( NullArgumentException iae )
        {
            // Expected.
        }

        try
        {
            Binary64 xval[] = buildD64(0.0, 1.0, 2.0, 3.0, 4.0);
            i.interpolate( xval, null );
            Assert.fail( "Failed to detect y null pointer" );
        }
        catch ( NullArgumentException iae )
        {
            // Expected.
        }

        // Data set arrays of different size.
        try {
            Binary64 xval[] = buildD64(0.0, 1.0);
            Binary64 yval[] = buildD64(0.0, 1.0, 2.0);
            i.interpolate(xval, yval);
            Assert.fail("Failed to detect data set array with different sizes.");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }
        // X values not sorted.
        try {
            Binary64 xval[] = buildD64(0.0, 1.0, 0.5);
            Binary64 yval[] = buildD64(0.0, 1.0, 2.0);
            i.interpolate(xval, yval);
            Assert.fail("Failed to detect unsorted arguments.");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }
        // Not enough data to interpolate.
        try {
            Binary64 xval[] = buildD64(0.0);
            Binary64 yval[] = buildD64(0.0);
            i.interpolate(xval, yval);
            Assert.fail("Failed to detect unsorted arguments.");
        } catch (MathIllegalArgumentException iae) {
            // Expected.
        }
    }

    /**
     * verifies that f(x[i]) = y[i] for i = 0..n-1 where n is common length.
     */
    protected void verifyInterpolation(UnivariateFunction f, double x[], double y[]) {
        for (int i = 0; i < x.length; i++) {
            Assert.assertEquals(y[i], f.value(x[i]), knotTolerance);
        }
    }

    /**
     * verifies that f(x[i]) = y[i] for i = 0..n-1 where n is common length.
     */
    protected <T extends CalculusFieldElement<T>> void verifyInterpolation(CalculusFieldUnivariateFunction<T> f,
                                                                           T x[], T y[]) {
        for (int i = 0; i < x.length; i++) {
            Assert.assertEquals( y[i].getReal(), f.value(x[i]).getReal(), knotTolerance);
        }
    }

    protected Binary64[] buildD64(double...c) {
        Binary64[] array = new Binary64[c.length];
        for (int i = 0; i < c.length; ++i) {
            array[i] = new Binary64(c[i]);
        }
        return array;
    }

    protected <T extends CalculusFieldElement<T>> void checkCoeffs(final double tolerance, final FieldPolynomialFunction<T> p,
                                                                   final double... ref) {
        final T[] c = p.getCoefficients();
        Assert.assertEquals(ref.length, c.length);
        for (int i = 0; i < ref.length; ++i) {
            Assert.assertEquals(ref[i], c[i].getReal(), tolerance);
        }
    }

}
