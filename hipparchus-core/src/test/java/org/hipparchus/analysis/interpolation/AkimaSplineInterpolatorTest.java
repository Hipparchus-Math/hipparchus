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
package org.hipparchus.analysis.interpolation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.hipparchus.Field;
import org.hipparchus.RealFieldElement;
import org.hipparchus.analysis.RealFieldUnivariateFunction;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.util.Decimal64;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.Precision;
import org.junit.Assert;
import org.junit.Test;

public class AkimaSplineInterpolatorTest {

    @Test
    public void testIllegalArguments()
    {
        // Data set arrays of different size.
        UnivariateInterpolator i = new AkimaSplineInterpolator();

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

        try
        {
            double xval[] = { 0.0, 1.0, 2.0, 3.0 };
            double yval[] = { 0.0, 1.0, 2.0, 3.0 };
            i.interpolate( xval, yval );
            Assert.fail( "Failed to detect insufficient data" );
        }
        catch ( MathIllegalArgumentException iae )
        {
            // Expected.
        }

        try
        {
            double xval[] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
            double yval[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
            i.interpolate( xval, yval );
            Assert.fail( "Failed to detect data set array with different sizes." );
        }
        catch ( MathIllegalArgumentException iae )
        {
            // Expected.
        }

        // X values not sorted.
        try
        {
            double xval[] = { 0.0, 1.0, 0.5, 7.0, 3.5, 2.2, 8.0 };
            double yval[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
            i.interpolate( xval, yval );
            Assert.fail( "Failed to detect unsorted arguments." );
        }
        catch ( MathIllegalArgumentException iae )
        {
            // Expected.
        }
    }

    @Test
    public void testIllegalArgumentsD64()
    {
        // Data set arrays of different size.
        FieldUnivariateInterpolator i = new AkimaSplineInterpolator();

        try
        {
            Decimal64 yval[] = buildD64(0.0, 1.0, 2.0, 3.0, 4.0);
            i.interpolate( null, yval );
            Assert.fail( "Failed to detect x null pointer" );
        }
        catch ( NullArgumentException iae )
        {
            // Expected.
        }

        try
        {
            Decimal64 xval[] = buildD64(0.0, 1.0, 2.0, 3.0, 4.0);
            i.interpolate( xval, null );
            Assert.fail( "Failed to detect y null pointer" );
        }
        catch ( NullArgumentException iae )
        {
            // Expected.
        }

        try
        {
            Decimal64 xval[] = buildD64(0.0, 1.0, 2.0, 3.0);
            Decimal64 yval[] = buildD64(0.0, 1.0, 2.0, 3.0);
            i.interpolate( xval, yval );
            Assert.fail( "Failed to detect insufficient data" );
        }
        catch ( MathIllegalArgumentException iae )
        {
            // Expected.
        }

        try
        {
            Decimal64 xval[] = buildD64(0.0, 1.0, 2.0, 3.0, 4.0);
            Decimal64 yval[] = buildD64(0.0, 1.0, 2.0, 3.0, 4.0, 5.0);
            i.interpolate( xval, yval );
            Assert.fail( "Failed to detect data set array with different sizes." );
        }
        catch ( MathIllegalArgumentException iae )
        {
            // Expected.
        }

        // X values not sorted.
        try
        {
            Decimal64 xval[] = buildD64(0.0, 1.0, 0.5, 7.0, 3.5, 2.2, 8.0);
            Decimal64 yval[] = buildD64(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
            i.interpolate( xval, yval );
            Assert.fail( "Failed to detect unsorted arguments." );
        }
        catch ( MathIllegalArgumentException iae )
        {
            // Expected.
        }
    }

    /*
     * Interpolate a straight line. <p> y = 2 x - 5 <p> Tolerances determined by performing same calculation using
     * Math.NET over ten runs of 100 random number draws for the same function over the same span with the same number
     * of elements
     */
    @Test
    public void testInterpolateLine()
    {
        final int numberOfElements = 10;
        final double minimumX = -10;
        final double maximumX = 10;
        final int numberOfSamples = 100;
        final double interpolationTolerance = 1e-15;
        final double maxTolerance = 1e-15;

        UnivariateFunction f = x -> 2 * x - 5;

        testInterpolation( minimumX, maximumX, numberOfElements, numberOfSamples, f, interpolationTolerance,
                           maxTolerance );
    }

    /*
     * Interpolate a straight line. <p> y = 2 x - 5 <p> Tolerances determined by performing same calculation using
     * Math.NET over ten runs of 100 random number draws for the same function over the same span with the same number
     * of elements
     */
    @Test
    public void testInterpolateLineD64()
    {
        final int numberOfElements = 10;
        final Decimal64 minimumX = new Decimal64(-10);
        final Decimal64 maximumX = new Decimal64(10);
        final int numberOfSamples = 100;
        final double interpolationTolerance = 1e-15;
        final double maxTolerance = 1e-15;

        RealFieldUnivariateFunction<Decimal64> f = x -> x.multiply(2).subtract(5);

        testInterpolation( minimumX, maximumX, numberOfElements, numberOfSamples, f, interpolationTolerance,
                           maxTolerance );
    }

    /*
     * Interpolate a straight line. <p> y = 3 x<sup>2</sup> - 5 x + 7 <p> Tolerances determined by performing same
     * calculation using Math.NET over ten runs of 100 random number draws for the same function over the same span with
     * the same number of elements
     */

    @Test
    public void testInterpolateParabola()
    {
        final int numberOfElements = 10;
        final double minimumX = -10;
        final double maximumX = 10;
        final int numberOfSamples = 100;
        final double interpolationTolerance = 7e-15;
        final double maxTolerance = 6e-14;

        UnivariateFunction f = x -> ( 3 * x * x ) - ( 5 * x ) + 7;

        testInterpolation( minimumX, maximumX, numberOfElements, numberOfSamples, f, interpolationTolerance,
                           maxTolerance );
    }

    /*
     * Interpolate a straight line. <p> y = 3 x<sup>2</sup> - 5 x + 7 <p> Tolerances determined by performing same
     * calculation using Math.NET over ten runs of 100 random number draws for the same function over the same span with
     * the same number of elements
     */

    @Test
    public void testInterpolateParabolaD64()
    {
        final int numberOfElements = 10;
        final Decimal64 minimumX = new Decimal64(-10);
        final Decimal64 maximumX = new Decimal64(10);
        final int numberOfSamples = 100;
        final double interpolationTolerance = 7e-15;
        final double maxTolerance = 6e-14;

        RealFieldUnivariateFunction<Decimal64> f = x -> x.multiply(x).multiply(3).
                                                        subtract(x.multiply(5)).
                                                        add(7);

        testInterpolation( minimumX, maximumX, numberOfElements, numberOfSamples, f, interpolationTolerance,
                           maxTolerance );
    }

    /*
     * Interpolate a straight line. <p> y = 3 x<sup>3</sup> - 0.5 x<sup>2</sup> + x - 1 <p> Tolerances determined by
     * performing same calculation using Math.NET over ten runs of 100 random number draws for the same function over
     * the same span with the same number of elements
     */
    @Test
    public void testInterpolateCubic()
    {
        final int numberOfElements = 10;
        final double minimumX = -3;
        final double maximumX = 3;
        final int numberOfSamples = 100;
        final double interpolationTolerance = 0.37;
        final double maxTolerance = 3.8;

        UnivariateFunction f = x -> ( 3 * x * x * x ) - ( 0.5 * x * x ) + ( 1 * x ) - 1;

        testInterpolation( minimumX, maximumX, numberOfElements, numberOfSamples, f, interpolationTolerance,
                           maxTolerance );
    }

    /*
     * Interpolate a straight line. <p> y = 3 x<sup>3</sup> - 0.5 x<sup>2</sup> + x - 1 <p> Tolerances determined by
     * performing same calculation using Math.NET over ten runs of 100 random number draws for the same function over
     * the same span with the same number of elements
     */
    @Test
    public void testInterpolateCubic64()
    {
        final int numberOfElements = 10;
        final Decimal64 minimumX = new Decimal64(-3);
        final Decimal64 maximumX = new Decimal64(3);
        final int numberOfSamples = 100;
        final double interpolationTolerance = 0.37;
        final double maxTolerance = 3.8;

        RealFieldUnivariateFunction<Decimal64> f = x -> x.multiply(x).multiply(x).multiply(3).
                                                        subtract(x.multiply(x).multiply(0.5)).
                                                        add(x).
                                                        subtract(1);
        
        testInterpolation( minimumX, maximumX, numberOfElements, numberOfSamples, f, interpolationTolerance,
                           maxTolerance );
    }

    private void testInterpolation( double minimumX, double maximumX, int numberOfElements, int numberOfSamples,
                                    UnivariateFunction f, double tolerance, double maxTolerance )
    {
        double expected;
        double actual;
        double currentX;
        final double delta = ( maximumX - minimumX ) / ( (double) numberOfElements );
        double xValues[] = new double[numberOfElements];
        double yValues[] = new double[numberOfElements];

        for ( int i = 0; i < numberOfElements; i++ )
        {
            xValues[i] = minimumX + delta * (double) i;
            yValues[i] = f.value( xValues[i] );
        }

        UnivariateFunction interpolation = new AkimaSplineInterpolator().interpolate( xValues, yValues );

        for ( int i = 0; i < numberOfElements; i++ )
        {
            currentX = xValues[i];
            expected = f.value( currentX );
            actual = interpolation.value( currentX );
            assertTrue( Precision.equals( expected, actual ) );
        }

        final RandomDataGenerator randomDataGenerator = new RandomDataGenerator(1234567L);

        double sumError = 0;
        for ( int i = 0; i < numberOfSamples; i++ )
        {
            currentX = randomDataGenerator.nextUniform(xValues[0], xValues[xValues.length - 1]);
            expected = f.value( currentX );
            actual = interpolation.value( currentX );
            sumError += FastMath.abs( actual - expected );
            assertEquals( expected, actual, maxTolerance );
        }

        assertEquals( 0.0, ( sumError / (double) numberOfSamples ), tolerance );
    }

    private <T extends RealFieldElement<T>> void testInterpolation(T minimumX, T maximumX,
                                                                   int numberOfElements, int numberOfSamples,
                                                                   RealFieldUnivariateFunction<T> f,
                                                                   double tolerance, double maxTolerance)
    {
        final Field<T> field = minimumX.getField();
        T expected;
        T actual;
        T currentX;
        final T delta = maximumX.subtract(minimumX).divide(numberOfElements);
        T xValues[] = MathArrays.buildArray(field, numberOfElements);
        T yValues[] = MathArrays.buildArray(field, numberOfElements);

        for ( int i = 0; i < numberOfElements; i++ )
        {
            xValues[i] = minimumX.add(delta.multiply(i));
            yValues[i] = f.value(xValues[i]);
        }

        RealFieldUnivariateFunction<T> interpolation =
                        new AkimaSplineInterpolator().interpolate( xValues, yValues );

        for ( int i = 0; i < numberOfElements; i++ )
        {
            currentX = xValues[i];
            expected = f.value( currentX );
            actual = interpolation.value( currentX );
            assertTrue(Precision.equals(expected.getReal(), actual.getReal()));
        }

        final RandomDataGenerator randomDataGenerator = new RandomDataGenerator(1234567L);

        double sumError = 0;
        for ( int i = 0; i < numberOfSamples; i++ )
        {
            final double r = randomDataGenerator.nextUniform(xValues[0].getReal(),
                                                             xValues[xValues.length - 1].getReal());
            currentX = field.getZero().add(r);
            expected = f.value( currentX );
            actual = interpolation.value( currentX );
            sumError += FastMath.abs(actual.subtract(expected)).getReal();
            assertEquals(expected.getReal(), actual.getReal(), maxTolerance);
        }

        assertEquals(0.0, sumError / numberOfSamples, tolerance);
    }

    private Decimal64[] buildD64(double...c) {
        Decimal64[] array = new Decimal64[c.length];
        for (int i = 0; i < c.length; ++i) {
            array[i] = new Decimal64(c[i]);
        }
        return array;
    }

}
