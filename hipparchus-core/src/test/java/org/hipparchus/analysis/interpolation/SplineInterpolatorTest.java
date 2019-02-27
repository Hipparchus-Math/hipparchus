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

import org.hipparchus.UnitTestUtils;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.polynomials.PolynomialFunction;
import org.hipparchus.analysis.polynomials.PolynomialSplineFunction;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test the SplineInterpolator.
 *
 */
public class SplineInterpolatorTest extends UnivariateInterpolatorAbstractTest {

    protected UnivariateInterpolator buildDoubleInterpolator() {
        return new SplineInterpolator();
    }

    protected FieldUnivariateInterpolator buildFieldInterpolator() {
        return new SplineInterpolator();
    }

    @Test
    public void testInterpolateSin() {
        double sineCoefficientTolerance = 1e-6;
        double sineInterpolationTolerance = 0.0043;
        double x[] =
            {
                0.0,
                FastMath.PI / 6d,
                FastMath.PI / 2d,
                5d * FastMath.PI / 6d,
                FastMath.PI,
                7d * FastMath.PI / 6d,
                3d * FastMath.PI / 2d,
                11d * FastMath.PI / 6d,
                2.d * FastMath.PI };
        double y[] = { 0d, 0.5d, 1d, 0.5d, 0d, -0.5d, -1d, -0.5d, 0d };
        UnivariateInterpolator i = buildDoubleInterpolator();
        UnivariateFunction f = i.interpolate(x, y);
        verifyInterpolation(f, x, y);
        verifyConsistency((PolynomialSplineFunction) f, x);

        /* Check coefficients against values computed using R (version 1.8.1, Red Hat Linux 9)
         *
         * To replicate in R:
         *     x[1] <- 0
         *     x[2] <- pi / 6, etc, same for y[] (could use y <- scan() for y values)
         *     g <- splinefun(x, y, "natural")
         *     splinecoef <- eval(expression(z), envir = environment(g))
         *     print(splinecoef)
         */
        PolynomialFunction polynomials[] = ((PolynomialSplineFunction) f).getPolynomials();
        double target[] = {y[0], 1.002676d, 0d, -0.17415829d};
        UnitTestUtils.assertEquals(polynomials[0].getCoefficients(), target, sineCoefficientTolerance);
        target = new double[]{y[1], 8.594367e-01, -2.735672e-01, -0.08707914};
        UnitTestUtils.assertEquals(polynomials[1].getCoefficients(), target, sineCoefficientTolerance);
        target = new double[]{y[2], 1.471804e-17,-5.471344e-01, 0.08707914};
        UnitTestUtils.assertEquals(polynomials[2].getCoefficients(), target, sineCoefficientTolerance);
        target = new double[]{y[3], -8.594367e-01, -2.735672e-01, 0.17415829};
        UnitTestUtils.assertEquals(polynomials[3].getCoefficients(), target, sineCoefficientTolerance);
        target = new double[]{y[4], -1.002676, 6.548562e-17, 0.17415829};
        UnitTestUtils.assertEquals(polynomials[4].getCoefficients(), target, sineCoefficientTolerance);
        target = new double[]{y[5], -8.594367e-01, 2.735672e-01, 0.08707914};
        UnitTestUtils.assertEquals(polynomials[5].getCoefficients(), target, sineCoefficientTolerance);
        target = new double[]{y[6], 3.466465e-16, 5.471344e-01, -0.08707914};
        UnitTestUtils.assertEquals(polynomials[6].getCoefficients(), target, sineCoefficientTolerance);
        target = new double[]{y[7], 8.594367e-01, 2.735672e-01, -0.17415829};
        UnitTestUtils.assertEquals(polynomials[7].getCoefficients(), target, sineCoefficientTolerance);

        //Check interpolation
        Assert.assertEquals(FastMath.sqrt(2d) / 2d,f.value(FastMath.PI/4d),sineInterpolationTolerance);
        Assert.assertEquals(FastMath.sqrt(2d) / 2d,f.value(3d*FastMath.PI/4d),sineInterpolationTolerance);
    }

    /**
     * Verifies that interpolating polynomials satisfy consistency requirement:
     *    adjacent polynomials must agree through two derivatives at knot points
     */
    protected void verifyConsistency(PolynomialSplineFunction f, double x[])
        {
        PolynomialFunction polynomials[] = f.getPolynomials();
        for (int i = 1; i < x.length - 2; i++) {
            // evaluate polynomials and derivatives at x[i + 1]
            Assert.assertEquals(polynomials[i].value(x[i +1] - x[i]), polynomials[i + 1].value(0), 0.1);
            Assert.assertEquals(polynomials[i].polynomialDerivative().value(x[i +1] - x[i]),
                                polynomials[i + 1].polynomialDerivative().value(0), 0.5);
            Assert.assertEquals(polynomials[i].polynomialDerivative().polynomialDerivative().value(x[i +1] - x[i]),
                                polynomials[i + 1].polynomialDerivative().polynomialDerivative().value(0), 0.5);
        }
    }

}
