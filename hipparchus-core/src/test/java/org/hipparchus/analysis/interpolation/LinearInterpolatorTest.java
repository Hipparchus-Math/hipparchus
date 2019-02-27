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
import org.hipparchus.analysis.RealFieldUnivariateFunction;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.polynomials.FieldPolynomialFunction;
import org.hipparchus.analysis.polynomials.FieldPolynomialSplineFunction;
import org.hipparchus.analysis.polynomials.PolynomialFunction;
import org.hipparchus.analysis.polynomials.PolynomialSplineFunction;
import org.hipparchus.util.Decimal64;
import org.junit.Test;

/**
 * Test the LinearInterpolator.
 */
public class LinearInterpolatorTest extends UnivariateInterpolatorAbstractTest {

    protected UnivariateInterpolator buildDoubleInterpolator() {
        return new LinearInterpolator();
    }

    protected FieldUnivariateInterpolator buildFieldInterpolator() {
        return new LinearInterpolator();
    }

    @Test
    public void testInterpolateLinear() {
        double x[] = { 0.0, 0.5, 1.0 };
        double y[] = { 0.0, 0.5, 0.0 };
        UnivariateInterpolator i = buildDoubleInterpolator();
        UnivariateFunction f = i.interpolate(x, y);
        verifyInterpolation(f, x, y);

        // Verify coefficients using analytical values
        PolynomialFunction polynomials[] = ((PolynomialSplineFunction) f).getPolynomials();
        double target[] = {y[0], 1d};
        UnitTestUtils.assertEquals(polynomials[0].getCoefficients(), target, coefficientTolerance);
        target = new double[]{y[1], -1d};
        UnitTestUtils.assertEquals(polynomials[1].getCoefficients(), target, coefficientTolerance);
    }

    @Test
    public void testInterpolateLinearD64() {
        Decimal64 x[] = buildD64(0.0, 0.5, 1.0);
        Decimal64 y[] = buildD64(0.0, 0.5, 0.0);
        FieldUnivariateInterpolator i = buildFieldInterpolator();
        RealFieldUnivariateFunction<Decimal64> f = i.interpolate(x, y);
        verifyInterpolation(f, x, y);

        // Verify coefficients using analytical values
        FieldPolynomialFunction<Decimal64> polynomials[] = ((FieldPolynomialSplineFunction<Decimal64>) f).getPolynomials();
        checkCoeffs(coefficientTolerance, polynomials[0], y[0].getReal(), +1.0);
        checkCoeffs(coefficientTolerance, polynomials[1], y[1].getReal(), -1.0);
    }

}
