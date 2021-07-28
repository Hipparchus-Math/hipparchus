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
package org.hipparchus.analysis.integration.gauss;

import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.util.Decimal64;
import org.hipparchus.util.Decimal64Field;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test of the {@link HermiteRuleFactory}.
 *
 */
public class FieldHermiteTest {
    private static final FieldGaussIntegratorFactory<Decimal64> factory = new FieldGaussIntegratorFactory<>(Decimal64Field.getInstance());

    @Test
    public void testNormalDistribution() {
        final Decimal64 oneOverSqrtPi = new Decimal64(1 / FastMath.sqrt(Math.PI));

        // By defintion, Gauss-Hermite quadrature readily provides the
        // integral of the normal distribution density.
        final int numPoints = 1;

        // Change of variable:
        //   y = (x - mu) / (sqrt(2) *  sigma)
        // such that the integrand
        //   N(x, mu, sigma)
        // is transformed to
        //   f(y) * exp(-y^2)
        final CalculusFieldUnivariateFunction<Decimal64> f = y -> oneOverSqrtPi;

        final FieldGaussIntegrator<Decimal64> integrator = factory.hermite(numPoints);
        final double result = integrator.integrate(f).getReal();
        final double expected = 1;
        Assert.assertEquals(expected, result, Math.ulp(expected));
    }

    @Test
    public void testNormalMean() {
        final Decimal64 sqrtTwo = new Decimal64(FastMath.sqrt(2));
        final Decimal64 oneOverSqrtPi = new Decimal64(1 / FastMath.sqrt(Math.PI));

        final Decimal64 mu = new Decimal64(12345.6789);
        final Decimal64 sigma = new Decimal64(987.654321);
        final int numPoints = 5;

        // Change of variable:
        //   y = (x - mu) / (sqrt(2) *  sigma)
        // such that the integrand
        //   x * N(x, mu, sigma)
        // is transformed to
        //   f(y) * exp(-y^2)
        final CalculusFieldUnivariateFunction<Decimal64> f =
                        y ->  oneOverSqrtPi.multiply(sqrtTwo.multiply(sigma).multiply(y).add(mu));

        final FieldGaussIntegrator<Decimal64> integrator = factory.hermite(numPoints);
        final double result = integrator.integrate(f).getReal();
        final double expected = mu.getReal();
        Assert.assertEquals(expected, result, 3 * Math.ulp(expected));
    }

    @Test
    public void testNormalVariance() {
        final Decimal64 twoOverSqrtPi = new Decimal64(2 / FastMath.sqrt(Math.PI));

        final Decimal64 sigma = new Decimal64(987.654321);
        final Decimal64 sigma2 = sigma.multiply(sigma);
        final int numPoints = 5;

        // Change of variable:
        //   y = (x - mu) / (sqrt(2) *  sigma)
        // such that the integrand
        //   (x - mu)^2 * N(x, mu, sigma)
        // is transformed to
        //   f(y) * exp(-y^2)
        final CalculusFieldUnivariateFunction<Decimal64> f =
                        y -> twoOverSqrtPi.multiply(sigma2).multiply(y).multiply(y);

        final FieldGaussIntegrator<Decimal64> integrator = factory.hermite(numPoints);
        final double result = integrator.integrate(f).getReal();
        final double expected = sigma2.getReal();
        Assert.assertEquals(expected, result, 10 * Math.ulp(expected));
    }
}
