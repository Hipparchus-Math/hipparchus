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
package org.hipparchus.analysis.integration.gauss;

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test of the {@link HermiteRuleFactory}.
 *
 */
class HermiteTest {
    private static final GaussIntegratorFactory factory = new GaussIntegratorFactory();

    @Test
    void testNormalDistribution() {
        final double oneOverSqrtPi = 1 / FastMath.sqrt(Math.PI);

        // By defintion, Gauss-Hermite quadrature readily provides the
        // integral of the normal distribution density.
        final int numPoints = 1;

        // Change of variable:
        //   y = (x - mu) / (sqrt(2) *  sigma)
        // such that the integrand
        //   N(x, mu, sigma)
        // is transformed to
        //   f(y) * exp(-y^2)
        final UnivariateFunction f = new UnivariateFunction() {
                public double value(double y) {
                    return oneOverSqrtPi; // Constant function.
                }
            };

        final GaussIntegrator integrator = factory.hermite(numPoints);
        final double result = integrator.integrate(f);
        final double expected = 1;
        assertEquals(expected, result, FastMath.ulp(expected));
    }

    @Test
    void testNormalMean() {
        final double sqrtTwo = FastMath.sqrt(2);
        final double oneOverSqrtPi = 1 / FastMath.sqrt(Math.PI);

        final double mu = 12345.6789;
        final double sigma = 987.654321;
        final int numPoints = 6;

        // Change of variable:
        //   y = (x - mu) / (sqrt(2) *  sigma)
        // such that the integrand
        //   x * N(x, mu, sigma)
        // is transformed to
        //   f(y) * exp(-y^2)
        final UnivariateFunction f = new UnivariateFunction() {
                public double value(double y) {
                    return oneOverSqrtPi * (sqrtTwo * sigma * y + mu);
                }
            };

        final GaussIntegrator integrator = factory.hermite(numPoints);
        final double result = integrator.integrate(f);
        final double expected = mu;
        assertEquals(expected, result, 5 * FastMath.ulp(expected));
    }

    @Test
    void testNormalVariance() {
        final double twoOverSqrtPi = 2 / FastMath.sqrt(Math.PI);

        final double sigma = 987.654321;
        final double sigma2 = sigma * sigma;
        final int numPoints = 5;

        // Change of variable:
        //   y = (x - mu) / (sqrt(2) *  sigma)
        // such that the integrand
        //   (x - mu)^2 * N(x, mu, sigma)
        // is transformed to
        //   f(y) * exp(-y^2)
        final UnivariateFunction f = new UnivariateFunction() {
                public double value(double y) {
                    return twoOverSqrtPi * sigma2 * y * y;
                }
            };

        final GaussIntegrator integrator = factory.hermite(numPoints);
        final double result = integrator.integrate(f);
        final double expected = sigma2;
        assertEquals(expected, result, 10 * FastMath.ulp(expected));
    }
}
