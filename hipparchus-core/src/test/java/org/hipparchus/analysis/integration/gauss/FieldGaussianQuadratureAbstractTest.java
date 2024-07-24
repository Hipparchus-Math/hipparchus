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

import org.hipparchus.util.Binary64;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;

/**
 * Base class for standard testing of Gaussian quadrature rules,
 * which are exact for polynomials up to a certain degree. In this test, each
 * monomial in turn is tested against the specified quadrature rule.
 *
 */
public abstract class FieldGaussianQuadratureAbstractTest {

    /**
     * Returns the expected value of the integral of the specified monomial.
     * The integration is carried out on the natural interval of the quadrature
     * rule under test.
     *
     * @param n Degree of the monomial.
     * @return the expected value of the integral of x<sup>n</sup>.
     */
    public abstract double getExpectedValue(final int n);

    /**
     * Checks that the value of the integral of each monomial
     * <code>x<sup>0</sup>, ... , x<sup>p</sup></code>
     * returned by the quadrature rule under test conforms with the expected
     * value. Here {@code p} denotes the degree of the highest polynomial for
     * which exactness is to be expected.
     */
    public void testAllMonomials(FieldGaussIntegrator<Binary64> integrator,
                                 int maxDegree, double eps, double numUlps) {
        for (int n = 0; n <= maxDegree; n++) {
            final double expected = getExpectedValue(n);

            final int p = n;
            final double actual = integrator.integrate(x -> FastMath.pow(x, p))
                            .getReal();

            // System.out.println(n + "/" + maxDegree + " " + integrator.getNumberOfPoints()
            //                    + " " + expected + " " + actual + " " + Math.ulp(expected));
            if (expected == 0) {
                Assertions.assertEquals(expected, actual, eps,
                                        "while integrating monomial x**" + n + " with a " + integrator.getNumberOfPoints() + "-point quadrature rule");
            } else {
                double err = FastMath.abs(actual - expected) / Math.ulp(
                                expected);
                Assertions.assertEquals(expected, actual,
                                        Math.ulp(expected) * numUlps,
                                        "while integrating monomial x**" + n + " with a " + +integrator.getNumberOfPoints() + "-point quadrature rule, " + " error was " + err + " ulps");
            }
        }
    }
}
