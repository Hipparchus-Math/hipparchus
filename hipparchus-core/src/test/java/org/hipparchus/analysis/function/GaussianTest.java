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

package org.hipparchus.analysis.function;

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.analysis.differentiation.UnivariateDifferentiableFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Test for class {@link Gaussian}.
 */
class GaussianTest {
    private final double EPS = Math.ulp(1d);

    @Test
    void testPreconditions() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new Gaussian(1, 2, -1);
        });
    }

    @Test
    void testSomeValues() {
        final UnivariateFunction f = new Gaussian();

        assertEquals(1 / FastMath.sqrt(2 * Math.PI), f.value(0), EPS);
    }

    @Test
    void testLargeArguments() {
        final UnivariateFunction f = new Gaussian();

        assertEquals(0, f.value(Double.NEGATIVE_INFINITY), 0);
        assertEquals(0, f.value(-Double.MAX_VALUE), 0);
        assertEquals(0, f.value(-1e2), 0);
        assertEquals(0, f.value(1e2), 0);
        assertEquals(0, f.value(Double.MAX_VALUE), 0);
        assertEquals(0, f.value(Double.POSITIVE_INFINITY), 0);
    }

    @Test
    void testDerivatives() {
        final UnivariateDifferentiableFunction gaussian = new Gaussian(2.0, 0.9, 3.0);
        final DerivativeStructure dsX = new DSFactory(1, 4).variable(0, 1.1);
        final DerivativeStructure dsY = gaussian.value(dsX);
        assertEquals( 1.9955604901712128349,   dsY.getValue(),              EPS);
        assertEquals(-0.044345788670471396332, dsY.getPartialDerivative(1), EPS);
        assertEquals(-0.22074348138190206174,  dsY.getPartialDerivative(2), EPS);
        assertEquals( 0.014760030401924800557, dsY.getPartialDerivative(3), EPS);
        assertEquals( 0.073253159785035691678, dsY.getPartialDerivative(4), EPS);
    }

    @Test
    void testDerivativeLargeArguments() {
        final Gaussian f = new Gaussian(0, 1e-50);

        DSFactory factory = new DSFactory(1, 1);
        assertEquals(0, f.value(factory.variable(0, Double.NEGATIVE_INFINITY)).getPartialDerivative(1), 0);
        assertEquals(0, f.value(factory.variable(0, -Double.MAX_VALUE)).getPartialDerivative(1), 0);
        assertEquals(0, f.value(factory.variable(0, -1e50)).getPartialDerivative(1), 0);
        assertEquals(0, f.value(factory.variable(0, -1e2)).getPartialDerivative(1), 0);
        assertEquals(0, f.value(factory.variable(0, 1e2)).getPartialDerivative(1), 0);
        assertEquals(0, f.value(factory.variable(0, 1e50)).getPartialDerivative(1), 0);
        assertEquals(0, f.value(factory.variable(0, Double.MAX_VALUE)).getPartialDerivative(1), 0);
        assertEquals(0, f.value(factory.variable(0, Double.POSITIVE_INFINITY)).getPartialDerivative(1), 0);
    }

    @Test
    void testDerivativesNaN() {
        final Gaussian f = new Gaussian(0, 1e-50);
        final DerivativeStructure fx = f.value(new DSFactory(1, 5).variable(0, Double.NaN));
        for (int i = 0; i <= fx.getOrder(); ++i) {
            assertTrue(Double.isNaN(fx.getPartialDerivative(i)));
        }
    }

    @Test
    void testParametricUsage1() {
        assertThrows(NullArgumentException.class, () -> {
            final Gaussian.Parametric g = new Gaussian.Parametric();
            g.value(0, null);
        });
    }

    @Test
    void testParametricUsage2() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Gaussian.Parametric g = new Gaussian.Parametric();
            g.value(0, new double[]{0});
        });
    }

    @Test
    void testParametricUsage3() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Gaussian.Parametric g = new Gaussian.Parametric();
            g.value(0, new double[]{0, 1, 0});
        });
    }

    @Test
    void testParametricUsage4() {
        assertThrows(NullArgumentException.class, () -> {
            final Gaussian.Parametric g = new Gaussian.Parametric();
            g.gradient(0, null);
        });
    }

    @Test
    void testParametricUsage5() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Gaussian.Parametric g = new Gaussian.Parametric();
            g.gradient(0, new double[]{0});
        });
    }

    @Test
    void testParametricUsage6() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Gaussian.Parametric g = new Gaussian.Parametric();
            g.gradient(0, new double[]{0, 1, 0});
        });
    }

    @Test
    void testParametricValue() {
        final double norm = 2;
        final double mean = 3;
        final double sigma = 4;
        final Gaussian f = new Gaussian(norm, mean, sigma);

        final Gaussian.Parametric g = new Gaussian.Parametric();
        assertEquals(f.value(-1), g.value(-1, new double[] {norm, mean, sigma}), 0);
        assertEquals(f.value(0), g.value(0, new double[] {norm, mean, sigma}), 0);
        assertEquals(f.value(2), g.value(2, new double[] {norm, mean, sigma}), 0);
    }

    @Test
    void testParametricGradient() {
        final double norm = 2;
        final double mean = 3;
        final double sigma = 4;
        final Gaussian.Parametric f = new Gaussian.Parametric();

        final double x = 1;
        final double[] grad = f.gradient(1, new double[] {norm, mean, sigma});
        final double diff = x - mean;
        final double n = FastMath.exp(-diff * diff / (2 * sigma * sigma));
        assertEquals(n, grad[0], EPS);
        final double m = norm * n * diff / (sigma * sigma);
        assertEquals(m, grad[1], EPS);
        final double s = m * diff / sigma;
        assertEquals(s, grad[2], EPS);
    }
}
