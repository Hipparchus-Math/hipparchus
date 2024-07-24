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

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.function.Expm1;
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;


/**
 * Test case for Divided Difference interpolator.
 * <p>
 * The error of polynomial interpolation is
 *     f(z) - p(z) = f^(n)(zeta) * (z-x[0])(z-x[1])...(z-x[n-1]) / n!
 * where f^(n) is the n-th derivative of the approximated function and
 * zeta is some point in the interval determined by x[] and z.
 * <p>
 * Since zeta is unknown, f^(n)(zeta) cannot be calculated. But we can bound
 * it and use the absolute value upper bound for estimates. For reference,
 * see <b>Introduction to Numerical Analysis</b>, ISBN 038795452X, chapter 2.
 *
 */
public final class DividedDifferenceInterpolatorTest {

    /**
     * Test of interpolator for the sine function.
     * <p>
     * |sin^(n)(zeta)| &lt;= 1.0, zeta in [0, 2*PI]
     */
    @Test
    public void testSinFunction() {
        UnivariateFunction f = new Sin();
        UnivariateInterpolator interpolator = new DividedDifferenceInterpolator();
        double[] x;
        double[] y;
        double z;
        double expected;
        double result;
        double tolerance;

        // 6 interpolating points on interval [0, 2*PI]
        int n = 6;
        double min = 0.0, max = 2 * FastMath.PI;
        x = new double[n];
        y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = min + i * (max - min) / n;
            y[i] = f.value(x[i]);
        }
        double derivativebound = 1.0;
        UnivariateFunction p = interpolator.interpolate(x, y);

        z = FastMath.PI / 4; expected = f.value(z); result = p.value(z);
        tolerance = FastMath.abs(derivativebound * partialerror(x, z));
        Assertions.assertEquals(expected, result, tolerance);

        z = FastMath.PI * 1.5; expected = f.value(z); result = p.value(z);
        tolerance = FastMath.abs(derivativebound * partialerror(x, z));
        Assertions.assertEquals(expected, result, tolerance);
    }

    /**
     * Test of interpolator for the exponential function.
     * <p>
     * |expm1^(n)(zeta)| &lt;= e, zeta in [-1, 1]
     */
    @Test
    public void testExpm1Function() {
        UnivariateFunction f = new Expm1();
        UnivariateInterpolator interpolator = new DividedDifferenceInterpolator();
        double[] x;
        double[] y;
        double z;
        double expected;
        double result;
        double tolerance;

        // 5 interpolating points on interval [-1, 1]
        int n = 5;
        double min = -1.0, max = 1.0;
        x = new double[n];
        y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = min + i * (max - min) / n;
            y[i] = f.value(x[i]);
        }
        double derivativebound = FastMath.E;
        UnivariateFunction p = interpolator.interpolate(x, y);

        z = 0.0; expected = f.value(z); result = p.value(z);
        tolerance = FastMath.abs(derivativebound * partialerror(x, z));
        Assertions.assertEquals(expected, result, tolerance);

        z = 0.5; expected = f.value(z); result = p.value(z);
        tolerance = FastMath.abs(derivativebound * partialerror(x, z));
        Assertions.assertEquals(expected, result, tolerance);

        z = -0.5; expected = f.value(z); result = p.value(z);
        tolerance = FastMath.abs(derivativebound * partialerror(x, z));
        Assertions.assertEquals(expected, result, tolerance);
    }

    /**
     * Test of parameters for the interpolator.
     */
    @Test
    public void testParameters() {
        UnivariateInterpolator interpolator = new DividedDifferenceInterpolator();

        try {
            // bad abscissas array
            double[] x = { 1.0, 2.0, 2.0, 4.0 };
            double[] y = { 0.0, 4.0, 4.0, 2.5 };
            UnivariateFunction p = interpolator.interpolate(x, y);
            p.value(0.0);
            Assertions.fail("Expecting MathIllegalArgumentException - bad abscissas array");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    /**
     * Returns the partial error term (z-x[0])(z-x[1])...(z-x[n-1])/n!
     */
    protected double partialerror(double[] x, double z) throws
        IllegalArgumentException {

        if (x.length < 1) {
            throw new IllegalArgumentException
                ("Interpolation array cannot be empty.");
        }
        double out = 1;
        for (int i = 0; i < x.length; i++) {
            out *= (z - x[i]) / (i + 1);
        }
        return out;
    }
}
