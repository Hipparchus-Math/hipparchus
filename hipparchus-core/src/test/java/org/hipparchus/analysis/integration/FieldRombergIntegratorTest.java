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
package org.hipparchus.analysis.integration;

import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;


/**
 * Test case for Romberg integrator.
 * <p>
 * Romberg algorithm is very fast for good behavior integrand. Test runs
 * show that for a default relative accuracy of 1E-6, it generally takes
 * takes less than 5 iterations for the integral to converge.
 *
 */
public final class FieldRombergIntegratorTest {

    /**
     * Test of integrator for the sine function.
     */
    @Test
    public void testSinFunction() {
        FieldUnivariateIntegrator<Binary64> integrator = new FieldRombergIntegrator<>(Binary64Field.getInstance());

        Binary64 min = new Binary64(0);
        Binary64 max = new Binary64(FastMath.PI);
        double expected = 2;
        double tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        double result = integrator.integrate(100, x -> x.sin(), min, max).getReal();
        Assert.assertTrue(integrator.getEvaluations() < 50);
        Assert.assertTrue(integrator.getIterations()  < 10);
        Assert.assertEquals(expected, result, tolerance);

        min = new Binary64(-FastMath.PI/3);
        max = new Binary64(0);
        expected = -0.5;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(100, x -> x.sin(), min, max).getReal();
        Assert.assertTrue(integrator.getEvaluations() < 50);
        Assert.assertTrue(integrator.getIterations()  < 10);
        Assert.assertEquals(expected, result, tolerance);
    }

    /**
     * Test of integrator for the quintic function.
     */
    @Test
    public void testQuinticFunction() {
        CalculusFieldUnivariateFunction<Binary64> f =
                        t -> t.subtract(1).multiply(t.subtract(0.5)).multiply(t).multiply(t.add(0.5)).multiply(t.add(1));
        FieldUnivariateIntegrator<Binary64> integrator = new FieldRombergIntegrator<>(Binary64Field.getInstance());

        Binary64 min = new Binary64(0);
        Binary64 max = new Binary64(1);
        double expected = -1.0 / 48;
        double tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        double result = integrator.integrate(100, f, min, max).getReal();
        Assert.assertTrue(integrator.getEvaluations() < 10);
        Assert.assertTrue(integrator.getIterations()  < 5);
        Assert.assertEquals(expected, result, tolerance);

        min = new Binary64(0);
        max = new Binary64(0.5);
        expected = 11.0 / 768;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(100, f, min, max).getReal();
        Assert.assertTrue(integrator.getEvaluations() < 10);
        Assert.assertTrue(integrator.getIterations()  < 5);
        Assert.assertEquals(expected, result, tolerance);

        min = new Binary64(-1);
        max = new Binary64(4);
        expected = 2048 / 3.0 - 78 + 1.0 / 48;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(100, f, min, max).getReal();
        Assert.assertTrue(integrator.getEvaluations() < 10);
        Assert.assertTrue(integrator.getIterations()  < 5);
        Assert.assertEquals(expected, result, tolerance);
    }

    /**
     * Test of parameters for the integrator.
     */
    @Test
    public void testParameters() {

        try {
            // bad interval
            new FieldRombergIntegrator<>(Binary64Field.getInstance()).integrate(1000, x -> x.sin(),
                                                                                 new Binary64(1), new Binary64(-1));
            Assert.fail("Expecting MathIllegalArgumentException - bad interval");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            // bad iteration limits
            new FieldRombergIntegrator<>(Binary64Field.getInstance(), 5, 4);
            Assert.fail("Expecting MathIllegalArgumentException - bad iteration limits");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            // bad iteration limits
            new FieldRombergIntegrator<>(Binary64Field.getInstance(), 10, 50);
            Assert.fail("Expecting MathIllegalArgumentException - bad iteration limits");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }
}
