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

import org.hipparchus.analysis.QuinticFunction;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;


/**
 * Test case for trapezoid integrator.
 * <p>
 * Test runs show that for a default relative accuracy of 1E-6, it
 * generally takes 10 to 15 iterations for the integral to converge.
 *
 */
final class TrapezoidIntegratorTest {

    /**
     * Test of integrator for the sine function.
     */
    @Test
    void testSinFunction() {
        UnivariateFunction f = new Sin();
        UnivariateIntegrator integrator = new TrapezoidIntegrator();
        double min, max, expected, result, tolerance;

        min = 0; max = FastMath.PI; expected = 2;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(10000, f, min, max);
        assertTrue(integrator.getEvaluations() < 2500);
        assertTrue(integrator.getIterations()  < 15);
        assertEquals(expected, result, tolerance);

        min = -FastMath.PI/3; max = 0; expected = -0.5;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(10000, f, min, max);
        assertTrue(integrator.getEvaluations() < 2500);
        assertTrue(integrator.getIterations()  < 15);
        assertEquals(expected, result, tolerance);
    }

    /**
     * Test of integrator for the quintic function.
     */
    @Test
    void testQuinticFunction() {
        UnivariateFunction f = new QuinticFunction();
        UnivariateIntegrator integrator = new TrapezoidIntegrator();
        double min, max, expected, result, tolerance;

        min = 0; max = 1; expected = -1.0/48;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(10000, f, min, max);
        assertTrue(integrator.getEvaluations() < 5000);
        assertTrue(integrator.getIterations()  < 15);
        assertEquals(expected, result, tolerance);

        min = 0; max = 0.5; expected = 11.0/768;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(10000, f, min, max);
        assertTrue(integrator.getEvaluations() < 2500);
        assertTrue(integrator.getIterations()  < 15);
        assertEquals(expected, result, tolerance);

        min = -1; max = 4; expected = 2048/3.0 - 78 + 1.0/48;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(10000, f, min, max);
        assertTrue(integrator.getEvaluations() < 5000);
        assertTrue(integrator.getIterations()  < 15);
        assertEquals(expected, result, tolerance);

    }

    /**
     * Test of parameters for the integrator.
     */
    @Test
    void testParameters() {
        UnivariateFunction f = new Sin();

        try {
            // bad interval
            new TrapezoidIntegrator().integrate(1000, f, 1, -1);
            fail("Expecting MathIllegalArgumentException - bad interval");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            // bad iteration limits
            new TrapezoidIntegrator(5, 4);
            fail("Expecting MathIllegalArgumentException - bad iteration limits");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            // bad iteration limits
            new TrapezoidIntegrator(10,99);
            fail("Expecting MathIllegalArgumentException - bad iteration limits");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }
}
