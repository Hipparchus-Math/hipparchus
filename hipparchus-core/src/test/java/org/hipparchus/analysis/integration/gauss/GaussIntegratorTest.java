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
import org.hipparchus.analysis.function.Constant;
import org.hipparchus.util.Pair;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test for {@link GaussIntegrator} class.
 *
 */
class GaussIntegratorTest {
    @Test
    void testGetWeights() {
        final double[] points = { 0, 1.2, 3.4 };
        final double[] weights = { 9.8, 7.6, 5.4 };

        final GaussIntegrator integrator
            = new GaussIntegrator(new Pair<double[], double[]>(points, weights));

        assertEquals(weights.length, integrator.getNumberOfPoints());

        for (int i = 0; i < integrator.getNumberOfPoints(); i++) {
            assertEquals(weights[i], integrator.getWeight(i), 0d);
        }
    }

    @Test
    void testGetPoints() {
        final double[] points = { 0, 1.2, 3.4 };
        final double[] weights = { 9.8, 7.6, 5.4 };

        final GaussIntegrator integrator
            = new GaussIntegrator(new Pair<double[], double[]>(points, weights));

        assertEquals(points.length, integrator.getNumberOfPoints());

        for (int i = 0; i < integrator.getNumberOfPoints(); i++) {
            assertEquals(points[i], integrator.getPoint(i), 0d);
        }
    }

    @Test
    void testIntegrate() {
        final double[] points = { 0, 1, 2, 3, 4, 5 };
        final double[] weights = { 1, 1, 1, 1, 1, 1 };

        final GaussIntegrator integrator
            = new GaussIntegrator(new Pair<double[], double[]>(points, weights));

        final double val = 123.456;
        final UnivariateFunction c = new Constant(val);

        final double s = integrator.integrate(c);
        assertEquals(points.length * val, s, 0d);
    }
}
