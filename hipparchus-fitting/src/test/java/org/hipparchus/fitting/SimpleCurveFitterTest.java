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
package org.hipparchus.fitting;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.analysis.ParametricUnivariateFunction;
import org.hipparchus.analysis.polynomials.PolynomialFunction;
import org.hipparchus.random.RandomDataGenerator;
import org.junit.jupiter.api.Test;

import java.util.Random;

/**
 * Test for class {@link SimpleCurveFitter}.
 */
class SimpleCurveFitterTest {
    @Test
    void testPolynomialFit() {
        final Random randomizer = new Random(53882150042L);
        final RandomDataGenerator randomDataGenerator = new RandomDataGenerator(64925784252L);

        final double[] coeff = { 12.9, -3.4, 2.1 }; // 12.9 - 3.4 x + 2.1 x^2
        final PolynomialFunction f = new PolynomialFunction(coeff);

        // Collect data from a known polynomial.
        final WeightedObservedPoints obs = new WeightedObservedPoints();
        for (int i = 0; i < 100; i++) {
            final double x = randomDataGenerator.nextUniform(-100, 100);
            obs.add(x, f.value(x) + 0.1 * randomizer.nextGaussian());
        }

        final ParametricUnivariateFunction function = new PolynomialFunction.Parametric();
        // Start fit from initial guesses that are far from the optimal values.
        final SimpleCurveFitter fitter
            = SimpleCurveFitter.create(function,
                                       new double[] { -1e20, 3e15, -5e25 });
        final double[] best = fitter.fit(obs.toList());

        UnitTestUtils.customAssertEquals("best != coeff", coeff, best, 2e-2);
    }
}
