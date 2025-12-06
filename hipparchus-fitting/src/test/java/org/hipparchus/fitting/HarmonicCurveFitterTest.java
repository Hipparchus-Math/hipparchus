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

import org.hipparchus.analysis.function.HarmonicOscillator;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class HarmonicCurveFitterTest {
    /**
     * Zero points is not enough observed points.
     */
    @Test
    void testPreconditions1() {
        assertThrows(MathIllegalArgumentException.class, () -> HarmonicCurveFitter.create().fit(new WeightedObservedPoints().toList()));
    }

    @Test
    void testNoError() {
        final double a = 0.2;
        final double w = 3.4;
        final double p = 4.1;
        final HarmonicOscillator f = new HarmonicOscillator(a, w, p);

        final WeightedObservedPoints points = new WeightedObservedPoints();
        for (double x = 0.0; x < 1.3; x += 0.01) {
            points.add(1, x, f.value(x));
        }

        final HarmonicCurveFitter fitter = HarmonicCurveFitter.create();
        final double[] fitted = fitter.fit(points.toList());
        assertEquals(a, fitted[0], 1.0e-13);
        assertEquals(w, fitted[1], 1.0e-13);
        assertEquals(p, MathUtils.normalizeAngle(fitted[2], p), 1e-13);

        final HarmonicOscillator ff = new HarmonicOscillator(fitted[0], fitted[1], fitted[2]);
        for (double x = -1.0; x < 1.0; x += 0.01) {
            assertTrue(FastMath.abs(f.value(x) - ff.value(x)) < 1e-13);
        }
    }

    @Test
    void test1PercentError() {
        final Random randomizer = new Random(64925784252L);
        final double a = 0.2;
        final double w = 3.4;
        final double p = 4.1;
        final HarmonicOscillator f = new HarmonicOscillator(a, w, p);

        final WeightedObservedPoints points = new WeightedObservedPoints();
        for (double x = 0.0; x < 10.0; x += 0.1) {
            points.add(1, x, f.value(x) + 0.01 * randomizer.nextGaussian());
        }

        final HarmonicCurveFitter fitter = HarmonicCurveFitter.create();
        final double[] fitted = fitter.fit(points.toList());
        assertEquals(a, fitted[0], 7.6e-4);
        assertEquals(w, fitted[1], 2.7e-3);
        assertEquals(p, MathUtils.normalizeAngle(fitted[2], p), 1.3e-2);
    }

    @Test
    void testTinyVariationsData() {
        final Random randomizer = new Random(64925784252L);

        final WeightedObservedPoints points = new WeightedObservedPoints();
        for (double x = 0.0; x < 10.0; x += 0.1) {
            points.add(1, x, 1e-7 * randomizer.nextGaussian());
        }

        final HarmonicCurveFitter fitter = HarmonicCurveFitter.create();
        fitter.fit(points.toList());

        // This test serves to cover the part of the code of "guessAOmega"
        // when the algorithm using integrals fails.
    }

    @Test
    void testInitialGuess() {
        final Random randomizer = new Random(45314242L);
        final double a = 0.2;
        final double w = 3.4;
        final double p = 4.1;
        final HarmonicOscillator f = new HarmonicOscillator(a, w, p);

        final WeightedObservedPoints points = new WeightedObservedPoints();
        for (double x = 0.0; x < 10.0; x += 0.1) {
            points.add(1, x, f.value(x) + 0.01 * randomizer.nextGaussian());
        }

        final HarmonicCurveFitter fitter = HarmonicCurveFitter.create()
            .withStartPoint(new double[] { 0.15, 3.6, 4.5 });
        final double[] fitted = fitter.fit(points.toList());
        assertEquals(a, fitted[0], 1.2e-3);
        assertEquals(w, fitted[1], 3.3e-3);
        assertEquals(p, MathUtils.normalizeAngle(fitted[2], p), 1.7e-2);
    }

    @Test
    void testUnsorted() {
        Random randomizer = new Random(64925784252L);
        final double a = 0.2;
        final double w = 3.4;
        final double p = 4.1;
        final HarmonicOscillator f = new HarmonicOscillator(a, w, p);

        // Build a regularly spaced array of measurements.
        final int size = 100;
        final double[] xTab = new double[size];
        final double[] yTab = new double[size];
        for (int i = 0; i < size; i++) {
            xTab[i] = 0.1 * i;
            yTab[i] = f.value(xTab[i]) + 0.01 * randomizer.nextGaussian();
        }

        // shake it
        for (int i = 0; i < size; i++) {
            int i1 = randomizer.nextInt(size);
            int i2 = randomizer.nextInt(size);
            double xTmp = xTab[i1];
            double yTmp = yTab[i1];
            xTab[i1] = xTab[i2];
            yTab[i1] = yTab[i2];
            xTab[i2] = xTmp;
            yTab[i2] = yTmp;
        }

        // Pass it to the fitter.
        final WeightedObservedPoints points = new WeightedObservedPoints();
        for (int i = 0; i < size; ++i) {
            points.add(1, xTab[i], yTab[i]);
        }

        final HarmonicCurveFitter fitter = HarmonicCurveFitter.create();
        final double[] fitted = fitter.fit(points.toList());
        assertEquals(a, fitted[0], 7.6e-4);
        assertEquals(w, fitted[1], 3.5e-3);
        assertEquals(p, MathUtils.normalizeAngle(fitted[2], p), 1.5e-2);
    }

    @Test
    void testMath844() {
        assertThrows(MathIllegalStateException.class, () -> {
            final double[] y = {0, 1, 2, 3, 2, 1,
                0, -1, -2, -3, -2, -1,
                0, 1, 2, 3, 2, 1,
                0, -1, -2, -3, -2, -1,
                0, 1, 2, 3, 2, 1, 0};
            final List<WeightedObservedPoint> points = new ArrayList<WeightedObservedPoint>();
            for (int i = 0; i < y.length; i++) {
                points.add(new WeightedObservedPoint(1, i, y[i]));
            }

            // The guesser fails because the function is far from an harmonic
            // function: It is a triangular periodic function with amplitude 3
            // and period 12, and all sample points are taken at integer abscissae
            // so function values all belong to the integer subset {-3, -2, -1, 0,
            // 1, 2, 3}.
            new HarmonicCurveFitter.ParameterGuesser(points);
        });
    }
}
