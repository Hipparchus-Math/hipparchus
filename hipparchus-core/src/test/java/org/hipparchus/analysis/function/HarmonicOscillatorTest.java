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
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Test for class {@link HarmonicOscillator}.
 */
public class HarmonicOscillatorTest {
    private final double EPS = Math.ulp(1d);

    @Test
    public void testSomeValues() {
        final double a = -1.2;
        final double w = 0.34;
        final double p = 5.6;
        final UnivariateFunction f = new HarmonicOscillator(a, w, p);

        final double d = 0.12345;
        for (int i = 0; i < 10; i++) {
            final double v = i * d;
            Assertions.assertEquals(a * FastMath.cos(w * v + p), f.value(v), 0);
        }
    }

    @Test
    public void testDerivative() {
        final double a = -1.2;
        final double w = 0.34;
        final double p = 5.6;
        final HarmonicOscillator f = new HarmonicOscillator(a, w, p);

        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            final DSFactory factory = new DSFactory(1, maxOrder);
            final double d = 0.12345;
            for (int i = 0; i < 10; i++) {
                final double v = i * d;
                final DerivativeStructure h = f.value(factory.variable(0, v));
                for (int k = 0; k <= maxOrder; ++k) {
                    final double trigo;
                    switch (k % 4) {
                        case 0:
                            trigo = +FastMath.cos(w * v + p);
                            break;
                        case 1:
                            trigo = -FastMath.sin(w * v + p);
                            break;
                        case 2:
                            trigo = -FastMath.cos(w * v + p);
                            break;
                        default:
                            trigo = +FastMath.sin(w * v + p);
                            break;
                    }
                    Assertions.assertEquals(a * FastMath.pow(w, k) * trigo,
                                        h.getPartialDerivative(k),
                                        Precision.EPSILON);
                }
            }
        }
    }

    @Test
    public void testParametricUsage1() {
        assertThrows(NullArgumentException.class, () -> {
            final HarmonicOscillator.Parametric g = new HarmonicOscillator.Parametric();
            g.value(0, null);
        });
    }

    @Test
    public void testParametricUsage2() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final HarmonicOscillator.Parametric g = new HarmonicOscillator.Parametric();
            g.value(0, new double[]{0});
        });
    }

    @Test
    public void testParametricUsage3() {
        assertThrows(NullArgumentException.class, () -> {
            final HarmonicOscillator.Parametric g = new HarmonicOscillator.Parametric();
            g.gradient(0, null);
        });
    }

    @Test
    public void testParametricUsage4() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final HarmonicOscillator.Parametric g = new HarmonicOscillator.Parametric();
            g.gradient(0, new double[]{0});
        });
    }

    @Test
    public void testParametricValue() {
        final double amplitude = 2;
        final double omega = 3;
        final double phase = 4;
        final HarmonicOscillator f = new HarmonicOscillator(amplitude, omega, phase);

        final HarmonicOscillator.Parametric g = new HarmonicOscillator.Parametric();
        Assertions.assertEquals(f.value(-1), g.value(-1, new double[] {amplitude, omega, phase}), 0);
        Assertions.assertEquals(f.value(0), g.value(0, new double[] {amplitude, omega, phase}), 0);
        Assertions.assertEquals(f.value(2), g.value(2, new double[] {amplitude, omega, phase}), 0);
    }

    @Test
    public void testParametricGradient() {
        final double amplitude = 2;
        final double omega = 3;
        final double phase = 4;
        final HarmonicOscillator.Parametric f = new HarmonicOscillator.Parametric();

        final double x = 1;
        final double[] grad = f.gradient(1, new double[] {amplitude, omega, phase});
        final double xTimesOmegaPlusPhase = omega * x + phase;
        final double a = FastMath.cos(xTimesOmegaPlusPhase);
        Assertions.assertEquals(a, grad[0], EPS);
        final double w = -amplitude * x * FastMath.sin(xTimesOmegaPlusPhase);
        Assertions.assertEquals(w, grad[1], EPS);
        final double p = -amplitude * FastMath.sin(xTimesOmegaPlusPhase);
        Assertions.assertEquals(p, grad[2], EPS);
    }
}
