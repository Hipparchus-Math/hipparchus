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

import org.hipparchus.analysis.FunctionUtils;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.analysis.differentiation.UnivariateDifferentiableFunction;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Test for class {@link Logit}.
 */
public class LogitTest {
    private final double EPS = Math.ulp(1d);

    @Test
    public void testPreconditions1() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final double lo = -1;
            final double hi = 2;
            final UnivariateFunction f = new Logit(lo, hi);

            f.value(lo - 1);
        });
    }

    @Test
    public void testPreconditions2() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final double lo = -1;
            final double hi = 2;
            final UnivariateFunction f = new Logit(lo, hi);

            f.value(hi + 1);
        });
    }

    @Test
    public void testSomeValues() {
        final double lo = 1;
        final double hi = 2;
        final UnivariateFunction f = new Logit(lo, hi);

        Assertions.assertEquals(Double.NEGATIVE_INFINITY, f.value(1), EPS);
        Assertions.assertEquals(Double.POSITIVE_INFINITY, f.value(2), EPS);
        Assertions.assertEquals(0, f.value(1.5), EPS);
    }

    @Test
    public void testDerivative() {
        final double lo = 1;
        final double hi = 2;
        final Logit f = new Logit(lo, hi);
        final DerivativeStructure f15 = f.value(new DSFactory(1, 1).variable(0, 1.5));

        Assertions.assertEquals(4, f15.getPartialDerivative(1), EPS);
    }

    @Test
    public void testDerivativeLargeArguments() {
        final Logit f = new Logit(1, 2);

        DSFactory factory = new DSFactory(1, 1);
        for (double arg : new double[] {
            Double.NEGATIVE_INFINITY, -Double.MAX_VALUE, -1e155, 1e155, Double.MAX_VALUE, Double.POSITIVE_INFINITY
            }) {
            try {
                f.value(factory.variable(0, arg));
                Assertions.fail("an exception should have been thrown");
            } catch (MathIllegalArgumentException ore) {
                // expected
            } catch (Exception e) {
                Assertions.fail("wrong exception caught: " + e.getMessage());
            }
        }
    }

    @Test
    public void testDerivativesHighOrder() {
        DerivativeStructure l = new Logit(1, 3).value(new DSFactory(1, 5).variable(0, 1.2));
        Assertions.assertEquals(-2.1972245773362193828, l.getPartialDerivative(0), 1.0e-16);
        Assertions.assertEquals(5.5555555555555555555,  l.getPartialDerivative(1), 9.0e-16);
        Assertions.assertEquals(-24.691358024691358025, l.getPartialDerivative(2), 2.0e-14);
        Assertions.assertEquals(250.34293552812071331,  l.getPartialDerivative(3), 2.0e-13);
        Assertions.assertEquals(-3749.4284407864654778, l.getPartialDerivative(4), 4.0e-12);
        Assertions.assertEquals(75001.270131585632282,  l.getPartialDerivative(5), 8.0e-11);
    }

    @Test
    public void testParametricUsage1() {
        assertThrows(NullArgumentException.class, () -> {
            final Logit.Parametric g = new Logit.Parametric();
            g.value(0, null);
        });
    }

    @Test
    public void testParametricUsage2() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Logit.Parametric g = new Logit.Parametric();
            g.value(0, new double[]{0});
        });
    }

    @Test
    public void testParametricUsage3() {
        assertThrows(NullArgumentException.class, () -> {
            final Logit.Parametric g = new Logit.Parametric();
            g.gradient(0, null);
        });
    }

    @Test
    public void testParametricUsage4() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Logit.Parametric g = new Logit.Parametric();
            g.gradient(0, new double[]{0});
        });
    }

    @Test
    public void testParametricUsage5() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Logit.Parametric g = new Logit.Parametric();
            g.value(-1, new double[]{0, 1});
        });
    }

    @Test
    public void testParametricUsage6() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Logit.Parametric g = new Logit.Parametric();
            g.value(2, new double[]{0, 1});
        });
    }

    @Test
    public void testParametricValue() {
        final double lo = 2;
        final double hi = 3;
        final Logit f = new Logit(lo, hi);

        final Logit.Parametric g = new Logit.Parametric();
        Assertions.assertEquals(f.value(2), g.value(2, new double[] {lo, hi}), 0);
        Assertions.assertEquals(f.value(2.34567), g.value(2.34567, new double[] {lo, hi}), 0);
        Assertions.assertEquals(f.value(3), g.value(3, new double[] {lo, hi}), 0);
    }

    @Test
    public void testValueWithInverseFunction() {
        final double lo = 2;
        final double hi = 3;
        final Logit f = new Logit(lo, hi);
        final Sigmoid g = new Sigmoid(lo, hi);
        RandomGenerator random = new Well1024a(0x49914cdd9f0b8db5l);
        final UnivariateDifferentiableFunction id = FunctionUtils.compose((UnivariateDifferentiableFunction) g,
                                                                (UnivariateDifferentiableFunction) f);

        DSFactory factory = new DSFactory(1, 1);
        for (int i = 0; i < 10; i++) {
            final double x = lo + random.nextDouble() * (hi - lo);
            Assertions.assertEquals(x, id.value(factory.variable(0, x)).getValue(), EPS);
        }

        Assertions.assertEquals(lo, id.value(factory.variable(0, lo)).getValue(), EPS);
        Assertions.assertEquals(hi, id.value(factory.variable(0, hi)).getValue(), EPS);
    }

    @Test
    public void testDerivativesWithInverseFunction() {
        double[] epsilon = new double[] { 1.0e-20, 4.0e-16, 3.0e-15, 2.0e-11, 3.0e-9, 1.0e-6 };
        final double lo = 2;
        final double hi = 3;
        final Logit f = new Logit(lo, hi);
        final Sigmoid g = new Sigmoid(lo, hi);
        RandomGenerator random = new Well1024a(0x96885e9c1f81cea5l);
        final UnivariateDifferentiableFunction id =
                FunctionUtils.compose((UnivariateDifferentiableFunction) g, (UnivariateDifferentiableFunction) f);
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            DSFactory factory = new DSFactory(1, maxOrder);
            double max = 0;
            for (int i = 0; i < 10; i++) {
                final double x = lo + random.nextDouble() * (hi - lo);
                final DerivativeStructure dsX = factory.variable(0, x);
                max = FastMath.max(max, FastMath.abs(dsX.getPartialDerivative(maxOrder) -
                                                     id.value(dsX).getPartialDerivative(maxOrder)));
                Assertions.assertEquals(dsX.getPartialDerivative(maxOrder),
                                    id.value(dsX).getPartialDerivative(maxOrder),
                                    epsilon[maxOrder]);
            }

            // each function evaluates correctly near boundaries,
            // but combination leads to NaN as some intermediate point is infinite
            final DerivativeStructure dsLo = factory.variable(0, lo);
            if (maxOrder == 0) {
                Assertions.assertTrue(Double.isInfinite(f.value(dsLo).getPartialDerivative(maxOrder)));
                Assertions.assertEquals(lo, id.value(dsLo).getPartialDerivative(maxOrder), epsilon[maxOrder]);
            } else if (maxOrder == 1) {
                Assertions.assertTrue(Double.isInfinite(f.value(dsLo).getPartialDerivative(maxOrder)));
                Assertions.assertTrue(Double.isNaN(id.value(dsLo).getPartialDerivative(maxOrder)));
            } else {
                Assertions.assertTrue(Double.isNaN(f.value(dsLo).getPartialDerivative(maxOrder)));
                Assertions.assertTrue(Double.isNaN(id.value(dsLo).getPartialDerivative(maxOrder)));
            }

            final DerivativeStructure dsHi = factory.variable(0, hi);
            if (maxOrder == 0) {
                Assertions.assertTrue(Double.isInfinite(f.value(dsHi).getPartialDerivative(maxOrder)));
                Assertions.assertEquals(hi, id.value(dsHi).getPartialDerivative(maxOrder), epsilon[maxOrder]);
            } else if (maxOrder == 1) {
                Assertions.assertTrue(Double.isInfinite(f.value(dsHi).getPartialDerivative(maxOrder)));
                Assertions.assertTrue(Double.isNaN(id.value(dsHi).getPartialDerivative(maxOrder)));
            } else {
                Assertions.assertTrue(Double.isNaN(f.value(dsHi).getPartialDerivative(maxOrder)));
                Assertions.assertTrue(Double.isNaN(id.value(dsHi).getPartialDerivative(maxOrder)));
            }

        }
    }
}
