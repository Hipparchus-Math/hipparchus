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
package org.hipparchus.optim.univariate;

import org.hipparchus.analysis.QuinticFunction;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.optim.MaxEval;
import org.hipparchus.optim.nonlinear.scalar.GoalType;
import org.hipparchus.random.JDKRandomGenerator;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

class MultiStartUnivariateOptimizerTest {
    @Test
    void testMissingMaxEval() {
        assertThrows(MathIllegalStateException.class, () -> {
            UnivariateOptimizer underlying = new BrentOptimizer(1e-10, 1e-14);
            JDKRandomGenerator g = new JDKRandomGenerator();
            g.setSeed(44428400075l);
            MultiStartUnivariateOptimizer optimizer = new MultiStartUnivariateOptimizer(underlying, 10, g);
            optimizer.optimize(new UnivariateObjectiveFunction(new Sin()),
                GoalType.MINIMIZE,
                new SearchInterval(-1, 1));
        });
    }

    @Test
    void testMissingSearchInterval() {
        assertThrows(MathIllegalStateException.class, () -> {
            UnivariateOptimizer underlying = new BrentOptimizer(1e-10, 1e-14);
            JDKRandomGenerator g = new JDKRandomGenerator();
            g.setSeed(44428400075l);
            MultiStartUnivariateOptimizer optimizer = new MultiStartUnivariateOptimizer(underlying, 10, g);
            optimizer.optimize(new MaxEval(300),
                new UnivariateObjectiveFunction(new Sin()),
                GoalType.MINIMIZE);
        });
    }

    @Test
    void testSinMin() {
        UnivariateFunction f = new Sin();
        UnivariateOptimizer underlying = new BrentOptimizer(1e-10, 1e-14);
        JDKRandomGenerator g = new JDKRandomGenerator();
        g.setSeed(44428400075l);
        MultiStartUnivariateOptimizer optimizer = new MultiStartUnivariateOptimizer(underlying, 10, g);
        optimizer.optimize(new MaxEval(300),
                           new UnivariateObjectiveFunction(f),
                           GoalType.MINIMIZE,
                           new SearchInterval(-100.0, 100.0));
        UnivariatePointValuePair[] optima = optimizer.getOptima();
        for (int i = 1; i < optima.length; ++i) {
            double d = (optima[i].getPoint() - optima[i-1].getPoint()) / (2 * FastMath.PI);
            assertTrue(FastMath.abs(d - FastMath.rint(d)) < 1.0e-8);
            assertEquals(-1.0, f.value(optima[i].getPoint()), 1.0e-10);
            assertEquals(f.value(optima[i].getPoint()), optima[i].getValue(), 1.0e-10);
        }
        assertTrue(optimizer.getEvaluations() > 200);
        assertTrue(optimizer.getEvaluations() < 300);
    }

    @Test
    void testQuinticMin() {
        // The quintic function has zeros at 0, +-0.5 and +-1.
        // The function has extrema (first derivative is zero) at 0.27195613 and 0.82221643,
        UnivariateFunction f = new QuinticFunction();
        UnivariateOptimizer underlying = new BrentOptimizer(1e-9, 1e-14);
        JDKRandomGenerator g = new JDKRandomGenerator();
        g.setSeed(4312000053L);
        MultiStartUnivariateOptimizer optimizer = new MultiStartUnivariateOptimizer(underlying, 5, g);

        UnivariatePointValuePair optimum
            = optimizer.optimize(new MaxEval(300),
                                 new UnivariateObjectiveFunction(f),
                                 GoalType.MINIMIZE,
                                 new SearchInterval(-0.3, -0.2));
        assertEquals(-0.27195613, optimum.getPoint(), 1e-9);
        assertEquals(-0.0443342695, optimum.getValue(), 1e-9);

        UnivariatePointValuePair[] optima = optimizer.getOptima();
        for (int i = 0; i < optima.length; ++i) {
            assertEquals(f.value(optima[i].getPoint()), optima[i].getValue(), 1e-9);
        }
        assertTrue(optimizer.getEvaluations() >= 50);
        assertTrue(optimizer.getEvaluations() <= 100);
    }

    @Test
    void testBadFunction() {
        UnivariateFunction f = new UnivariateFunction() {
                public double value(double x) {
                    if (x < 0) {
                        throw new LocalException();
                    }
                    return 0;
                }
            };
        UnivariateOptimizer underlying = new BrentOptimizer(1e-9, 1e-14);
        JDKRandomGenerator g = new JDKRandomGenerator();
        g.setSeed(4312000053L);
        MultiStartUnivariateOptimizer optimizer = new MultiStartUnivariateOptimizer(underlying, 5, g);

        try {
            optimizer.optimize(new MaxEval(300),
                               new UnivariateObjectiveFunction(f),
                               GoalType.MINIMIZE,
                               new SearchInterval(-0.3, -0.2));
            fail();
        } catch (LocalException e) {
            // Expected.
        }

        // Ensure that the exception was thrown because no optimum was found.
        assertNull(optimizer.getOptima()[0]);
    }

    private static class LocalException extends RuntimeException {
        private static final long serialVersionUID = 1194682757034350629L;
    }

}
