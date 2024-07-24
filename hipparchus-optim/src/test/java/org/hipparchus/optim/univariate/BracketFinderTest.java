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

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.optim.nonlinear.scalar.GoalType;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Test for {@link BracketFinder}.
 */
public class BracketFinderTest {

    @Test
    public void testCubicMin() {
        final BracketFinder bFind = new BracketFinder();
        final UnivariateFunction func = new UnivariateFunction() {
                public double value(double x) {
                    if (x < -2) {
                        return value(-2);
                    }
                    else  {
                        return (x - 1) * (x + 2) * (x + 3);
                    }
                }
            };

        bFind.search(func, GoalType.MINIMIZE, -2 , -1);
        final double tol = 1e-15;
        // Comparing with results computed in Python.
        Assertions.assertEquals(-2, bFind.getLo(), tol);
        Assertions.assertEquals(-1, bFind.getMid(), tol);
        Assertions.assertEquals(0.61803399999999997, bFind.getHi(), tol);
    }

    @Test
    public void testCubicMax() {
        final BracketFinder bFind = new BracketFinder();
        final UnivariateFunction func = new UnivariateFunction() {
                public double value(double x) {
                    if (x < -2) {
                        return value(-2);
                    }
                    else  {
                        return -(x - 1) * (x + 2) * (x + 3);
                    }
                }
            };

        bFind.search(func, GoalType.MAXIMIZE, -2 , -1);
        final double tol = 1e-15;
        Assertions.assertEquals(-2, bFind.getLo(), tol);
        Assertions.assertEquals(-1, bFind.getMid(), tol);
        Assertions.assertEquals(0.61803399999999997, bFind.getHi(), tol);
    }

    @Test
    public void testMinimumIsOnIntervalBoundary() {
        final UnivariateFunction func = new UnivariateFunction() {
                public double value(double x) {
                    return x * x;
                }
            };

        final BracketFinder bFind = new BracketFinder();

        bFind.search(func, GoalType.MINIMIZE, 0, 1);
        Assertions.assertTrue(bFind.getLo() <= 0);
        Assertions.assertTrue(0 <= bFind.getHi());

        bFind.search(func, GoalType.MINIMIZE, -1, 0);
        Assertions.assertTrue(bFind.getLo() <= 0);
        Assertions.assertTrue(0 <= bFind.getHi());
    }

    @Test
    public void testIntervalBoundsOrdering() {
        final UnivariateFunction func = new UnivariateFunction() {
                public double value(double x) {
                    return x * x;
                }
            };

        final BracketFinder bFind = new BracketFinder();

        bFind.search(func, GoalType.MINIMIZE, -1, 1);
        Assertions.assertTrue(bFind.getLo() <= 0);
        Assertions.assertTrue(0 <= bFind.getHi());

        bFind.search(func, GoalType.MINIMIZE, 1, -1);
        Assertions.assertTrue(bFind.getLo() <= 0);
        Assertions.assertTrue(0 <= bFind.getHi());

        bFind.search(func, GoalType.MINIMIZE, 1, 2);
        Assertions.assertTrue(bFind.getLo() <= 0);
        Assertions.assertTrue(0 <= bFind.getHi());

        bFind.search(func, GoalType.MINIMIZE, 2, 1);
        Assertions.assertTrue(bFind.getLo() <= 0);
        Assertions.assertTrue(0 <= bFind.getHi());
    }
}
