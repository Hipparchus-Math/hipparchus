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
package org.hipparchus.util;

import org.hipparchus.exception.MathIllegalStateException;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Tests for ContinuedFraction.
 */
class ContinuedFractionTest {

    @Test
    void testGoldenRatio() throws Exception {
        ContinuedFraction cf = new ContinuedFraction() {

            @Override
            public double getA(int n, double x) {
                return 1.0;
            }

            @Override
            public double getB(int n, double x) {
                return 1.0;
            }
        };

        double gr = cf.evaluate(0.0, 10e-9);
        assertEquals(1.61803399, gr, 10e-9);
    }

    @Test
    void testNonConvergentContinuedFraction() {
        assertThrows(MathIllegalStateException.class, () -> {
            ContinuedFraction cf = new ContinuedFraction() {

                @Override
                public double getA(int n, double x) {
                    return 1.0;
                }

                @Override
                public double getB(int n, double x) {
                    return 1.0;
                }

            };

            cf.evaluate(0.0, 10e-9, 10);
        });
    }

    @Test
    void testInfinityDivergence() {
        assertThrows(MathIllegalStateException.class, () -> {
            ContinuedFraction cf = new ContinuedFraction() {

                @Override
                public double getA(int n, double x) {
                    return 1. / n;
                }

                @Override
                public double getB(int n, double x) {
                    return 1.0;
                }

            };

            cf.evaluate(1);
        });
    }
}
