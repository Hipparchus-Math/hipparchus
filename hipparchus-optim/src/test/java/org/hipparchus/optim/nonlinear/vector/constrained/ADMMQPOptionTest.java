/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.optim.nonlinear.vector.constrained;


import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class ADMMQPOptionTest {

    @Test
    void testSettersDouble() {
        // GIVEN
        final ADMMQPOption option = new ADMMQPOption();
        // WHEN
        final double expectedValue = 1.;
        option.setAlpha(expectedValue);
        option.setEps(expectedValue);
        option.setEpsInfeasible(expectedValue);
        option.setRhoMax(expectedValue);
        option.setRhoMin(expectedValue);
        option.setRhoMax(expectedValue);
        option.setSigma(expectedValue);
        // THEN
        final double tolerance = 0.;
        assertEquals(expectedValue, option.getAlpha(), tolerance);
        assertEquals(expectedValue, option.getEps(), tolerance);
        assertEquals(expectedValue, option.getEpsInfeasible(), tolerance);
        assertEquals(expectedValue, option.getRhoMax(), tolerance);
        assertEquals(expectedValue, option.getRhoMin(), tolerance);
        assertEquals(expectedValue, option.getSigma(), tolerance);
    }

    @Test
    void testSettersBoolean() {
        // GIVEN
        final ADMMQPOption option = new ADMMQPOption();
        // WHEN
        final boolean expectedValue = true;
        option.setPolishing(expectedValue);
        option.setUpdateRho(expectedValue);
        option.setScaling(expectedValue);
        // THEN
        assertEquals(expectedValue, option.isPolishing());
        assertEquals(expectedValue, option.updateRho());
        assertEquals(expectedValue, option.isScaling());
    }

    @Test
    void testSettersIntegers() {
        // GIVEN
        final ADMMQPOption option = new ADMMQPOption();
        // WHEN
        final int expectedValue = 10;
        option.setMaxRhoIteration(expectedValue);
        option.setPolishingIteration(expectedValue);
        option.setScaleMaxIteration(expectedValue);
        // THEN
        assertEquals(expectedValue, option.getMaxRhoIteration());
        assertEquals(expectedValue, option.getPolishIteration());
        assertEquals(expectedValue, option.getScaleMaxIteration());
    }

}
