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
package org.hipparchus.distribution.continuous;

import org.hipparchus.util.Precision;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Test cases for NakagamiDistribution.
 */
public class NakagamiDistributionTest extends RealDistributionAbstractTest {

    @Test
    void testParameters() {
        NakagamiDistribution d = makeDistribution();
        assertEquals(0.5, d.getShape(), Precision.EPSILON);
        assertEquals(1, d.getScale(), Precision.EPSILON);
    }

    @Test
    void testSupport() {
        NakagamiDistribution d = makeDistribution();
        assertEquals(0, d.getSupportLowerBound(), Precision.EPSILON);
        assertTrue(Double.isInfinite(d.getSupportUpperBound()));
        assertTrue(d.isSupportConnected());
    }

    @Override
    public NakagamiDistribution makeDistribution() {
        return new NakagamiDistribution(0.5, 1);
    }

    @Override
    public double[] makeCumulativeTestPoints() {
        return new double[] {
            0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2
        };
    }

    @Override
    public double[] makeDensityTestValues() {
        return new double[] {
            0.0000000, 0.7820854, 0.7365403, 0.6664492, 0.5793831, 0.4839414,
            0.3883721, 0.2994549, 0.2218417, 0.1579003, 0.1079819
        };
    }

    @Override
    public double[] makeCumulativeTestValues() {
        return new double[] {
            0.0000000, 0.1585194, 0.3108435, 0.4514938, 0.5762892, 0.6826895,
            0.7698607, 0.8384867, 0.8904014, 0.9281394, 0.9544997
        };
    }

}

