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
 * Test cases for LogisticsDistribution.
 */
public class LogisticsDistributionTest extends RealDistributionAbstractTest {

    @Test
    void testParameters() {
        LogisticDistribution d = makeDistribution();
        assertEquals(2, d.getLocation(), Precision.EPSILON);
        assertEquals(5, d.getScale(), Precision.EPSILON);
    }

    @Test
    void testSupport() {
        LogisticDistribution d = makeDistribution();
        assertTrue(Double.isInfinite(d.getSupportLowerBound()));
        assertTrue(Double.isInfinite(d.getSupportUpperBound()));
        assertTrue(d.isSupportConnected());
    }

    @Override
    public LogisticDistribution makeDistribution() {
        return new LogisticDistribution(2, 5);
    }

    @Override
    public double[] makeCumulativeTestPoints() {
        return new double[] {
            -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5
        };
    }

    @Override
    public double[] makeDensityTestValues() {
        return new double[] {
            0.03173698, 0.03557889, 0.03932239, 0.04278194, 0.04575685, 0.04805215,
            0.04950331, 0.05000000, 0.04950331, 0.04805215, 0.04575685
        };
    }

    @Override
    public double[] makeCumulativeTestValues() {
        return new double[] {
            0.1978161, 0.2314752, 0.2689414, 0.3100255, 0.3543437, 0.4013123,
            0.4501660, 0.5000000, 0.5498340, 0.5986877, 0.6456563
        };
    }

}

