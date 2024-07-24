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
package org.hipparchus.stat.descriptive.moment;

import org.hipparchus.stat.descriptive.StorelessUnivariateStatisticAbstractTest;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Test cases for the {@link Kurtosis} class.
 */
public class KurtosisTest extends StorelessUnivariateStatisticAbstractTest{

    @Override
    public Kurtosis getUnivariateStatistic() {
        return new Kurtosis();
    }

    @Override
    public double expectedValue() {
        return this.kurt;
    }

    /**
     * Make sure Double.NaN is returned iff n &lt; 4
     */
    @Test
    void testNaN() {
        Kurtosis kurt = getUnivariateStatistic();
        assertTrue(Double.isNaN(kurt.getResult()));
        kurt.increment(1d);
        assertTrue(Double.isNaN(kurt.getResult()));
        kurt.increment(1d);
        assertTrue(Double.isNaN(kurt.getResult()));
        kurt.increment(1d);
        assertTrue(Double.isNaN(kurt.getResult()));
        kurt.increment(1d);
        assertFalse(Double.isNaN(kurt.getResult()));
    }

}
