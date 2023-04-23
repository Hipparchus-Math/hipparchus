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

import static org.junit.Assert.assertTrue;

import org.hipparchus.stat.descriptive.StorelessUnivariateStatisticAbstractTest;
import org.junit.Test;

/**
 * Test cases for the {@link FirstMoment} class.
 */
public class FirstMomentTest extends StorelessUnivariateStatisticAbstractTest {

    @Override
    public FirstMoment getUnivariateStatistic() {
        return new FirstMoment();
    }

    @Override
    public double expectedValue() {
        return this.mean;
    }

    @Test
    public void testSpecialValues() {
        final FirstMoment mean = new FirstMoment();

        mean.clear();
        mean.increment(Double.POSITIVE_INFINITY);
        mean.increment(1d);
        assertTrue(Double.isNaN(mean.getResult()));

        mean.clear();
        mean.increment(Double.POSITIVE_INFINITY);
        mean.increment(-1d);
        assertTrue(Double.isNaN(mean.getResult()));

        mean.clear();
        mean.increment(Double.NEGATIVE_INFINITY);
        mean.increment(1d);
        assertTrue(Double.isNaN(mean.getResult()));

        mean.clear();
        mean.increment(Double.NEGATIVE_INFINITY);
        mean.increment(-1d);
        assertTrue(Double.isNaN(mean.getResult()));

        mean.clear();
        mean.increment(Double.POSITIVE_INFINITY);
        mean.increment(Double.POSITIVE_INFINITY);
        assertTrue(Double.isNaN(mean.getResult()));

        mean.clear();
        mean.increment(Double.NEGATIVE_INFINITY);
        mean.increment(Double.NEGATIVE_INFINITY);
        assertTrue(Double.isNaN(mean.getResult()));

        mean.clear();
        mean.increment(Double.POSITIVE_INFINITY);
        mean.increment(Double.NEGATIVE_INFINITY);
        assertTrue(Double.isNaN(mean.getResult()));

        mean.clear();
        mean.increment(Double.NEGATIVE_INFINITY);
        mean.increment(Double.POSITIVE_INFINITY);
        assertTrue(Double.isNaN(mean.getResult()));

        mean.clear();
        mean.increment(Double.NaN);
        mean.increment(Double.POSITIVE_INFINITY);
        assertTrue(Double.isNaN(mean.getResult()));

        mean.clear();
        mean.increment(Double.NaN);
        mean.increment(Double.NEGATIVE_INFINITY);
        assertTrue(Double.isNaN(mean.getResult()));

        mean.clear();
        mean.increment(Double.NaN);
        mean.increment(0d);
        assertTrue(Double.isNaN(mean.getResult()));
    }
}
