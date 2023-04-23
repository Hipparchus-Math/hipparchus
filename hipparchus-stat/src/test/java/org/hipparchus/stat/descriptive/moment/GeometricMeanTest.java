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
import org.junit.Assert;
import org.junit.Test;

/**
 * Test cases for the {@link GeometricMean} class.
 */
public class GeometricMeanTest extends StorelessUnivariateStatisticAbstractTest{

    @Override
    public GeometricMean getUnivariateStatistic() {
        return new GeometricMean();
    }

    @Override
    public double expectedValue() {
        return this.geoMean;
    }

    @Test
    public void testSpecialValues() {
        GeometricMean mean = getUnivariateStatistic();
        // empty
        Assert.assertTrue(Double.isNaN(mean.getResult()));

        // finite data
        mean.increment(1d);
        Assert.assertFalse(Double.isNaN(mean.getResult()));

        // add 0 -- makes log sum blow to minus infinity, should make 0
        mean.increment(0d);
        Assert.assertEquals(0d, mean.getResult(), 0);

        // add positive infinity - note the minus infinity above
        mean.increment(Double.POSITIVE_INFINITY);
        Assert.assertTrue(Double.isNaN(mean.getResult()));

        // clear
        mean.clear();
        Assert.assertTrue(Double.isNaN(mean.getResult()));

        // positive infinity by itself
        mean.increment(Double.POSITIVE_INFINITY);
        Assert.assertEquals(Double.POSITIVE_INFINITY, mean.getResult(), 0);

        // negative value -- should make NaN
        mean.increment(-2d);
        Assert.assertTrue(Double.isNaN(mean.getResult()));
    }

}
