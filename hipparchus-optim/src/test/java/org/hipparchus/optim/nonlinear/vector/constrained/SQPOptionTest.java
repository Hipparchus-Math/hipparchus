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


import org.junit.Assert;
import org.junit.Test;

public class SQPOptionTest {

    @Test
    public void testSettersDouble() {
        // GIVEN
        final SQPOption option = new SQPOption();
        // WHEN
        final double expectedValue = 1.;
        option.setB(expectedValue);
        option.setEps(expectedValue);
        option.setMu(expectedValue);
        option.setRhoCons(expectedValue);
        option.setSigmaMax(expectedValue);
        // THEN
        final double tolerance = 0.;
        Assert.assertEquals(expectedValue, option.getB(), tolerance);
        Assert.assertEquals(expectedValue, option.getEps(), tolerance);
        Assert.assertEquals(expectedValue, option.getMu(), tolerance);
        Assert.assertEquals(expectedValue, option.getRhoCons(), tolerance);
        Assert.assertEquals(expectedValue, option.getSigmaMax(), tolerance);
    }

    @Test
    public void testSettersBoolean() {
        // GIVEN
        final SQPOption option = new SQPOption();
        // WHEN
        final boolean expectedValue = true;
        option.setUseFunHessian(expectedValue);
        // THEN
        Assert.assertEquals(expectedValue, option.getUseFunHessian());
    }

    @Test
    public void testSettersIntegers() {
        // GIVEN
        final SQPOption option = new SQPOption();
        // WHEN
        final int expectedValue = 10;
        option.setConvCriteria(expectedValue);
        option.setMaxLineSearchIteration(expectedValue);
        option.setQpMaxLoop(expectedValue);
        // THEN
        Assert.assertEquals(expectedValue, option.getConvCriteria());
        Assert.assertEquals(expectedValue, option.getMaxLineSearchIteration());
        Assert.assertEquals(expectedValue, option.getQpMaxLoop());
    }

}
