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

import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.junit.Assert;
import org.junit.Test;

public class TwiceDifferentiableFunctionTest {

    @Test
    public void testValue() {
        // GIVEN
        final double[] x = new double[] { 2. };
        final TestTwiceDifferentiableFunction function = new TestTwiceDifferentiableFunction();
        // WHEN
        final double actualValue = function.value(x);
        // THEN
        final double expectedValueElement = x[0];
        Assert.assertEquals(expectedValueElement, actualValue, 0);
    }

    @Test
    public void testGradient() {
        // GIVEN
        final double[] x = new double[] { 2. };
        final TestTwiceDifferentiableFunction function = new TestTwiceDifferentiableFunction();
        // WHEN
        final RealVector actualGradient = function.gradient(x);
        // THEN
        final double expectedValueElement = 1.;
        Assert.assertEquals(expectedValueElement, actualGradient.getEntry(0), 0);
    }

    @Test
    public void testHessian() {
        // GIVEN
        final double[] x = new double[] { 2. };
        final TestTwiceDifferentiableFunction function = new TestTwiceDifferentiableFunction();
        // WHEN
        final RealMatrix actualGradient = function.hessian(x);
        // THEN
        final double expectedValueElement = 0.;
        Assert.assertEquals(expectedValueElement, actualGradient.getEntry(0, 0), 0);
    }

    private static class TestTwiceDifferentiableFunction extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return 1;
        }

        @Override
        public double value(RealVector x) {
            return x.getEntry(0);
        }

        @Override
        public RealVector gradient(RealVector x) {
            return MatrixUtils.createRealVector(new double[] { 1. });
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            return MatrixUtils.createRealMatrix(1, 1);
        }

    }

}
