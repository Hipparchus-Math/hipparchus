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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.junit.Assert;
import org.junit.Test;

public class AbstractSQPOptimizerTest {

    @Test
    public void testLagrangianGradX() {
        // GIVEN
        final TestSQPOptimizer testSQPOptimizer = new TestSQPOptimizer();
        final TwiceDifferentiableFunction multivariateFunction = new RosenbrockFunction();
        final ObjectiveFunction objectiveFunction = new ObjectiveFunction(multivariateFunction);
        testSQPOptimizer.parseOptimizationData(objectiveFunction);
        final RealVector expectedVector = MatrixUtils.createRealVector(new double[] { 1 });
        // WHEN
        final RealVector actualVector = testSQPOptimizer.lagrangianGradX(expectedVector, null,null, null);
        // THEN
        for (int i = 0; i < expectedVector.getDimension(); i++) {
            Assert.assertEquals(expectedVector.getEntry(i), actualVector.getEntry(i), 0);
        }
    }

    @Test
    public void testParseOptimizationData() {
        // GIVEN
        final TestSQPOptimizer testSQPOptimizer = new TestSQPOptimizer();
        final SQPOption expectedOptions = new SQPOption();
        final TwiceDifferentiableFunction multivariateFunction = new RosenbrockFunction();
        final ObjectiveFunction objectiveFunction = new ObjectiveFunction(multivariateFunction);
        // WHEN
        testSQPOptimizer.parseOptimizationData(objectiveFunction, expectedOptions);
        // THEN
        Assert.assertEquals(expectedOptions, testSQPOptimizer.getSettings());
    }

    @Test
    public void testParseOptimizationDataException() {
        // GIVEN
        final TestSQPOptimizer testSQPOptimizer = new TestSQPOptimizer();
        final EqualityConstraint equalityConstraint = new TestEqualityConstraints(100000);
        final TwiceDifferentiableFunction multivariateFunction = new RosenbrockFunction();
        final ObjectiveFunction objectiveFunction = new ObjectiveFunction(multivariateFunction);
        // WHEN
        try {
            testSQPOptimizer.parseOptimizationData(objectiveFunction, equalityConstraint);
            Assert.fail();
        } catch (final MathIllegalArgumentException exception) {
            Assert.assertEquals("rank of constraints must be lesser than domain dimension, but 100,000 >= 2",
                    exception.getMessage());
        }
    }

    @Test
    public void testParseOptimizationDataException2() {
        // GIVEN
        final TestSQPOptimizer testSQPOptimizer = new TestSQPOptimizer();
        final EqualityConstraint equalityConstraint = new TestEqualityConstraints(0);
        final TwiceDifferentiableFunction multivariateFunction = new RosenbrockFunction();
        final ObjectiveFunction objectiveFunction = new ObjectiveFunction(multivariateFunction);
        // WHEN
        try {
            testSQPOptimizer.parseOptimizationData(objectiveFunction, equalityConstraint);
            Assert.fail();
        } catch (final MathIllegalArgumentException exception) {
            Assert.assertEquals(LocalizedCoreFormats.ZERO_NOT_ALLOWED.getSourceString(),
                    exception.getMessage());
        }
    }

    private static class TestSQPOptimizer extends AbstractSQPOptimizer {

        TestSQPOptimizer() {
            super();
        }

        @Override
        protected LagrangeSolution doOptimize() {
            return null;
        }
    }

    private static class TestEqualityConstraints extends EqualityConstraint {

        private final int dimension;

        TestEqualityConstraints(final int dimension) {
            super(MatrixUtils.createRealVector(new double[dimension]));
            this.dimension = dimension;
        }

        @Override
        public int dim() {
            return dimension;
        }

        @Override
        public RealVector value(RealVector x) {
            return null;
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            return null;
        }
    }

}
