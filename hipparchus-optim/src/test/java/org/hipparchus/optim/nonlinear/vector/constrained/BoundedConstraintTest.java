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
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class BoundedConstraintTest {

    @Test
    void testConstructor() {
        // GIVEN
        final RealVector lowerBound = MatrixUtils.createRealVector(new double[] { -1. });
        final RealVector upperBound = MatrixUtils.createRealVector(new double[] { 1. });
        // WHEN
        final BoundedConstraint boundedConstraint = new TestBoundedConstraint(lowerBound, upperBound);
        // THEN
        assertEquals(lowerBound.getEntry(0), boundedConstraint.getLowerBound().getEntry(0), 0.);
        assertEquals(upperBound.getEntry(0), boundedConstraint.getUpperBound().getEntry(0), 0.);
    }

    @Test
    void testConstructorNullLowerBound() {
        // GIVEN
        final RealVector upperBound = MatrixUtils.createRealVector(new double[1]);
        // WHEN
        final BoundedConstraint boundedConstraint = new TestBoundedConstraint(null, upperBound);
        // THEN
        assertEquals(boundedConstraint.getLowerBound().getDimension(),
                boundedConstraint.getUpperBound().getDimension());
    }

    @Test
    void testConstructorNullUpperBound() {
        // GIVEN
        final RealVector lowerBound = MatrixUtils.createRealVector(new double[1]);
        // WHEN
        final BoundedConstraint boundedConstraint = new TestBoundedConstraint(lowerBound, null);
        // THEN
        assertEquals(boundedConstraint.getLowerBound().getDimension(),
                boundedConstraint.getUpperBound().getDimension());
    }

    private static class TestBoundedConstraint extends BoundedConstraint {

        TestBoundedConstraint(RealVector lower, RealVector upper) {
            super(lower, upper);
        }

        @Override
        public int dim() {
            return 0;
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
