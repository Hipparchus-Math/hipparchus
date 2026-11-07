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

import org.hipparchus.optim.OptimizationData;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class SQPOptimizerGMTest extends AbstractTestAbstractSQPOptimizerTest {

    protected ConstraintOptimizer buildOptimizer() {
        return new SQPOptimizerGM();
    }

    @Test
    void testWithEqualityConstraintsOnly() {
        final QuadraticFunction q = new QuadraticFunction(new double[][] { { 1.0, 0.0 }, { 0.0, 1.0 } },
                new double[] { 1.0, 0.0 },
                0.0);

        final LinearEqualityConstraint eqc = new LinearEqualityConstraint(new double[][] { { 1.0, 0.0 } },
                new double[] { 1.0 });

        final ConstraintOptimizer optimizer = buildOptimizer();
        final OptimizationData[] data = new OptimizationData[] { new ObjectiveFunction(q), eqc };
        final LagrangeSolution    solution  = optimizer.optimize(data);

        assertEquals(1.5, solution.getValue(), 1.e-4);
    }

    @Test
    void testWithInequalityConstraintsOnly() {
        final QuadraticFunction q = new QuadraticFunction(new double[][] { { 1.0, 0.0 }, { 0.0, 1.0 } },
                new double[] { 1.0, 0.0 },
                0.0);

        final LinearInequalityConstraint eqc = new LinearInequalityConstraint(new double[][] { { 1.0, 0.0 } },
                new double[] { 1.0 });

        final ConstraintOptimizer optimizer = buildOptimizer();
        final OptimizationData[] data = new OptimizationData[] { new ObjectiveFunction(q), eqc };
        final LagrangeSolution    solution  = optimizer.optimize(data);

        assertEquals(1.5, solution.getValue(), 1.e-4);
    }

}
