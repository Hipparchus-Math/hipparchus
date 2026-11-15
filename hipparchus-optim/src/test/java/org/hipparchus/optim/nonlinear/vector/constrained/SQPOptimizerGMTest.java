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
import org.hipparchus.optim.InitialGuess;
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

    // The tests below only add coverage to this deprecated class, they will be removed with it in 5.0.

    @Test
    void testFunctionHessian() {
        final SQPOptimizerGM optimizer = new SQPOptimizerGM();
        optimizer.getSettings().setUseFunHessian(true);
        final LagrangeSolution solution = optimizer.optimize(createOptimizationData());
        // observed errors: 1.9e-6 on the solution, 3.3e-6 on the value
        assertEquals(0.0,
                     MatrixUtils.createRealVector(new double[] { 1, 1 }).subtract(solution.getX()).getL1Norm(),
                     5.0e-6);
        assertEquals(8., solution.getValue(), 5.0e-6);
    }

    @Test
    void testConvergenceCriterion0() {
        final SQPOptimizerGM optimizer = new SQPOptimizerGM();
        optimizer.getSettings().setConvCriteria(0);
        final LagrangeSolution solution = optimizer.optimize(createOptimizationData());
        // observed errors: 1.0e-10 on the solution, 9.7e-10 on the value
        assertEquals(0.0,
                     MatrixUtils.createRealVector(new double[] { 1, 1 }).subtract(solution.getX()).getL1Norm(),
                     1.0e-9);
        assertEquals(8., solution.getValue(), 2.0e-9);
    }

    private OptimizationData[] createOptimizationData() {
        final QuadraticFunction q = new QuadraticFunction(new double[][] { { 4.0, -2.0 }, { -2.0, 4.0 } },
                new double[] { 6.0, 0.0 },
                0.0);
        // y = 1
        final LinearEqualityConstraint eqc = new LinearEqualityConstraint(new double[][] { { 0.0, 1.0 } },
                new double[] { 1.0 });
        // x > 0, y > 0, x + y > 2
        final LinearInequalityConstraint ineqc = new LinearInequalityConstraint(new double[][] { { 1.0, 0.0 }, { 0.0, 1.0 }, { 1.0, 1.0 } },
                new double[] { 0.0, 0.0, 2.0 });
        return new OptimizationData[] {
            new ObjectiveFunction(q), eqc, ineqc, new InitialGuess(new double[] { 3.5, 3.5 })
        };
    }

}
