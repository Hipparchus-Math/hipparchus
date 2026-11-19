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
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.OptimizationData;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public abstract class AbstractConstrainedOptimizerTest {

    /** Build the optimizer.
     * @return built optimizer
     */
    protected abstract ConstraintOptimizer buildOptimizer();

    /** Test one problem.
     * @param expectedSolution expected solution
     * @param solutionTolerance tolerance on solution (L₁ norm)
     * @param expectedMultipliers expected multipliers
     * @param multipliersTolerance tolerance on multipliers  (L₁ norm)
     * @param expectedValue expected objective function value
     * @param valueTolerance tolerance on objective function value
     * @param objectiveFunction objective function
     * @param initialGuess initial guess (may be null)
     * @param constraints contraints
     */
    protected void doTestProblem(final double[] expectedSolution,
                                 final double solutionTolerance,
                                 final double[] expectedMultipliers,
                                 final double multipliersTolerance,
                                 final double expectedValue,
                                 final double valueTolerance,
                                 final ObjectiveFunction objectiveFunction,
                                 final double[] initialGuess,
                                 final Constraint... constraints) {

        // find optimum solution
        final ConstraintOptimizer optimizer = buildOptimizer();
        final OptimizationData[] data = new OptimizationData[constraints.length + (initialGuess == null ? 2 : 3)];
        data[0] = objectiveFunction;
        data[1] = new MatrixDecompositionTolerance(2.0e-11);
        System.arraycopy(constraints, 0, data, 2, constraints.length);
        if (initialGuess != null) {
            data[data.length - 1] = new InitialGuess(initialGuess);
        }
        final LagrangeSolution    solution  = optimizer.optimize(data);

        // check result
        assertEquals(0.0,
                     MatrixUtils.createRealVector(expectedSolution).subtract(solution.getX()).getL1Norm(),
                     solutionTolerance);
        assertEquals(0.0,
                     MatrixUtils.createRealVector(expectedMultipliers).subtract(solution.getLambda()).getL1Norm(),
                     multipliersTolerance);
        assertEquals(expectedValue, solution.getValue(), valueTolerance);

        // check neighboring points either violate constraints or have worst objective function
        for (int i = 0; i < expectedSolution.length; ++i) {

            final RealVector plusShift = MatrixUtils.createRealVector(expectedSolution);
            plusShift.addToEntry(i, 2 * solutionTolerance);
            boolean plusIsFeasible = true;
            for (final Constraint constraint : constraints) {
                plusIsFeasible &= constraint.overshoot(constraint.value(plusShift)) <= 0;
            }
            if (plusIsFeasible) {
                // the plusShift point fulfills all constraints,
                // so it must have worst objective function than the expected optimum
                assertTrue(objectiveFunction.getObjectiveFunction().value(plusShift.toArray()) > expectedValue);
            }
            
            final RealVector minusShift = MatrixUtils.createRealVector(expectedSolution);
            minusShift.addToEntry(i, -2 * solutionTolerance);
            boolean minusIsFeasible = true;
            for (final Constraint constraint : constraints) {
                minusIsFeasible &= constraint.overshoot(constraint.value(minusShift)) <= 0;
            }
            if (minusIsFeasible) {
                // the minusShift point fulfills all constraints,
                // so it must have worst objective function than the expected optimum
                assertTrue(objectiveFunction.getObjectiveFunction().value(minusShift.toArray()) > expectedValue);
            }
            
        }

    }

}
