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
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

@Deprecated
class SQPOptimizerSTest extends AbstractTestAbstractSQPOptimizerTest {

    protected ConstraintOptimizer buildOptimizer() {
        return new SQPOptimizerS();
    }

    @Test
    void test2() {
        QuadraticFunction q = new QuadraticFunction(new double[][] { { 6.0, 2.0 }, { 2.0, 8.0 } },
                                                    new double[] { 5.0, 1.0 },
                                                    0.0);

        // constraint x + 2y = 1
        LinearEqualityConstraint eqc = new LinearEqualityConstraint(new double[][] { { 1.0, 2.0 } },
                                                                    new double[] { 1.0 });

        // x > 0, y > 0
        LinearInequalityConstraint ineqc = new LinearInequalityConstraint(new double[][] { { 1.0, 0.0 }, { 0.0, 1.0 } },
                                                                          new double[] { 0.0, 0.0 });


        doTestProblem(new double[] {  0, 0.5 },     1.9e-15,
                      new double[] { 2.5, 3.5, 0 }, 7.8e-5,
                      1.5, 4.0e-15,
                      new ObjectiveFunction(q),
                      new double[] { 10.5, 10.5 },
                      eqc, ineqc);

    }

    @Test
    void testHockShittkowski71() {
        doTestProblem(new double[] { 1, 4.74293167, 3.82123882, 1.37939596 }, 2.2e-4,
                      new double[] { -0.16145839, 0.55229016, 1.08782965, 0, 0, 0, 0, 0, 0, 0 }, 5.6e-5,
                      17.01401698, 1.4e-7,
                      new ObjectiveFunction(new HockSchittkowskiFunction71()),
                      new double[] { 1, 5, 5, 1 },
                      new HockSchittkowskiConstraintInequality71(),
                      new HockSchittkowskiConstraintEquality71());
    }

    @Test
    @Disabled // this test fails completely
    void testHockShittkowski72() {
        doTestProblem(new double[] { 193.12529425, 180.14766487, 184.58883790, 168.82104861 }, 1.5,
                      new double[] { 7693.73706410, 41453.54250351 }, 1.1e-8,
                      727.68284564, 1.0e-8,
                      new ObjectiveFunction(new HockSchittkowskiFunction72()),
                      new double[] { 1, 5, 5, 1 },
                      new HockSchittkowskiConstraintInequality72());
    }

    @Test
    @Disabled // this test fails completely
    void testHockShittkowski77() {
        doTestProblem(new double[] { 1.16617194, 1.18211086, 1.38025671, 1.50603641, 0.61092012 }, 1.4e-8,
                      new double[] { 0.08553981, 0.03187858 }, 1.0e-8,
                      0.24150486, 1.0e-8,
                      new ObjectiveFunction(new HockSchittkowskiFunction77()),
                      new double[] { 2, 2, 2, 3, 3 },
                      new HockSchittkowskiConstraintEquality77());
    }

    @Test
    void testHockShittkowski78() {
        doTestProblem(new double[] { -1.71714365, 1.59570987, 1.82724583, 0.76364341, 0.76364341 }, 2.2e-7,
                      new double[] { -0.74445225, 0.70358075, -0.09680628 }, 9.0e-6,
                      -2.91970350, 1.0e-8,
                      new ObjectiveFunction(new HockSchittkowskiFunction78()),
                      new double[] { -2, 1.5, 2, 1, 1 },
                      new HockSchittkowskiConstraintEquality78());
    }

    @Test
    void testRosenbrock() {
        doTestProblem(new double[] { 1, 1 }, 9.6e-7,
                      new double[] { 0, 0, 0, 0, 0}, 1.5e-10,
                      0.0, 1.6e-13,
                      new ObjectiveFunction(new RosenbrockFunction()),
                      new double[] { 2, 2 },
                      new RosenbrookConstraint(new double[]{ -2, -1.5, -1.5, -1.5, -1.5 }));
    }

    @Test
    void testLowMaxLineSearchAndConvergenceCriterion0() {
        // GIVEN
        final OptimizationData[] data = createOptimizationData();
        final SQPOption option = new SQPOption();
        option.setMaxLineSearchIteration(2);
        option.setConvCriteria(0);
        data[data.length - 1] = option;

        // WHEN
        final SQPOptimizerS optimizer = new SQPOptimizerS();
        optimizer.parseOptimizationData(data);
        final LagrangeSolution    solution  = optimizer.optimize(data);

        // THEN
        final double[] expectedSolution = new double[] { 1, 1 };
        assertEquals(0.0,
                MatrixUtils.createRealVector(expectedSolution).subtract(solution.getX()).getL1Norm(), 2.5e-5);
        assertEquals(8., solution.getValue(), 2e-4);
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

        final OptimizationData[] constraints = new OptimizationData[] { eqc, ineqc };
        final OptimizationData[] data = new OptimizationData[constraints.length + 3];
        final ObjectiveFunction objectiveFunction =  new ObjectiveFunction(q);
        data[0] = objectiveFunction;
        System.arraycopy(constraints, 0, data, 1, constraints.length);
        final double[] initialGuess = new double[] { -3.5, 3.5 };
        data[data.length - 2] = new InitialGuess(initialGuess);
        return data;
    }

}
