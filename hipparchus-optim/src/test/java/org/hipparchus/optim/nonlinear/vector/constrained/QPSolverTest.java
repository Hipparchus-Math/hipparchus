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
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;


public class QPSolverTest {

    @Test
    public void testGoldfarbIdnaniExampleNoInverseL() {
        doTestGoldfarbIdnaniExample(false);
    }

    @Test
    public void testGoldfarbIdnaniExampleInverseL() {
        doTestGoldfarbIdnaniExample(true);
    }

    private void doTestGoldfarbIdnaniExample(final boolean inverseL) {

        QuadraticFunction q = new QuadraticFunction(new double[][] { { 4.0, -2.0 }, { -2.0, 4.0 } },
                                                    new double[] { 6.0, 0.0 },
                                                    0.0);

        // y = 1
        LinearEqualityConstraint eqc = new LinearEqualityConstraint(new double[][] { { 1.0, 1.0 } },
                                                                    new double[] { 3.0 });


        // x > 0, y > 0, x + y > 2
        LinearInequalityConstraint ineqc = new LinearInequalityConstraint(new double[][] { { 1.0, 0.0 }, { 0.0, 1.0 }, { 1.0, 1.0 } },
                                                                          new double[] { 0.0, 0.0, 2.0 });
       

       

       

        
        QPDualActiveSolver solver = new QPDualActiveSolver();

        final double s3 = FastMath.sqrt(3.0);
        LagrangeSolution solution = inverseL ?
                                    solver.optimize(new ObjectiveFunction(q), eqc, ineqc,
                                                    new InverseCholesky(MatrixUtils.createRealMatrix(new double[][] {
                                                            { 0.5, 0.0 }, { 0.5 / s3 , 1.0 / s3 }
                                                    })),
                                                    new MatrixDecompositionTolerance(1.0e-13)) :
                                    solver.optimize(new ObjectiveFunction(q), eqc, ineqc);

        RealVector x = solution.getX();

        double expectedObj = solution.getValue();

        assertArrayEquals(new double[]{ 1.0, 2.0 }, x.toArray(), 1e-8);
        assertEquals(12.0, expectedObj, 1e-8);
    }

    @Test
    public void testGoldfarbIdnaniExampleInverseCholeskyMode() {
        // same problem as above, but the quadratic term is given directly as the
        // inverse of the Cholesky factor of [[4, -2], [-2, 4]]
        final double s3 = FastMath.sqrt(3.0);
        QuadraticFunction q = new QuadraticFunction(new double[][] { { 0.5, 0.0 }, { 0.5 / s3, 1.0 / s3 } },
                                                    new double[] { 6.0, 0.0 },
                                                    0.0);
        LinearEqualityConstraint eqc = new LinearEqualityConstraint(new double[][] { { 1.0, 1.0 } },
                                                                    new double[] { 3.0 });
        LinearInequalityConstraint ineqc = new LinearInequalityConstraint(new double[][] { { 1.0, 0.0 }, { 0.0, 1.0 }, { 1.0, 1.0 } },
                                                                          new double[] { 0.0, 0.0, 2.0 });
        QPDualActiveSolverOption option = new QPDualActiveSolverOption();
        option.setMatrixMode(QPMatrixMode.INVCHOLESKY);

        LagrangeSolution solution = new QPDualActiveSolver().optimize(option, new ObjectiveFunction(q), eqc, ineqc);

        assertArrayEquals(new double[]{ 1.0, 2.0 }, solution.getX().toArray(), 1e-8);
        assertEquals(12.0, solution.getValue(), 1e-8);
    }

    @Test
    public void testNonSymmetricHessian() {
        QuadraticFunction q = new QuadraticFunction(new double[][] { { 4.0, -2.0 }, { 3.0, 4.0 } },
                                                    new double[] { 6.0, 0.0 },
                                                    0.0);
        LinearInequalityConstraint ineqc = new LinearInequalityConstraint(new double[][] { { 1.0, 0.0 } },
                                                                          new double[] { 0.0 });

        LagrangeSolution solution = new QPDualActiveSolver().optimize(new ObjectiveFunction(q), ineqc);

        assertEquals(0, solution.getX().getDimension());
        assertEquals(QPDualActiveSolver.ERROR_CHOLESKY_DECOMPOSITION, solution.getLambda().getEntry(0), 0.0);
    }
}
