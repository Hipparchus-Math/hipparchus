/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements. See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 *
 * https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS343Test {

    private static final int DIM = 3;
    private static final double C_FACTOR = 0.0201e-6; // .201D-1 * .1D-6

    static final class HS343Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }

        // F(X) = -0.0201e-6 * X1^4 * X2 * X3^2
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            return -C_FACTOR * Math.pow(x1, 4) * x2 * (x3 * x3);
        }

        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double x1_4 = Math.pow(x1, 4);
            double x1_3 = x1_4 / x1;
            double x3_2 = x3 * x3;

            // GF(1) = -(4 * 0.0201e-6) * X1^3 * X2 * X3^2
            double g1 = -4.0 * C_FACTOR * x1_3 * x2 * x3_2;

            // GF(2) = -0.0201e-6 * X1^4 * X3^2
            double g2 = -C_FACTOR * x1_4 * x3_2;

            // GF(3) = -(2 * 0.0201e-6) * X1^4 * X2 * X3
            double g3 = -2.0 * C_FACTOR * x1_4 * x2 * x3;

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            // Hessian matrix is highly non-linear and not needed for a standard SQP test.
            throw new UnsupportedOperationException("Hessian matrix is not implemented for this test case.");
        }
    }

    static final class HS343Ineq extends InequalityConstraint {

        HS343Ineq() { super(new ArrayRealVector(new double[2])); }

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            // G(1) = 675 - X1^2 * X2 >= 0
            double g1 = 675.0 - (x1 * x1 * x2);

            // G(2) = 0.419 - 1e-6 * X1^2 * X3^2 >= 0
            double g2 = 0.419 - 1.0e-6 * (x1 * x1 * x3 * x3);

            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double x1_2 = x1 * x1;
            double x3_2 = x3 * x3;

            double[][] J = new double[2][DIM];

            // G1: dG1/dX1 = -2*X1*X2, dG1/dX2 = -X1^2, dG1/dX3 = 0
            J[0][0] = -2.0 * x1 * x2;
            J[0][1] = -x1_2;
            J[0][2] = 0.0;

            // G2: dG2/dX1 = -2e-6*X1*X3^2, dG2/dX2 = 0, dG2/dX3 = -2e-6*X1^2*X3
            J[1][0] = -2.0e-6 * x1 * x3_2;
            J[1][1] = 0.0;
            J[1][2] = -2.0e-6 * x1_2 * x3;

            return MatrixUtils.createRealMatrix(J);
        }
    }

    private static double[] start() {
        return new double[]{22.3, 0.5, 125.0};
    }

    @Test
    public void testHS343() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();

        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // Box constraints: 0 <= X1 <= 36, 0 <= X2 <= 5, 0 <= X3 <= 125
        SimpleBounds bounds = new SimpleBounds(
            new double[]{0.0, 0.0, 0.0},
            new double[]{36.0, 5.0, 125.0}
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS343Obj()),
                new HS343Ineq(),
                bounds
        );

        double f = sol.getValue();
        final double fExpected = -5.6847825;

        assertEquals(fExpected, f, 1.0e-5 * (Math.abs(fExpected) + 1.0), "objective mismatch");


    }
}
