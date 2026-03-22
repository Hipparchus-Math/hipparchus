/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
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

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS340Test {

    private static final int DIM = 3;

    static final class HS340Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }

        // F(X) = -X1*X2*X3
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            return -x1 * x2 * x3;
        }

        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            // GF(1) = -X2*X3
            // GF(2) = -X1*X3
            // GF(3) = -X1*X2
            return new ArrayRealVector(new double[]{-x2 * x3, -x1 * x3, -x1 * x2}, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double[][] H = new double[DIM][DIM];
            // H12 = -X3
            H[0][1] = -x3;
            H[1][0] = -x3;

            // H13 = -X2
            H[0][2] = -x2;
            H[2][0] = -x2;

            // H23 = -X1
            H[1][2] = -x1;
            H[2][1] = -x1;

            return MatrixUtils.createRealMatrix(H);
        }
    }

    static final class HS340Ineq extends InequalityConstraint {

        HS340Ineq() { super(new ArrayRealVector(new double[1])); }

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            // G(1) = 1.8 - X1 - 2*X2 - 2*X3 >= 0
            double g1 = 1.8 - x1 - 2.0 * x2 - 2.0 * x3;

            return new ArrayRealVector(new double[]{g1}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double[][] J = new double[1][DIM];

            // dG1/dX: -1.0, -2.0, -2.0
            J[0][0] = -1.0;
            J[0][1] = -2.0;
            J[0][2] = -2.0;

            return MatrixUtils.createRealMatrix(J);
        }
    }

    private static double[] start() {
        return new double[]{1.0, 1.0, 1.0};
    }

    @Test
    public void testHS340() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();

        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // Box constraints: X1 <= 1.0, X1, X2, X3 >= 0.0
        SimpleBounds bounds = new SimpleBounds(
            new double[]{0.0, 0.0, 0.0},
            new double[]{1.0, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS340Obj()),
                new HS340Ineq(),
                bounds
        );

        double f = sol.getValue();
        final double fExpected = -0.054;

        assertEquals(fExpected, f, 1.0e-5 * (Math.abs(fExpected) + 1.0), "objective mismatch");

    }
}
