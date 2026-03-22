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

import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;

public class HS345Test {

    private static final int DIM = 3;
    private static final double CONST_TERM = 4.0 + Math.sqrt(18.0); // 4 + sqrt(18)

    static final class HS345Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }

        // F(X) = (X1-1)^2 + (X1-X2)^2 + (X2-X3)^4
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double d1 = x1 - 1.0;
            double d2 = x1 - x2;
            double d3 = x2 - x3;

            return d1 * d1 + d2 * d2 + Math.pow(d3, 4);
        }

        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double d3_3 = Math.pow(x2 - x3, 3);

            // GF(1) = 2*(X1-1) + 2*(X1-X2)
            double g1 = 2.0 * (x1 - 1.0) + 2.0 * (x1 - x2);

            // GF(2) = -2*(X1-X2) + 4*(X2-X3)^3
            double g2 = -2.0 * (x1 - x2) + 4.0 * d3_3;

            // GF(3) = -4*(X2-X3)^3
            double g3 = -4.0 * d3_3;

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            // Hessian is identical to HS344's objective function
            double d3_2 = Math.pow(x.getEntry(1) - x.getEntry(2), 2);

            double[][] H = new double[DIM][DIM];

            H[0][0] = 4.0;
            H[0][1] = -2.0;
            H[1][0] = -2.0;
            H[1][1] = 2.0 + 12.0 * d3_2;
            H[1][2] = -12.0 * d3_2;
            H[2][1] = -12.0 * d3_2;
            H[2][2] = 12.0 * d3_2;

            return MatrixUtils.createRealMatrix(H);
        }
    }

    static final class HS345Ineq extends InequalityConstraint {

        HS345Ineq() { super(new ArrayRealVector(new double[1])); }

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            // The Fortran constraint is: X1*(1 + X2^2) + X3^4 - 4 - sqrt(18) <= 0
            // We convert to G >= 0: 4 + sqrt(18) - (X1*(1 + X2^2) + X3^4) >= 0
            double g1 = CONST_TERM - (x1 * (1.0 + x2 * x2) + Math.pow(x3, 4));

            return new ArrayRealVector(new double[]{g1}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double x2_2 = x2 * x2;
            double x3_3 = Math.pow(x3, 3);

            double[][] J = new double[1][DIM];

            // dG1/dX1 = -(1 + X2^2)
            J[0][0] = -(1.0 + x2_2);

            // dG1/dX2 = -2*X1*X2
            J[0][1] = -2.0 * x1 * x2;

            // dG1/dX3 = -4*X3^3
            J[0][2] = -4.0 * x3_3;

            return MatrixUtils.createRealMatrix(J);
        }
    }

    private static double[] start() {
        return new double[]{0.0, 0.0, 0.0};
    }

    @Test
    public void testHS345() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();

        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }



        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS345Obj()),
                new HS345Ineq()

        );

        double f = sol.getValue();
        final double fExpected = 0.032568200;
        final double tolerance = 1.0e-5 * (Math.abs(fExpected) + 1.0);

    // Assert: The found value 'f' must be close to OR better (less than) the reference 'fExpected'.
    // We check if f <= fExpected + tolerance
      assertTrue(f <= fExpected + tolerance, "The found objective value (" + f +
        ") is significantly worse than the expected value (" + fExpected + ")");


    }
}
