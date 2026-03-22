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
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS330Test {

    private static final int DIM = 2;
    private static final double MIN_BOUND = 0.0001;

    static final class HS330Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // F(X) = 0.044*X1^3/X2^2 + 1/X1 + 0.0592*X1/X2^3
            return 0.044 * Math.pow(x1, 3) / Math.pow(x2, 2) + 1.0 / x1 + 0.0592 * x1 / Math.pow(x2, 3);
        }
        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // dF/dX1 = 0.132*X1^2/X2^2 - X1^(-2) + 0.0592/X2^3
            double g1 = 0.132 * x1 * x1 / (x2 * x2) - 1.0 / (x1 * x1) + 0.0592 / (x2 * x2 * x2);

            // dF/dX2 = -0.088*X1^3/X2^3 - 0.1776*X1/X2^4
            double g2 = -0.088 * Math.pow(x1, 3) / Math.pow(x2, 3) - 0.1776 * x1 / Math.pow(x2, 4);

            return new ArrayRealVector(new double[]{g1, g2}, false);
        }
        @Override public RealMatrix hessian(RealVector x) {
            // Hessian not provided in the original subroutine
            throw new UnsupportedOperationException("Hessian matrix is not implemented for this test case.");
        }
    }

    static final class HS330Ineq extends InequalityConstraint {

        HS330Ineq() { super(new ArrayRealVector(new double[1])); }

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // G(1) = 1 - 8.62*X2^3/X1 >= 0
            double g1 = 1.0 - 8.62 * Math.pow(x2, 3) / x1;

            return new ArrayRealVector(new double[]{g1}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            double[][] J = new double[1][DIM];

            // dG1/dX1 = 8.62*X2^3/X1^2
            J[0][0] = 8.62 * Math.pow(x2, 3) / (x1 * x1);

            // dG1/dX2 = -25.86*X2^2/X1
            J[0][1] = -25.86 * (x2 * x2) / x1;

            return MatrixUtils.createRealMatrix(J);
        }
    }

    private static double[] start() {
        return new double[]{2.5, 2.5};
    }

    @Test
    public void testHS330() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();

        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // Box constraints: 0.0001 <= X1, X2 <= 5.0
        SimpleBounds bounds = new SimpleBounds(
            new double[]{0.0001, Double.NEGATIVE_INFINITY},
            new double[]{Double.POSITIVE_INFINITY, 5.0}
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS330Obj()),
                new HS330Ineq(),
                bounds
        );

        double f = sol.getValue();
        final double fExpected = 1.6205833;

        assertEquals(fExpected, f, 1.0e-5 * (Math.abs(fExpected) + 1.0), "objective mismatch");


    }
}
