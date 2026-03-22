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

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS334Test {

    private static final int DIM = 3;
    private static final int DATA_SIZE = 15;

    // Dependent variable data Y (15 points)
    private static final double[] Y_DATA = {
        0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 0.58, 0.73, 0.96, 1.34, 2.1, 4.39
    };

    static final class HS334Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }

        // F(X) = sum( ( Y(I) - (X1 + I / (X2*VI + X3*WI)) )^2 )
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double fx = 0.0;

            for (int i = 0; i < DATA_SIZE; i++) {
                double ui = i + 1;             // I
                double vi = 16.0 - ui;         // VI
                double wi = Math.min(ui, vi);  // WI

                double denominator = x2 * vi + x3 * wi;
                // Residual F(I)
                double fi = Y_DATA[i] - (x1 + ui / denominator);
                fx += fi * fi;
            }
            return fx;
        }

        @Override public RealVector gradient(RealVector x) {
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double g1 = 0.0;
            double g2 = 0.0;
            double g3 = 0.0;

            for (int i = 0; i < DATA_SIZE; i++) {
                double ui = i + 1;
                double vi = 16.0 - ui;
                double wi = Math.min(ui, vi);

                double denominator = x2 * vi + x3 * wi;
                double denominator_sq = denominator * denominator;

                // Residual F(I)
                double fi = Y_DATA[i] - (x.getEntry(0) + ui / denominator);

                // DF(I,1) = -1.0
                double dfi_dx1 = -1.0;

                // DF(I,2) = (I * VI) / (X2*VI + X3*WI)^2
                // d/dX2 (-(X1 + I/D)) = - ( -I/D^2 * dD/dX2 ) = I*VI / D^2
                double dfi_dx2 = (ui * vi) / denominator_sq;

                // DF(I,3) = (I * WI) / (X2*VI + X3*WI)^2
                double dfi_dx3 = (ui * wi) / denominator_sq;

                // Gradient GF(j) = 2 * sum(F(I) * DF(I, j))
                g1 += dfi_dx1 * fi * 2.0;
                g2 += dfi_dx2 * fi * 2.0;
                g3 += dfi_dx3 * fi * 2.0;
            }
            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            // Hessian matrix is not implemented for this test case.
            throw new UnsupportedOperationException("Hessian matrix is not implemented for this test case.");
        }
    }

    private static double[] start() {
        return new double[]{1.0, 1.0, 1.0};
    }

    @Test
    public void testHS334() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();

        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }


        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS334Obj())

        );

        double f = sol.getValue();
        final double fExpected = 0.0082148773;

        assertEquals(fExpected, f, 1.0e-6 * (Math.abs(fExpected) + 1.0), "objective mismatch");


    }
}
