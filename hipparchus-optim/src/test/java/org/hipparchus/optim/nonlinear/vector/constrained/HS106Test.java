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

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

/** HS TP106 (Schittkowski): 8 variabili, 3 diseguaglianze lineari + 3 nonlineari. */
public class HS106Test {

    // Bounds finali dal Fortran:
    // x1 in [100, 1e4], x2 in [1e3, 1e4], x3 in [1e3, 1e4],
    // x4..x8 in [10, 1e3]
    private static final double[] LB = {
        100.0, 1000.0, 1000.0, 10.0, 10.0, 10.0, 10.0, 10.0
    };
    private static final double[] UB = {
        1.0e4, 1.0e4, 1.0e4, 1.0e3, 1.0e3, 1.0e3, 1.0e3, 1.0e3
    };

    /** f(x) = x1 + x2 + x3. */
    private static class TP106Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 8; }
        @Override public double value(RealVector X) {
            return X.getEntry(0) + X.getEntry(1) + X.getEntry(2);
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    /** 6 diseguaglianze in forma g(x) >= 0 (MODE=4 del Fortran). */
    private static class TP106Ineq extends InequalityConstraint {
        TP106Ineq() { super(new ArrayRealVector(new double[]{0,0,0,0,0,0})); }
        @Override public int dim() { return 8; }

        @Override public RealVector value(RealVector X) {
            final double x1 = X.getEntry(0);
            final double x2 = X.getEntry(1);
            final double x3 = X.getEntry(2);
            final double x4 = X.getEntry(3);
            final double x5 = X.getEntry(4);
            final double x6 = X.getEntry(5);
            final double x7 = X.getEntry(6);
            final double x8 = X.getEntry(7);

            // Lineari:
            final double g1 = -2.5e-3*(x4 + x6) + 1.0;
            final double g2 = -2.5e-3*(x5 + x7 - x4) + 1.0;
            final double g3 = -0.01*(x8 - x5) + 1.0;

            // Non lineari:
            final double g4 = -833.33252*x4 - 100.0*x1 + 8.3333333e4 + x1*x6;
            final double g5 = -1.25e3*x5 - x2*x4 + 1.25e3*x4 + x2*x7;
            final double g6 = -1.25e6 - x3*x5 + 2.5e3*x5 + x3*x8;

            return new ArrayRealVector(new double[]{ g1, g2, g3, g4, g5, g6 });
        }

        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testHS106() {
        // Guess dal Fortran
        final double[] x0 = {
            5000.0, 5000.0, 5000.0, 200.0, 350.0, 150.0, 225.0, 425.0
        };

        final InitialGuess guess = new InitialGuess(x0);
        final SimpleBounds bounds = new SimpleBounds(LB, UB);

        final SQPOptimizerS2 opt = new SQPOptimizerS2();
        opt.setDebugPrinter(System.out::println);

        final LagrangeSolution sol = opt.optimize(
            guess,
            new ObjectiveFunction(new TP106Obj()),
            new TP106Ineq(),
            bounds
        );

        // FEX (Fortran): 0.70492480D+04 = 7049.248
        final double expected = 7049.248;
        assertEquals(expected, sol.getValue(), 1e-3);
    }
}
