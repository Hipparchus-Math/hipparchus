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

/** HS TP116. 13 variabili, 15 vincoli di disuguaglianza (nessuna uguaglianza). */
public class HS116Test {

    // ---- Bounds dal MODE=1 --------------------------------------------------
    // Per i=1..10: XL(i)=0.1, ma XL(4)=1e-4 override; XL(9)=500; XL(11)=1; XL(12)=XL(13)=1e-4
    private static final double[] LB = {
        0.1,      // x1
        0.1,      // x2
        0.1,      // x3
        1e-4,     // x4
        0.1,      // x5
        0.1,      // x6
        0.1,      // x7
        0.1,      // x8
        500.0,    // x9
        0.1,      // x10
        1.0,      // x11
        1e-4,     // x12
        1e-4      // x13
    };

    // XU(1..3)=1; XU(7..9)=1000; XU(11..13)=150; XU(4)=0.1; XU(5)=0.9; XU(6)=0.9; XU(10)=500
    private static final double[] UB = {
        1.0,      // x1
        1.0,      // x2
        1.0,      // x3
        0.1,      // x4
        0.9,      // x5
        0.9,      // x6
        1000.0,   // x7
        1000.0,   // x8
        1000.0,   // x9
        500.0,    // x10
        150.0,    // x11
        150.0,    // x12
        150.0     // x13
    };

    // ---- Obiettivo: f(x) = x11 + x12 + x13 ---------------------------------
    private static class TP116Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 13; }
        @Override public double value(RealVector X) {
            return X.getEntry(10) + X.getEntry(11) + X.getEntry(12);
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    // ---- 15 disuguaglianze in forma g(x) >= 0 (ordine identico al Fortran) --
    private static class TP116Ineq extends InequalityConstraint {
        TP116Ineq() { super(new ArrayRealVector(new double[15])); }
        @Override public int dim() { return 13; }
        @Override public RealVector value(RealVector X) {
            final double x1=X.getEntry(0),  x2=X.getEntry(1),  x3=X.getEntry(2),
                         x4=X.getEntry(3),  x5=X.getEntry(4),  x6=X.getEntry(5),
                         x7=X.getEntry(6),  x8=X.getEntry(7),  x9=X.getEntry(8),
                         x10=X.getEntry(9), x11=X.getEntry(10), x12=X.getEntry(11),
                         x13=X.getEntry(12);

            final double[] g = new double[15];
            int k=0;

            // G(1)..G(5)
            g[k++] = x3 - x2;
            g[k++] = x2 - x1;
            g[k++] = 1.0 - 2.0e-3 * (x7 - x8);
            g[k++] = x11 + x12 + x13 - 50.0;
            g[k++] = 250.0 - x11 - x12 - x13;

            // G(6)..G(8)
            g[k++] = x13 - 1.262626 * x10 + 1.231059 * x3 * x10;
            g[k++] = x5 - 0.03475 * x2 - 0.975 * x2 * x5 + 9.75e-3 * x2 * x2;
            g[k++] = x6 - 0.03475 * x3 - 0.975 * x3 * x6 + 9.75e-3 * x3 * x3;

            // G(9)..G(12)
            g[k++] = x5 * x7 - x1 * x8 - x4 * x7 + x4 * x8;
            g[k++] = -2.0e-3 * (x2 * x9 + x5 * x8 - x1 * x8 - x6 * x9) - x6 - x5 + 1.0;
            g[k++] = x2 * x9 - x3 * x10 - x6 * x9 - 500.0 * (x2 - x6) + x2 * x10;
            g[k++] = x2 - 0.9 - 2.0e-3 * (x2 * x10 - x3 * x10);

            // G(13)..G(15)
            g[k++] = x4 - 0.03475 * x1 - 0.975 * x1 * x4 + 9.75e-3 * x1 * x1;
            g[k++] = x11 - 1.262626 * x8 + 1.231059 * x1 * x8;
            g[k++] = x12 - 1.262626 * x9 + 1.231059 * x2 * x9;

            return new ArrayRealVector(g);
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    // ---- Test ---------------------------------------------------------------
    @Test
    public void testHS116() {
        final InitialGuess guess = new InitialGuess(new double[]{
            0.5,   0.8,   0.9,   0.1,   0.14,  0.5,   489.0, 80.0, 650.0,
            450.0, 150.0, 150.0, 150.0
        });

        final SimpleBounds bounds = new SimpleBounds(LB, UB);

        final SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        optimizer.setDebugPrinter(System.out::println);

        final LagrangeSolution sol = optimizer.optimize(
            guess,
            new ObjectiveFunction(new TP116Obj()),
            new TP116Ineq(),
            bounds
        );

        // FEX = 0.975884089805D+02
        assertEquals(97.5884089805, sol.getValue(), 1e-2);
    }
}
