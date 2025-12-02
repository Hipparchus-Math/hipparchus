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

/** HS TP71 (Schittkowski). f = x1*x4*(x1+x2+x3) + x3; 1 ineq (≥), 1 eq. */
public class HS071Test {

    // Bounds (MODE=1): 1 <= xi <= 5
    private static final double[] LB = { 1.0, 1.0, 1.0, 1.0 };
    private static final double[] UB = { 5.0, 5.0, 5.0, 5.0 };

    /** Obiettivo di TP71. */
    private static class HS071Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 4; }
        @Override public double value(RealVector X) {
            final double x1 = X.getEntry(0);
            final double x2 = X.getEntry(1);
            final double x3 = X.getEntry(2);
            final double x4 = X.getEntry(3);
            return x1 * x4 * (x1 + x2 + x3) + x3;
        }
        @Override public RealVector gradient(RealVector x) {
            throw new UnsupportedOperationException("Analytical gradient not provided");
        }
        @Override public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Analytical Hessian not provided");
        }
    }

    /** Vincolo di disuguaglianza (forma ≥): (x1*x2*x3*x4 - 25)/25 ≥ 0. */
    private static class HS071Ineq extends InequalityConstraint {
        HS071Ineq() { super(new ArrayRealVector(new double[] { 0.0 })); }
        @Override public int dim() { return 4; }
        @Override public RealVector value(RealVector X) {
            final double x1 = X.getEntry(0);
            final double x2 = X.getEntry(1);
            final double x3 = X.getEntry(2);
            final double x4 = X.getEntry(3);
            final double g = (x1 * x2 * x3 * x4 - 25.0) ; // ≥ 0
            return new ArrayRealVector(new double[] { g });
        }
        @Override public RealMatrix jacobian(RealVector x) {
            throw new UnsupportedOperationException("Analytical Jacobian not provided");
        }
    }

    /** Vincolo di uguaglianza: (sum(xi^2) - 40) = 0. */
    private static class HS071Eq extends EqualityConstraint {
        HS071Eq() { super(new ArrayRealVector(new double[] { 0.0 })); }
        @Override public int dim() { return 4; }
        @Override public RealVector value(RealVector X) {
            final double x1 = X.getEntry(0);
            final double x2 = X.getEntry(1);
            final double x3 = X.getEntry(2);
            final double x4 = X.getEntry(3);
            final double h = (x1*x1 + x2*x2 + x3*x3 + x4*x4 - 40.0) ; // = 0
            return new ArrayRealVector(new double[] { h });
        }
        @Override public RealMatrix jacobian(RealVector x) {
            throw new UnsupportedOperationException("Analytical Jacobian not provided");
        }
    }

    @Test
    public void testHS071() {
        final InitialGuess guess = new InitialGuess(new double[] { 1.0, 5.0, 5.0, 1.0 });
        final SimpleBounds bounds = new SimpleBounds(LB, UB);

        final SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        optimizer.setDebugPrinter(System.out::println);

        final LagrangeSolution sol = optimizer.optimize(
            guess,
            new ObjectiveFunction(new HS071Obj()),
            new HS071Eq(),     // equality
            new HS071Ineq(),   // inequality (≥)
            bounds
        );

        // FEX dal Fortran: 0.170140172895D+02
        final double expected = 17.0140172895;
        assertEquals(expected, sol.getValue(), 1e-6);
    }
}
