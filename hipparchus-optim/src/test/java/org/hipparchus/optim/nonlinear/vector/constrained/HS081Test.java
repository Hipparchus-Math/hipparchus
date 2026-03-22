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
import org.hipparchus.util.FastMath;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

/** HS TP81. 5 vars, 3 nonlinear equalities, no inequalities. */
public class HS081Test {

    // Bounds: x1..x2 in [-2.3, 2.3]; x3..x5 in [-3.2, 3.2]
    private static final double[] LB = { -2.3, -2.3, -3.2, -3.2, -3.2 };
    private static final double[] UB = {  2.3,  2.3,  3.2,  3.2,  3.2 };

    /** f(x) = exp(x1*x2*x3*x4*x5) - 0.5*(x1^3 + x2^3 + 1)^2. */
    private static class TP81Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 5; }
        @Override public double value(RealVector X) {
            final double x1 = X.getEntry(0);
            final double x2 = X.getEntry(1);
            final double x3 = X.getEntry(2);
            final double x4 = X.getEntry(3);
            final double x5 = X.getEntry(4);
            final double prod = x1 * x2 * x3 * x4 * x5;
            final double v1 = x1*x1*x1 + x2*x2*x2 + 1.0;
            return FastMath.exp(prod) - 0.5 * v1 * v1;
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    /** Equalities (MODE=4): h(x) = 0 with h1, h2, h3. */
    private static class TP81Eq extends EqualityConstraint {
        TP81Eq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0 })); }
        @Override public int dim() { return 5; }
        @Override public RealVector value(RealVector X) {
            final double x1=X.getEntry(0), x2=X.getEntry(1), x3=X.getEntry(2),
                         x4=X.getEntry(3), x5=X.getEntry(4);
            final double h1 = x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 - 10.0;
            final double h2 = x2*x3 - 5.0*x4*x5;
            final double h3 = x1*x1*x1 + x2*x2*x2 + 1.0;
            return new ArrayRealVector(new double[]{ h1, h2, h3 });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testHS081() {
        final InitialGuess guess = new InitialGuess(new double[]{ -2.0,  2.0,  2.0, -1.0, -1.0 });
        final SimpleBounds bounds = new SimpleBounds(LB, UB);

        final SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        optimizer.setDebugPrinter(System.out::println);

        final LagrangeSolution sol = optimizer.optimize(
                guess,
                new ObjectiveFunction(new TP81Obj()),
                new TP81Eq(),
                bounds
        );

        // FEX (TP81): 0.539498477749D-01
        assertEquals(0.0539498477749, sol.getValue(), 1e-4);
    }
}
