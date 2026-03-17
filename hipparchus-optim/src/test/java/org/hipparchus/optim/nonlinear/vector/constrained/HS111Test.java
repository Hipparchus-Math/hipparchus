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

/** HS TP111 (Schittkowski). 10 variabili, 3 vincoli di uguaglianza. */
public class HS111Test {

    private static final double[] C = {
        -6.089, -17.164, -34.054, -5.914, -24.721,
        -14.986, -24.1,  -10.708, -26.662, -22.179
    };

    private static final double[] LB = fill(-100.0, 10);
    private static final double[] UB = fill( 100.0, 10);

    private static double[] fill(double v, int n) {
        double[] a = new double[n];
        for (int i = 0; i < n; i++) a[i] = v;
        return a;
    }

    /** f(x) = sum_i exp(x_i) * (C_i + x_i - log(T)), con T = sum_j exp(x_j). */
    private static class TP111Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 10; }

        @Override public double value(RealVector X) {
            final int n = 10;
            double T = 0.0;
            for (int i = 0; i < n; i++) {
                T += FastMath.exp(X.getEntry(i));
            }
            final double logT = FastMath.log(T);

            double S = 0.0;
            for (int i = 0; i < n; i++) {
                final double xi = X.getEntry(i);
                S += FastMath.exp(xi) * (C[i] + xi - logT);
            }
            return S;
        }

        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    /** h(x) = [h1, h2, h3] = 0 accorpate. */
    private static class TP111EqAll extends EqualityConstraint {
        TP111EqAll() { super(new ArrayRealVector(new double[] { 0.0, 0.0, 0.0 })); }
        @Override public int dim() { return 10; }

        @Override public RealVector value(RealVector X) {
            final double x1=X.getEntry(0), x2=X.getEntry(1), x3=X.getEntry(2),
                         x4=X.getEntry(3), x5=X.getEntry(4), x6=X.getEntry(5),
                         x7=X.getEntry(6), x8=X.getEntry(7), x9=X.getEntry(8),
                         x10=X.getEntry(9);

            final double h1 = FastMath.exp(x1)
                            + 2.0 * FastMath.exp(x2)
                            + 2.0 * FastMath.exp(x3)
                            + FastMath.exp(x6)
                            + FastMath.exp(x10)
                            - 2.0;

            final double h2 = FastMath.exp(x4)
                            + 2.0 * FastMath.exp(x5)
                            + FastMath.exp(x6)
                            + FastMath.exp(x7)
                            - 1.0;

            final double h3 = FastMath.exp(x3)
                            + FastMath.exp(x7)
                            + FastMath.exp(x8)
                            + 2.0 * FastMath.exp(x9)
                            + FastMath.exp(x10)
                            - 1.0;

            return new ArrayRealVector(new double[] { h1, h2, h3 });
        }

        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testTP111() {
        final double[] x0 = new double[10];
        for (int i = 0; i < 10; i++) x0[i] = -2.3;

        final InitialGuess guess = new InitialGuess(x0);
        final SimpleBounds bounds = new SimpleBounds(LB, UB);
        SQPOption sqpOption = new SQPOption();
        sqpOption.setMaxLineSearchIteration(50);
        sqpOption.setEps(1e-7);
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        final LagrangeSolution sol = optimizer.optimize(
            guess,
            new ObjectiveFunction(new TP111Obj()),
            new TP111EqAll(),
            bounds
        );

        final double expected = -47.7610902637; // FEX Fortran
        HSProblemTestUtils.assertExpectedObjective(expected, sol);
    }
}
