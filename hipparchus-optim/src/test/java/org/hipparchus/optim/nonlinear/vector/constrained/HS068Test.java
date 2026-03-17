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
import org.hipparchus.special.Erf;
import org.hipparchus.util.FastMath;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS068Test {

    // Fortran TP68 (KN1=1): A(1)=1e-4, B(1)=1, D(1)=1, Z(1)=24
    private static final double a = 1.0e-4;
    private static final double b = 1.0;
    private static final double d = 1.0;      // D(1)
    private static final double n = 24.0;     // Z(1)

    // Bounds from MODE=1
    private static final double[] lb = { 1.0e-4, 0.0, 0.0, 1.0e-4 };
    private static final double[] ub = { 1.0,    100.0, 2.0, 2.0 };

    /** Φ(z): CDF normale standard, equivalente a MDNORD in Fortran. */
    private static double phi(double z) {
        return 0.5 * (1.0 + Erf.erf(z / FastMath.sqrt(2.0)));
    }

    /** f(x) come in TP68, con clamp su x1 e uso di expm1 per stabilità. */
    private static class HS068Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 4; }

        @Override public double value(RealVector X) {
            final double[] x = X.toArray();
            // X1 = min(max(1e-8, x[0]), 10) original
//            final double x1 = FastMath.min(10.0, FastMath.max(1.0e-8, x[0]));
            // X1 without clamp
            final double x1 = x[0];
            final double x3 = x[2];
            final double x4 = x[3];

            // v = exp(x1) - 1
            final double v = FastMath.expm1(x1);

            // FX = (a*n - x4 * (b*v - x3) / (v + x4)) / x1
            return (a * n - x4 * (b * v - x3) / (v + x4)) / x1;
        }

        @Override public RealVector gradient(RealVector x) {
            throw new UnsupportedOperationException("Analytical gradient not provided");
        }

        @Override public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Analytical Hessian not provided");
        }
    }

    /** g(x) = 0 (2 vincoli) come in TP68 usando Φ al posto di MDNORD. */
    private static class HS068Eq extends EqualityConstraint {
        HS068Eq() {
            // target = [0, 0]
            super(new ArrayRealVector(new double[] { 0.0, 0.0 }));
        }

        @Override public RealVector value(RealVector X) {
            final double[] x = X.toArray();
            final double x2 = x[1];
            final double x3 = x[2];
            final double x4 = x[3];
            final double rtN = FastMath.sqrt(n);

            // g1 = x3 - 2*Phi(-x2)
            final double g1 = x3 - 2.0 * phi(-x2);

            // g2 = x4 - Phi(-x2 + sqrt(n)) - Phi(-x2 - sqrt(n))
            final double g2 = x4 - phi(-x2 + rtN) - phi(-x2 - rtN);

            return new ArrayRealVector(new double[] { g1, g2 });
        }

        @Override public RealMatrix jacobian(RealVector x) {
            throw new UnsupportedOperationException("Analytical Jacobian not provided");
        }

        @Override public int dim() { return 4; }
    }

    @Test
    public void testHS068() {
        final InitialGuess guess = new InitialGuess(new double[] { 1, 1, 1, 1 });
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        final double expected = -0.920425020704; // FEX dal Fortran
        final SimpleBounds bounds = new SimpleBounds(lb, ub);

        final LagrangeSolution sol =
            optimizer.optimize(guess, new ObjectiveFunction(new HS068Obj()), new HS068Eq(), bounds);

        HSProblemTestUtils.assertExpectedObjective(expected, sol);
    }
}
