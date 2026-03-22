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
import org.junit.jupiter.api.Test;

import static java.lang.Math.*;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS110Test {

    // Reference solution (all variables equal) and objective
    private static final double[] X_REF = {
            9.35025654733, 9.35025654733, 9.35025654733, 9.35025654733, 9.35025654733,
            9.35025654733, 9.35025654733, 9.35025654733, 9.35025654733, 9.35025654733
    };
    private static final double F_REF = -45.7784697153;

    // Bounds from the Fortran setup
    private static final double LO = 2.001;
    private static final double UP = 9.999;

    /** Objective f(x) and gradient; piecewise identical to TP110. */
    private static final class HS110Objective extends TwiceDifferentiableFunction {
        @Override public int dim() { return 10; }

        /** Product sign-power term S = sign(T) * |T|^0.2 with T=∏x_i. */
        private static double S(final RealVector x) {
            double T = 1.0;
            for (int i = 0; i < x.getDimension(); i++) {
                T *= x.getEntry(i);
            }
            return copySign(pow(abs(T), 0.2), T);
        }

        private static boolean inPenaltyBranch(final RealVector x) {
            for (int i = 0; i < x.getDimension(); i++) {
                final double xi = x.getEntry(i);
                if (xi <= 2.0 || xi >= 10.0) {
                    return true;
                }
            }
            return false;
        }

        @Override public double value(final RealVector x) {
            final double s = S(x);

            // Penalty branch when any xi is outside (2, 10) (Fortran uses <=2 or >=10)
            if (inPenaltyBranch(x)) {
                double sum = 0.0;
                for (int i = 0; i < x.getDimension(); i++) {
                    final double d = x.getEntry(i) - 5.0;
                    sum += d * d;
                }
                return sum + 1.0e3 - 45.8;
            }

            // Regular branch: U - S, where U = Σ [ ln(xi-2)^2 + ln(10-xi)^2 ]
            double U = 0.0;
            for (int i = 0; i < x.getDimension(); i++) {
                final double xi = x.getEntry(i);
                U += pow(log(xi - 2.0), 2) + pow(log(10.0 - xi), 2);
            }
            return U - s;
        }

        @Override public RealVector gradient(final RealVector x) {
            final int n = x.getDimension();
            final double[] g = new double[n];

            // Penalty branch gradient: 2*(xi - 5)
            if (inPenaltyBranch(x)) {
                for (int i = 0; i < n; i++) {
                    g[i] = 2.0 * (x.getEntry(i) - 5.0);
                }
                return new ArrayRealVector(g, false);
            }

            // Regular branch: dU/dxi - dS/dxi, with dS/dxi = 0.2*S/xi  (all xi>0 here)
            final double s = S(x);
            for (int i = 0; i < n; i++) {
                final double xi = x.getEntry(i);
                final double term =
                        2.0 * (log(xi - 2.0) / (xi - 2.0) - log(10.0 - xi) / (10.0 - xi))
                      - 0.2 * s / xi;
                g[i] = term;
            }
            return new ArrayRealVector(g, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    private static LagrangeSolution solve() {
        final double[] x0 = {9,9,9,9,9,9,9,9,9,9};
        final double[] lo = {LO,LO,LO,LO,LO,LO,LO,LO,LO,LO};
        final double[] up = {UP,UP,UP,UP,UP,UP,UP,UP,UP,UP};

        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        optimizer.setDebugPrinter(System.out::println);

        return optimizer.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS110Objective()),
                new SimpleBounds(lo, up)
        );
    }

    @Test
    public void testHS110() {
        final LagrangeSolution sol = solve();

        // Objective value
        assertEquals(F_REF, sol.getValue(), 1e-6 * (abs(F_REF) + 1.0), "objective mismatch");



    }
}
