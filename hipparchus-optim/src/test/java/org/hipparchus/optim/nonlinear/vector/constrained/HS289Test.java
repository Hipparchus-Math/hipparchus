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

import org.hipparchus.linear.*;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS289Test {

    static final class HS289Obj extends TwiceDifferentiableFunction {
        private final int n;
        HS289Obj(int n) { this.n = n; }
        @Override public int dim() { return n; }

        @Override public double value(RealVector x) {
            double s2 = 0.0;
            for (int i = 0; i < n; i++) s2 += x.getEntry(i)*x.getEntry(i);
            return 1.0 - Math.exp(-s2 / 60.0);
        }

        @Override public RealVector gradient(RealVector x) {
            double s2 = 0.0;
            for (int i = 0; i < n; i++) s2 += x.getEntry(i)*x.getEntry(i);
            double fx = 1.0 - Math.exp(-s2 / 60.0);
            double coeff = (fx - 1.0) * (-2.0 / 60.0); // = (1 - fx) / 30
            double[] g = new double[n];
            for (int i = 0; i < n; i++) g[i] = coeff * x.getEntry(i);
            return new ArrayRealVector(g, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    private static double[] tp289Init(int n) {
        double[] x0 = new double[n];
        for (int i = 0; i < n; i++) {
            int I = i + 1;                         // Fortran 1-based
            double sign = (I % 2 == 0) ? 1.0 : -1.0; // (-1)^I
            x0[i] = sign * (1.0 + I / 30.0);
        }
        return x0;
    }

    @Test
    public void testHS289() {
        final int n = 30;
        double[] x0 = tp289Init(n);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        opt.setDebugPrinter(System.out::println);

        LagrangeSolution sol = opt.optimize(
            new InitialGuess(x0),
            new ObjectiveFunction(new HS289Obj(n))
        );

        double f = sol.getValue();
        double expected = 0.0; // FEX
        assertEquals(expected, f, 1.0e-6 * (Math.abs(expected) + 1.0), "objective mismatch");
    }
}
