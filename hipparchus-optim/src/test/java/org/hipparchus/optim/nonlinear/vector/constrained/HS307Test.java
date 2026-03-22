/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
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


public class HS307Test {

    static final class HS307Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }

        @Override
        public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            double fx = 0.0;
            for (int i = 1; i <= 10; i++) {
                double wi = i;
                double xi1 = wi * x1;
                if (xi1 > 20.0) xi1 = 0.0;            // Fortran: IF (XI1.GT.20) XI1=0
                double xi2 = wi * x2;
                if (xi2 > 20.0) xi2 = 0.0;            // Fortran: IF (XI2.GT.20) XI2=0
                double Fi = 2.0 + 2.0 * wi - Math.exp(xi1) - Math.exp(xi2);
                fx += Fi * Fi;
            }
            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            double g1 = 0.0, g2 = 0.0;
            for (int i = 1; i <= 10; i++) {
                double wi = i;
                double xi1 = wi * x1;
                if (xi1 > 20.0) xi1 = 0.0;            // same guard used for exp in gradient
                double xi2 = wi * x2;
                if (xi2 > 20.0) xi2 = 0.0;
                double e1 = Math.exp(xi1);
                double e2 = Math.exp(xi2);
                double Fi = 2.0 + 2.0 * wi - e1 - e2;

                // DF/ Dx1 = - wi * exp(xi1); DF/Dx2 = - wi * exp(xi2)
                double dFdx1 = -wi * e1;
                double dFdx2 = -wi * e2;

                // grad of sum Fi^2 is 2 * Fi * grad(Fi)
                g1 += 2.0 * Fi * dFdx1;
                g2 += 2.0 * Fi * dFdx2;
            }
            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    private LagrangeSolution solve() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();
        opt.setDebugPrinter(System.out::println); // richiesto

        final double[] x0 = {0.3, 0.4};
        final double[] lo = {0.0, 0.0};
        final double[] up = {1.0e10, 1.0e10};

        return opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS307Obj()),
                new SimpleBounds(lo, up)
        );
    }

    @Test
    public void testHS307() {
        final double fExpected = 0.12436e3;
        LagrangeSolution sol = solve();
        double f = sol.getValue();
        assertEquals(fExpected, f, 1.0e-2 * (Math.abs(fExpected) + 1.0), "objective mismatch");
    }
}
