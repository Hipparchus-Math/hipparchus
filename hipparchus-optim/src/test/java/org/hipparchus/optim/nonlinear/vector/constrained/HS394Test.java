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

public class HS394Test {


    private static class TP394Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 20; }

        @Override public double value(RealVector X) {
            double fx = 0.0;
            for (int i = 0; i < 20; i++) {
                final double xi = X.getEntry(i);
                final double ii = i + 1.0;
                fx += ii * (xi*xi + xi*xi*xi*xi);
            }
            return fx;
        }

        @Override public RealVector gradient(RealVector X) {
            final double[] g = new double[20];
            for (int i = 0; i < 20; i++) {
                final double xi = X.getEntry(i);
                final double ii = i + 1.0;
                g[i] = ii * (2.0*xi + 4.0*xi*xi*xi);
            }
            return new ArrayRealVector(g);
        }

        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }


    private static class TP394Eq extends EqualityConstraint {
        TP394Eq() { super(new ArrayRealVector(new double[] { 0.0 })); }

        @Override public int dim() { return 20; }

        @Override public RealVector value(RealVector X) {
            double s = 0.0;
            for (int i = 0; i < 20; i++) {
                final double xi = X.getEntry(i);
                s += xi*xi;
            }
            return new ArrayRealVector(new double[]{ s - 1.0 });
        }

        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testTP394() {

        final double[] x0 = new double[20];
        for (int i = 0; i < 20; i++) x0[i] = 2.0;

        // Nessun bound
        final double[] lb = new double[20];
        final double[] ub = new double[20];
        for (int i = 0; i < 20; i++) {
            lb[i] = Double.NEGATIVE_INFINITY;
            ub[i] = Double.POSITIVE_INFINITY;
        }

        final SQPOptimizerS2 opt = new SQPOptimizerS2();
         if (Boolean.getBoolean("hipparchus.debug.sqp")) {
          opt.setDebugPrinter(System.out::println);
          }

        final LagrangeSolution sol = opt.optimize(
            new InitialGuess(x0),
            new ObjectiveFunction(new TP394Obj()),
            new TP394Eq(),
            new SimpleBounds(lb, ub)
        );

       
        assertEquals(1.9166667, sol.getValue(), 1e-4);
    }
}
