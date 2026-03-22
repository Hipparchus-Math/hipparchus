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


public class HS293Test {


    private static final int N = 50;


    private static class TP293Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return N; }

        @Override public double value(RealVector X) {
            double s = 0.0;
            for (int i = 0; i < N; i++) {
                double xi = X.getEntry(i);
                s += (i + 1) * xi * xi;
            }
            return s * s;
        }

        @Override public RealVector gradient(RealVector X) {
            double s = 0.0;
            for (int i = 0; i < N; i++) s += (i + 1) * X.getEntry(i) * X.getEntry(i);
            double[] g = new double[N];
            for (int i = 0; i < N; i++) g[i] = 4.0 * (i + 1) * s * X.getEntry(i);
            return new ArrayRealVector(g);
        }

        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testHS293() {

        double[] x0 = new double[N];
        for (int i = 0; i < N; i++) x0[i] = 1.0;


        final InitialGuess guess = new InitialGuess(x0);


        final SQPOptimizerS2 opt = new SQPOptimizerS2();
        opt.setDebugPrinter(System.out::println);

        final LagrangeSolution sol = opt.optimize(
            guess,
            new ObjectiveFunction(new TP293Obj())
        );

        assertEquals(0.0, sol.getValue(), 1e-3);
    }
}
