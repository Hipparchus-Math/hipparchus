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
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;


public class HS311Test {

    static final class HS311Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return 2;
        }

        @Override
        public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            double t1 = x1 * x1 + x2 - 11.0;
            double t2 = x1 + x2 * x2 - 7.0;
            return t1 * t1 + t2 * t2;
        }

        @Override
        public RealVector gradient(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);


            double t1 = x1 * x1 + x2 - 11.0;
            double t2 = x1 + x2 * x2 - 7.0;


            double g1 = 4.0 * x1 * t1 + 2.0 * t2;
            double g2 = 2.0 * t1 + 4.0 * x2 * t2;

            return new ArrayRealVector(new double[] { g1, g2 }, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided for HS311");
        }
    }

    private LagrangeSolution solve() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();

        opt.setDebugPrinter(System.out::println);
        double[] x0 = { 1.0, 1.0 };


        return opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS311Obj())
        );
    }

    @Test
    public void testHS311() {

        final double fExpected = 0.0;
        LagrangeSolution sol = solve();
        double f = sol.getValue();
        assertEquals(fExpected, f, 1.0e-6 * (Math.abs(fExpected) + 1.0),
                     "objective mismatch for HS311");
    }
}
