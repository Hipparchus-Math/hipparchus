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
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;


public class HS312Test {

    
    static final class HS312Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return 2;
        }

        @Override
        public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            
            double A = x1 * x1 + 12.0 * x2 - 10.0;
            double B = 49.0 * (x1 * x1 + x2 * x2) + 84.0 * x1 + 2324.0 * x2 - 681.0;

            return A * A + B * B;
        }

        @Override
        public RealVector gradient(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            double A = x1 * x1 + 12.0 * x2 - 10.0;
            double B = 49.0 * (x1 * x1 + x2 * x2) + 84.0 * x1 + 2324.0 * x2 - 681.0;
            double g1 = 2.0 * (2.0 * x1 * A + B * (98.0 * x1 + 84.0));
            double g2 = 2.0 * (12.0 * A + B * (98.0 * x2 + 2324.0));

            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided for HS312");
        }
    }

    private LagrangeSolution solve() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();
       
        opt.setDebugPrinter(System.out::println);

       
        double[] x0 = {10.0, 10.0};

      
        return opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS312Obj())
        );
    }

    @Test
    public void testHS312() {
       
        final double fExpected = 0.0;
        LagrangeSolution sol = solve();
        double f = sol.getValue();
        assertEquals(fExpected, f, 1.0e-6 * (Math.abs(fExpected) + 1.0));
    }
}

