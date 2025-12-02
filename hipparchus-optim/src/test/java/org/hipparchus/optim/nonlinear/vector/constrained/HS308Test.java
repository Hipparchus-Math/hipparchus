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
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;


public class HS308Test {

    static final class HS308Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }

        @Override
        public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            double F1 = x1*x1 + x2*x2 + x1*x2;
            double F2 = Math.sin(x1);
            double F3 = Math.cos(x2);
            return F1*F1 + F2*F2 + F3*F3;
        }

        @Override
        public RealVector gradient(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            double F1 = x1*x1 + x2*x2 + x1*x2;
            double F2 = Math.sin(x1);
            double F3 = Math.cos(x2);

            // DF entries exactly as in the Fortran block
            double DF11 = 2.0*x1 + x2;
            double DF12 = 2.0*x2 + x1;
            double DF21 = 2.0*Math.sin(x1)*Math.cos(x1);
            double DF22 = 0.0;
            double DF31 = 0.0;
            double DF32 = -2.0*Math.cos(x2)*Math.sin(x2);

            double g1 = 2.0*F1*DF11 + 2.0*F2*DF21 + 2.0*F3*DF31;
            double g2 = 2.0*F1*DF12 + 2.0*F2*DF22 + 2.0*F3*DF32;

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

        // Start: X(1)=3.0, X(2)=1.0 
        double[] x0 = {3.0, 1.0};

        
        return opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS308Obj())
        );
    }

    @Test
    public void testHS308() {
        
        final double fExpected = 0.77319906;
        LagrangeSolution sol = solve();
        double f = sol.getValue();
        assertEquals(fExpected, f, 1.0e-6 * (Math.abs(fExpected) + 1.0), "objective mismatch");
    }
}
