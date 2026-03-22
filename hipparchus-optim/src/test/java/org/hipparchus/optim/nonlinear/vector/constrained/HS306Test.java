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
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS306Test {

    
    static final class HS306Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }

        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0), x2 = x.getEntry(1);
            double A  = Math.exp(-(x1 + x2));
            double B  = 2.0 * x1 * x1 + 3.0 * x2 * x2;
            return -A * B;
        }

        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0), x2 = x.getEntry(1);
            double A  = Math.exp(-(x1 + x2));
            double B  = 2.0 * x1 * x1 + 3.0 * x2 * x2;
            double g1 =  A * (B - 4.0 * x1); // d/dx1
            double g2 =  A * (B - 6.0 * x2); // d/dx2
            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override public org.hipparchus.linear.RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    private static double[] start() { return new double[]{1.0, 1.0}; } // X(1)=1, X(2)=1

//    @Test
//    public void testHS306() {
//        SQPOptimizerS2 opt = new SQPOptimizerS2();
//        opt.setDebugPrinter(System.out::println);
//        SQPOption sqpOption=new SQPOption();
//        sqpOption.setGradientMode(GradientMode.FORWARD);
//        LagrangeSolution sol = opt.optimize(
//                sqpOption,
//               new InitialGuess(start()),
//                new ObjectiveFunction(new HS306Obj())
//        );
//
//        double f = sol.getValue();
//        double fExpected = -1.1036; 
//        assertEquals(fExpected, f, 1.0e-6 * (Math.abs(fExpected) + 1.0), "objective mismatch");
//    }
}
