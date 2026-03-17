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
import org.hipparchus.util.FastMath;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS332Test {

    /** Objective function for HS332. */
    private static class HS332Obj extends TwiceDifferentiableFunction {
        private static final double PI = 3.1415926535d;
        @Override public int dim() { return 2; }
        @Override public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double PIM = PI / 36.0; // from Fortran: PI/.36D+1
            double fx = 0.0;
            for (int i = 1; i <= 100; i++) {
                double TR = PI * (1.0/3.0 + ( (i - 1.0) / 180.0 )); // PI*((1/3)+((i-1)/180))
                double A  = FastMath.log(TR);
                double B  = FastMath.sin(TR);
                double C  = FastMath.cos(TR);
                double XXX = (A + x2) * B + x1 * C;
                double YYY = (A + x2) * C - x1 * B;
                fx += PIM * (XXX*XXX + YYY*YYY);
            }
            return fx;
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    
    private static class HS332Ineq extends InequalityConstraint {
        private static final double PI = 3.1415926535d;
        HS332Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double PIM = 180.0 / PI; // from Fortran: .18D+3/PI
            double pbig = -360.0;          // -3.6D+3 in Fortran, any small start

            for (int i = 1; i <= 100; i++) {
                double TR = PI * (1.0/3.0 + ( (i - 1.0) / 180.0 ));
                double A  = 1.0 / TR - x1;
                double B  = FastMath.log(TR) + x2;
                double pang = PIM * FastMath.atan(FastMath.abs(A / B));
                if (pang > pbig) { pbig = pang; }
            }
            double g1 = 30.0 - pbig; // ≥ 0
            double g2 = pbig + 30.0; // ≥ 0
            return new ArrayRealVector(new double[]{ g1, g2 });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 2; }
    }

    @Test
    public void testHS332Bounds() {
      
        InitialGuess guess = new InitialGuess(new double[]{ 0.75, 0.75 });

        // Optimizer instance
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        // Bounds handled separately: 0 ≤ x_i ≤ 1.5
        SimpleBounds bounds = new SimpleBounds(
                new double[]{ 0.0, 0.0 },
                new double[]{ 1.5, 1.5 }
        );

        
       
        double expectedF = 11.495015;

        LagrangeSolution sol = optimizer.optimize(
                guess,
                new ObjectiveFunction(new HS332Obj()),
                new HS332Ineq(),
                bounds
        );
        
        
        HSProblemTestUtils.assertExpectedObjective(expectedF, sol);
    }
}
