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


public class HS355Test {

    
    private static final double[] LB = { 0.1, 0.1, 0.0, 0.0 };
    private static final double BIG = Double.POSITIVE_INFINITY;
    private static final double[] UB = { BIG, BIG, BIG, BIG };

    
    private static final class Rs {
        final double r1, r2, r3, r4;
        Rs(double x1, double x2, double x3, double x4) {
            final double h0 = x2 * x4;
            r1 = 11.0 - (x1 + x2 - x3) * x4;
            r2 = x1 + 10.0 * x2 - x3 + x4 + h0 * (x3 - x1);
            r3 = 11.0 - (4.0 * x1 + 4.0 * x2 - x3) * x4;
            r4 = 2.0 * x1 + 20.0 * x2 - 0.5 * x3 + 2.0 * x4
               + 2.0 * h0 * (x3 - 4.0 * x1);
        }
    }

    /** f(x) = R1^2 + R2^2. */
    private static class TP355Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 4; }
        @Override public double value(RealVector X) {
            final double x1 = X.getEntry(0), x2 = X.getEntry(1),
                         x3 = X.getEntry(2), x4 = X.getEntry(3);
            final Rs r = new Rs(x1, x2, x3, x4);
            return r.r1 * r.r1 + r.r2 * r.r2;
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    
    private static class TP355Eq extends EqualityConstraint {
        TP355Eq() { super(new ArrayRealVector(new double[]{ 0.0 })); }
        @Override public int dim() { return 4; }
        @Override public RealVector value(RealVector X) {
            final double x1 = X.getEntry(0), x2 = X.getEntry(1),
                         x3 = X.getEntry(2), x4 = X.getEntry(3);
            final Rs r = new Rs(x1, x2, x3, x4);
            final double g = (r.r1*r.r1 + r.r2*r.r2) - (r.r3*r.r3 + r.r4*r.r4);
            return new ArrayRealVector(new double[]{ g });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testHS355() {
        
        // TP355 MODE=1: X(1)=0.1, X(2)=0.1, X(3)=0.0, X(4)=0.0
        final InitialGuess guess = new InitialGuess(new double[]{ 0.1, 0.1, 0.0, 0.0 });
        final SimpleBounds bounds = new SimpleBounds(LB, UB);

        final SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        SQPOption option=new SQPOption();
        option.setGradientMode(GradientMode.CENTRAL);
        final LagrangeSolution sol = opt.optimize(
                option,
            guess,
            new ObjectiveFunction(new TP355Obj()),
            new TP355Eq(),
           bounds
        );
        
        final double expected = 69.675463;
       
         HSProblemTestUtils.assertBetterObjective(expected, sol);
        
       
    }
}