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
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS058Test {

    private static class HS058Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }
        @Override public double value(RealVector x) {
            return ((100.0 * FastMath.pow((x.getEntry(1) - FastMath.pow(x.getEntry(0), 2)), 2)) + FastMath.pow((1 - x.getEntry(0)), 2));
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS058Ineq extends InequalityConstraint {
        HS058Ineq() { super(new ArrayRealVector(new double[]{  0.0, 0.0,0.0,0.0,0.0  })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{
            FastMath.pow(x.getEntry(1), 2) - x.getEntry(0),                      // x2^2 - x1 ≥ 0
            FastMath.pow(x.getEntry(0), 2) - x.getEntry(1),                      // x1^2 - x2 ≥ 0
            FastMath.pow(x.getEntry(0), 2) + FastMath.pow(x.getEntry(1), 2) - 1, // x1^2 + x2^2 - 1 ≥ 0
            x.getEntry(0) - (-2.0),                                              // x1 ≥ -2
            0.5 - x.getEntry(0)                                                  // x1 ≤ 0.5
        });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 2; }
    }

    @Test
    public void testHS058() {
         SQPOption sqpOption=new SQPOption();
        InitialGuess guess = new InitialGuess(new double[]{-2.0, 1.0});
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        double val = 3.19033354957;
        LagrangeSolution sol = optimizer.optimize(sqpOption,guess, new ObjectiveFunction(new HS058Obj()), new HS058Ineq());
        HSProblemTestUtils.assertExpectedObjective(val, sol);
        
    }
}
