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
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;


public class HS005Test {

    private static class HS005Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }
        @Override public double value(RealVector x) {
            return ((((FastMath.sin((x.getEntry(0) + x.getEntry(1))) + FastMath.pow((x.getEntry(0) - x.getEntry(1)), 2)) - (1.5 * x.getEntry(0))) + (2.5 * x.getEntry(1))) + 1);
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS005Ineq extends InequalityConstraint {
        HS005Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (x.getEntry(0)) - (-1.5), (4) - (x.getEntry(0)), (x.getEntry(1)) - (-3), (3) - (x.getEntry(1)) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 2; }
    }

    @Test
    public void testHS005() {
        InitialGuess guess = new InitialGuess(new double[]{0, 0});
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        double val = -1.91322207;
        LagrangeSolution sol = optimizer.optimize(guess, new ObjectiveFunction(new HS005Obj()), new HS005Ineq());
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
    
    @Test
    public void testHS005Bounds() {
        InitialGuess guess = new InitialGuess(new double[]{0, 0});
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        double val = -1.91322207;
         SimpleBounds bounds=new SimpleBounds(new double[]{-1.5,-3},
                                             new double[]{4.0,3.0});
        
        LagrangeSolution sol = optimizer.optimize(guess, new ObjectiveFunction(new HS005Obj()), bounds);
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
