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
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS100Test {

    private static class HS100Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 7; }
        @Override public double value(RealVector x) {
            return  (((((((((FastMath.pow((x.getEntry(0) - 10), 2) + (5 * FastMath.pow((x.getEntry(1) - 12), 2))) + FastMath.pow(x.getEntry(2), 
4)) + (3 * FastMath.pow((x.getEntry(3) - 11), 2))) + (10 * FastMath.pow(x.getEntry(4), 6))) + (7 * FastMath.pow(x.getEntry(5), 2))) + FastMath.pow(x.getEntry(6), 4)) - ((4 * x.getEntry(5)) * x.getEntry(6))) - (10 * x.getEntry(5))) - (8 * x.getEntry(6)));
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS100Ineq extends InequalityConstraint {
        HS100Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (127) - ((((((2 * FastMath.pow(x.getEntry(0), 2)) + (3 * FastMath.pow(x.getEntry(1), 4))) + x.getEntry(2)) + (4 * FastMath.pow(x.getEntry(3), 2))) + (5 * x.getEntry(4)))), (282) - ((((((7 * x.getEntry(0)) + (3 * x.getEntry(1))) + (10 * FastMath.pow(x.getEntry(2), 2))) + x.getEntry(3)) - x.getEntry(4))), (196) - (((((23 * x.getEntry(0)) + FastMath.pow(x.getEntry(1), 2)) + (6 * FastMath.pow(x.getEntry(5), 2))) - (8 * x.getEntry(6)))), ((((((((-4) * FastMath.pow(x.getEntry(0), 2)) - FastMath.pow(x.getEntry(1), 2)) + ((3 * x.getEntry(0)) * x.getEntry(1))) - (2 * FastMath.pow(x.getEntry(2), 2))) - (5 * x.getEntry(5))) + (11 * x.getEntry(6)))) - (0) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 7; }
    }

    @Test
    public void testHS100() {
        InitialGuess guess = new InitialGuess(new double[]{1, 2, 0, 4, 0, 1, 1});
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        double val = 680.6300573;
        LagrangeSolution sol = optimizer.optimize(guess, new ObjectiveFunction(new HS100Obj()), new HS100Ineq());
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
