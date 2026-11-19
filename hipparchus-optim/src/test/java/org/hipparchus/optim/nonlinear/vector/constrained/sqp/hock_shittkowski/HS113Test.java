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

public class HS113Test {

    private static class HS113Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 10; }
        @Override public double value(RealVector x) {
            return (((((((((((((FastMath.pow(x.getEntry(0), 2) + FastMath.pow(x.getEntry(1), 2)) + (x.getEntry(0) * x.getEntry(1))) - (14 * x.getEntry(0))) - (16 * x.getEntry(1))) + FastMath.pow((x.getEntry(2) - 10), 2)) + (4 * FastMath.pow((x.getEntry(3) - 5), 2))) + FastMath.pow((x.getEntry(4) - 3), 2)) + (2 * FastMath.pow((x.getEntry(5) - 1), 2))) + (5 * FastMath.pow(x.getEntry(6), 2))) + (7 * FastMath.pow((x.getEntry(7) - 11), 2))) + (2 * FastMath.pow((x.getEntry(8) - 10), 2))) + FastMath.pow((x.getEntry(9) - 7), 2)) + 45);
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS113Ineq extends InequalityConstraint {
        HS113Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (((((105 - (4 * x.getEntry(0))) - (5 * x.getEntry(1))) + (3 * x.getEntry(6))) - (9 * x.getEntry(7)))) - (0), ((((((-10) * x.getEntry(0)) + (8 * x.getEntry(1))) + (17 * x.getEntry(6))) - (2 * x.getEntry(7)))) - (0), ((((((8 * x.getEntry(0)) - (2 * x.getEntry(1))) - (5 * x.getEntry(8))) + (2 * x.getEntry(9))) + 12)) - (0), (((((((-3) * FastMath.pow((x.getEntry(0) - 2), 2)) - (4 * FastMath.pow((x.getEntry(1) - 3), 2))) - (2 * FastMath.pow(x.getEntry(2), 2))) + (7 * x.getEntry(3))) + 120)) - (0), (((((((-5) * FastMath.pow(x.getEntry(0), 2)) - (8 * x.getEntry(1))) - FastMath.pow((x.getEntry(2) - 6), 2)) + (2 * x.getEntry(3))) + 40)) - (0), (((((((-0.5) * FastMath.pow((x.getEntry(0) - 8), 2)) - (2 * FastMath.pow((x.getEntry(1) - 4), 2))) - (3 * FastMath.pow(x.getEntry(4), 2))) + x.getEntry(5)) + 30)) - (0), ((((((-FastMath.pow(x.getEntry(0), 2)) - (2 * FastMath.pow((x.getEntry(1) - 2), 2))) + ((2 * x.getEntry(0)) * x.getEntry(1))) - (14 * x.getEntry(4))) + (6 * x.getEntry(5)))) - (0), (((((3 * x.getEntry(0)) - (6 * x.getEntry(1))) - (12 * FastMath.pow((x.getEntry(8) - 8), 2))) + (7 * x.getEntry(9)))) - (0) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 10; }
    }

    @Test
    public void testHS113() {
        InitialGuess guess = new InitialGuess(new double[]{2, 3, 5, 5, 1, 2, 7, 3, 6, 10});
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        double val = 24.3062091;
        LagrangeSolution sol = optimizer.optimize(guess, new ObjectiveFunction(new HS113Obj()), new HS113Ineq());
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
