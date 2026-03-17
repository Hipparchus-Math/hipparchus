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

public class HS104Test {

    private static class HS104Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 8; }
        @Override public double value(RealVector x) {
            return ((((((0.4 * FastMath.pow(x.getEntry(0), 0.67)) * FastMath.pow(x.getEntry(6), (-0.67))) + ((0.4 * FastMath.pow(x.getEntry(1), 0.67)) * FastMath.pow(x.getEntry(7), (-0.67)))) + 10) - x.getEntry(0)) - x.getEntry(1));
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS104Ineq extends InequalityConstraint {
        HS104Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (((1 - ((0.0588 * x.getEntry(4)) * x.getEntry(6))) - (0.1 * x.getEntry(0)))) - (0), ((((1 - ((0.0588 * x.getEntry(5)) * x.getEntry(7))) - (0.1 * x.getEntry(0))) - (0.1 * x.getEntry(1)))) - (0), ((((1 - ((4 * x.getEntry(2)) / x.getEntry(4))) - (2 / (FastMath.pow(x.getEntry(2), 0.71) * x.getEntry(4)))) - ((0.0588 * x.getEntry(6)) / FastMath.pow(x.getEntry(2), 1.3)))) - (0), ((((1 - ((4 * x.getEntry(3)) / x.getEntry(5))) - (2 / (FastMath.pow(x.getEntry(3), 0.71) * x.getEntry(5)))) - ((0.0588 * x.getEntry(7)) / FastMath.pow(x.getEntry(3), 1.3)))) - (0), (((((((0.4 * FastMath.pow(x.getEntry(0), 0.67)) * FastMath.pow(x.getEntry(6), (-0.67))) + ((0.4 * FastMath.pow(x.getEntry(1), 0.67)) * FastMath.pow(x.getEntry(7), (-0.67)))) + 10) - x.getEntry(0)) - x.getEntry(1))) - (0.1), (4.2) - (((((((0.4 * FastMath.pow(x.getEntry(0), 0.67)) * FastMath.pow(x.getEntry(6), (-0.67))) + ((0.4 * FastMath.pow(x.getEntry(1), 0.67)) * FastMath.pow(x.getEntry(7), (-0.67)))) + 10) - x.getEntry(0)) - x.getEntry(1))), (10) - (x.getEntry(0)), (10) - (x.getEntry(1)), (10) - (x.getEntry(2)), (10) - (x.getEntry(3)), (10) - (x.getEntry(4)), (10) - (x.getEntry(5)), (10) - (x.getEntry(6)), (10) - (x.getEntry(7)) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 8; }
    }

    @Test
    public void testHS104() {
        InitialGuess guess = new InitialGuess(new double[]{6, 3, 0.4, 0.2, 6, 6, 1, 0.5});
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        double val = 3.9511634396;
        LagrangeSolution sol = optimizer.optimize(guess, new ObjectiveFunction(new HS104Obj()), new HS104Ineq());
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
