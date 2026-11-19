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
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS064Test {

    private static class HS064Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 3; }
        @Override public double value(RealVector x) {
            return ((((((5.0 * x.getEntry(0)) + (50000.0 / x.getEntry(0))) + (20.0 * x.getEntry(1))) + (72000.0 / x.getEntry(1))) + (10.0 * x.getEntry(2))) + (144000.0 / x.getEntry(2)));
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS064Ineq extends InequalityConstraint {
        HS064Ineq() { super(new ArrayRealVector(new double[]{ 0.0})); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (1.0) - ((((4.0 / x.getEntry(0)) + (32.0 / x.getEntry(1))) + (120.0 / x.getEntry(2)))) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 3; }
    }

    @Test
    public void testHS064() {
        SQPOption sqpOption=new SQPOption();
        sqpOption.setGradientMode(GradientMode.FORWARD);

 final double[] LB = { 1e-5, 1e-5, 1e-5 };
  final double[] UB = { Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        final SimpleBounds bounds = new SimpleBounds(LB, UB);
        InitialGuess guess = new InitialGuess(new double[]{1.0 ,1.0, 1.0});
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        double val = 6299.842428;
        LagrangeSolution sol = optimizer.optimize(sqpOption,bounds,guess, new ObjectiveFunction(new HS064Obj()), new HS064Ineq());
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
