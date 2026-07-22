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
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS017Test {

    private static class HS017Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }
        @Override public double value(RealVector x) {
            return ((100 * FastMath.pow((x.getEntry(1) - FastMath.pow(x.getEntry(0), 2)), 2)) + FastMath.pow((1 - x.getEntry(0)), 2));
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS017Ineq extends InequalityConstraint {
        HS017Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0, 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (((-x.getEntry(0)) + FastMath.pow(x.getEntry(1), 2))) - (0), ((FastMath.pow(x.getEntry(0), 2) - x.getEntry(1))) - (0), (x.getEntry(0)) - (((-1.0) / 2)), ((1.0 / 2)) - (x.getEntry(0)), (1) - (x.getEntry(1)) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 2; }
    }
    
    private static class HS017IneqNoBounds extends InequalityConstraint {
        HS017IneqNoBounds() { super(new ArrayRealVector(new double[]{ 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (((-x.getEntry(0)) + FastMath.pow(x.getEntry(1), 2))) - (0), ((FastMath.pow(x.getEntry(0), 2) - x.getEntry(1))) - (0) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 2; }
    }

    
    @Test
    public void testHS017Bound() {
        InitialGuess guess = new InitialGuess(new double[]{-2, 1});
        SimpleBounds bounds=new SimpleBounds(new double[]{-2.0,Double.NEGATIVE_INFINITY},
                                             new double[]{0.5,1.0});
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        SQPOption sqpOption=new SQPOption();
        sqpOption.setGradientMode(GradientMode.CENTRAL);
        double val = 1.0;
        LagrangeSolution sol = optimizer.optimize(sqpOption,guess, new ObjectiveFunction(new HS017Obj()), new HS017IneqNoBounds(),bounds);
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
