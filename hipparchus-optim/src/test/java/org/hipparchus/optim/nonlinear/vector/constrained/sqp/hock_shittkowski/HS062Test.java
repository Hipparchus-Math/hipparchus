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
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS062Test {

    private static class HS062Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 3; }
        @Override public double value(RealVector x) {
            return ((-32.174) * (((255.0 * FastMath.log(((((x.getEntry(0) + x.getEntry(1)) + x.getEntry(2)) + 0.03) / ((((0.09 * x.getEntry(0)) + x.getEntry(1)) + x.getEntry(2)) + 0.03)))) + (280.0 * FastMath.log((((x.getEntry(1) + x.getEntry(2)) + 0.03) / (((0.07 * x.getEntry(1)) + x.getEntry(2)) + 0.03))))) + (290.0 * FastMath.log(((x.getEntry(2) + 0.03) / ((0.13 * x.getEntry(2)) + 0.03))))));
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS062Eq extends EqualityConstraint {
        HS062Eq() { super(new ArrayRealVector(new double[]{ 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (((x.getEntry(0) + x.getEntry(1)) + x.getEntry(2))) - (1.0) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 3; }
    }

    private static class HS062Ineq extends InequalityConstraint {
        HS062Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (x.getEntry(0)) - (0.0), (1.0) - (x.getEntry(0)), (x.getEntry(1)) - (0.0), (1.0) - (x.getEntry(1)), (x.getEntry(2)) - (0.0), (1.0) - (x.getEntry(2)) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 3; }
    }

    @Test
    public void testHS062() {
        InitialGuess guess = new InitialGuess(new double[]{0.7, 0.2, 0.1});
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        double val = -26272.51448;
        SQPOption sqpOption=new SQPOption();
       
        LagrangeSolution sol = optimizer.optimize(sqpOption,guess, new ObjectiveFunction(new HS062Obj()), new HS062Eq(), new HS062Ineq());
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
