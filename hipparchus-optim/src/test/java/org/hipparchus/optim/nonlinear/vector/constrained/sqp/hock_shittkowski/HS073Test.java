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

public class HS073Test {

    private static class HS073Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 4; }
        @Override public double value(RealVector x) {
            return ((((24.55 * x.getEntry(0)) + (26.75 * x.getEntry(1))) + (39 * x.getEntry(2))) + (40.5 * x.getEntry(3)));
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS073Eq extends EqualityConstraint {
        HS073Eq() { super(new ArrayRealVector(new double[]{ 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ ((((x.getEntry(0) + x.getEntry(1)) + x.getEntry(2)) + x.getEntry(3))) - (1.0) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 4; }
    }

    private static class HS073Ineq extends InequalityConstraint {
        HS073Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0 ,0.0,0.0,0.0,0.0})); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (((((2.3 * x.getEntry(0)) + (5.6 * x.getEntry(1))) + (11.1 * x.getEntry(2))) + (1.3 * x.getEntry(3)))) - (5), (((((12 * x.getEntry(0)) + (11.9 * x.getEntry(1))) + (41.8 * x.getEntry(2))) + (52.1 * x.getEntry(3)))) - ((21 + (1.645 * FastMath.sqrt(((((0.28 * FastMath.pow(x.getEntry(0), 2)) + (0.19 * FastMath.pow(x.getEntry(1), 2))) + (20.5 * FastMath.pow(x.getEntry(2), 2))) + (0.62 * FastMath.pow(x.getEntry(3), 2))))))),x.getEntry(0),x.getEntry(1),x.getEntry(2),x.getEntry(3) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 4; }
    }

    @Test
    public void testHS073() {
        InitialGuess guess = new InitialGuess(new double[]{1.0, 1.0, 1.0, 1.0});
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        double val = 29.894378;
        SQPOption sqpOption=new SQPOption();
        sqpOption.setEps(1.0e-9);
        LagrangeSolution sol = optimizer.optimize(sqpOption,guess, new ObjectiveFunction(new HS073Obj()), new HS073Eq(), new HS073Ineq());
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
