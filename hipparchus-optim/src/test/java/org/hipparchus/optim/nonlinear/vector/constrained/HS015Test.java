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

public class HS015Test {

    private static class HS015Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }
        @Override public double value(RealVector x) {
            return ((100 * FastMath.pow((x.getEntry(1) - FastMath.pow(x.getEntry(0), 2)), 2)) + FastMath.pow((1 - x.getEntry(0)), 2));
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS015Ineq extends InequalityConstraint {
        HS015Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ ((x.getEntry(0) * x.getEntry(1))) - (1), ((x.getEntry(0) + FastMath.pow(x.getEntry(1), 2))) - (0), ((1.0 / 2.0)) - (x.getEntry(0)) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 2; }
    }

    private static class HS015IneqNoBounds extends InequalityConstraint {
        HS015IneqNoBounds() { super(new ArrayRealVector(new double[]{ 0.0, 0.0})); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ ((x.getEntry(0) * x.getEntry(1))) - (1), ((x.getEntry(0) + FastMath.pow(x.getEntry(1), 2))) - (0) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 2; }
    }

    @Test
    public void testHS015() {
        InitialGuess guess = new InitialGuess(new double[]{-2, 1});
        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            optimizer.setDebugPrinter(System.out::println);
        }
        double val = 306.5;
        LagrangeSolution sol = optimizer.optimize(guess, new ObjectiveFunction(new HS015Obj()), new HS015Ineq());
        assertEquals(val, sol.getValue(), 1e-1);
    }

    @Test
    public void testHS015Bound() {
        InitialGuess guess = new InitialGuess(new double[]{-2, 0.4});
        SimpleBounds bounds=new SimpleBounds(new double[]{Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY},
                                             new double[]{0.5,Double.POSITIVE_INFINITY});
        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            optimizer.setDebugPrinter(System.out::println);
        }
        double val = 306.5;
        LagrangeSolution sol = optimizer.optimize(guess, new ObjectiveFunction(new HS015Obj()), new HS015IneqNoBounds(),bounds);
        assertEquals(val, sol.getValue(), 1e-3);
    }
}
