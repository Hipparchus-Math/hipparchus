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

public class HS034Test {

    private static class HS034Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 3; }
        @Override public double value(RealVector x) {
            return (-x.getEntry(0));
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS034Ineq extends InequalityConstraint {
        HS034Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (x.getEntry(1)) - (FastMath.exp(x.getEntry(0))), (x.getEntry(2)) - (FastMath.exp(x.getEntry(1))), (100.0) - (x.getEntry(0)), (100.0) - (x.getEntry(1)), (10.0) - (x.getEntry(2)), (x.getEntry(0)) - (0), (x.getEntry(1)) - (0), (x.getEntry(2)) - (0) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 3; }
    }

    private static class HS034IneqNoBound extends InequalityConstraint {
        HS034IneqNoBound() { super(new ArrayRealVector(new double[]{ 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (x.getEntry(1)) - (FastMath.exp(x.getEntry(0))), (x.getEntry(2)) - (FastMath.exp(x.getEntry(1)))   });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 3; }
    }

//    @Test
//    public void testHS034() {
//        InitialGuess guess = new InitialGuess(new double[]{0, 1.05, 2.9});
//        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
//        double val = (-FastMath.log(FastMath.log(10)));
//        LagrangeSolution sol = optimizer.optimize(guess, new ObjectiveFunction(new HS034Obj()), new HS034Ineq());
//        assertEquals(val, sol.getValue(), 1e-6);
//    }
//
     @Test
    public void testHS034Bound() {
        InitialGuess guess = new InitialGuess(new double[]{0, 1.05, 2.9});
        SimpleBounds bounds = new SimpleBounds(
                new double[]{ 0.0, 0.0,0.0 },
                new double[]{ 100.0,100.0,10.0 }
        );
        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            optimizer.setDebugPrinter(System.out::println);
        }
        double val = (-FastMath.log(FastMath.log(10)));
        LagrangeSolution sol = optimizer.optimize(bounds,guess, new ObjectiveFunction(new HS034Obj()), new HS034Ineq());
        assertEquals(val, sol.getValue(), 1e-6);
    }
}
