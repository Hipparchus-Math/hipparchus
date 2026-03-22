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

public class HS042Test {

    private static class HS042Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 4; }
        @Override public double value(RealVector x) {
            return (((FastMath.pow((x.getEntry(0) - 1), 2) + FastMath.pow((x.getEntry(1) - 2), 2)) + FastMath.pow((x.getEntry(2) - 3), 2)) + FastMath.pow((x.getEntry(3) - 4), 2));
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS042Eq extends EqualityConstraint {
        HS042Eq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (x.getEntry(0)) - (2), ((FastMath.pow(x.getEntry(2), 2) + FastMath.pow(x.getEntry(3), 2))) - (2) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 4; }
    }

    @Test
    public void testHS042() {
        InitialGuess guess = new InitialGuess(new double[]{1, 1, 1, 1});
        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            optimizer.setDebugPrinter(System.out::println);
        }
        double val = (28 - (10 * FastMath.sqrt(2)));
        LagrangeSolution sol = optimizer.optimize(guess, new ObjectiveFunction(new HS042Obj()), new HS042Eq());
        assertEquals(val, sol.getValue(), 1e-6);
    }
}
