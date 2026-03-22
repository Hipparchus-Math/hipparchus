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

public class HS020Test {

    private static class HS020Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }
        @Override public double value(RealVector x) {
            return ((100.0 * FastMath.pow((x.getEntry(1) - FastMath.pow(x.getEntry(0), 2)), 2)) + FastMath.pow((1 - x.getEntry(0)), 2));
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS020Ineq extends InequalityConstraint {
        HS020Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0, 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ ((x.getEntry(0) + FastMath.pow(x.getEntry(1), 2))) - (0), ((FastMath.pow(x.getEntry(0), 2) + x.getEntry(1))) - (0), ((FastMath.pow(x.getEntry(0), 2) + FastMath.pow(x.getEntry(1), 2.0))) - (1.0), (x.getEntry(0)) - (((-1.0) / 2.0)), ((1.0 / 2.0)) - (x.getEntry(0)) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 2; }
    }
    
    private static class HS020Ineq1 extends InequalityConstraint {
        HS020Ineq1() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0})); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ ((x.getEntry(0) + FastMath.pow(x.getEntry(1), 2))) - (0), ((FastMath.pow(x.getEntry(0), 2) + x.getEntry(1))) - (0), ((FastMath.pow(x.getEntry(0), 2) + FastMath.pow(x.getEntry(1), 2.0))) - (1.0) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 2; }
    }

    @Test
    public void testHS020() {
        SQPOption sqpOption=new SQPOption();
        sqpOption.setMaxLineSearchIteration(20);
        InitialGuess guess = new InitialGuess(new double[]{0.1, 1.0});
        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            optimizer.setDebugPrinter(System.out::println);
        }
        double val = 81.5-25*FastMath.sqrt(3.0);
        LagrangeSolution sol = optimizer.optimize(sqpOption,guess, new ObjectiveFunction(new HS020Obj()), new HS020Ineq());
        assertEquals(val, sol.getValue(), 1e-3);
    }
    double[]lb=new double[]{-1.0/2.0,Double.NEGATIVE_INFINITY};
    double[]ub=new double[]{1.0/2.0,Double.POSITIVE_INFINITY};
    
    @Test
    public void testHS020Bound() {
        
        InitialGuess guess = new InitialGuess(new double[]{0.1, 1.0});
        SQPOption sqpOption=new SQPOption();
        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            optimizer.setDebugPrinter(System.out::println);
        }
        double val = 81.5-25*FastMath.sqrt(3.0);  //(81.5D0-25.D0*DSQRT(3.D0))
        LagrangeSolution sol = optimizer.optimize(sqpOption,guess, new ObjectiveFunction(new HS020Obj()), new HS020Ineq1(),new SimpleBounds(lb, ub));
        assertEquals(val, sol.getValue(), 1e-3);
    }
}
