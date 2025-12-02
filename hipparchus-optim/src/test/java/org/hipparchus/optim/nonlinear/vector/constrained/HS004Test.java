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


public class HS004Test {

    private static class HS004Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }
        @Override public double value(RealVector x) {
            return ((FastMath.pow((x.getEntry(0) + 1), 3) / 3) + x.getEntry(1));
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS004Ineq extends InequalityConstraint {
        HS004Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (x.getEntry(0)) - (1), (x.getEntry(1)) - (0) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 2; }
    }
     

    @Test
    public void testHS004() {
        InitialGuess guess = new InitialGuess(new double[]{1.125, 0.125});
        SQPOption sqpOption=new SQPOption();
        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            optimizer.setDebugPrinter(System.out::println);
        }
        double val = (8.0 / 3.0);
        LagrangeSolution sol = optimizer.optimize(sqpOption,guess, new ObjectiveFunction(new HS004Obj()), new HS004Ineq());
        assertEquals(val, sol.getValue(), sqpOption.getEps()*10.0*(1.0+val));
    }
    
    @Test
    public void testHS004Bounds() {
        InitialGuess guess = new InitialGuess(new double[]{1.125, 0.125});
        SQPOption sqpOption = new SQPOption();
        sqpOption.setMaxLineSearchIteration(50);
        sqpOption.setB(0.5);
        sqpOption.setMu(1.0e-4);
        sqpOption.setEps(1e-7);
        
        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            optimizer.setDebugPrinter(System.out::println);
        }
        optimizer.setDebugPrinter(System.out::println);
        double val = (8.0 / 3.0);
        SimpleBounds bounds=new SimpleBounds(new double[]{1.0,0.0},
                                             new double[]{Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY});
        LagrangeSolution sol = optimizer.optimize(sqpOption, guess,new ObjectiveFunction(new HS004Obj()), bounds);
        assertEquals(val, sol.getValue(), sqpOption.getEps()*10.0*(1.0+val));
    }
}
