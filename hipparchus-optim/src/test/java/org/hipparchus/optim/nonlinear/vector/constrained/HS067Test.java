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
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

public class HS067Test {

    private static class HS067Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 10; }
        @Override public double value(RealVector x) {
            return (-((((((0.063 * x.getEntry(3)) * x.getEntry(6)) - (5.04 * x.getEntry(0))) - (3.36 * x.getEntry(4))) - (0.035 * x.getEntry(1))) - (10 * x.getEntry(2))));
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS067Eq extends EqualityConstraint {
        HS067Eq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (x.getEntry(4)) - (((1.22 * x.getEntry(3)) - x.getEntry(0))), (x.getEntry(7)) - (((x.getEntry(1) + x.getEntry(4)) / x.getEntry(0))), (x.getEntry(3)) - (((0.01 * x.getEntry(0)) * ((112 + (13.167 * x.getEntry(7))) - (0.6667 * FastMath.pow(x.getEntry(7), 2))))), (x.getEntry(6)) - ((((86.35 + (1.098 * x.getEntry(7))) - (0.038 * FastMath.pow(x.getEntry(7), 2))) + (0.325 * (x.getEntry(5) - 89)))), (x.getEntry(9)) - (((3 * x.getEntry(6)) - 133)), (x.getEntry(8)) - ((35.82 - (0.222 * x.getEntry(9)))), (x.getEntry(5)) - (((98000 * x.getEntry(2)) / ((x.getEntry(3) * x.getEntry(8)) + (1000 * x.getEntry(2))))) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 10; }
    }
    
     // Bounds from MODE=1
    private static final double[] LB = { 1e-5, 1e-5, 1e-5, 0,0,85,90,3,0.01,145 };
    private static final double[] UB = { 2000,16000,120, 5000,2000,93,95,12,4,162 };
    private static class HS067Ineq extends InequalityConstraint {
        HS067Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0,0,0,0,0,0,0,0,0,0,0,0,0,0})); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{ (x.getEntry(0)) - (1e-5), (2000) - (x.getEntry(0)), 
                                                     (x.getEntry(1)) - (1e-5), (16000) - (x.getEntry(1)), 
                                                     (x.getEntry(2)) - (1e-5), (120) - (x.getEntry(2)),
                                                     (x.getEntry(3)) - (0), (5000) - (x.getEntry(3)), 
                                                     (x.getEntry(4)) - (0), (2000) - (x.getEntry(4)), 
                                                     (x.getEntry(5)) - (85), (93) - (x.getEntry(5)), 
                                                     (x.getEntry(6)) - (90), (95) - (x.getEntry(6)), 
                                                     (x.getEntry(7)) - (3), (12) - (x.getEntry(7)),
                                                      (x.getEntry(8)) - (0.01), (4) - (x.getEntry(8)),
                                                       (x.getEntry(9)) - (145), (162) - (x.getEntry(9))});
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 10; }
    }
//    @Disabled
//    @Test
//    public void testHS067() {
//        InitialGuess guess = new InitialGuess(new double[]{1745, 12000, 110,0,0,0,0,0,0,0});
//        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
//       
//        double val = -1162.02698006;
//        LagrangeSolution sol = optimizer.optimize(guess, new ObjectiveFunction(new HS067Obj()), new HS067Eq(), new HS067Ineq());
//        HSProblemTestUtils.assertExpectedObjective(val, sol);
//    }
    
    @Test
    public void testHS067Bounds() {
        
        final SimpleBounds bounds = new SimpleBounds(LB, UB);
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
       
        double val = -1162.02698006;
        LagrangeSolution sol = optimizer.optimize( 
                new ObjectiveFunction(new HS067Obj()),
                new HS067Eq(), 
                bounds);
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
