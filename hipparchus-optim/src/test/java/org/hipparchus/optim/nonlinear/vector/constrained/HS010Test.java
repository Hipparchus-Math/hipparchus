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

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;


public class HS010Test {

    private static class HS010Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }
        @Override public double value(RealVector x) {
            return (x.getEntry(0) - x.getEntry(1));
        }
        @Override public RealVector gradient(RealVector x) { 
           return  new ArrayRealVector(new double[]{1.0,-1.0}); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS010Ineq extends InequalityConstraint {
        HS010Ineq() { super(new ArrayRealVector(new double[]{ 0.0 })); }
        @Override public RealVector value(RealVector x) {
            double x1=x.getEntry(0);
            double x2=x.getEntry(1);
            return new ArrayRealVector(new double[]{(-3.0*x1*x1+2.0*x1*x2-x2*x2+1.0 )});
        }
        @Override public RealMatrix jacobian(RealVector x) { 
            double x1=x.getEntry(0);
            double x2=x.getEntry(1);
            double g1=-6.0*x1+2.0*x2;      
            double g2=2.0*(x1-x2);  
            RealMatrix J = new Array2DRowRealMatrix(1, 2);
            J.setEntry(0, 0, g1);
            J.setEntry(0, 1, g2);
            return  J;}
        
        @Override public int dim() { return 2; }
    }

    @Test
    public void testHS010() {
        InitialGuess guess = new InitialGuess(new double[]{-10.0, 10.0});
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        SQPOption option=new SQPOption();
        option.setGradientMode(GradientMode.EXTERNAL);
        double val = -1.0;
        LagrangeSolution sol = optimizer.optimize(option,guess, new ObjectiveFunction(new HS010Obj()), new HS010Ineq());
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
