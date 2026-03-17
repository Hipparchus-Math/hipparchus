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
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS013Test {

    private static class HS013Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }
        @Override public double value(RealVector x) {
            double x1=x.getEntry(0);
             double x2=x.getEntry(1);
            return (x1-2.0)*(x1-2.0) + x2*x2;
        }
        @Override public RealVector gradient(RealVector x) {
             double x1=x.getEntry(0);
             double x2=x.getEntry(1);
             double g1=2.0*(x1-2.0);
             double g2=2.0*x2;
            return new ArrayRealVector(new double[]{g1,g2})
                    
                    ; }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS013Ineq extends InequalityConstraint {
        HS013Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0 })); }
        @Override public RealVector value(RealVector x) {
            return new ArrayRealVector(new double[]{((FastMath.pow((1.0 - x.getEntry(0)), 3.0)) - (x.getEntry(1))), (x.getEntry(0)) - (0.0), (x.getEntry(1)) - (0.0) });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 2; }
    }
    
    private static class HS013IneqNoBound extends InequalityConstraint {
        HS013IneqNoBound() { super(new ArrayRealVector(new double[]{ 0.0})); }
        @Override public RealVector value(RealVector x) {
             double x1=x.getEntry(0);
             double x2=x.getEntry(1);
            return new ArrayRealVector(new double[]{(1.0-x1)*(1.0-x1)*(1.0-x1) - x2});
        }
        @Override public RealMatrix jacobian(RealVector x) { 
             final double x1 = x.getEntry(0);
            final double  x2 = x.getEntry(1);
            final RealMatrix J = new Array2DRowRealMatrix(1, 2);
            double g1=-3.0*(1.0-x1)*(1.0-x1);
            double g2=-1.0;
            J.setEntry(0, 0, g1);
            J.setEntry(0, 1, g2);
           return J; }
        @Override public int dim() { return 2; }
    }

   
    
     @Test
    public void testHS013() {
         SQPOption sqpOption=new SQPOption();
          
        sqpOption.setGradientMode(GradientMode.EXTERNAL);
        
        
        InitialGuess guess = new InitialGuess(new double[]{0.01, 0.01});
        SimpleBounds bounds=new SimpleBounds(new double[]{0.0,0.0},
                                             new double[]{Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY});
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        double val = 1.0;
        LagrangeSolution sol = optimizer.optimize(sqpOption,guess, new ObjectiveFunction(new HS013Obj()), new HS013IneqNoBound (),bounds);
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
