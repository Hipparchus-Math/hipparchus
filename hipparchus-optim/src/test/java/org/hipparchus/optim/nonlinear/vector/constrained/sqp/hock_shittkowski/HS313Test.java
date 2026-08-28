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
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS313Test {

    
    static final class HS313Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }

        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            return 0.001 * Math.pow(x1 - 3.0, 2) - (x2 - x1) + Math.exp(20.0 * (x2 - x1));
        }

        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);         
            double xh = 20.0 * Math.exp(20.0 * (x2 - x1));
            double g1 = 1.0 + 0.0002 * (x1 - 3.0) - xh;
            double g2 = xh - 1.0;
            
            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override public org.hipparchus.linear.RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    private static double[] start() { 
       
        return new double[]{0.0, -1.0}; 
    }

    @Test
    public void testHS313() {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer() ;

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS313Obj())
        );

        double f = sol.getValue();
       
        final double fExpected = 0.199786;
        
      
        assertEquals(fExpected, f, 1.0e-5, "objective mismatch");
        
        
       
    }
}
