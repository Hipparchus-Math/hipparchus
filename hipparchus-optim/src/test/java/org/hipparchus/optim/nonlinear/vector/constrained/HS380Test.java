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
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;


public class HS380Test {

   
    private static final double[] A = {
        -0.00133172, -0.002270927, -0.00248546,
        -4.67, -4.671973, -0.00814, -0.008092,
        -0.005, -0.000909, -0.00088, -0.00119
    };

  
    private static final double[] C = {
        5.367373e-2,  2.1863746e-2, 9.7733533e-2, 6.6940803e-3,
        1.0e-6,       1.0e-5,       1.0e-6,       1.0e-10,     1.0e-8,   1.0e-2,
        1.0e-4,       1.0898645e-1, 1.6108052e-4, 1.0e-23,     1.9304541e-6, 1.0e-3,
        1.0e-6,       1.0e-5,       1.0e-6,       1.0e-9,      1.0e-9,   1.0e-3,
        1.0e-3,       1.0898645e-1, 1.6108052e-5, 1.0e-23,     1.9304541e-8, 1.0e-5,
        1.1184059e-4, 1.0e-4
    };

   

    private static final int N = 12;
    private static final double[] LB, UB;

    static {
        LB = new double[N];
        UB = new double[N];
        for (int i = 0; i < N; i++) {
            LB[i] = 0.1;
            UB[i] = 100.0;
        }
    }

    
    private static class TP380Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return N; }

        @Override public double value(RealVector X) {
            double fx = 1.0; // 0.1D+1
            for (int i = 0; i < 11; i++) {
                double xi = X.getEntry(i);
                if (xi < 1e-14) xi = 1e-14;
                fx *= FastMath.pow(xi, A[i]);
            }
            return 1.0e5 * fx;
        }

        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    
    private static class TP380Ineq extends InequalityConstraint {
        TP380Ineq() { super(new ArrayRealVector(new double[]{0,0,0})); }

        @Override public int dim() { return N; }

        @Override public RealVector value(RealVector X) {
            final double x1 = X.getEntry(0),  x2 = X.getEntry(1),  x3 = X.getEntry(2),
                         x4 = X.getEntry(3),  x5 = X.getEntry(4),  x6 = X.getEntry(5),
                         x7 = X.getEntry(6),  x8 = X.getEntry(7),  x9 = X.getEntry(8),
                         x10= X.getEntry(9),  x11= X.getEntry(10), x12= X.getEntry(11);

            final double g1 = 1.0
                    - C[0]*x1 - C[1]*x2 - C[2]*x3 - C[3]*x4*x5;

            final double g2 = 1.0
                    - C[4]*x1 - C[5]*x2 - C[6]*x3
                    - C[7]*x4*x12
                    - C[8]*x5 / x12
                    - C[9]*x6 / x12
                    - C[10]*x7 * x12
                    - C[11]*x4 * x5
                    - C[12]*x2 * x5 / x12
                    - C[13]*x2 * x4 * x5
                    - C[14]*(x2 / x4) * x5 / (x12*x12)
                    - C[15]*x10 / x12;

            final double g3 = 1.0
                    - C[16]*x1 - C[17]*x2 - C[18]*x3
                    - C[19]*x4 - C[20]*x5 - C[21]*x6
                    - C[22]*x8
                    - C[23]*x4 * x5
                    - C[24]*x2 * x5
                    - C[25]*x2 * x4 * x5
                    - C[26]*x2 * x5 / x4
                    - C[27]*x9
                    - C[28]*x1 * x9
                    - C[29]*x11;

            return new ArrayRealVector(new double[]{ g1, g2, g3 });
        }

        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    

    @Test
    public void testHS380() {
        final double[] x0 = new double[N];
        for (int i = 0; i < N; i++) x0[i] = 4.0; // 0.4D+1

        final InitialGuess guess = new InitialGuess(x0);
        final SimpleBounds bounds = new SimpleBounds(LB, UB);
         
        final SQPOptimizerS2 opt = new SQPOptimizerS2();
         if (Boolean.getBoolean("hipparchus.debug.sqp")) {
          opt.setDebugPrinter(System.out::println);
          }

        final LagrangeSolution sol = opt.optimize(
            guess,
            new ObjectiveFunction(new TP380Obj()),
            new TP380Ineq(),
            bounds
            
        );

        
        final double expected = 3.1682215;
        assertEquals(expected, sol.getValue(), 1e-2);
    }
}
