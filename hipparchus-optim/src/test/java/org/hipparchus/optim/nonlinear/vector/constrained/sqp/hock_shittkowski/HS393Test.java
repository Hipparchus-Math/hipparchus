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
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS393Test {

  
    private static class TP393Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 48; }

        @Override public double value(RealVector Xv) {
            final double[] x = Xv.toArray();
            double E = 0.0;

            
            for (int i = 0; i < 12; i++) {
                final double c = 1.0 - x[i];
                E += 10.0 * c * c;
            }
            
            for (int i = 24; i <= 35; i++) {
                final double c = x[i] - 1.0;
                E += 1000.0 * (0.1 + 2.0 * c * (c + Math.sqrt(0.1 + c * c))) / 4.0;
            }
           
            for (int i = 36; i <= 41; i++) {
                final double c = x[i] - 1.0;
                E += 2000.0 * (0.1 + 2.0 * c * (c + Math.sqrt(0.1 + c * c))) / 4.0;
            }
            
            for (int i = 42; i <= 47; i++) {
                E += 100.0 * x[i];
            }
            return E / 1000.0;
        }

        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    
    private static class TP393Ineq extends InequalityConstraint {
        TP393Ineq() { super(new ArrayRealVector(new double[]{ 0.0 })); } // g(x) >= 0
        @Override public int dim() { return 48; }

        @Override public RealVector value(RealVector Xv) {
            return new ArrayRealVector(new double[]{ phiTP393B(Xv.toArray()) });
        }

        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    
    private static double phiTP393B(double[] X) {
        final double[] A = {
            0.9, 0.8, 1.1, 1.0, 0.7, 1.1, 1.0, 1.0, 1.1,
            0.9, 0.8, 1.2, 0.9, 1.2, 1.2, 1.0, 1.0, 0.9
        };
        final double[] U = new double[18];

        
        for (int i = 1; i <= 6; i++) {
            final int K1 = i + 24;  // -> X[24..29]
            final int K2 = i + 42;  // -> X[42..47]
            final int K3 = i + 12;  // -> X[12..17]
            final double alp = X[K1-1]*X[K1-1] * A[i-1] * 2.0*X[K2-1]/(1.0+X[K2-1]) * X[K3-1];
            final double xi  = X[i-1];
            U[i-1] = xi*xi / (xi + alp);
        }

        
        for (int i = 7; i <= 12; i++) {
            final int K1 = i + 24;  // -> X[30..35]
            final int K2 = i + 36;  // -> X[42..47] (sì, condivisi)
            final int K3 = i + 12;  // -> X[18..23]
            final double alp = X[K1-1]*X[K1-1] * A[i-1] * 2.0*X[K2-1]/(1.0+X[K2-1]) * X[K3-1];
            final double sum = X[i-1] + U[i-7];
            U[i-1] = (sum*sum) / (sum + alp);
        }

        
        for (int i = 13; i <= 15; i++) {
            final int K1 = 2*(i-10)+1;  // usa U[K1], U[K1+1]
            final int K2 = i + 24;      // -> X[36..38]
            final double alp = X[K2-1]*X[K2-1] * A[i-1];
            final double sum = U[K1-1] + U[K1];
            U[i-1] = (sum*sum) / (sum + alp);
        }

        
        for (int i = 16; i <= 18; i++) {
            final int K2 = i + 24;      // -> X[39..41]
            final double alp = X[K2-1]*X[K2-1] * A[i-1];
            final double sum = U[i-4];
            U[i-1] = (sum*sum) / (sum + alp);
        }

        final double R = U[15] + U[16] + U[17];
        return 1.5 - R;
    }

    
    private static class TP393Eq extends EqualityConstraint {
        TP393Eq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0 })); }
        @Override public int dim() { return 48; }

        @Override public RealVector value(RealVector Xv) {
            final double[] x = Xv.toArray();
            double s1 = 0.0, s2 = 0.0;
            for (int i = 0; i < 12; i++) s1 += x[i];
            for (int i = 12; i < 24; i++) s2 += x[i]; 
            return new ArrayRealVector(new double[]{ 12.0 - s1, 12.0 - s2 });
        }

        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testTP393() {
        
        final double[] x0 = new double[48];
        for (int i = 0; i < 24; i++) x0[i] = 1.0;
        for (int i = 24; i < 30; i++) x0[i] = 1.3;
        for (int i = 30; i < 48; i++) x0[i] = 1.0;

        
        final double[] lb = new double[48];
        final double[] ub = new double[48];
        for (int i = 0; i < 48; i++) {
            lb[i] = 0.002;
            ub[i] = (i < 24) ? 2.0 : Double.POSITIVE_INFINITY;
        }

        final InitialGuess guess = new InitialGuess(x0);
        final SimpleBounds bounds = new SimpleBounds(lb, ub);

        final SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        final LagrangeSolution sol = opt.optimize(
            guess,
            new ObjectiveFunction(new TP393Obj()),
            new TP393Ineq(),
            new TP393Eq(),
            bounds
        );

        final double expected = 0.86337998; 
        HSProblemTestUtils.assertExpectedObjective(expected, sol);
    }
}
