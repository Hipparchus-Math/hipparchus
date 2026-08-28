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

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS395Test {

    private static class HS395Obj extends TwiceDifferentiableFunction {
        @Override
        public int dim() { return 50; }

        @Override
        public double value(RealVector x) {
            double sum = 0.0;
            //FX=FX+DBLE(I)*(X(I)**2+X(I)**4)
            for (int i = 0; i < x.getDimension(); i++) {
                double xi = x.getEntry(i);
                double idx = 1.0*i + 1.0; // 1-based index
                sum += idx * (xi * xi + FastMath.pow(xi, 4));
            }
            return sum;
        }

       
        @Override
        public RealVector gradient(RealVector x) {
            final double[] g = new double[50];
            for (int i = 0; i < x.getDimension(); i++) {
                final double xi = x.getEntry(i);
                final double idx = i + 1.0;
                g[i] = idx * (2.0 * xi + 4.0 * xi * xi * xi);
            }
            return new ArrayRealVector(g, false);
        }

        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS395Eq extends EqualityConstraint {
        HS395Eq() {
            super(new ArrayRealVector(new double[]{0.0}));
        }

        @Override
        public RealVector value(RealVector x) {
            double sum = 0.0;
            for (int i = 0; i < x.getDimension(); i++) {
                sum += FastMath.pow(x.getEntry(i), 2);
            }
            return new ArrayRealVector(new double[]{sum - 1.0});
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            final double[][] j = new double[1][50];
            for (int i = 0; i < x.getDimension(); i++) {
                j[0][i] = 2.0 * x.getEntry(i);
            }
            return new Array2DRowRealMatrix(j, false);
        }

        @Override public int dim() { return 50; }
    }

    @Test
    public void testHS395() {
        double[] start = new double[50];
        for (int i = 0; i < 50; i++) {
            start[i] = 2.0;
        }
        
        InitialGuess guess = new InitialGuess(start);
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        SQPOption option=new SQPOption();
        option.setGradientMode(GradientMode.FORWARD);
        double val = 1.9166668;

        LagrangeSolution sol = optimizer.optimize(
                option,
           guess,
            new ObjectiveFunction(new HS395Obj()),
            new HS395Eq()
        );

        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
  
}
