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

import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS299Test {

    private static class HS299Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return 100;
        }

        @Override
        public double value(RealVector x) {
            int n = dim();
            double fx = 0.0;
            for (int i = 0; i < n - 1; i++) {
                double xi = x.getEntry(i);
                double xip1 = x.getEntry(i + 1);
                double fi = 10.0 * (xip1 - xi * xi);
                double gi = 1.0 - xi;
                fx += fi * fi + gi * gi;
            }
            return fx * 1.0e-4;
        }

        @Override
        public RealVector gradient(RealVector x) {
            throw new UnsupportedOperationException();
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            return null;
        }
    }

    @Test
    public void testHS299() {
        SQPOption sqpOption = new SQPOption();
        sqpOption.setGradientMode(GradientMode.CENTRAL);
        double[] start = new double[100];
        for (int i = 0; i < 100; i++) {
            start[i] = -1.2;
        }
        for (int i = 1; i < 100; i += 2) {
            start[i] = 1.0;
        }

        InitialGuess guess = new InitialGuess(start);
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        

        double val = 0.0;
        LagrangeSolution sol = optimizer.optimize(sqpOption, guess, new ObjectiveFunction(new HS299Obj()));

        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
