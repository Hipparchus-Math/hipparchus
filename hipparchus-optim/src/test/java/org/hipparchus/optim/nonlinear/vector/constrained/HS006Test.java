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
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.MaxIter;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;


public class HS006Test {

    private static class HS006Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }
        @Override public double value(RealVector x) {
            final double x0     = x.getEntry(0);
            final double oMx0   = 1 - x0;
            return oMx0 * oMx0;
        }
        @Override public RealVector gradient(RealVector x) {
            final double x0     = x.getEntry(0);
            return MatrixUtils.createRealVector(new double[] {
                    2 * (x0 - 1), 0
            });
        }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static class HS006Eq extends EqualityConstraint {
        HS006Eq() { super(new ArrayRealVector(new double[]{ 0.0 })); }
        @Override public RealVector value(RealVector x) {
            final double x0     = x.getEntry(0);
            final double x1     = x.getEntry(1);
            return new ArrayRealVector(new double[] { (10 * (x1 - x0 * x0)) - (0) });
        }
        @Override public RealMatrix jacobian(RealVector x) {
            final double x0 = x.getEntry(0);
            return MatrixUtils.createRealMatrix(new double[][] {
                    { -20 * x0, 10 }
            });
        }
        @Override public int dim() { return 2; }
    }

    @Test
    public void testHS006ExternalGradient() {
        doTestHS006(GradientMode.EXTERNAL);
    }

    @Test
    public void testHS006ForwardGradient() {
        doTestHS006(GradientMode.FORWARD);
    }

    @Test
    public void testHS006CentralGradient() {
        doTestHS006(GradientMode.CENTRAL);
    }

    private void doTestHS006(final GradientMode gradientMode) {
        SQPOption sqpOption = new SQPOption();

        sqpOption.setGradientMode(gradientMode);
        InitialGuess guess = new InitialGuess(new double[]{-1.2, 1});
        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            optimizer.setDebugPrinter(System.out::println);
        }
        double val = 0.0;
        LagrangeSolution sol = optimizer.optimize(sqpOption, guess, new MaxIter(100),
                                                  new ObjectiveFunction(new HS006Obj()), new HS006Eq());
        assertEquals(val, sol.getValue(), sqpOption.getEps()*10.0*(1.0+val));
    }
}
