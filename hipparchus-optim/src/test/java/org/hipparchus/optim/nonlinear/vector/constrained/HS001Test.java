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

import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS001Test {

    private static class HS001Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }
        @Override public double value(final RealVector x) {

            final double x0     = x.getEntry(0);
            final double x1     = x.getEntry(1);
            final double x1Mx02 = x1 - x0  * x0;
            final double oMx0   = 1 - x0;
            return 100 * x1Mx02 * x1Mx02 + oMx0 * oMx0;
        }
        @Override public RealVector gradient(final RealVector x) {
            final double x0     = x.getEntry(0);
            final double x1     = x.getEntry(1);
            final double x1Mx02 = x1 - x0  * x0;
            final double a      = 200 * x1Mx02;
            return MatrixUtils.createRealVector(new double[] {
                    2 * (x0 * (1 - a) - 1),
                    a
            });
        }
        @Override public RealMatrix hessian(final RealVector x) {
            throw new UnsupportedOperationException();
        }
    }
    
    @Test
    public void testHS001ExternalGradient() {
        doTestHS001(GradientMode.EXTERNAL);
    }

    @Test
    public void testHS001ForwardGradient() {
        doTestHS001(GradientMode.FORWARD);
    }

    @Test
    public void testHS001CentralGradient() {
        doTestHS001(GradientMode.CENTRAL);
    }

    private void doTestHS001(final GradientMode gradientMode) {
        SQPOption sqpOption=new SQPOption();
        sqpOption.setGradientMode(gradientMode);
        InitialGuess guess = new InitialGuess(new double[]{-2, 1});
        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            optimizer.setDebugPrinter(System.out::println);
        }
        double val = 0.0;
        LagrangeSolution sol = optimizer.optimize(sqpOption, guess, new ObjectiveFunction(new HS001Obj()));
        
        assertEquals(val, sol.getValue(), sqpOption.getEps()*10.0*(1.0+val));
    }
}
