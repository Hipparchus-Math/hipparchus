/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 * https://www.apache.org/licenses/LICENSE-2.0
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
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.MatrixUtils;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS319Test {

    private static final int DIM = 2;

    static final class HS319Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }
        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            return Math.pow(x1 - 20.0, 2) + Math.pow(x2 + 20.0, 2);
        }
        @Override public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double g1 = 2.0 * x1 - 40.0;
            double g2 = 2.0 * x2 + 40.0;
            return new ArrayRealVector(new double[]{g1, g2}, false);
        }
        @Override public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    static final class HS319Eq extends EqualityConstraint {

        HS319Eq() { super(new ArrayRealVector(new double[1])); }

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            // G(1)=0.01*X(1)^2 + X(2)^2/16 - 1 = 0
            double g1 = 0.01 * x1 * x1 + x2 * x2 / 16.0 - 1.0;
            return new ArrayRealVector(new double[]{g1}, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double[][] J = new double[1][DIM];

            // GG(1,1)=0.02*X(1), GG(1,2)=2*X(2)/16
            J[0][0] = 0.02 * x1;
            J[0][1] = 2.0 * x2 / 16.0;

            return MatrixUtils.createRealMatrix(J);
        }
    }

    private static double[] start() {
        return new double[]{0.0, 0.0};
    }

    @Test
    public void testHS319() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();

        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS319Obj()),
                new HS319Eq()
        );

        double f = sol.getValue();
        final double fExpected = 452.4044;

        assertEquals(fExpected, f, 1.0e-6 * (Math.abs(fExpected) + 1.0), "objective mismatch");

    }
}
