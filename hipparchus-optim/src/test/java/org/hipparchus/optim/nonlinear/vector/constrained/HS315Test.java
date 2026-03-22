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
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.linear.RealMatrix;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;


public class HS315Test {


    private static final int DIM = 2;


    static final class HS315Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return DIM; }

        @Override public double value(RealVector x) {

            return -x.getEntry(1);
        }

        @Override public RealVector gradient(RealVector x) {

            return new ArrayRealVector(new double[]{0.0, -1.0}, false);
        }

        @Override public RealMatrix hessian(RealVector x) {

            throw new UnsupportedOperationException("Hessian non fornita");
        }
    }


    static final class HS315Ineq extends InequalityConstraint {


        public HS315Ineq() { super(new ArrayRealVector(new double[3])); }

        @Override public int dim() { return DIM; }

        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double[] g = new double[3];

            // G1: 1 - 2*X2 + X1 >= 0
            g[0] = 1.0 - 2.0 * x2 + x1;

            // G2: X1^2 + X2^2 >= 0 (Vincolo irrilevante/sempre soddisfatto)
            g[1] = x1 * x1 + x2 * x2;

            // G3: 1 - X1^2 - X2^2 >= 0 (Unità cerchio)
            g[2] = 1.0 - x1 * x1 - x2 * x2;

            return new ArrayRealVector(g, false);
        }

        @Override public RealMatrix jacobian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double[][] J = new double[3][2];

            // G1 Jacobiana
            J[0][0] = 1.0;
            J[0][1] = -2.0;

            // G2 Jacobiana
            J[1][0] = 2.0 * x1;
            J[1][1] = 2.0 * x2;

            // G3 Jacobiana
            J[2][0] = -2.0 * x1;
            J[2][1] = -2.0 * x2;

            return org.hipparchus.linear.MatrixUtils.createRealMatrix(J);
        }
    }



    private static double[] start() {

        return new double[]{-0.1, -0.9};
    }

    @Test
    public void testHS315() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();
        opt.setDebugPrinter(System.out::println);

        final double fExpected = -0.8;

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                new ObjectiveFunction(new HS315Obj()),
                new HS315Ineq()
        );

        double f = sol.getValue();
        assertEquals(fExpected, f, 1.0e-3, "Objective mismatch");


    }
}
